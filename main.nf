nextflow.enable.dsl=2

//Process to take in ONT sequencing log, find all individual barcode files, concatenate and rename with isolate name at desired location

process CONCATENATE_FASTQ {
    tag "$meta.id"
    publishDir "${params.outdir}/combined_fastqs", mode: 'copy'

    input:
    tuple val(meta), val(isolate), val(barcode), val(parent_dir), val(run_id), path(fastq_files)

    output:
    path("${run_id}_${barcode}_${isolate}_combined.fastq.gz"), emit: combined_fastq

    script:
    """
    #Concatenate all fastq files for the run_id/barcode/isolate combination 
    cat ${fastq_files} > ${run_id}_${barcode}_${isolate}_combined.fastq.gz
    """
}

//runs seqkit stats on the combined fastq files to generate report used to predict genome coverage and filter for downstream processing.
process SEQKIT_STATS {
    publishDir "${params.outdir}/seqkit_stats", mode: 'copy'

    conda "${moduleDir}/envs/seqkit.yml"
    
    container {
    workflow.containerEngine == 'singularity' && !(task.ext?.singularity_pull_docker_container ?: false) ?
        'oras://community.wave.seqera.io/library/seqkit:2.12.0--ec0d76090cceee7c' :
        'community.wave.seqera.io/library/seqkit:2.12.0--5e60eb93d3a59212'
    }

    input:
    path(combined_fastq)

    output:
    path("seqkit_stats_report.tsv"), emit: seqkit_stats_report

    script:
    """
    #Run seqkit stats on the combined fastq file and save the output
    seqkit stats -T ${combined_fastq} > seqkit_stats_report.tsv
    """
}



process PROCESS_STATS_WITH_R {

    label "process_single"
    publishDir "${params.outdir}/processed_stats", mode: 'copy'

    conda "${moduleDir}/envs/r_processing.yml"
    container { 
        workflow.containerEngine == 'singularity' && !(task.ext.singularity_pull_docker_container ?: false) ?
            'oras://community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-stringr:9418abefbed22a4d' :
            'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-stringr:5897cbde71ea29a5'
    }

    input:
    path(seqkit_stats_report)

    output:
    path("processed_seqkit_stats.csv")
    path("filtered_seqkit_stats.csv")
    path("bacass_prep_sample_sheet.csv")
    
    script:
    """
    cat > run_stats.R <<'RS'
    #!/usr/bin/env Rscript
    
    library(dplyr)
    library(data.table)
    library(stringr)
    
    # Read data
    stats <- data.table::fread("${seqkit_stats_report}")
    combine_ont <- as.logical("${params.combine_ont}")
    read_type <- as.character("${params.read_type}")
    fastq_path <- "${params.outdir}/combined_fastqs/"

    # Process - using str_match with capture group
    
    if (combine_ont) {
        cat("Calculating combined coverage by combining runs for isolates...\\n")
        stats <- stats |>
        dplyr::mutate(
                genome_size   = ${params.genome_size},
                predicted_cov = sum_len / genome_size,
                # Extract isolate name from filename using capture group
                isolate = stringr::str_match(file, "_barcode[0-9]+_([^_]+)_combined")[,2]
            ) |>
            dplyr::group_by(isolate) |>
            dplyr::mutate(combined_cov = sum(predicted_cov)) |>
            dplyr::ungroup()
        
        filtered_stats <- stats |>
            dplyr::filter(combined_cov >= ${params.cov_threshold})
    } else {
        cat("Calculating coverage for individual runs without combining...\\n")
        stats <- stats |>
            dplyr::mutate(
                genome_size   = ${params.genome_size},
                predicted_cov = sum_len / genome_size,
                # Extract isolate name from filename using capture group
                isolate = stringr::str_match(file, "_barcode[0-9]+_([^_]+)_combined")[,2]
            )
        filtered_stats <- stats |>
            dplyr::filter(predicted_cov >= ${params.cov_threshold})
    }
    write.csv(stats, "processed_seqkit_stats.csv", row.names = FALSE)
    write.csv(filtered_stats, "filtered_seqkit_stats.csv", row.names = FALSE)

    #Now, make versions of the bacass sample sheet depending on read_type parameter
    if (read_type == "long") {
        sample_sheet <- filtered_stats |>
            dplyr::mutate(
                R1 = NA, 
                R2 = NA, 
                LongFastQ = paste0(fastq_path, "/", file),
                GenomeSize = ${params.genome_size}
            ) |>
            dplyr::rename(ID = isolate) |>
            dplyr::select(ID, R1, R2, LongFastQ, GenomeSize)
    }

    write.csv(sample_sheet, "bacass_prep_sample_sheet.csv", row.names = FALSE)

    RS
    
    chmod +x run_stats.R
    ./run_stats.R
    """
}


workflow {
    main:
    // Validate read type parameter selction.
    if (!(params.read_type in ['long', 'short', 'hybrid'])) {
        error "Invalid read_type: ${params.read_type}. Must be one of: long, short, hybrid"
    }

    //create a channel for inputs from the CSV log file
    ont_ch = params.read_type in ["long", "hybrid"] 
        ? channel.fromPath(params.ont_log, checkIfExists: true)
                          .splitCsv(sep: ",", header: true)
                          .map { row -> 
                              // Extract the relevant fields from the CSV row
                              def isolate = row['Isolate']
                              def barcode = row['Barcode']
                              def parent_dir = row['InputDir']
                              def run_id = row['DataFolder']
                              def fastq_files = file("${parent_dir}/${barcode}/*fastq.gz")
                              def meta = [id: isolate, barcode: barcode, run_id: run_id]
                              // Fixed: tuple order matches process input signature
                              tuple(meta, isolate, barcode, parent_dir, run_id, fastq_files)
                          }
                          .view { meta, isolate, barcode, parent_dir, run_id, fastq_files -> 
                              "Run ID: $run_id, Isolate: $isolate, Barcode: $barcode, Number of FASTQ Files: ${fastq_files.size()}" 
                          }
        : channel.empty() // If read_type is not long or hybrid, create an empty channel

    illumina_ch = params.read_type in ["short", "hybrid"]
        ? channel.fromPath(params.illumina_log, checkIfExists: true)
                          .splitCsv(sep: ",", header: true)
                          .map { row -> 
                              // Extract the relevant fields from the CSV row
                              def isolate = row['Isolate']
                              // Fixed: get paths from CSV columns
                              def read_1 = file(row['R1'])
                              def read_2 = file(row['R2'])
                              def meta = [id: isolate]
                              tuple(meta, isolate, read_1, read_2)
                          }               
        : channel.empty() // If read_type is not short or hybrid, create an empty channel   

    if (params.read_type in ["long", "hybrid"]) {
        //feed the channel into the concatenate process
        CONCATENATE_FASTQ(ont_ch)
        ont_combined_fastqs = CONCATENATE_FASTQ.out.combined_fastq.collect()

        SEQKIT_STATS(ont_combined_fastqs)

        // insert process here with R script to calculate predicted coverage, optionally combine isolate
        PROCESS_STATS_WITH_R(SEQKIT_STATS.out.seqkit_stats_report)
    }
    


}