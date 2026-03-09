nextflow.enable.dsl=2

//Process to take in ONT sequencing log, find all individual barcode files, concatenate and rename with isolate name at desired location

process STAGE_FASTQS {
    tag "$meta.id"
    publishDir "${params.outdir}/processed_fastqs", mode: 'copy'

    input:
    tuple val(meta), val(isolate), val(barcode), val(parent_dir), val(run_id), path(fastq_files) 

    output:
    //tuple val(meta), val(isolate), val(barcode), val(run_id),
        path("${run_id}_${barcode}_${isolate}.reads.fastq.gz"), emit: processed_fastq

    script:
    // Build a space-separated, deterministically ordered list of filenames
    def list = fastq_files.collect { it.getName() }.join(' ')
    """
    if [ ${fastq_files.size()} -eq 1 ]; then
        # Single file: just link to a consistent output name for downstream steps
        ln -s ${list} ${run_id}_${barcode}_${isolate}.reads.fastq.gz
    else
        # Multiple files: byte-wise concat of gzip members (fast, no recompression)
        cat ${list} > ${run_id}_${barcode}_${isolate}.reads.fastq.gz
    fi
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
        path(reads) //tuple val(meta), val(isolate), val(barcode), val(run_id), 

    output:
        // one report per merged/symlinked FASTQ tuple val(meta), val(isolate), val(barcode), val(run_id),
        path("combined.seqkit.stats.tsv"), emit: combined_seqkit_stats

    """
    seqkit stats -T ${reads} > combined.seqkit.stats.tsv
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
    path(seqkit_report_file)

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
    #get all files in the seqkit stats output results directory
    stats <- data.table::fread("${seqkit_report_file}")
    print(stats)
    print(summary(stats))
    read_type <- as.character("${params.read_type}")
    fastq_path <- "${params.outdir}/combined_fastqs/"

    cat("Calculating combined coverage with and without combining runs for isolates...\\n")
    
    stats <- stats |>
        dplyr::mutate(
                genome_size   = ${params.genome_size},
                predicted_cov = sum_len / genome_size,
                # Extract isolate name from filename using capture group
                isolate = stringr::str_match(file, "_barcode[0-9]+_([^_]+).reads")[,2]
                isolate = stringr::str_remove(isolate, ".")
            ) |>
            dplyr::group_by(isolate) |>
            dplyr::mutate(combined_cov = sum(predicted_cov)) |>
            dplyr::ungroup()
    
    print(stats)
        
    comb_filtered_stats <- stats |>
        dplyr::filter(combined_cov >= ${params.cov_threshold}) |>
        dplyr::arrange(isolate, desc(combined_cov))
    
    ind_filtered_stats <- stats |>
        dplyr::filter(predicted_cov >= ${params.cov_threshold}) |>
        dplyr::arrange(isolate, desc(predicted_cov))
    }
    write.csv(stats, "processed_seqkit_stats.csv", row.names = FALSE)
    write.csv(comb_filtered_stats, "combined_filtered_seqkit_stats.csv", row.names = FALSE)
    write.csv(ind_filtered_stats, "individual_filtered_seqkit_stats.csv", row.names = FALSE)

    #Now, make versions of the bacass sample sheet depending on read_type parameter
    if (read_type == "long") {
        ind_sample_sheet <- ind_filtered_stats |>
            dplyr::mutate(
                R1 = NA, 
                R2 = NA, 
                LongFastQ = paste0(fastq_path, file),
                Fast5 = NA,
                GenomeSize = ${params.genome_size}
            ) |>
            dplyr::rename(ID = isolate) |>
            dplyr::select(ID, R1, R2, LongFastQ, GenomeSize)

        comb_sample_sheet <- comb_filtered_stats |>
            dplyr::mutate(
                R1 = NA, 
                R2 = NA, 
                LongFastQ = paste0(fastq_path, file),
                Fast5 = NA,
                GenomeSize = ${params.genome_size}
            ) |>
            dplyr::rename(ID = isolate) |>
            dplyr::select(ID, R1, R2, LongFastQ, GenomeSize)

        write.csv(ind_sample_sheet, "bacass_individual_run_filt_${params.cov_threshold}X_thresh_ont_only.csv", row.names = FALSE)
        write.csv(comb_sample_sheet, "bacass_combined_runs_filt_${params.cov_threshold}X_thresh_ont_only.csv", row.names = FALSE)

    }

    if (read_type == "hybrid") {
        sr_log <- data.table::fread("${params.sr_log}")
        sr_log <- sr_log |>
            dplyr::mutate(ID = Isolate, R1 = R1, R2 = R2) |>
            dplyr::select(ID, R1, R2)

        ind_sample_sheet <- ind_filtered_stats |>
            dplyr::mutate(
                LongFastQ = paste0(fastq_path, file),
                Fast5 = NA
                GenomeSize = ${params.genome_size}
            ) |>
            dplyr::rename(ID = isolate) |>
            dplyr::select(ID, LongFastQ, Fast5, GenomeSize) |>
            dplyr::left_join(sr_log, by = "ID") |>
            dplyr::select(ID, R1, R2, LongFastQ, Fast5, GenomeSize)

        comb_sample_sheet <- comb_filtered_stats |>
            dplyr::mutate(
                LongFastQ = paste0(fastq_path, file),
                Fast5 = NA,
                GenomeSize = ${params.genome_size}
            ) |>
            dplyr::rename(ID = isolate) |>
            dplyr::select(ID, LongFastQ, Fast5, GenomeSize) |>
            dplyr::left_join(sr_log, by = "ID") |>
            dplyr::select(ID, R1, R2, LongFastQ, Fast5, GenomeSize)

        write.csv(ind_sample_sheet, "bacass_individual_run_filt_${params.cov_threshold}X_thresh_hybrid.csv", row.names = FALSE)
        write.csv(comb_sample_sheet, "bacass_combined_runs_filt_${params.cov_threshold}X_thresh_hybrid.csv", row.names = FALSE)

    }
    
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

    //create a channel for inputs from the ONT CSV log file
    ont_ch = params.read_type in ["long", "hybrid"] 
        ? channel.fromPath(params.ont_log, checkIfExists: true)
                .splitCsv(sep: ",", header: true)
                .map { row ->
                    def isolate    = row['Isolate']
                    def barcode    = row['Barcode']
                    def inputPath  = row['Input']        // can be a dir OR a file
                    def run_id     = row['DataFolder']
                    def inputFile = file(inputPath)
                    // Detect file vs directory
                    def fastq_files
                    if (inputFile.isDirectory()) {
                        fastq_files = inputFile.listFiles().findAll { it.name.endsWith(".fastq.gz") }
                    }
                    else if (inputFile.isFile() && inputFile.name.endsWith(".fastq.gz")) {
                        fastq_files = [ inputFile ] as List
                    }
                    else {
                        throw new IllegalArgumentException(
                            "InputDir must be a fastq.gz file or directory containing fastq.gz files: $inputPath"
                        )
                    }

                    def meta = [
                        id: isolate,
                        barcode: barcode,
                        run_id: run_id
                    ]

                    tuple(meta, isolate, barcode, inputPath, run_id, fastq_files)
                }
                .view { meta, isolate, barcode, parent_path, run_id, fastq_files ->
    "Run ID: ${run_id}, Isolate: ${isolate}, Barcode: ${barcode}, FASTQ Count: ${fastq_files.size()}"
                }
        : channel.empty() // If read_type is not long or hybrid, create an empty channel

    //create a short read channel. If type is short or hybrid, will pull from the input 'sr_log', otherwise creates empty channel.
    sr_ch = params.read_type in ["short", "hybrid"]
        ? channel.fromPath(params.sr_log, checkIfExists: true)
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

    //Pipeline logic for long/hybrid. 
    if (params.read_type in ["long", "hybrid"]) {

        //feed the channel into the concatenate process
        STAGE_FASTQS(ont_ch)

        processed_fastqs = STAGE_FASTQS.out.processed_fastq.collect()

        SEQKIT_STATS(processed_fastqs)

        // insert process here with R script to calculate predicted coverage, optionally combine isolate
        PROCESS_STATS_WITH_R(SEQKIT_STATS.out.combined_seqkit_stats)
    }
}