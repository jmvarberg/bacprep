nextflow.enable.dsl=2
// Resolve params.outdir relative to the launch directory, normalize, and make absolute
def abs_outdir = file(params.outdir).toAbsolutePath().normalize().toString()


//Process to take in ONT sequencing log, find all individual barcode files, concatenate and rename with isolate name at desired location

process STAGE_FASTQS {
    tag "$meta.id"
    publishDir "${params.outdir}/staged_fastqs", mode: 'copy'

    input: 
    tuple val(meta), val(run_id), val(isolate), val(inputFile), val(type), path(fastq_files) 

    output:
    //tuple val(meta), val(isolate), val(barcode), val(run_id),
    path "RunID_${run_id}_${isolate}_${type}_*.fastq.gz", emit: processed_fastq

    script:
    // Normalize to a List<Path> and sort deterministically by name
    def files = (fastq_files instanceof List) ? fastq_files : [ fastq_files ]
    def sorted = files.sort { it.name }

    // Safely quoted names for the shell (single-quoted, with escapes), space-separated
    def quotedList = sorted.collect { it.name.inspect() }.join(' ')

    // Basename without .fastq or .fastq.gz for the single-file case
    def basename = { String n -> n.replaceFirst(/\.fastq(\.gz)?$/, '') }
    def orig = (sorted.size() == 1) ? basename(sorted[0].name) : 'combined'

    """
    set -euo pipefail

    if [ ${sorted.size()} -eq 1 ]; then
        # Single file: just link to a consistent output name for downstream steps
        ln -s ${quotedList} "RunID_${run_id}_${isolate}_${type}_${orig}.fastq.gz"
    else
        # Multiple files: gzip members can be concatenated byte-wise safely
        cat ${quotedList} > "RunID_${run_id}_${isolate}_${type}_${orig}.fastq.gz"
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
        path("seqkit_stats_report.tsv"), emit: combined_seqkit_stats

    """
    seqkit stats -T ${reads} > seqkit_stats_report.tsv
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
    sr_log <- "${params.sr_log}"
    print(paste0("sr_log: ", sr_log))

    #Check if the short read log (sr_log) param was given for short and hybrid sample sheets. If so, validate path and read in. 
    if(!is.null(sr_log) || sr_log != "null" || sr_log == "") {
        if(!file.exists(sr_log)) {
            stop("Path to sr_log file doesn't exist. Check and re-submit.")
        } else {
            sr_log <- data.table::fread(sr_log)
            #validate that columns needed exist
            req_sr_cols = c("Isolate", "R1", "R2")
            missing = setdiff(req_sr_cols, colnames(sr_log))
            if (length(missing) > 0) {
                stop(paste0("Short Read log missing required columns: ", missing))
            }
        }
    }

    # Process - using str_match with capture group
    
    if (combine_ont) {
        cat("Calculating combined coverage by combining runs for isolates...\\n")
        stats <- stats |>
        dplyr::mutate(
                genome_size   = ${params.genome_size},
                predicted_cov = sum_len / genome_size,
                # Extract isolate name from filename using capture group
                Isolate = stringr::str_match(file, "_barcode[0-9]+_([^_]+)_combined")[,2]
            ) |>
            dplyr::group_by(Isolate) |>
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
                Isolate = stringr::str_match(file, "_barcode[0-9]+_([^_]+)_combined")[,2]
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
                LongFastQ = paste0(fastq_path, file),
                Fast5 = NA,
                GenomeSize = ${params.genome_size}
            ) |>
            dplyr::rename(ID = Isolate) |>
            dplyr::select(ID, R1, R2, LongFastQ, Fast5, GenomeSize)
    } else if (read_type == "hybrid") {
        sample_sheet = filtered_stats |>
            dplyr::left_join(sr_log, by = "Isolate") |>
            dplyr::mutate(
                LongFastQ = paste0(fastq_path, file),
                Fast5 = NA,
                GenomeSize = ${params.genome_size}) |>
            dplyr::rename(ID = Isolate) |>
            dplyr::select(ID, R1, R2, LongFastQ, Fast5, GenomeSize)
    }

    write.csv(sample_sheet, "bacass_prep_sample_sheet.csv", row.names = FALSE)

    RS
    
    chmod +x run_stats.R
    ./run_stats.R
    """
}

process CREATE_SHINY_INPUT_R {

    label "process_single"
    publishDir "${params.outdir}/shiny", mode: 'copy'

    conda "${moduleDir}/envs/r_processing.yml"
    container { 
        workflow.containerEngine == 'singularity' && !(task.ext.singularity_pull_docker_container ?: false) ?
            'oras://community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-stringr:9418abefbed22a4d' :
            'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-stringr:5897cbde71ea29a5'
    }

    input:
    path(seqkit_stats_report)
    path(input_log)

    output:
    path("modified_seqkit_stats.csv")
    path("bacprep_log_shiny_input.csv")
    
    script:
    """
    cat > run_stats.R <<'RS'
    #!/usr/bin/env Rscript
    
    
    suppressPackageStartupMessages({
        library(dplyr)
        library(data.table)
        library(stringr)
    })

    
    # Read data
    stats <- data.table::fread("${seqkit_stats_report}")
    input_log <- data.table::fread("${input_log}")
    fastq_path <- "${abs_outdir}/staged_fastqs"

    #Add seqkit stats output to the input log file

    #Step 1: prepare seqkit stats for merging
    stats_mod <- stats |>
        dplyr::mutate(RunID = as.integer(stringr::str_match(file, "RunID_([0-9]+)")[,2]),
                      path = file.path(fastq_path, file))

    write.csv(stats_mod, "modified_seqkit_stats.csv", row.names = F)
    
    #merge seqkit stats with input log using RunID
    input_with_stats <- input_log |>
        dplyr::left_join(stats_mod, by = "RunID") |>
        dplyr::mutate(
            genome_size   = ${params.genome_size},
            predicted_cov = sum_len / genome_size)
    
    write.csv(input_with_stats, "bacprep_log_shiny_input.csv", row.names = FALSE)

    RS
    
    chmod +x run_stats.R
    ./run_stats.R
    """
}

workflow {
    main:

    //create a channel for inputs from the input CSV log file
    input_log_ch = Channel.fromPath(params.input, checkIfExists: true)

    input_ch = Channel.fromPath(params.input, checkIfExists: true)
                .splitCsv(sep: ",", header: true)
                .map { row ->
                    def run_id     = row['RunID']
                    def isolate    = row['Isolate']
                    def data       = row['Data']  // can be a dir OR a file
                    def type       = row['Type']       
                    def inputFile = file(data)
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
                        run_id: run_id,
                        type: type
                    ]

                    tuple(meta, run_id, isolate, inputFile, type, fastq_files)
                }
                .view { meta, run_id, isolate, inputFile, type, fastq_files ->
    "Run ID: ${run_id}, Isolate: ${isolate}, Type: ${type}, FASTQ Count: ${fastq_files.size()}"
                }

    //Pipeline logic 

        def ch_fastqs //make unified channel to handle with or without staging approaches

        if(!params.skip_staging) { //if staging not skipped, then stage all fastq files in the input log, collect into ch_fastqs
        //feed the channel into the concatenate process
        STAGE_FASTQS(input_ch)
        ch_fastqs = STAGE_FASTQS.out.processed_fastq.collect()
        } else {
            //--------------------------------------------
            // Skip staging; user provides a directory with ready FASTQs
            //--------------------------------------------
            if (!params.staged_dir) //check/require path to dir with fastqs if not staged
                error "skip_staging=true but --staged_dir was not provided. Resubmit with a --staged_dir as a valid path to directory with fastq files."

            def staged = file(params.staged_dir)
            if (!staged.exists()) //make sure that path to the staged_dir exists
                error "Directory provided as param '--staged_dir' does not exist: ${staged}"

            // Build a channel of paths from pre-staged directory.
            
            // Use a plain interpolated glob string; no resolve()
            def ch_staged_paths = Channel.fromPath(
                "${staged.toString()}/*.fastq.gz",
                checkIfExists: true
            ).ifEmpty { error "No FASTQ(.fastq.gz) files found in staged_dir: ${staged}" }

            ch_fastqs = ch_staged_paths.collect()
        }
            
    SEQKIT_STATS(ch_fastqs)
    CREATE_SHINY_INPUT_R(SEQKIT_STATS.out.combined_seqkit_stats, input_log_ch)
    // insert process here with R script to calculate predicted coverage, optionally combine isolate
    //PROCESS_STATS_WITH_R(SEQKIT_STATS.out.combined_seqkit_stats)
}