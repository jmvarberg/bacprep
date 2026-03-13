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
}