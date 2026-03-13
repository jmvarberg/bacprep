# README for bacass_prep workflow

# Required Input

__Sequencing Log:__ A CSV file sample sheet capturing information for all sequencing runs associated with the project. Required columns include:

* __RunID__ - unique numeric identifier
* __Isolate__ - character: sample/isolate name ('EC001')
* __Data__ - path to fastq.gz/fq.gz files for sequencing data for this sample. Can either be path to a file, or to a directory containing multiple files. If points to a directory, then will combine all files in the directory into a single file for downstream processing.
* __Type__ - ONT (long read) or SR (illumina/short read).

Additional metadata columns can be included in the log for record keeping purposes but are not necessary for the pipeline to run.


# Required Parameters:

## Input Parameters
* input = path to CSV input log file (--intput)
* skip_staging = (true/false), logical to skip staging/copying the fastq files into the local directory. Use if you already have all of the fastq files prepared for analysis and don't require any merging/concatenating.
* staged_dir = path to directory containing the fastq files that are ready for processing. Required if skip_staging = true.

## Output Parameters
* outdir = name of output directory to publish results to (default = 'results')

## Analysis Parameters
* genome_size = 5100000 //genome size in bp, used for coverage calculation
