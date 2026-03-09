# README for bacass_prep workflow

Required Inputs

ONT Log: This sample sheet captures information about runs for ONT long read samples. Required columns include:

* Isolate - character: sample/isolate name ('EC001')
* Barcode - character: barcode value for the isolate from the run ('barcode01')
* DataFolder - character: name of parent folder for the run - assumes that all outputs are in subdirectory 'fastq_pass' in this folder. This is only used for naming/tracking isolate-barcode-run results.
* InputDir - path: full filepath to the 'fastq_pass' directory that contains sequencing results. Assumption is that output files are organized in subdirectories by barcode, with multiple .fastq.gz files per barcode. 


Required Parameters:

## Input Parameters
* ont_log = log file containing information on ONT sequencing runs (default = ''). Required for read_type = 'long' or 'hybrid'
* sr_log = log file containing information on short-read Illumina sequencing runs (default = ''). Required for read_type 'short' or 'hybrid'

## Output Parameters
* outdir = name of output directory to publish results to (default = 'results')

## Help
* help = logical, show help or not (default = false)

## Analysis Parameters
* genome_size = 5100000 //genome size in bp, used for coverage calculation
* cov_threshold = 50 //minimum predicted coverage to include in output sample sheet for bacass
* read_type = type of reads being processed to generate bacass samples sheet. Accepted values: 'long', 'short', 'hybrid'). If 'hybrid', must provide both ONT and Illumina log files. Only isolates with both ONT and Illumina reads will be included in output sample sheet for bacass.
* combine_ont = logical, whether to combine separate ONT runs of same isolate into one combined fastq file for bacass input (default = false)


Development To Do:

* Don't want to have to make multiple copies of the fastq files for each type of analysis (combined vs. not, for example). Should be an option to point to an existsing location if fastqs are already processed OR if they don't need to be concatenated. 

Or, if we just change the 'combine_ont' flag, then maybe it just re-runs the R process? 
Yes, this is correct. Once you've concatenated them, changing the other parameters for the R processing doesn't require re-running the cat or seqkit stats processes.

* Maybe just get rid of the option for combining? It doesn't cost anything here except for the R script. Can just have two versions of the bacass sample sheet out, one filtered on per sample, one filtered on per isolate (combined) coverage values.

* Need to build in functionality for creating bacass sample sheets for short and hybrid options.
* Will need to have R function that reads in the short read log, left joins with the ONT log after filtering, and outputs correct R1/R2 paths.
* Need to think about where the fastq files should be stored? What would make sense with respect to where bacass will be run? Think within the context of bacass prep > bacass > annotyper workflow.


DEV Branch Refactor:

The reason we need this workflow is to accomplish the following:

1) Pull ONT sequencing results from many runs/locations into one location and combine fastqs into one file per sample/run.
2) Run seqkit stats on each isolates fastq file to calculate predicted coverage. We want to pre-filter to only keep strains with sufficient predicted coverage based on a user-defined threshold.
3) Allow for samples to be retained if combining ONT runs from multiple sequencing runs allows it to pass the threshold.
4) Prepare a sample sheet for bacass that only includes ONT samples that should pass the minimum sequence coverage threshold and optionally combine the ONT samples with Illumina short read samples to use for polishing of assemblies with Polypolish.

To accomplish this, the user must provide:

1) ONT sample log with run, barcode, and isolate information and a path to the location where the fastq files are located.
2) Desired minimum predicted coverage value (defaults to 50x)
3) Logical input on whether samples should be allowed to be combined to reach the threshold?
4) Logical input on whether there are short read samples to be added to the bacass sample sheet output?
5) If Yes to short read, a CSV file input with run, barcode and isolate information. Isolate ID must match exactly the isolate information in the ONT long read table.

To Do/Changes:

1) Changed ONT channel specification to handle both paths to directories with multiple fastq files to be concatenated, as well as to individual FASTQ files.
