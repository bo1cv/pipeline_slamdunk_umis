# pipeline_slamdunk_umis

## The pipeline performs the following:
   * deduplicate reads with umi_tools
   * Performs SLAM-seq analysis with SLAM-DUNK
   Docs: https://t-neumann.github.io/slamdunk/docs.html#docstart
   * Determine mRNA half-life based on Tcounts output from SLAM-DUNK


## Inputs needed      
   1. Inputs Reads
   Input reads must be the lone present in the cwd with suffix ".fastq.gz"
   Input reads names should have the following structure:
   trimmed-cellline-timepointinhours-R#(replicate number)
   example: trimmed-cho-3h-R1.fastq.gz
   2. A genome file
   3. A bed file of 3'UTRs


## The pipeline outputs the different files/directories
   * sample_description.tsv file with sample names and time points asked by
       slamdunk_all
   * xxx_processed.fastq files : ouputs of the umi tools extract
   * map directory : contains outputs from slamdunk map
       "xxx_slamdunk_mapped.bam" and  .log
   * filter directory : contains output from
       - slamdunk filter "xxx_slamdunk_mapped_filtered.bam" with index .bai
       and .log
       - umi_tools dedup "xxx_slamdunk_mapped_filtered_dedup.bam" with
       index .bai
       - slamdunk alleyoop summary "xxx_slamdunk_mapped_filtered_summary.tsv"
       and .log, "xxx_slamdunk_mapped_filtered_summary_PCA.txt"
   * snp directory : contains outputs from slamdunk snp
       "xxx_slamdunk_mapped_filtered_snp.vcf" and .log
   * count directory : contains outputs from
       - slamdunk count "xxx_slamdunk_mapped_filtered_dedup_tcount.tsv"
       and .log with Conversion rates for each transcript (used to
           determine half-life),
           "xxx_slamdunk_mapped_filtered_dedup_tcount_mins.bedgraph" and
           "xxx_slamdunk_mapped_filtered_dedup_tcount_plus.bedgraph"
       - cgat combine_tables  "xxx-aggConvRate.tsv" and "xxx-agg-ReadsCPM.tsv"
   * halflife directories contains output of Rscripts
       - ConvRate_processed.tsv : Conversion Rates table of each sample and
           replicate, normalized to time-point 0h and background substracted
           when no4su samples are provided
       - halflife_unfiltered.csv : selected halflife between the 3 models
           half-life were not filtered meaning abberant half-life are present
           OR here only model 1 hal-lifes
       - model_fitting_summary.txt : summary of half-life below 0h or above
           the max time-point
           OR NOT if only model 1
       - models_halflife_decay_aic.tsv: result of all models if 3 models,
           nothing if only model 1
       - bootstrap_halflife.tsv : bootstrapped half-lifes (iteration = 1000)
       - halflife_percentileCIs.tsv : Pencentile CI calculated on bootstrapped
           half-lifes, wth other stats like mean and median
       - halflife_filtered.tsv : filtered half-lifes for h-l that were out
       of the bs CI, < 0h or > 24h
       - halflife_filtered.log : number of transcripts filtered

## Requirements

On top of the default CGAT setup, the pipeline requires the following
* Software:
    - python (v3.8.12 with pysam v0.17.0 when built)
    - slamdunk (v0.4.3 when built)
    - samtools (v1.12 when built)
    - umi_tools (v1.0.1 when built)
    - meme (v5.3.0 when built)
* R modules:
   - optparse
   - matrixStats
   - tidyverse
   - stringr
   - foreach
   - doParallel
* Python modules
   - pyteiser (pip install) https://github.com/goodarzilab/pyteiser

## Configuration
The pipeline requires a configured :file: `pipeline.yml` file.

Make a directory with your project name.
Configure the pipeline with `python [path_to_repo]/pipeline_slamdunk_umis.py config`.
A pipeline.log and pipeline.yml file(s) will be added to your new directory.
Modify the pipeline.yml according to your project (specify annotation database and directory, database for uploading the outputs; specify options for Salmon quantification).

## Pipeline use
Run the pipeline with `python [path_to_repo]/pipeline_slamdunk_umis.py make full -v5`.

For running the pipeline on a large set of samples, submit the pipeline onto the cluster (sharc), using a submit_pipeline custom script.

## P.S
/!\/!\Only model 1 included in pipeline/!\/!\
The Rscript directory contains 4 scripts to calculate h-l from conversion rates.
all rep aboceCPM: model(s) are given T>C Conversion rates of transcripts
with a CPM above threshold in all replicates of each time point.
2rep aboveCPM: model(s) were given model(s) are given T>C Conversion rates of
transcripts with a CPM above threshold in at least 2 replicates of each time points

3model: Rscript with the 3 models: these scripts are not included in the pipeline,
but I kept the code in case.

model1: Rscript with only model 1.
