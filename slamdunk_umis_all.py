"""===========================
Pipeline slamdunk_all.py
===========================

Overview
========

This pipeline computes conversion Rate from fastq SLAM-seq data

/!\ This piepeine is still prone to modifications

files :file:``pipeline.yml` and :file:`conf.py`.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/slamdunk_all.py config

Input files
-----------

Input reads must be the lone present in the cwd with suffix ".fastq.gz"
Input reads names should have the following structure:
trimmed-cellline-timepointinhours-R#(replicate number)
example: trimmed-cho-3h-R1.fastq.gz

A genome file

A bed file of 3'UTRs

The pipeline configuration file pipeline.yml.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:
    - python (v3.8.12 with pysam v0.17.0 when built)
    - slamdunk (v0.4.3 when built)
    - samtools (v1.12 when built)
    - umi_tools (v1.0.1 when built)
    - meme (v5.3.0 when built)

Required R modules:
   - optparse
   - matrixStats
   - tidyverse
   - stringr

Pipeline output
===============

The pipeline outputs the different files/directories
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
    * halflife directory contains output of Rscript slampy_halflife_bed
        - ConvRate_processed.tsv : Conversion Rates table of each sample and
            replicate, normalized to time-point 0h and background substracted
            when no4su samples are provided
        - halflife_unfiltered.csv : selected halflife between the 3 models
            half-life were not filtered meaning abberant half-life are present
        - model_fitting_summary.txt : summary of half-life below 0h or above
            the max time-point
        - halflife_filtered.tsv : filtered half-life (to be determined)
            half life betweem -1 and 0 are replace by 0
            haf life above 50h are filtered out
        - bootstrap_halflife.tsv : bootstrapped half-lifes (iteration = 100)
        - .......

See each function for an explanation of each job ran
See slamdunk documentation for more details on slamdunk and alleyoop functions:
https://t-neumann.github.io/slamdunk/docs.html


Code
====

"""

import sys
import os
import sqlite3
import csv
import re
import glob

from cgatcore import pipeline as P
import cgat.GTF as GTF
import cgatcore.iotools as IOTools
from ruffus import *

#Load config file options
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])
########################

@transform("*.fastq.gz",
           regex("(.+).fastq.gz"),
           r"\1_processed.fastq")
def umi_extract(infile, outfile):
    '''Moves the UMI from the read to the read name'''
    n_threads = PARAMS["number_of_threads"]
    log_file = P.snip(outfile, ".fastq")+".log"
    job_memory="8G"
    statement = '''
    umi_tools extract --stdin=%(infile)s --bc-pattern=NNNNNN --log=%(log_file)s --stdout %(outfile)s
    '''
    P.run(statement,
          job_memory="8G")

def make_conf(outfile, list_samples, pattern):
    '''Creates the slamdunk sample file.
    The sample file has four columns:
    1. Filename
    2. Sample name. This is currently the filename with the .fastq removed
    3. is "chase"
    4. time point in minutes. This is retreived from the file anme assuming trimmed-[A-Za-z]-(.+)h-R[0-9]_processed.fastq'''

    with open(outfile, 'w', newline='') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        for i in list_samples:
            samples_grouping = re.search(pattern,i)
            if samples_grouping is None:
                tsv_output.writerow((i, i[:-6], "chase", 0))
            else:
                tsv_output.writerow((i, i[:-6], "chase", float(samples_grouping.groups()[0]) * 60))

@originate("sample_description.tsv")
def make_config_file(outfile):
    '''Takes the reads input names and outputs them to a tsv of filenames for slamdunk dunks'''
    list_samples = glob.glob('*_{}'.format("*processed.fastq"))
    pattern = "trimmed-[A-Za-z0-9]+-(.+)h-R[0-9]_processed.fastq"
    make_conf(outfile, list_samples, pattern)


@follows(mkdir("map"))
@split([make_config_file,umi_extract],
       ["map/{}_slamdunk_mapped.bam".format(P.snip(sample, ".fastq"))
       for sample in glob.glob('*processed.fastq')])
def slamdunk_map(infiles, outfiles):
    '''slamdunk map dunk'''
    infiles = infiles[0]
    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
    bed_file = PARAMS["bed_file_dir"]
    n_threads = PARAMS["number_of_threads"]
    five_trimming = PARAMS["5_trimming"]
    multi_ali = PARAMS["max_alignments_per_read"]
    out_dir = os.getcwd()+"/map"
    job_memory = PARAMS["memory"]
    statement = '''
    slamdunk map -r %(genome_file)s
                 -o %(out_dir)s
                 -5 %(five_trimming)s
                 -n %(multi_ali)s
                 -t %(n_threads)s
                 --skip-sam %(infiles)s
    '''
    P.run(statement,
          job_memory= job_memory,
          job_threads=n_threads)

@follows(mkdir("filter"))
@transform(slamdunk_map,
           regex("map/(.+)_slamdunk_mapped.bam"),
           r"filter/\1_slamdunk_mapped_filtered.bam")
def slamdunk_filter(infiles, outfiles):
    '''slamdunk filter dunk'''
    bed_file = PARAMS["bed_file_dir"]
    n_threads = 1
    mismatch = PARAMS["nm_mismatch"]
    out_dir = os.getcwd()+"/filter"
    job_memory="8G"
    statement = '''
    slamdunk filter -b %(bed_file)s
                    -o %(out_dir)s
                    -t %(n_threads)s
                    -nm %(mismatch)s
                    %(infiles)s
    '''
    P.run(statement,
          job_memory="8G",
          job_threads=n_threads)

@transform(slamdunk_filter,
           regex("filter/(.+).bam"),
           r"filter/\1_dedup.bam")
def umi_dedup(infile, outfile):
    '''Deduplicates the filtered BAM files using the UMI and position'''
    log_file = P.snip(outfile, ".bam")
    job_memory="4G"
    statement = '''
    umi_tools dedup -I %(infile)s
                    -S %(outfile)s
                    -L %(log_file)s
    '''
    P.run(statement,
          job_memory="4G")

@follows(mkdir("snp"))
@transform(umi_dedup,
           regex("filter/(.+).bam"),
           r"snp/\1_snp.vcf")
def slamdunk_snp(infiles, outfiles):
    '''slamdunk snp dunk'''
    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
    n_threads = 1
    out_dir = os.getcwd()+"/snp"
    job_memory="4G"
    statement = '''
    slamdunk snp -r %(genome_file)s
                 -o %(out_dir)s
                 -t %(n_threads)s
                 %(infiles)s
    '''
    P.run(statement,
          job_memory="4G",
          job_threads=n_threads)

@transform(umi_dedup,
           regex("filter/(.+)_dedup.bam"),
           r"filter/\1_dedup.bam.bai")
def dedup_indexing(infiles, outfiles):
    '''Index dedup.bam files'''
    statement = '''
    samtools index %(infiles)s
    '''
    P.run(statement)

@follows(slamdunk_snp,dedup_indexing, mkdir("count"))
@transform(umi_dedup,
           regex("filter/(.+)_dedup.bam"),
           r"count/\1_dedup_tcount.tsv")
def slamdunk_count(infiles, outfiles):
    '''slamdunk count dunk'''
    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
    n_threads = 1
    bed_file = PARAMS["bed_file_dir"]
    rl = PARAMS["max_read_length"]
    out_dir = os.getcwd()+"/count"
    snp_dir = str.join(os.getcwd(), "/snp")
    job_memory="4G"
    statement = '''
    slamdunk count -r %(genome_file)s
                   -b %(bed_file)s
                   -o %(out_dir)s
                   -s %(snp_dir)s
                   -l %(rl)s
                   -t %(n_threads)s
                   %(infiles)s
    '''
    P.run(statement,
          job_memory="4G",
          job_threads=n_threads)

@collate(slamdunk_count,
         regex("count/(.+)-(R.).+tsv"),
         r"count/\1-aggConvRate.tsv")
def merge_ConversionRate(infiles, outfile):
    '''Aggregate Conversion Rate of replicates into a single table'''
    infiles = str.join(" ", infiles)
    statement='''cgat combine_tables %(infiles)s
                    -S %(outfile)s -c 1,2,3,4,5,6
                    -k 7
                    --regex-filename='.+-(R[0-9]).+tsv'
                    --use-file-prefix'''
    P.run(statement)


@collate(slamdunk_count,
         regex("count/(.+)-(R.).+tsv"),
         r"count/\1-aggReadsCPM.tsv")
def merge_CPM(infiles, outfile):
    '''Aggregate ReadsCPM of replicates'''
    infiles = str.join(" ", infiles)
    statement='''cgat combine_tables
                     %(infiles)s
                     -S %(outfile)s
                     -c 1,2,3,4,5,6
                     -k 8
                     --regex-filename='.+-(R[0-9]).+tsv'
                     --use-file-prefix'''
    P.run(statement)

@follows(dedup_indexing,mkdir("rates"))
@transform(umi_dedup,
           regex("filter/(.+).bam"),
           r"rates/\1_overallrates.csv")
def slamdunk_qc_rates(infiles,outfiles):
    '''Computes the overall conversion rates in your reads and plots them
    as a barplot'''
    genome_file = os.path.abspath(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
    bed_file = PARAMS["bed_file_dir"]
    outfiles = os.path.dirname(os.path.abspath(outfiles))
    n_threads = 1
    job_memory="4G"
    statement = '''
    alleyoop rates -o %(outfiles)s
                   -r %(genome_file)s
                   -t %(n_threads)s
                   %(infiles)s
    '''
    P.run(statement,
          job_memory="2G",
          job_threads=n_threads)


@follows(dedup_indexing, mkdir("utrrates"))
@transform(umi_dedup,
           regex("filter/(.+).bam"),
           r"utrrates/\1_mutationrates_utr.csv")
def slamdunk_qc_utrrates(infiles,outfiles):
    '''Checks the individual conversion rates per 3’UTR and plots them as
    boxplots over the entire realm of 3’UTRs'''
    genome_file = os.path.abspath(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
    bed_file = PARAMS["bed_file_dir"]
    outfiles = os.path.dirname(os.path.abspath(outfiles))
    n_threads = 1
    job_memory="2G"
    statement = '''
    alleyoop utrrates -o %(outfiles)s
                      -r %(genome_file)s
                      -t %(n_threads)s
                      -b %(bed_file)s
                      -l 50
                      %(infiles)s
    '''
    P.run(statement,
          job_memory="2G",
          job_threads=n_threads)


@transform(umi_dedup,
            suffix(".bam"),
            r"_summary.tsv")
def slamdunk_summary(infile,outfile):
    '''Basic statistics of the mapping process'''
    statement = '''
    alleyoop summary -o %(outfile)s
                     -t /count
                     %(infile)s
    '''
    P.run(statement)

@follows(mkdir("halflife"))
@merge(slamdunk_count,
       r"halflife/halflife_filtered.tsv")
def run_curvefit(infile, outfile):
    '''From slamdunk count files:
    - computes half lifes
    - boostrap half-life
    (more coming)'''
    cdir = os.getcwd()+"/count"
    threshold = PARAMS["cpm_threshold"]
    bg_sub = PARAMS["background"]
    out_dir = os.getcwd()+"/halflife"
    script_path = os.path.join((os.path.dirname(__file__)),
                               "Rscripts",
                               "slampy_halflife_bed.R")
    statement = '''Rscript %(script_path)s
                         --count-directory %(cdir)s
                         --output-directory %(out_dir)s
                         --cpm--threshold %(threshold)s
                         --background %(bg_sub)s'''
    P.run(statement,
          job_memory="8G")


# @follows(mkdir("meme.dir"), run_curvefit)
# @transform("count/*.bed",
#            regex("count/(.+).bed"),
#            r"meme.dir/\1.fasta")
# def getfasta(infile, outfile):
#     '''From most and less stable bed file, get sequences and generate fasta for streme'''
#     genome_file = os.path.abspath(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
#     statement = '''
#     cat %(infile)s | cgat bed2fasta -L bed2fasta.log --genome-file %(genome_file)s > %(outfile)s
#     '''
#     P.run(statement)

# @follows(getfasta)
# @transform("meme.dir/*.fasta",
#            regex("meme.dir/(.+).fasta"),
#            r"meme.dir/\1.dir")
# def streme(infile, outfile):
#     '''From most and less stable bed file, get sequences and generate fasta for streme'''
#     genome_file = os.path.abspath(os.path.join(PARAMS["genome_dir"], PARAMS["genome"] + ".fa"))
#     statement = '''
#     streme --p %(infile)s --minw 8 --maxw 15 --oc %(outfile)s
#     '''
#     P.run(statement,
#     job_memory="4G",
#     job_threads=2)

@follows(merge_ConversionRate,merge_CPM,slamdunk_qc_rates,slamdunk_qc_utrrates,slamdunk_summary,run_curvefit)
def full():
    '''Later alligator'''
    pass

P.main()
