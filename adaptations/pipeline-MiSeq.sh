#!/bin/bash
###############################################################################
# QIIME 2 pipeline for 23S–5S amplicon data (paired-end MiSeq reads)
#
# This script:
# 1. Runs QIIME 2 inside a Singularity container
# 2. Imports paired-end FASTQ files using a manifest
# 3. Performs denoising with DADA2
# 4. Generates summary visualizations
# 5. Trains a custom Naive Bayes classifier on a homemade 23S–5S reference
# 6. Assigns taxonomy and produces taxonomic bar plots
#
# Author: Camille JACQUELINE
# Environment: Singularity + QIIME 2
###############################################################################

###############################################################################
# Environment setup
###############################################################################

# Move to directory containing the Singularity image
cd /srv/scratch/singularity/

# Start an interactive QIIME 2 shell
singularity shell qiime2.sif

# Define a writable temporary directory for QIIME 2
export TMPDIR=/yourdirectory/

###############################################################################
# Project directory
###############################################################################

# Set working directory containing FASTQ files
cd /yourdirectory/

###############################################################################
# Manifest file creation (paired-end, Phred33)
#
# NOTE:
# - This uses the *legacy* comma-separated manifest format.
# - Newer QIIME 2 versions recommend tab-separated manifests.
# - Sample IDs are extracted from filenames using "_" as delimiter.
###############################################################################

# Generate forward-read entries
ls *_R1.fastq.gz | \
awk -v P="$(pwd)" -F'_' 'BEGIN{OFS=","}{print $1, P"/"$0,"forward"}' \
> forward_manifest_temp

# Generate reverse-read entries
ls *_R2.fastq.gz | \
awk -v P="$(pwd)" -F'_' 'BEGIN{OFS=","}{print $1, P"/"$0,"reverse"}' \
> reverse_manifest_temp

# Combine header and both read directions into final manifest
cat <(echo sample-id,absolute-filepath,direction) \
    forward_manifest_temp \
    reverse_manifest_temp \
> paired-33-manifest

###############################################################################
# Import paired-end FASTQ files into QIIME 2
###############################################################################

qiime tools import \
  --input-path paired-33-manifest \
  --output-path paired-end-demux.qza \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33

# Generate demultiplexing summary (quality plots)
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv

###############################################################################
# Primer trimming
###############################################################################
# Trim first set of primers 
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences paired-end-demux.qza \
  --p-front-f TNNNNNNNNNNNNNNNNNNC \
  --p-front-r GNNNNNNNNNNNNNNNNNT \
  --p-error-rate 0 \
  --p-cores 48 \
  --o-trimmed-sequences trimmed23s5s.qza \
  --p-discard-untrimmed \
  --verbose

qiime demux summarize \
  --i-data trimmed23s5s.qza \
  --o-visualization trimmed23s5s_summary.qzv

###############################################################################
# DADA2 denoising (23S–5S amplicon)
#
# Parameters:
# - Reads truncated at 148 bp (forward and reverse), length of the illumina reads
# - High number of training reads for error model
# - Multi-threaded execution
###############################################################################

qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed23s5s.qza \
--p-trunc-len-f 280 \
--p-trunc-len-r 280 \
--p-min-fold-parent-over-abundance 2 \
--p-n-reads-learn 500000 \
--p-n-threads 50 \
--o-representative-sequences rep-seqs-dada2-23s5s.qza \
--o-table table-dada2-23s5s.qza \
--o-denoising-stats stats-dada2-23s5s.qza \
--verbose

###############################################################################
# DADA2 statistics and feature table summaries
###############################################################################

qiime metadata tabulate \
  --m-input-file stats-dada2-23s5s.qza \
  --o-visualization stats-dada2-23s5s.qzv

# Rename outputs for consistency
mv rep-seqs-dada2-23s5s.qza rep-seqs-23s5s.qza
mv table-dada2-23s5s.qza table-23s5s.qza

# Summarize feature table
qiime feature-table summarize \
  --i-table table-23s5s.qza \
  --o-visualization table-23s5s.qzv

# Visualize representative sequences
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-23s5s.qza \
  --o-visualization rep-seqs-23s5s.qzv

###Following steps are similar to main pipeline
