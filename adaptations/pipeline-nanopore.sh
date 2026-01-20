#!/bin/bash
###############################################################################
# QIIME 2 pipeline for 23S–5S amplicon data (Nanopore data)
#
# This script:
# 1. Runs QIIME 2 inside a Singularity container
# 2. Imports nanopore files using a manifest
# 3. Performs denoising with DADA2
# 4. Generates summary visualizations
# 5. Assigns taxonomy and produces taxonomic bar plots
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
# Manifest file creation 
#
# NOTE:
# - This uses the *legacy* comma-separated manifest format.
# - Newer QIIME 2 versions recommend tab-separated manifests.
###############################################################################

# Generate forward-read entries
ls *.fastq.gz | awk -v P="$(pwd)" 'BEGIN{OFS="\t"; print "sample-id","absolute-filepath","direction"} {
    split($0, a, "_");
    print a[1], P"/"$0, "forward"
}' > nanopore-manifest.tsv

###############################################################################
# Import Nanopore data files into QIIME 2
###############################################################################

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path nanopore-manifest.tsv \
  --output-path nanopore-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Generate demultiplexing summary (quality plots)
qiime demux summarize \
  --i-data nanopore-demux.qza \
  --o-visualization nanopore-demux.qzv

###############################################################################
# Subsampling nanopore data and filtering according to q-score
###############################################################################

qiime demux subsample-single \
  --i-sequences nanopore-demux.qza \
  --p-fraction 0.1 \
  --o-subsampled-sequences nanopore-demux-subset.qza

qiime quality-filter q-score \
  --i-demux nanopore-demux-subset.qza \
  --o-filtered-sequences filtered.qza \
  --o-filter-stats filter-stats.qza

qiime metadata tabulate \
  --m-input-file filter-stats.qza \
  --o-visualization filter-stats.qzv

###############################################################################
# Taxonomic assignment using BLAST-like alignment (VSEARCH)
#
# VSEARCH consensus taxonomy assignment
# Parameters:
# - Percent identity threshold set to 97%
###############################################################################

qiime vsearch dereplicate-sequences \
  --i-sequences filtered.qza \
  --o-dereplicated-table table.qza \
  --o-dereplicated-sequences rep-seqs.qza 

###############################################################################
# Visualization of VSEARCH-based taxonomic assignment
###############################################################################

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

###############################################################################
# Taxonomic assignment using Naive Bayesien Classifier 
#
# Parameters:
# - Percent identity threshold set to 70% by default
###############################################################################
qiime feature-classifier classify-sklearn \
  --p-n-jobs 48 \
  --i-classifier 23s5sdb_V1.3.0_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy-23s5s.qza

###############################################################################
# Visualization of NB taxonomic assignment
###############################################################################
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy-23s5s.qza \
  --o-visualization rartaxa-bar-plot-23s5s.qzv
