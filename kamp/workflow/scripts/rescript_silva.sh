#!/bin/bash

OUT_DIR=$1
JOBS=$2

# Get silva ref
fo_rna=$OUT_DIR/silva-138.1-ssu-nr99-rna-seqs.qza
fo_tax=$OUT_DIR/silva-138.1-ssu-nr99-tax.qza
if [ ! -f $fo_rna ] || [ ! -f $fo_tax ]; then
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --p-include-species-labels \
    --p-no-rank-propagation \
    --o-silva-sequences $fo_rna \
    --o-silva-taxonomy $fo_tax 
fi

fo_dna=$OUT_DIR/silva-138.1-ssu-nr99-seqs.qza
if [ ! -f $fo_dna ]; then
qiime rescript reverse-transcribe \
    --i-rna-sequences $fo_rna \
    --o-dna-sequences $fo_dna
fi
rm -f $fo_rna

# Culling low-quality sequences with cull-seqs
fo_cull=$OUT_DIR/silva-138.1-ssu-nr99-seqs-cleaned.qza
if [ ! -f $fo_cull ]; then
qiime rescript cull-seqs \
    --i-sequences $fo_dna \
    --o-clean-sequences $fo_cull \
    --p-n-jobs $JOBS
fi
rm -f $fo_dna

# Filtering sequences by length and taxonomy
fo_filt=$OUT_DIR/silva-138.1-ssu-nr99-seqs-filt.qza
fo_discard=$OUT_DIR/silva-138.1-ssu-nr99-seqs-discard.qza
if [ ! -f $fo_filt ]; then
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences $fo_cull \
    --i-taxonomy $fo_tax \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs $fo_filt \
    --o-discarded-seqs $fo_discard
fi
rm -f $fo_cull $fo_discard

# Dereplication of sequences and taxonomy
fo_drep=$OUT_DIR/silva_seqs.qza
fo_tax_drep=$OUT_DIR/silva_tax.qza
if [ ! -f $fo_drep ] || [ ! -f $fo_tax_drep ]; then
qiime rescript dereplicate \
    --i-sequences $fo_filt \
    --i-taxa $fo_tax \
    --p-perc-identity 1.0 \
    --p-rank-handles 'silva' \
    --p-mode 'uniq' \
    --o-dereplicated-sequences $fo_drep \
    --o-dereplicated-taxa $fo_tax_drep \
    --p-threads $JOBS
fi
rm -f $fo_filt $fo_tax