#!/bin/sh
#fqs="laca/workflow/resources/data"
#database="test/Database"
#wd="test"
fqs=$1
database=$2
wd=$3

laca init -f $fqs -d $database -w $wd --kmerbin --fqs-min 0 \
--cluster NGSpeciesID --NGSpeciesID2 --cluster isONcorCon --cluster kmerCon \
--cluster clustCon --cluster umiCon

# JSON object sintax
# https://stackoverflow.com/questions/66461259/snakemake-specifying-configuration-parameter-from-nested-config-yaml-file-on-the
laca run all -w $wd \
--config kraken2="'{prebuilt: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20220607.tar.gz}'" \
--config classifier="kraken2"
