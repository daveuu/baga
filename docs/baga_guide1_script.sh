#!/usr/bin/env bash

set -e

# update this to your email for using NCBI Entrez!
EMAIL=my.email.address@me.com

if [ "$EMAIL" = "my.email.address@me.com" ];
then
    echo "Please edit this script and update the EMAIL variable to your email address for use with NCBI Entrez"
    exit 1
fi

# to select a more recent baga version (commit), browse here:
# https://github.com/daveuu/baga/commits/master
# (click the clipboard icon next to the commit you want)

BAGAVERSION=a54b6de7e16ac37dd89f401339b5a3f24d294113

BAGA_CLI=baga/baga_cli.py
NUM_CPUS=7
GENOME=NC_011770.1
GENOME_LEN=6601757
READSNAME=Liverpool
GATK=~/GenomeAnalysisTK.jar
# if your system version is java 1.7, use:
# JRE17=java 
# else download java runtime environment
JRE17=~/Downloads/jre1.7.0_80/bin/java

## some preliminary checks
if [ -f $GATK ];
then
   echo "$GATK exists: using it."
else
   echo "$GATK does not exist: please download at least version 3.3 from https://www.broadinstitute.org/gatk/download/ or see https://baga.readthedocs.org/en/latest/guide1/#get-the-genome-analysis-toolkit-gatk for more information"
   exit 1
fi

if [ -f $JRE17 ];
then
   echo "$JRE17 exists: using it."
else
   echo "$JRE17 does not exist: see https://baga.readthedocs.org/en/latest/guide1/#get-the-genome-analysis-toolkit-gatk for more information"
   exit 1
fi

if [ -f phylogeny_sample_names.txt ];
then
   echo "phylogeny_sample_names.txt exists: using it."
else
   wget https://baga.readthedocs.org/en/latest/phylogeny_sample_names.txt
fi


## Get BAGA
{
git clone https://github.com/daveuu/baga.git
cd baga
git checkout $BAGAVERSION
cd ..
} || {
cd baga
git fetch origin
git reset --hard origin/master
git checkout $BAGAVERSION
cd ..
}


## 1. Collect the data

$BAGA_CLI Dependencies \
--checkgetfor CollectData

$BAGA_CLI CollectData \
--reads_group_name $READSNAME \
--reads_download ERR953517 ERR953518 ERR953481 ERR953487 ERR953492 ERR953511 ERR953519 ERR953520 ERR953521 ERR953522 ERR953523 ERR953524 ERR953525 ERR953526 ERR953527 ERR953528

$BAGA_CLI CollectData \
--genomes $GENOME \
--email $EMAIL


## 2. Prepare the short read data

$BAGA_CLI Dependencies \
--checkgetfor PrepareReads

$BAGA_CLI PrepareReads \
--reads_name $READSNAME \
--subsample_to_cov 80 $GENOME_LEN

$BAGA_CLI PrepareReads \
--reads_name $READSNAME \
--adaptors --trim \
--max_cpus $NUM_CPUS

# optionally remove some intermediate files
#baga/baga_cli.py PrepareReads --reads_name Liverpool --delete_intermediates


## 3. Align short reads

$BAGA_CLI Dependencies \
--checkgetfor AlignReads

$BAGA_CLI AlignReads \
--reads_name $READSNAME \
--genome_name $GENOME \
--align --deduplicate

$BAGA_CLI AlignReads \
--reads_name $READSNAME \
--genome_name $GENOME \
--indelrealign \
--GATK_jar_path $GATK --JRE_1_7_path $JRE17



## 4. BAGA filter: rearrangements

$BAGA_CLI Dependencies \
--checkgetfor Structure

$BAGA_CLI Structure \
--reads_name $READSNAME \
--genome_name $GENOME \
--ratio_threshold 0.4 \
--check

$BAGA_CLI Structure \
--reads_name $READSNAME \
--genome_name $GENOME \
--plot


$BAGA_CLI Structure \
--reads_name $READSNAME \
--genome_name $GENOME \
--include_samples ERR953521 \
--plot_range 4391000 4401500


## 5. BAGA filter: repeats

$BAGA_CLI Dependencies \
--checkgetfor Repeats

$BAGA_CLI Repeats \
--genome_name $GENOME \
--find

$BAGA_CLI Repeats \
--genome_name $GENOME \
--plot



## 6. Call variants

$BAGA_CLI CallVariants \
--reads_name $READSNAME \
--genome_name $GENOME \
--calleach \
--calljoint \
--hardfilter \
--recalibrate \
--max_cpus 7 --max_memory 8 \
--GATK_jar_path $GATK --JRE_1_7_path $JRE17


$BAGA_CLI CallVariants \
--reads_name $READSNAME \
--genome_name $GENOME \
--calleach \
--calljoint \
--hardfilter \
--max_cpus 7 --max_memory 8 \
--GATK_jar_path $GATK --JRE_1_7_path $JRE17



## 7. Apply BAGA filters

$BAGA_CLI Dependencies \
--checkgetfor ComparativeAnalysis

$BAGA_CLI FilterVariants \
--reads_name $READSNAME \
--genome_name $GENOME \
--filters genome_repeats rearrangements


## 8. Multiple sequence alignment for comparative analyses

$BAGA_CLI ComparativeAnalysis \
--reads_name $READSNAME \
--genome_name $GENOME \
--build_MSA \
--include_invariants


## 9. Infer Phylogeny

$BAGA_CLI ComparativeAnalysis \
--infer_phylogeny \
--path_to_MSA $GENOME__$READSNAME_SNPs.phy \
--out_group ${GENOME:0:10}


## 10. Infer genetic imports

$BAGA_CLI ComparativeAnalysis \
--infer_recombination \
--path_to_MSA $GENOME__$READSNAME_SNPs.phy \
--path_to_tree $GENOME__$READSNAME_SNPs_rooted.phy_phyml_tree


## 11. Summarise: plot tree and generate tables

$BAGA_CLI ComparativeAnalysis \
--plot_phylogeny \
--plot_transfers \
--path_to_tree $GENOME__$READSNAME_SNPs_rooted.phy_phyml_tree \
--use_names "phylogeny_sample_names.txt" \
--genome_name $GENOME \
--out_group ${GENOME:0:10}

$BAGA_CLI SummariseVariants \
--genome_name $GENOME \
--reads_name $READSNAME \
--filters GATK genome_repeats rearrangements \
--cumulative

$BAGA_CLI SummariseVariants \
--simple  \
--genome_names $GENOME \
--vcfs_paths variants/$GENOME/$READSNAME__$GENOME_2_samples_hardfiltered_[SI]*s__F_genome_repeats__F_rearrangements.vcf

