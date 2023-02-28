#!/bin/bash
if [ $# -lt 1 ]; then
    echo "error.. need args"
    exit 1
fi

core=$1
out_path=$2
Log=$3

ATAC_R1=$4
ATAC_R2=$5
RNA_R1=$6
RNA_R2=$7

mkdir -p ${out_path}/fastqc

## Skip ATAC ##
if [ $4 != "None" ];then
## 1 ##
echo "Starting to process ATAC data:" >> ${Log}
echo "1.1 Total read pairs:" >> ${Log}
zcat ${ATAC_R2}| awk '{if (NR%4==2){print $0}}' - | wc -l >> ${Log}
echo "1.2 Proper barcode2:" >> ${Log}
zcat ${ATAC_R2}| awk '{if (NR%4==2){print $0}}' - | cut -c 1-15,24-38 | sort - | uniq -c | sed 's/^[ ]*//g' | sort -k1nr | head >> ${Log}
echo "1.3 Proper barcode1:" >> ${Log}
zcat ${ATAC_R2}| awk '{if (NR%4==2){print $0}}' - | cut -c 39-53,62-76 | sort - | uniq -c | sed 's/^[ ]*//g' | sort -k1nr | head >> ${Log}

# 2 or more linker 2
echo "1.4 Multiple barcode1:" >> ${Log}
zcat ${ATAC_R2}| grep "CCCATGATCGTCCGA[ATCG]\{22,24\}CCCATGATCGTCCGA" | wc -l >> ${Log}
echo "1.5 Multiple barcode2:" >> ${Log}
zcat ${ATAC_R2}| grep "ATCCACGTGCTTGAG[ATCG]\{22,24\}ATCCACGTGCTTGAG" | wc -l >> ${Log}
zcat ${ATAC_R2}| grep "ATCCACGTGCTTGAG[ATCG]\{22,24\}ATCCACGTGCTTGAG" > ${out_path}/1_ATAC_2ormore_linker2.read2.txt
# see: 2_check_multiple_barcode2.R

## 2 ##
skewer -f sanger -t $core -m tail -l 35 -x CTGTCTCTTATACACATCTCATGGCACGACTGCA -o ${out_path}/ATAC_skewer -z ${ATAC_R1}

## use FastQC to validate proper trimming and check overall sequence data quality ##
fastqc -o ${out_path}/fastqc -t $core ${out_path}/ATAC_skewer-trimmed.fastq.gz # -d ${out_path}/fastqc/tmp
fi


## Skip RNA ##
if [ -z ${RNA_R2} ];then exit; fi

## 1 ##
echo "Starting to process RNA data:" >> ${Log}
echo "1.1 Total read pairs:" >> ${Log}
zcat ${RNA_R2}| awk '{if (NR%4==2){print $0}}' - | wc -l >> ${Log}
echo "1.2 Proper barcode 2:" >> ${Log}
zcat ${RNA_R2}| awk '{if (NR%4==2){print $0}}' - | cut -c 1-15,24-38 | sort - | uniq -c | sed 's/^[ ]*//g' | sort -k1nr | head >> ${Log}
echo "1.3 Proper barcode1:" >> ${Log}
zcat ${RNA_R2}| awk '{if (NR%4==2){print $0}}' - | cut -c 39-53,62-76 | sort - | uniq -c | sed 's/^[ ]*//g' | sort -k1nr | head >> ${Log}

# 2 or more linker 2
echo "1.4 Multiple barcode1:" >> ${Log}
zcat ${RNA_R2}| grep "CCCATGATCGTCCGA[ATCG]\{22,24\}CCCATGATCGTCCGA" | wc -l >> ${Log}
echo "1.5 Multiple barcode2:" >> ${Log}
zcat ${RNA_R2}| grep "ATCCACGTGCTTGAG[ATCG]\{22,24\}ATCCACGTGCTTGAG" | wc -l >> ${Log}
zcat ${RNA_R2}| grep "ATCCACGTGCTTGAG[ATCG]\{22,24\}ATCCACGTGCTTGAG" >  ${out_path}/1_RNA_2ormore_linker2.read2.txt
# see: 1_check_multiple_barcode2.R

## 2 ##
skewer -f sanger -t $core -m tail -l 85 -x CTGTCTCTTATACACATCT -o ${out_path}/RNA_ME -z ${RNA_R1}
skewer -f sanger -t $core -m tail -l 85 -x CTGTCTCTTATACACATCT -o ${out_path}/RNA_ME2 -z ${out_path}/RNA_ME-trimmed.fastq.gz

zcat ${out_path}/RNA_ME2-trimmed.fastq.gz | grep -v -E "CTGTCTCTTATACACATCT|CATGGCACGACTGCA|TCGGACGATCATGGG|CAAGTATGCAGCGCG|CTCAAGCACGTGGAT|@|\+|F|,|:|#" -B 1 -A 2 - | grep -v "--" - |\
gzip -c - > ${out_path}/RNA_ME2-trimmed.filter.fastq.gz
zcat ${out_path}/RNA_ME2-trimmed.filter.fastq.gz | seqkit subseq -r 1:-51 | gzip - > ${out_path}/RNA_ME2-trimmed.filter.cut_3_50.fastq.gz
rm ${out_path}/RNA_ME-trimmed.fastq.gz ${out_path}/RNA_ME2-trimmed.fastq.gz ${out_path}/RNA_ME2-trimmed.filter.fastq.gz

## use FastQC to validate proper trimming and check overall sequence data quality ##
fastqc -o ${out_path}/fastqc -t $core ${out_path}/RNA_ME2-trimmed.filter.cut_3_50.fastq.gz
