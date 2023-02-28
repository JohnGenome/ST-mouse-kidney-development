#!/bin/bash
if [ $# -lt 1 ]; then
    echo "error.. need args"
    exit 1
fi

export PATH=/software/biosoft/software/python/anaconda3-python3-2018/bin:$PATH
export LD_LIBRARY_PATH=/software/biosoft/software/OpenBLAS/bin/lib:$LD_LIBRARY_PATH

ppn=$1
name=$2
Log=$3
out_path=$4
fq_file=$5
index_mm10=$6
gtf_mm10=$7
gene_model=$8

mkdir -p $out_path
cd $out_path
echo -e "\nStarting to analysis RNA data as bulk assay:" >> ${Log}
###Mapping
hisat2 -p $ppn --dta-cufflinks --rna-strandness F -x $index_mm10 -U $fq_file --un-gz Unmaped.fq.gz| \
samtools view -F 4 -q 10 -bS - | \
samtools sort - -o ${name}.sorted.bam # -@ $ppn

###Track
samtools index ${name}.sorted.bam
bamCoverage -b ${name}.sorted.bam --normalizeUsing RPKM -o ${name}.fpkm.bw -p ${ppn} --binSize 10

### cufflinks ####
mkdir 1_cufflinks
cufflinks -p $ppn -G $gtf_mm10 -o 1_cufflinks ${name}.sorted.bam # --library-type F

echo "4.1 Detected genes:" >> ${Log}
awk 'BEGIN{Nh=0;Nm=0;Nl=0}{if(NR>1){if($10>10){Nh=Nh+1}else if($10>1){Nm=Nm+1}else if($10>0.1){Nl=Nl+1}}}END{print "High: "Nh"\nMod: "Nm"\nLow: "Nl}' \
1_cufflinks/genes.fpkm_tracking >> ${Log}

### RSeQC
# http://rseqc.sourceforge.net/#read-distribution-py
echo "4.2 Read distribution:" >> ${Log}
read_distribution.py -i ${name}.sorted.bam -r $gene_model >> ${Log}
date >> ${Log}
rm Unmaped.fq.gz
