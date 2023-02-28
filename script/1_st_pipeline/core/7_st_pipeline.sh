#!/bin/bash
if [ $# -lt 1 ]; then
    echo "error.. need args"
    exit 1
fi

ppn=$1
mem=$2
name=$3
path=$4
FW=$5
RV=$6
MAP=$7
ANN=$8
CONT=$9
ID_File=${10}
gene_model=${11}
OUTPUT=1_stpipeline
TMP=1_stpipeline/tmp

cd $path
mkdir -p $TMP

st_pipeline_run.py \
  --output-folder $OUTPUT \
  --temp-folder $TMP \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --ids $ID_File \
  --ref-map $MAP \
  --ref-annotation $ANN \
  --expName $name \
  --htseq-no-ambiguous \
  --verbose \
  --threads $ppn \
  --log-file $OUTPUT/${name}_log.txt \
  --star-two-pass-mode \
  --no-clean-up \
  --contaminant-index $CONT \
  --disable-clipping \
  --min-length-qual-trimming 30 \
  --star-sort-mem-limit ${mem}000000000  \
  $FW $RV

convertEnsemblToNames.py \
  --annotation $ANN \
  --output $OUTPUT/${name}_stdata.symbol.tsv \
  $OUTPUT/${name}_stdata.tsv

st_qa-new.py --input-data $OUTPUT/${name}_stdata.tsv

###Track
samtools index $TMP/mapped.bam
bamCoverage -b $TMP/mapped.bam --normalizeUsing RPKM -o $TMP/mapped.fpkm.bw -p ${ppn} --binSize 10

### RSeQC
# http://rseqc.sourceforge.net/#read-distribution-py
read_distribution.py -i $TMP/mapped.bam -r $gene_model

# samtools
echo -e "\n## reads after remove rRNA ##"
samtools view -c $TMP/contaminated_clean.bam
echo -e "## reads after mapping ##"
samtools view -c $TMP/mapped.bam
