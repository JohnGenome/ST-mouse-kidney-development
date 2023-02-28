core=2
mem=55
name=RNA_E16_220704_Slice29
RNA_R1=/path/to/R1.fq.gz
RNA_R2=/path/to/R2.fq.gz
out_path=/path/to/out_path/

MAP=/path/to/STAR_index
ANN=/path/to/Mus_musculus.GRCm38.102.chr.cover_intron.gtf
CONT=/path/to/STAR_index
ID_File=/path/to/combine_barcode.round2round1_index1_index2.v3_big.txt
gene_model=/path/to/mm10_RefSeq.fromRSeQC.bed

core/5_extract_barcode.sh $name $RNA_R1 $RNA_R2 $out_path/fastq

R_path=/path/to/Rscript
zcat $out_path/fastq/${name}_2.extract.fq.gz | awk 'NR%4==2{print substr($1,1,16)}' > $out_path/fastq/${name}_2.extract.barcode
$R_path core/9_check_barcode.R -b $ID_File -p $out_path/fastq -n ${name}_2.extract.barcode
rm $out_path/fastq/${name}_2.extract.barcode

core/7_st_pipeline.sh $core $mem $name $out_path \
${out_path}/fastq/${name}_2.extract.fq.gz ${out_path}/fastq/${name}_1.extract.fq.gz $MAP $ANN $CONT $ID_File $gene_model

