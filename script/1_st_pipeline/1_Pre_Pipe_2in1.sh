core=2
mem=40
genome=mm10
Log=/path/to/.log
RNA_R1=/path/to/R1.fq.gz
RNA_R2=/path/to/R2.fq.gz
out_path=/path/to/out_path/

order_path=None
refe_path=core/barcode_total.100.txt
R_path=/path/to/Rscript

# 0
if [ $order_path == "None" ];then
order_path=core/barcode_add_order.n1_n2.v6_big_25to97_9to98.txt
else
mkdir $out_path
$R_path core/0_prepare_combine_barcode.R -1 $order_path -2 $refe_path -l $Log \
-o $out_path/combine_barcode.round2round1_index1_index2.reorder.txt
fi

# 1
core/1_check_linker.sh $core $out_path $Log $ATAC_R1 $ATAC_R2 $RNA_R1 $RNA_R2

# 2
$R_path core/2_check_multiple_barcode2.R -q $core -l $Log -d $out_path -i "1_RNA_2ormore_linker2.read2.txt" -b $refe_path -o $order_path

