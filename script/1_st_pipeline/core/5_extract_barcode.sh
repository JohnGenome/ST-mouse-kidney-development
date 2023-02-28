#!/bin/bash
if [ $# -lt 1 ]; then
    echo "error.. need args"
    exit 1
fi

name=$1
RNA_R1=$2
RNA_R2=$3
out_path=$4

mkdir -p ${out_path}
## 1 ##
umi_tools extract --extract-method=regex \
--bc-pattern2="^(?P<umi_1>.)(?P<discard_1>ATCGGCGTACGACT){s<=1}.{8}(?P<discard_2>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=1}.{8}(?P<discard_3>CCCATGATCGTCCGATGCAGTCGTGCCATGAGATGTGTATAAGAGACAG){s<=2,i<=1,d<=1}.{10}(?P<discard_4>.*)$" \
-I ${RNA_R1} -S ${out_path}/${name}_1.extract.fq.gz --read2-in=${RNA_R2} --read2-out=${out_path}/${name}_2.extract.fq.gz \
-L ${out_path}/extract.log

