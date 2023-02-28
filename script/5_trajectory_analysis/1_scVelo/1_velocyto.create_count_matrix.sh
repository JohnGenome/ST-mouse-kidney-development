core=4
samMem=8G
genoGTF=/path/to/Mus_musculus.GRCm38.102.chr.gtf
rmskGTF=/path/to/mm10_rmsk.gtf
bcFile=/path/to/core/barcode_total_r2r1.10000.txt

bamPath=/path/to/bam/
outPath=/path/to/outPath/
log=/path/to/.log

## [1] Check that BAM file has been sorted (`SO:coordinate`)
samtools view -H ${bamPath}/mapped.bam | grep @HD
# samtools sort ${bamPath}/mapped.bam -o ${bamPath}/mapped.sort.bam

## [2] Add CB/UB tag to BAM file
samtools view -h ${bamPath}/mapped.bam | sed 's/B0:Z/CB:Z/g' | sed 's/B3:Z/UB:Z/g' | samtools view -bS - -o ${bamPath}/mapped.CBflag.bam

## [3] Get cell-barcode-sorted BAM file
samtools sort -t CB -O BAM -o ${bamPath}/cellsorted_mapped.CBflag.bam -m ${samMem} -@ ${core} ${bamPath}/mapped.CBflag.bam

## [4] Run Velocyto
velocyto run -b ${bcFile} -o ${outPath} -m ${rmskGTF} ${bamPath}/mapped.CBflag.bam ${genoGTF} >${log} 2>&1

rm ${bamPath}/mapped.CBflag.bam ${bamPath}/cellsorted_mapped.CBflag.bam
