library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

mKid.E15_4 = mKid.merge[,mKid.merge$orig.ident=="E15_4"]

## 0 Variable gene
mKid.E15_4 = SCTransform(mKid.E15_4,assay = "RNA",variable.features.n = 10000)
var_gene.ls = VariableFeatures(mKid.E15_4[["SCT"]])

## 1 Expressed gene
corrected_umi.ls = apply(mKid.E15_4@assays[["SCT"]]@counts, 1, sum)
pdf("Plot_E15_4/1_filter_expressed_genes.pdf")
ggplot(data.frame("corrected_umi"=log2(corrected_umi.ls+1)),aes(x=corrected_umi)) + 
  geom_density() + geom_vline(xintercept=6,colour="blue",linetype="longdash") + theme_bw()
dev.off()
table(corrected_umi.ls > 64)
# FALSE  TRUE 
# 8289 10956
express_gene.ls = names(corrected_umi.ls)[corrected_umi.ls > 64]
table(var_gene.ls %in% express_gene.ls)
# FALSE  TRUE 
# 3774  6226

## 2 Expressed LR pairs
lr_pair.df = read.table("../core/mouse_ligand_receptors.txt",header = F,stringsAsFactors = F) # 1067
lr_pair.df = lr_pair.df[lr_pair.df$V1 %in% express_gene.ls & lr_pair.df$V2 %in% express_gene.ls,] # 268
lr_pair.df = lr_pair.df[lr_pair.df$V1 %in% var_gene.ls | lr_pair.df$V2 %in% var_gene.ls,] # 248
row.names(lr_pair.df) = NULL

## 3 Export
meta.df = mKid.E15_4@meta.data[,c("Ann_v5","coord","nA","nB")]
correct_count.m = as.matrix(mKid.E15_4@assays[["SCT"]]@counts[sort(unique(c(lr_pair.df$V1,lr_pair.df$V2))),])
base::save(list = c("correct_count.m","lr_pair.df","meta.df"),file = "Plot_E15_4/2_CCC_prepare.RData")





