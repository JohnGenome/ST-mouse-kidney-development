library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

## 0 Load data ##
count.m = read.table("RNA_Slice29_stdata.symbol.tsv")
count.m = as.matrix(count.m)
tmp = CreateSeuratObject(counts = t(count.m), project = "E16_1", assay = "RNA")
tmp$channel = "10um"
tmp$stage = "E16"
tmp = tmp[,tmp$nCount_RNA > 500]

tmp$coord = paste0("E16_1_",row.names(tmp@meta.data))
tmp$nA = as.integer(str_split_fixed(tmp$coord,"_|x",4)[,3])
tmp$nB = as.integer(str_split_fixed(tmp$coord,"_|x",4)[,4])
E16_1 = tmp

barcode.df = read.table("../core/combine_barcode.round2round1_index1_index2.v6_big_25to97_9to98.txt",stringsAsFactors = F)
row.names(barcode.df) = paste(barcode.df$V2,barcode.df$V3,sep = "x")
E16_1$barcodeBA = barcode.df[row.names(E16_1@meta.data),"V1"]
E16_1$barcodeBA = paste0("mapped_RCI4O:",E16_1$barcodeBA,"x") # change prefix !!!
rm(list = c("count.m","tmp","barcode.df"))

ann.df = read.csv("../74_Combine_E11_E12_E13_E15_E16_E18_Adult/Out/6_metadata.Annv5.csv",row.names = 1,stringsAsFactors = F)
ann.df = ann.df[ann.df$orig.ident == "E16_1",]

## 1 KNN k=8 ##
ct.ls = c(8,7,24,6,25,12,10,2,27,20,28,14,29)
E16.Npn = E16_1[,E16_1$coord %in% ann.df$coord[ann.df$Ann_v5 %in% ct.ls]]
E16.Npn$Ann_v5 = ann.df[E16.Npn$coord,"Ann_v5"]
Idents(E16.Npn) = factor(E16.Npn$Ann_v5,levels = 0:32)

E16.Npn <- SCTransform(E16.Npn, assay = "RNA", verbose = T,variable.features.n = 3000)
E16.Npn <- RunPCA(E16.Npn, assay = "SCT", verbose = FALSE)
E16.Npn <- FindNeighbors(E16.Npn, reduction = "pca", dims = 1:k, k.param = 8)
KNN.m = as.matrix(E16.Npn@graphs$SCT_nn)
row.names(KNN.m) = E16.Npn$barcodeBA
colnames(KNN.m) = E16.Npn$barcodeBA

## 2 Meta cell information ##
dup_num.ls = unlist(lapply(1:dim(KNN.m)[1], function(i){
  if(i==1) return(0)
  a = max(KNN.m[1:(i-1),] %*% KNN.m[i,])
  return(a)
}))
KNN.df = data.frame("barcodeBA"=E16.Npn$barcodeBA,
                    "coord"=E16.Npn$coord,
                    "Ann_v5"=E16.Npn$Ann_v5,
                    "Top_Ann"=NA,
                    "duplicate" = dup_num.ls,stringsAsFactors = F)
for (cti in ct.ls) {
  KNN.df[,paste0("numC",cti)] = KNN.m %*% (1*(E16.Npn$Ann_v5==cti))
  KNN.df$Top_Ann[KNN.df[,paste0("numC",cti)] >= 4] = cti
}
#   2   6   7   8  10  12  14  20  24  25  27  28  29 NA
# 269 127  83 315 539 184  33  27  88  34  16  39  94 17

# filter duplicates
KNN.df = KNN.df[(KNN.df$duplicate < 6) & (!is.na(KNN.df$Top_An)),]
KNN.m = KNN.m[KNN.df$barcodeBA,] # 1566 1865

## 3 UMAP embedding ##
mKid.merge.Npns$UMAP1 = mKid.merge.Npns@reductions$umap@cell.embeddings[,1]
mKid.merge.Npns$UMAP2 = mKid.merge.Npns@reductions$umap@cell.embeddings[,2]
KNN.df$UMAP1 = mKid.merge.Npns$UMAP1[KNN.df$coord]
KNN.df$UMAP2 = mKid.merge.Npns$UMAP2[KNN.df$coord]

write.csv(KNN.m,"Out/1_KNN_matrix.Npns.filtered.csv")
write.csv(KNN.df,"Out/1_KNN_metadata.Npns.filtered.csv")

