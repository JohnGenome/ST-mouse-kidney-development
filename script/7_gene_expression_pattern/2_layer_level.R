library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

mKid.merge.CD = mKid.merge[,mKid.merge$orig.ident=="E18_3" & mKid.merge$Ann_v5 %in% c("0","30")]

RPS.m = as.matrix(mKid.merge.CD@assays$SCT@data)
rowsum.ls = apply(RPS.m, 1, sum)
rowsum.ls = sort(rowsum.ls[rowsum.ls>600],decreasing = T)
RPS.m = RPS.m[names(rowsum.ls),] # order genes according to total UMI

# 2 merge layers
layer.m = Reduce(cbind,lapply(2:5, function(i){
  RPS.i = RPS.m[,mKid.merge.CD$layer==i]
  return(apply(RPS.i, 1, mean))
}))
colnames(layer.m) = 2:5
layer.m = layer.m[apply(layer.m,1,max) > 2.5*apply(layer.m,1,min),] # 144

m1 = row.names(layer.m[layer.m[,2]>layer.m[,1] &
                         layer.m[,3]>layer.m[,2] & 
                         layer.m[,4]>layer.m[,3],]) # 23
m2 = row.names(layer.m[layer.m[,2]>layer.m[,1] &
                         layer.m[,3]>layer.m[,2] & 
                         layer.m[,4]<layer.m[,3],]) # 21
m3 = row.names(layer.m[layer.m[,2]>layer.m[,1] &
                         layer.m[,3]<layer.m[,2] & 
                         layer.m[,4]<layer.m[,3],]) # 26
m4 = row.names(layer.m[layer.m[,2]<layer.m[,1] &
                         layer.m[,3]<layer.m[,2] & 
                         layer.m[,4]<layer.m[,3],]) # 24
layer.norm = t(scale(t(layer.m)))[c(m1,m2,m3,m4),]
pdf(paste0("Plot_CD/1_heatmap.manual_select.pdf"),width = 5,height = 14)
htmap = Heatmap(layer.norm,cluster_columns = F,cluster_rows = F,show_row_names = T)
print(htmap)
dev.off()
final_genes = data.frame("gene"=c(m1,m2,m3,m4),"layer"=c(rep(5,length(m1)),
                                                         rep(4,length(m2)),
                                                         rep(3,length(m3)),
                                                         rep(2,length(m4))))
