library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

RPS.m = as.matrix(mKid.merge.global@assays$SCT@data)
RPS.m = RPS.m[row.names(RPS.m) %in% tf.ls,] # 1118 tfs
rowsum.ls = apply(RPS.m, 1, sum)
rowsum.ls = sort(rowsum.ls[rowsum.ls>2000],decreasing = T)
RPS.m = RPS.m[names(rowsum.ls),] # order genes according to total UMI

# 2 merge bins
bin.m = Reduce(cbind,lapply(1:37, function(i){
  RPS.i = RPS.m[,mKid.merge.global$dep_bin==i]
  return(apply(RPS.i, 1, mean))
}))
colnames(bin.m) = 1:37
bin.m = bin.m[apply(bin.m,1,max) > 3*apply(bin.m,1,min),] # 167

# 3 k-means
bin.norm = t(scale(t(bin.m)))
kmeans_gradient <- function(k){
  km = kmeans(bin.norm,centers = k,iter.max = 20,nstart = 2,algorithm = "Hartigan-Wong")
  km.class = fitted(km,method=c("classes"))
  save(list = c("km","km.class"),file = paste0("Plot_global_bin/kms/km",k,".RData"))
  bin.norm.reorder = bin.norm[sort.int(km.class,index.return = T)$ix,]
  pdf(paste0("Plot_global_bin/kms/km",k,"_heatmap.pdf"),width = 20,height = 30)
  htmap = Heatmap(bin.norm.reorder,cluster_columns = F,cluster_rows = F,show_row_names = T,column_split = rep(c("L1","L2","L3","L4","L5"),c(2,10,7,10,8)),column_gap = unit(4,"mm"))
  print(htmap)
  dev.off()
}
lapply(seq(4,10,1),kmeans_gradient)
load("Plot_global_bin/kms/km6.RData")
# table(km.class)
# km.class
# 1  2  3  4  5  6
# 39 13 25 37 24 29

bin.norm.reorder = bin.norm[sort.int(km.class,index.return = T)$ix,]
pdf("Plot_global_bin/1_km6_heatmap.pdf",width = 20,height = 30)
htmap = Heatmap(bin.norm.reorder,cluster_columns = F,cluster_rows = F,show_row_names = T,column_split = rep(c("L1","L2","L3","L4","L5"),c(2,10,7,10,8)),column_gap = unit(4,"mm"))
print(htmap)
dev.off()

# 3a k-means sort
name_i = list("C1"=intersect(names(rowsum.ls),names(km.class)[km.class==4]),
              "C2"=intersect(names(rowsum.ls),names(km.class)[km.class==2]),
              "C3"=intersect(names(rowsum.ls),names(km.class)[km.class==6]),
              "C4"=intersect(names(rowsum.ls),names(km.class)[km.class==5]),
              "C5"=intersect(names(rowsum.ls),names(km.class)[km.class==3]),
              "C6"=intersect(names(rowsum.ls),names(km.class)[km.class==1]))
final_genes = data.frame("gene"=Reduce(c,name_i),"cluster"=rep(names(name_i),unlist(lapply(name_i,length))))

bin.norm.reorder = bin.norm[final_genes$gene,]
pdf("Plot_global_bin/1a_km6_heatmap.reorder.pdf",width = 18,height = 30)
htmap = Heatmap(bin.norm.reorder,cluster_columns = F,cluster_rows = F,show_row_names = T,
                row_split = final_genes$cluster,row_gap = unit(2,"mm"),
                column_split = rep(c("L1","L2","L3","L4","L5"),c(2,10,7,10,8)),column_gap = unit(4,"mm"))
print(htmap)
dev.off()

# 4 mannul order
bin.norm = t(scale(t(bin.m)))
bin.norm = bin.norm[order(unlist(lapply(1:dim(bin.norm)[1],
                                        function(i){
                                          median(order(bin.norm[i,],decreasing = T)[1:3])
                                        }))),]
pdf("Plot_global_axon/1_layer_ordered_heatmap.median_top3.pdf",width = 7,height = 10)
htmap = Heatmap(bin.norm,cluster_columns = F,cluster_rows = F,show_row_names = T,column_split = rep(c("L1","L2","L3","L4","L5"),c(2,10,7,10,8)),column_gap = unit(2,"mm"))
print(htmap)
dev.off()
