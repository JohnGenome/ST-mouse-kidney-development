# Ref:
# https://jdblischak.github.io/singleCellSeq/analysis/cell-cycle.html
# https://github.com/mahmoudibrahim/KidneyMap/blob/bd00adae97d49f15eb1733ec6505bbbd45dfa713/make_intergrated_maps/human_CD10positive.r
# https://github.com/mahmoudibrahim/KidneyMap/blob/bd00adae97d49f15eb1733ec6505bbbd45dfa713/demos/ecm_score_demo.r

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

# 2 cell cycle
flexible_normalization <- function(data_in,by_row=TRUE){
  if(by_row){
    row_mean <- apply(data_in,1,mean)
    row_sd   <- apply(data_in,1,sd)
    output <- data_in
    for(i in 1:dim(data_in)[1]){
      output[i,] <- (data_in[i,] - row_mean[i])/row_sd[i]
    }
  }
  if(!by_row){
    col_mean <- apply(data_in,2,mean)
    col_sd   <- apply(data_in,2,sd)
    output <- data_in
    for(i in 1:dim(data_in)[2]){
      output[,i] <- (data_in[,i] - col_mean[i])/col_sd[i]
    }
  }
  output
}

G1_S = read.table("/path/to/G1S")[[1]]
G2_M = read.table("/path/to/G2M")[[1]]
G0 = read.table("/path/to/G0")[[1]]
cycle_list = list(G1_S, G2_M, G0);rm(list = c("G1_S","G2_M","G0"))

log1pRPS.m = log1p(as.matrix(mKid.merge@assays$SCT@data))
phase_score = lapply(cycle_list,function(xx){
  log1pRPS.m <- log1pRPS.m[rownames(log1pRPS.m) %in% unlist(xx),]
  log1pRPS.m <- log1pRPS.m[apply(log1pRPS.m,1,sum)>50,]
  combined_matrix <- rbind(log1pRPS.m,average=apply(log1pRPS.m,2,mean)) # Filter genes inconsistent with average expression
  cor_matrix <- cor(t(combined_matrix))
  cor_vector <- cor_matrix[,dim(cor_matrix)[1]]
  log1pRPS.m_restricted <- log1pRPS.m[rownames(log1pRPS.m) %in% names(cor_vector[cor_vector >= 0.1]),]
  print(sum(cor_vector >= 0.1))
  print(sum(cor_vector < 0.1)) # ~20%
  apply(log1pRPS.m_restricted,2,mean)
})

phase_score = matrix(unlist(phase_score), ncol = length(cycle_list), byrow = FALSE)
colnames(phase_score) = c("G1_S","G2_M","G0")
phase_score_normed <- flexible_normalization(phase_score,by_row=FALSE)
phase_score_normed_normed <- flexible_normalization(phase_score_normed,by_row=TRUE)
cell_phase = apply(phase_score_normed_normed,1,function(x) colnames(phase_score_normed_normed)[which.max(x)])
table(cell_phase, Idents(mKid.merge))
# cell_phase    0    1    2    3    4    5    6    7    8    9   10   11   12
# G0   2181 1933 1203  688  862  922  615  188  187  961  413  516  598
# G1_S  629  618  805  680  681  583  668  497  618  420  775  483  486
# G2_M  437  274  635  833  805  470  546 1035  852  299  478  539  373
# 
# cell_phase   13   14   15   16   17   18   19   20   21   22   23
# G0    276  447  153  107  200  244  302  330  151   78   91
# G1_S  411  392  333  194  227  155  173  156  180  207   35
# G2_M  624  310  596  404  223  249  174  110  251  256   35

# 3 save
mKid.merge$G0 = phase_score[,"G0"]
mKid.merge$G1_S = phase_score[,"G1_S"]
mKid.merge$G2_M = phase_score[,"G2_M"]
mKid.merge$G0_norm = phase_score_normed_normed[,"G0"]
mKid.merge$G1_S_norm = phase_score_normed_normed[,"G1_S"]
mKid.merge$G2_M_norm = phase_score_normed_normed[,"G2_M"]
mKid.merge$cell_phase = cell_phase

# 4 Plot
pdf("Plot/5_cellcycle_score_UMAP.pdf")
DimPlot(mKid.merge, reduction = "umap", label = F,group.by = "cell_phase",shuffle = T)
FeaturePlot(mKid.merge, features = "G0_norm",min.cutoff='q1', max.cutoff='q99') + scale_colour_gradientn(colours = c("blue3","#F0F0F0","red3"))
FeaturePlot(mKid.merge, features = "G1_S_norm",min.cutoff='q1', max.cutoff='q99') + scale_colour_gradientn(colours = c("blue3","#F0F0F0","red3"))
FeaturePlot(mKid.merge, features = "G2_M_norm",min.cutoff='q1', max.cutoff='q99') + scale_colour_gradientn(colours = c("blue3","#F0F0F0","red3"))
VlnPlot(mKid.merge, features = c("G0_norm","G1_S_norm","G2_M_norm"),stack = T,flip = T) + NoLegend()
dev.off()

