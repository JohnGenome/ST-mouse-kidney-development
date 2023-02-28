# Ref:
# https://jdblischak.github.io/singleCellSeq/analysis/cell-cycle.html
# https://github.com/mahmoudibrahim/KidneyMap/blob/bd00adae97d49f15eb1733ec6505bbbd45dfa713/make_intergrated_maps/human_CD10positive.r
# https://github.com/mahmoudibrahim/KidneyMap/blob/bd00adae97d49f15eb1733ec6505bbbd45dfa713/demos/ecm_score_demo.r

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

# 2 ECM
ECM = read.table("/path/to/ECM")[[1]]
cycle_list = list(ECM);rm(list = c("ECM"))

log1pRPS.m = log1p(as.matrix(mKid.merge@assays$SCT@data))
phase_score = lapply(cycle_list,function(xx){
  log1pRPS.m <- log1pRPS.m[rownames(log1pRPS.m) %in% unlist(xx),]
  log1pRPS.m <- log1pRPS.m[apply(log1pRPS.m,1,sum)>50,]
  combined_matrix <- rbind(log1pRPS.m,average=apply(log1pRPS.m,2,mean)) # Filter genes inconsistent with average expression
  cor_matrix <- cor(t(combined_matrix))
  cor_vector <- cor_matrix[,dim(cor_matrix)[1]]
  log1pRPS.m_restricted <- log1pRPS.m[rownames(log1pRPS.m) %in% names(cor_vector[cor_vector >= 0.1]),]
  print(sum(cor_vector >= 0.1))
  print(sum(cor_vector < 0.1))
  apply(log1pRPS.m_restricted,2,mean)
})

# 3 save
mKid.merge$ECM = phase_score[[1]]

# 4 Plot
pdf("Plot_ECM/1_ECM_score_UMAP.pdf")
FeaturePlot(mKid.merge, features = "ECM",min.cutoff='q1', max.cutoff='q99') + scale_colour_gradientn(colours = c("blue3","#F0F0F0","red3"))
VlnPlot(mKid.merge, features = "ECM",sort = T) + NoLegend() + scale_fill_manual(values = color.ls5)
dev.off()

# 2D map
pdf("Plot_ECM/2_ECM_score_spatial_map.pdf",width = 6.3,height = 6.3)
maxi = 1.25
mKid.merge$ECM2 = pmin(mKid.merge$ECM,1.25)
p1 = ggplot(mKid.merge@meta.data[mKid.merge@meta.data$orig.ident=="E11_1",],aes(x=nB,y=nA,color=ECM2)) +
  geom_point(shape=15,size=1.5,show.legend=T) + theme_minimal() +
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB","E11_1")) +
  expand_limits(y=c(0,96),x=c(0,96)) +
  scale_colour_gradient(low = "#F0F0F0",high = "blue4",limits = c(0,maxi))
print(p1)
for (Samplei in c("E11_1","E11_2","E12_1","E12_2","E12_3","E12_4","E13_1","E13_2","E15_1","E15_2","E15_3","E15_4","E16_1","E18_1","E18_2","E18_3","E18_4")) {
  p1 = ggplot(mKid.merge@meta.data[mKid.merge@meta.data$orig.ident==Samplei,],aes(x=nB,y=nA,color=ECM2)) +
    geom_point(shape=15,size=1.5,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,96),x=c(0,96)) +
    scale_colour_gradient(low = "#F0F0F0",high = "blue4",limits = c(0,maxi))
  print(p1)
}
dev.off()




