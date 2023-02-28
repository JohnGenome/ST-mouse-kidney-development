library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

pseudotime.df = read.csv("Plot/4_monocle.CD_system.trajectory.pseudotime.csv",row.names = 1)
pseudotime.df = pseudotime.df[intersect(row.names(pseudotime.df),Cells(mKid.CDs)),,drop=FALSE]
mKid.CDs$pseudotime = -1
mKid.CDs@meta.data[row.names(pseudotime.df),"pseudotime"] = pseudotime.df$V1
mKid.CDs = mKid.CDs[,mKid.CDs$pseudotime > -1 & mKid.CDs$pseudotime < 27]

# 2
RPS.m = as.matrix(mKid.CDs@assays$SCT@data)
RPS.m = RPS.m[-grep("^mt\\.",row.names(RPS.m)),] # rm MT.* genes
RPS.m = RPS.m[apply(RPS.m,1,sum)>4000 | row.names(RPS.m) %in% c("Ret","Grip1","Slc14a2"),] # 2569 + Ret
RPS.m = RPS.m[apply(RPS.m,1,sd)>2 | row.names(RPS.m) %in% c("Ret","Grip1","Slc14a2"),] # 2391 + Ret
RPS.m = RPS.m[order(apply(RPS.m,1,function(i){
  k = sum(i>0)
  cutoff = sort(i,decreasing = T)[round(min(k/10,length(i)*0.03))]
  # cutoff = quantile(i,0.93) # top 296 cells
  return(sd(mKid.CDs$pseudotime[i>cutoff]))
})),] # rm unconcentrated genes
RPS.m = RPS.m[1:50,]
RPS.m = RPS.m[order(apply(RPS.m,1,function(i){
  k = sum(i>0)
  cutoff = sort(i,decreasing = T)[round(min(k/10,length(i)*0.03))]
  return(median(mKid.CDs$pseudotime[i>cutoff]))
})),]
RPS.m = RPS.m[,order(mKid.CDs$pseudotime)]

library("viridis")
library(RColorBrewer)
pdf("Plot/4a_monocle.CD_system.dynamic_genes.CDs_color.pdf",width = 16,height = 10)
Heatmap(t(scale(t(RPS.m))),cluster_columns = F,cluster_rows = F,show_row_names = T,show_column_names = F)
tmp=t(scale(t(RPS.m)));tmp[tmp>2]=2;tmp[tmp<-2]=-2;Heatmap(tmp,cluster_columns = F,cluster_rows = F,show_row_names = T,show_column_names = F)
Heatmap(t(scale(t(RPS.m))),cluster_columns = F,cluster_rows = F,show_row_names = T,show_column_names = F,
        col = circlize::colorRamp2(c(-2,0,2), c("blue","white","red")))
Heatmap(t(scale(t(RPS.m))),cluster_columns = F,cluster_rows = F,show_row_names = T,show_column_names = F,
        col = circlize::colorRamp2(seq(2,-2,length.out=9), brewer.pal(n=9,name="RdYlBu")))
Heatmap(t(scale(t(RPS.m))),cluster_columns = F,cluster_rows = F,show_row_names = T,show_column_names = F,
        col = circlize::colorRamp2(seq(-2,2,length.out=9), viridis(9)))
dev.off()

# 3
pdf("Plot/4b_monocle.CD_system.dynamic_genes.dot_plot.pdf",width = 8,height = 3)
for (genei in row.names(RPS.m)) { # c("Ret","Grip1","Slc14a2")
  plot.df = data.frame("pseudotime"=mKid.CDs$pseudotime,"Expr"=RPS.m[genei,mKid.CDs$coord],"cell_type"=mKid.CDs$Ann_v5)
  p1 = ggplot(plot.df,aes(x=pseudotime,y=Expr)) +
    geom_point(aes(color=cell_type),size=1) + scale_color_manual(values=c("15"="#C77CFF","0"="#7CAE00","30"="#006400")) +
    theme_bw() + xlab(genei) + geom_smooth(method = "loess",colour = "black",formula = "1.5*y ~ x")
  print(p1)
}
dev.off()


