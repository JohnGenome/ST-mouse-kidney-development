# ====================================================================
# 项目名称： 基于空间邻域的空间转录组功能配体-受体对推断软件
# 脚本名称： neighborhood_based_spatial_LR_inference.R
# 作者： 喻昊
# 日期： 2024-11-08
# 版本： 1.0
# 描述： 
#     ** 属于生物信息学软件，处理空间转录组文库数据 **
#     本软件基于图模型和基于空间邻域的置换检验方法，在空间转录组数据集
#     中实现空间邻近细胞类型对的检测、细胞类型对之间功能配体-受体对的
#     推断、计算结果的可视化。通过引入空间邻域的概念（包含细胞类型对及
#     其周边细胞）弥补了空间转录组数据稀疏性的缺点，可用于包括但不限于
#     DBiT-seq空间转录组数据的定量分析和深度挖掘。
# 依赖包：
#     - parallel
# ====================================================================

# -------------------- 1. 环境设置 --------------------

# 清空工作空间
rm(list = ls())

# 设置工作目录（可选）
# setwd("path/to/your/directory")

# 加载必要的包
library(ggplot2)
library(dplyr)

# -------------------- 2. 定义函数 --------------------

# 函数1：根据距离计算任意两种细胞类型之间距离小于2*10um的细胞对的数量
get_neighbor_count <- function(meta.df,dist.m,clus.ls,permut=FALSE){
  # 确定细胞类型种类
  sam.ls = unique(meta.df$orig.ident)
  # 如果permut==T，则打乱细胞类型标签，统计随机情况的结果
  if(permut){ 
    new_order = unlist(Reduce(c,lapply(sam.ls,function(sami) sample(which(meta.df$orig.ident==sami)))))
    meta.df$Ann = meta.df$Ann[new_order]
  }
  dist.df = Reduce(rbind,lapply(sam.ls, function(sami) data.frame("cti"=factor(rep(meta.df$Ann[meta.df$orig.ident==sami], times=sum(meta.df$orig.ident==sami)),levels = clus.ls), # cautious 
                                                                  "ctj"=factor(rep(meta.df$Ann[meta.df$orig.ident==sami], each=sum(meta.df$orig.ident==sami)),levels = clus.ls), # !!!
                                                                  "dist"=as.vector(dist.m[meta.df$orig.ident==sami,meta.df$orig.ident==sami]))))
  dist.df = dist.df[dist.df$dist < 2,] # 保留距离小于2*10um的细胞对
  neighbor.df = as.data.frame(table(dist.df$cti,dist.df$ctj))
  colnames(neighbor.df) = c("cti","ctj","count")
  return(neighbor.df)
}


# -------------------- 3. 数据导入与预处理 --------------------

# 导入数据
sample_i = "your_sample_name"
cut_off = 30 # 细胞数量低于阈值的细胞类型将被过滤
meta.df = read.csv("path/to/your/file",stringsAsFactors = F,row.names = 1)

# 数据清洗与预处理
cluster.ls = table(meta.df$Ann)
cluster.ls = names(cluster.ls)[cluster.ls > cut_off]
meta.df = meta.df[meta.df$Ann %in% cluster.ls,]

# -------------------- 4. Step1 图模型构建 --------------------

# 1 Prepare distance matrix
distA.m = meta.df$nA %*% matrix(1,nrow = 1,ncol = dim(meta.df)[1]) - matrix(1,nrow = dim(meta.df)[1],ncol = 1) %*% t(meta.df$nA)
distB.m = meta.df$nB %*% matrix(1,nrow = 1,ncol = dim(meta.df)[1]) - matrix(1,nrow = dim(meta.df)[1],ncol = 1) %*% t(meta.df$nB)
dist.m = sqrt(distA.m**2 + distB.m**2)
dist.m = dist.m + diag(100,nrow = dim(meta.df)[1])
rm(list = c("cut_off","distA.m","distB.m"))
print(paste0(sample_i,": Step 1 finished"))

# 2 Permutation test
times = 2000
neighbor.df = get_neighbor_count(meta.df,dist.m,cluster.ls)
system.time(random.m <- mclapply(1:times, function(i) get_neighbor_count(meta.df,dist.m,cluster.ls,permut=T)$count,mc.cores=47)) # cost much times!!!
# user    system   elapsed
# 41278.271 23987.661  2112.267
random.m = Reduce(cbind,random.m)
random.m = as.matrix(random.m)

neighbor.df$expect = apply(random.m,1,mean)
neighbor.df$pvalue = apply((random.m-neighbor.df$count)>0,1,sum)/times
neighbor.df$oe = neighbor.df$count/neighbor.df$expect

rm(list = c("dist.m","meta.df","random.m","times"))
base::save.image(file = "path/to/output/0_cell_type_neighbor_permutation.RData")
write.csv(neighbor.df,file = "path/to/output/1_save_neighboring.csv",row.names = F)

# 3 Visualization
for (pcut_off in c(0.01,0.05)) {
  # 3 Heatmap
  oe.m = matrix(0,length(cluster.ls),length(cluster.ls))
  row.names(oe.m) = levels.ls[cluster.ls]
  colnames(oe.m) = levels.ls[cluster.ls]
  tmp <- lapply(1:dim(neighbor.df)[1], function(i){
    if(neighbor.df$pvalue[i] < pcut_off){
      oe.m[neighbor.df$cti[i],neighbor.df$ctj[i]] <<- neighbor.df$oe[i]
    }
  })
  pdf(paste0("path/to/output/2_oe_dist_between_cell_type.",pcut_off,".pdf"),width = 7,height = 7)
  print(Heatmap(oe.m,cluster_rows=T,cluster_columns=T,rect_gp = gpar(col = "white", lwd = 1),
                col = colorRamp2(c(0, 2, 5), c("blue", "white", "red")),
                cell_fun = function(j, i, x, y, width, height, fill){
                  if(!is.na(oe.m[i, j])) if(oe.m[i, j] > 0) grid.text(sprintf("%.2f", oe.m[i, j]), x, y, gp = gpar(fontsize = 5))
                }))
  dev.off()
  # 4 Network
  neighbor.df2 = neighbor.df[!duplicated(neighbor.df$expect) | !duplicated(neighbor.df$count),]
  neighbor.df2 = neighbor.df2[neighbor.df2$cti != neighbor.df2$ctj,]
  neighbor.df2 = neighbor.df2[neighbor.df2$pvalue < pcut_off,]
  oe.graph = graph.data.frame(neighbor.df2,directed = F)
  pdf(paste0("path/to/output/3_neighboring_cell_type_network.",pcut_off,".pdf"))
  print(plot(oe.graph, edge.width = -log10(neighbor.df2$pvalue+0.00001),
             edge.color = colorRamp2(c(0, 4), c("grey","red"))(neighbor.df2$oe),
             vertex.size = 5,vertex.label.cex=.7,vertex.label.dist=1))
  print(plot(oe.graph, edge.width = -log10(neighbor.df2$pvalue+0.00001),
             edge.color = colorRamp2(c(0, 4), c("grey","red"))(neighbor.df2$oe),
             vertex.size = 5,vertex.label.cex=.7,vertex.label.dist=1))
  dev.off()
}

# -------------------- 5. Step2 筛选成对的核心差异表达基因 --------------------

LR.df = read.table("path/to/mouse_ligand_receptors.txt")
DEG.df = read.csv("path/to/7_de_genes.AnnV5.filterpixel.csv",stringsAsFactors = F)
neighbor.df = read.csv("path/to/output/1_save_neighboring.csv",stringsAsFactors = F)

## get core DEGs ##
#1 avg_log2FC >= 0.8
#2 Top2 cluster for each gene
#3 avg_log2FC >= 2/3*(Top1 avg_log2FC)
DEG.df$ID = paste(DEG.df$cluster,DEG.df$gene,sep = "_")
DEG.df$is.core = 0
coreDEG.df = DEG.df[DEG.df$avg_log2FC>=0.8,]
coreDEG.df$top2 = tapply(coreDEG.df$avg_log2FC,coreDEG.df$gene,
                         function(x) ifelse(length(x)>1,sort(x,decreasing = T)[2],sort(x,decreasing = T)[1]))[coreDEG.df$gene]
coreDEG.df$two.thirds.top1 = tapply(coreDEG.df$avg_log2FC,coreDEG.df$gene,function(x) 2/3*max(x))[coreDEG.df$gene]
coreDEG.df$is.core[(coreDEG.df$avg_log2FC>=coreDEG.df$top2) & (coreDEG.df$avg_log2FC>=coreDEG.df$two.thirds.top1)] = 1
DEG.df$is.core[DEG.df$ID %in% coreDEG.df$ID[coreDEG.df$is.core==1]] = 1 # 1875 core and 6088 non-core DEGs

write.csv(DEG.df,"path/to/output/7a_de_genes.AnnV5.filterpixel.iscore.csv",row.names = F)

## get core DEGs as L or R ##
DEG.df = DEG.df[DEG.df$gene %in% c(LR.df$V1,LR.df$V2),] # 650 (197 core and 453 non-core)

neighbor.df = neighbor.df[(neighbor.df$pvalue<0.01) & (neighbor.df$cti != neighbor.df$ctj),]
levels.ls = 0:32
names(levels.ls) = c("0_oCD","1_iMFib","2_aLoH","3_Ukn","4_UEndo","5_oMFib","6_PTprog","7_cNP","8_NP","9_mMfib",
                     "10_PT","11_Capil","12_Pod","13_oCFib","14_DT","15_UBT","16_CUkn","17_UFib","18_TC","19_LA",
                     "20_dLoH","21_USMC","22_Uprog","23_RBC","24_Dprog","25_Podprog","26_Mesa","27_Corin","28_preDT","29_MD",
                     "30_iCD","31_Utran","32_iCFib")
neighbor.df$cti = levels.ls[neighbor.df$cti]
neighbor.df$ctj = levels.ls[neighbor.df$ctj]

get_neighbor_partner <- function(id1,DEG.d,neighbor.d,LR.d){
  ci = DEG.d$cluster[DEG.d$ID==id1]
  gi = DEG.d$gene[DEG.d$ID==id1]
  cj = neighbor.d$ctj[neighbor.d$cti==ci]
  gj = c(LR.d$V2[LR.d$V1==gi],LR.d$V1[LR.d$V2==gi])
  subDEG = DEG.d[(DEG.d$cluster %in% cj) & (DEG.d$gene %in% gj),]
  typei = ifelse(gi %in% LR.d$V1,"L","R")
  df1 = data.frame("IDi"=id1,"Clusteri"=ci,"Genei"=gi,"log2FCi"=DEG.d$avg_log2FC[DEG.d$ID==id1],"corei"=1)
  df2 = data.frame("IDj"=subDEG$ID,"Clusterj"=subDEG$cluster,"Genej"=subDEG$gene,"log2FCj"=subDEG$avg_log2FC,"corej"=subDEG$is.core)
  if (dim(df2)[1]==0) df2 = data.frame("IDj"=NA,"Clusterj"=NA,"Genej"=NA,"log2FCj"=NA,"corej"=NA)
  df2 = df2[order(df2$log2FCj,decreasing = T),] # for each LR pair, select best matched cell type
  df2 = df2[!duplicated(df2$Genej),]
  if (typei=="L") {
    final.d = cbind(df1,df2)
  }else{
    final.d = cbind(df2,df1)
  }
  final.d$neigh_p = neighbor.d$pvalue[match(paste(final.d$Clusteri,final.d$Clusterj),
                                            paste(neighbor.d$cti,neighbor.d$ctj))]
  final.d$neigh_oe = neighbor.d$oe[match(paste(final.d$Clusteri,final.d$Clusterj),
                                         paste(neighbor.d$cti,neighbor.d$ctj))]
  colnames(final.d) = c("ID_L","Cluster_L","Gene_L","log2FC_L","core_L",
                        "ID_R","Cluster_R","Gene_R","log2FC_R","core_R","neigh_p","neigh_oe") # unify colnames for later rbind()
  return(final.d)
}
# DEG.df[DEG.df$is.core==1,c(6,7,8,9)][1:10,]
# get_neighbor_partner("8_Gdnf",DEG.df,neighbor.df,LR.df)
final.df = Reduce(rbind,lapply(DEG.df$ID[DEG.df$is.core==1],
                               get_neighbor_partner,DEG.d=DEG.df,neighbor.d=neighbor.df,LR.d=LR.df)) # 227
final.df = final.df[!(duplicated(final.df)),] # 204 (54 only L + 64 only R + 86 pairs)
write.csv(final.df,"path/to/output/1_high_confidence_LR_pairs.at_least_one_core_deg.neighbor_based.top1_celltype.csv",row.names = F)

# -------------------- 6. Step3 基于邻域的置换检验 --------------------

library(Seurat)
library(parallel)
library(circlize)
library(igraph)
library(ggplot2)
library(paletteer)
library("viridis")

# Functions
get_neighbor_cells <- function(meta.df,dist.m,cti,ctj,permut=FALSE){
  dist.df = data.frame("celli"=rep(meta.df$coord, times=dim(meta.df)[1]),
                       "cellj"=rep(meta.df$coord, each=dim(meta.df)[1]),
                       "cti"=rep(meta.df$Ann, times=dim(meta.df)[1]),
                       "ctj"=rep(meta.df$Ann, each=dim(meta.df)[1]),
                       "dist"=as.vector(dist.m),stringsAsFactors = F)
  dist.df = dist.df[dist.df$dist < 2 & dist.df$cti == cti & dist.df$ctj == ctj,]
  return(list(unique(dist.df$celli),unique(dist.df$cellj)))
}

get_neighborhood <- function(meta.df,dist.m,nei_cell,nbhood=FALSE,cti=NA,ctj=NA,permut=FALSE){
  if(!nbhood){
    nei_cell = c(nei_cell[[1]],nei_cell[[2]])
    dist.df = data.frame("celli"=rep(meta.df$coord, times=dim(meta.df)[1]),
                         "cellj"=rep(meta.df$coord, each=dim(meta.df)[1]),
                         "cti"=rep(meta.df$Ann, times=dim(meta.df)[1]),
                         "ctj"=rep(meta.df$Ann, each=dim(meta.df)[1]),
                         "dist"=as.vector(dist.m),stringsAsFactors = F)
    dist.df = dist.df[dist.df$dist < 2 & dist.df$cellj %in% nei_cell,]
    nei_cell= unique(dist.df$celli)
    return(nei_cell)
  }
  meta.df = meta.df[meta.df$coord %in% nei_cell,]
  sam.ls = unique(meta.df$orig.ident)
  if(permut){
    new_order = unlist(Reduce(c,lapply(sam.ls,function(sami) sample(which(meta.df$orig.ident==sami)))))
    meta.df$Ann = meta.df$Ann[new_order]
  }
  return(list(meta.df$coord[meta.df$Ann==cti],meta.df$coord[meta.df$Ann==ctj]))
}

get_LR_expr_score2 <- function(count.m,nei_cell,only1pair.df){
  LR_score1 = apply(count.m[only1pair.df$Gene_L,nei_cell[[1]],drop=FALSE], 1, mean)
  LR_score2 = apply(count.m[only1pair.df$Gene_R,nei_cell[[2]],drop=FALSE], 1, mean)
  names(LR_score1) = paste0(only1pair.df$Gene_L,"_",only1pair.df$Gene_R)
  names(LR_score2) = paste0(only1pair.df$Gene_L,"_",only1pair.df$Gene_R)
  return(list(LR_score1,LR_score2))
}

get_LR_expr <- function(count.m,cti,ctj,only1pair.df){
  LR_01.1 = apply(count.m[only1pair.df$Gene_L,meta.df$coord[meta.df$Ann==cti],drop=FALSE], 1, mean)
  LR_01.2 = apply(count.m[only1pair.df$Gene_R,meta.df$coord[meta.df$Ann==ctj],drop=FALSE], 1, mean)
  LR_01 = apply(cbind(LR_01.1,LR_01.2), 1, min) # 1*(LR_01.1 > 0.1 & LR_01.2 > 0.1)
  names(LR_01) = paste0(only1pair.df$Gene_L,"_",only1pair.df$Gene_R)
  return(LR_01)
}

do_permutation_test <- function(meta.df,dist.m,only1pair.df) {
  neighbor_i=c(only1pair.df$Cluster_L[1],only1pair.df$Cluster_R[1])
  name_i=paste0(only1pair.df$Cluster_L[1],"_to_",only1pair.df$Cluster_R[1])
  print(paste("Start:",name_i))
  if(is.na(neighbor_i[1])|is.na(neighbor_i[2])){
    only1pair.df$permu_p = NA
    only1pair.df$permu_oe = NA
    only1pair.df$permu_expr = NA
    return(only1pair.df)
  }
  ## 1 Get neighborhood ##
  # for direct CM-UBT pairs
  neighbor_i_cell = get_neighbor_cells(meta.df,dist.m,neighbor_i[1],neighbor_i[2])
  meta.df$interact = "Others"
  meta.df$interact[meta.df$coord %in% neighbor_i_cell[[1]]] = neighbor_i[1]
  meta.df$interact[meta.df$coord %in% neighbor_i_cell[[2]]] = neighbor_i[2]
  # for neighborhood of direct CM-UBT pairs
  neighborhood_i = get_neighborhood(meta.df,dist.m,neighbor_i_cell)
  meta.df$nbhood = "Others"
  meta.df$nbhood[meta.df$coord %in% neighborhood_i] = "Local"
  meta.df$nbhood[meta.df$coord %in% neighborhood_i & meta.df$Ann == neighbor_i[1]] = neighbor_i[1]
  meta.df$nbhood[meta.df$coord %in% neighborhood_i & meta.df$Ann == neighbor_i[2]] = neighbor_i[2]
  # for cell type label
  meta.df$shape = "background"
  meta.df$shape[meta.df$Ann==neighbor_i[1]] = neighbor_i[1]
  meta.df$shape[meta.df$Ann==neighbor_i[2]] = neighbor_i[2]
  pdf(paste0("path/to/output/",name_i,".1_spatial_map_overview.pdf"))
  color.ls = c("Others"="grey","Local"="green","L"="red","R"="blue")
  names(color.ls)[3:4] = neighbor_i
  shape.ls = c("background"=15,"L"=19,"R"=17)
  names(shape.ls)[2:3] = neighbor_i
  for (samplei in c("E15_4","E16_1","E18_3")) {
    p1 = ggplot(meta.df[meta.df$orig.ident==samplei,],aes(x=nB,y=nA,color=interact)) + 
      geom_point(aes(shape=shape),size=1.9,show.legend=F) + theme_void() + 
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
      expand_limits(y=c(0,96),x=c(0,96)) + 
      scale_colour_manual(values = color.ls) +
      scale_shape_manual(values = shape.ls)
    print(p1)
    p1 = ggplot(meta.df[meta.df$orig.ident==samplei,],aes(x=nB,y=nA,color=nbhood)) + 
      geom_point(aes(shape=shape),size=1.9,show.legend=F) + theme_void() + 
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
      expand_limits(y=c(0,96),x=c(0,96)) + 
      scale_colour_manual(values = color.ls) +
      scale_shape_manual(values = shape.ls)
    print(p1)
    p1 = ggplot(meta.df[meta.df$orig.ident==samplei,],aes(x=nB,y=nA,color=nbhood)) + 
      geom_point(aes(shape=shape),size=1.9) + theme_minimal() + 
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
      expand_limits(y=c(0,96),x=c(0,96)) + 
      scale_colour_manual(values = color.ls) +
      scale_shape_manual(values = shape.ls)
    print(p1)
  }
  dev.off()
  ## 2 permutation test ##
  neighbor_i_score = get_LR_expr_score2(correct_count.m,
                                        get_neighborhood(meta.df,dist.m,neighborhood_i,nbhood = T,
                                                         cti = neighbor_i[1],ctj = neighbor_i[2],permut = F),
                                        only1pair.df)
  neighbor_i_score.shuff = Reduce(function(i,j) return(list(cbind(i[[1]],j[[1]]),cbind(i[[2]],j[[2]]))),
                                  mclapply(1:5000,function(i){
                                    cell_k = get_neighborhood(meta.df,dist.m,neighborhood_i,nbhood = T,
                                                              cti = neighbor_i[1],ctj = neighbor_i[2],permut = T)
                                    return(get_LR_expr_score2(correct_count.m,cell_k,only1pair.df))
                                  },mc.cores = 4))
  neighbor_i_fc = (neighbor_i_score[[1]]/apply(neighbor_i_score.shuff[[1]],1,mean) +
                     neighbor_i_score[[2]]/apply(neighbor_i_score.shuff[[2]],1,mean))/2
  neighbor_i_pv = cbind(apply((neighbor_i_score.shuff[[1]]-neighbor_i_score[[1]])>0,1,sum)/5000,
                        apply((neighbor_i_score.shuff[[2]]-neighbor_i_score[[2]])>0,1,sum)/5000)
  neighbor_i_pv = apply(neighbor_i_pv, 1, max) # mean
  neighbor_i_expr = get_LR_expr(correct_count.m,neighbor_i[1],neighbor_i[2],only1pair.df)
  only1pair.df$permu_p = neighbor_i_pv
  only1pair.df$permu_oe = neighbor_i_fc
  only1pair.df$permu_expr = neighbor_i_expr
  ## 3 2D plot ##
  for (i in 1:dim(only1pair.df)[1]) {
    pdf(paste0("path/to/output/",name_i,".2_spatial_map_expression.",only1pair.df$Gene_L[i],"_",only1pair.df$Gene_R[i],".pdf"))
    meta.df$Li = correct_count.m[only1pair.df$Gene_L[i],meta.df$coord]
    meta.df$Ri = correct_count.m[only1pair.df$Gene_R[i],meta.df$coord]
    meta.df$Li = pmin(meta.df$Li,sort(meta.df$Li,decreasing = T)[30])
    meta.df$Ri = pmin(meta.df$Ri,sort(meta.df$Ri,decreasing = T)[30])
    p1 = ggplot(meta.df[meta.df$orig.ident=="E18_3",],aes(x=nB,y=nA,color=Li)) +
      geom_point(aes(shape=shape),size=1.9,show.legend=T) + theme_minimal() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",only1pair.df$Gene_L[i])) +
      expand_limits(y=c(0,96),x=c(0,96)) +
      scale_color_viridis(option = "B",begin = 0.1,end = 0.9) +
      # scale_colour_gradient(low = "grey",high = "blue4",limits = c(0,sort(meta.df$Li,decreasing = T)[30])) +
      scale_shape_manual(values = shape.ls)
    print(p1)
    p1 = ggplot(meta.df[meta.df$orig.ident=="E18_3",],aes(x=nB,y=nA,color=Ri)) +
      geom_point(aes(shape=shape),size=1.9,show.legend=T) + theme_minimal() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",only1pair.df$Gene_R[i])) +
      expand_limits(y=c(0,96),x=c(0,96)) +
      scale_color_viridis(option = "B",begin = 0.1,end = 0.9) +
      # scale_colour_gradient(low = "grey",high = "blue4",limits = c(0,sort(meta.df$Li,decreasing = T)[30])) +
      scale_shape_manual(values = shape.ls)
    print(p1)
    for (samplei in c("E15_4","E16_1","E18_3")) {
      p1 = ggplot(meta.df[meta.df$orig.ident==samplei,],aes(x=nB,y=nA,color=Li)) +
        geom_point(aes(shape=shape),size=1.9,show.legend=F) + theme_void() +
        scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
        scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",only1pair.df$Gene_L[i])) +
        expand_limits(y=c(0,96),x=c(0,96)) +
        scale_color_viridis(option = "B",begin = 0.1,end = 0.9) +
        # scale_colour_gradient(low = "grey",high = "blue4",limits = c(0,sort(meta.df$Li,decreasing = T)[30])) +
        scale_shape_manual(values = shape.ls)
      print(p1)
      p1 = ggplot(meta.df[meta.df$orig.ident==samplei,],aes(x=nB,y=nA,color=Ri)) +
        geom_point(aes(shape=shape),size=1.9,show.legend=F) + theme_void() +
        scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
        scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",only1pair.df$Gene_R[i])) +
        expand_limits(y=c(0,96),x=c(0,96)) +
        scale_color_viridis(option = "B",begin = 0.1,end = 0.9) +
        # scale_colour_gradient(low = "grey",high = "blue4",limits = c(0,sort(meta.df$Li,decreasing = T)[30])) +
        scale_shape_manual(values = shape.ls)
      print(p1)
    }
    dev.off()
  }
  print(name_i)
  return(only1pair.df)
}

## 0 Load data ##
final.df = read.csv("path/to/output/1_high_confidence_LR_pairs.at_least_one_core_deg.neighbor_based.top1_celltype.csv",stringsAsFactors = F)
# Seurat object
mKid.merge = readRDS("path/to/output/8_SCT_transform&Merged_data.only_embryo.depthnormed.Annv5.filtered.cellcycle.depth.RDS")
mKid.merge = mKid.merge[,mKid.merge$orig.ident %in% c("E15_4","E16_1","E18_3")]
meta.df = mKid.merge@meta.data[,c("orig.ident","Ann","coord","nA","nB")]
correct_count.m = as.matrix(mKid.merge@assays[["SCT"]]@data[sort(unique(c(final.df$Gene_L,final.df$Gene_R))),])

## 1 get distance between celli cellj ##
# meta.df = meta.df[order(meta.df$Ann),] # prepare work: shuffle cell coord while keep cell identity ; see also: get_neighbor_cells(permut=T)
distA.m = meta.df$nA %*% matrix(1,nrow = 1,ncol = dim(meta.df)[1]) - matrix(1,nrow = dim(meta.df)[1],ncol = 1) %*% t(meta.df$nA)
distB.m = meta.df$nB %*% matrix(1,nrow = 1,ncol = dim(meta.df)[1]) - matrix(1,nrow = dim(meta.df)[1],ncol = 1) %*% t(meta.df$nB)
dist.m = sqrt(distA.m**2 + distB.m**2)
dist.m = dist.m + diag(100,nrow = dim(meta.df)[1])
rm(list = c("distA.m","distB.m"))
meta.df$nA = meta.df$nA %% 200

## 2 Run ##
final_test.df = Reduce(rbind,lapply(split(final.df,f=paste(final.df$Cluster_L,final.df$Cluster_R)),
                                    do_permutation_test,meta.df=meta.df,dist.m=dist.m)) # about 2min for each turn
write.csv(final_test.df,"path/to/output/1a_high_confidence_LR_pairs.at_least_one_core_deg.neighbor_based.top1_celltype.permu_test.csv",row.names = F)

# do_permutation_test(meta.df=meta.df,dist.m=dist.m,only1pair.df=final.df[(final.df$Cluster_L==8) & (final.df$Cluster_R==15) 
#                                                                         & !(is.na(final.df$Cluster_L)) & !(is.na(final.df$Cluster_R)),])

final.df = read.csv("path/to/output/1a_high_confidence_LR_pairs.at_least_one_core_deg.neighbor_based.top1_celltype.permu_test.csv",stringsAsFactors = F)
load("path/to/CellChatDB.mouse.rda")
final.df$annotation = NA
final.df$annotation = CellChatDB.mouse$interaction$annotation[match(final.df$Gene_L,CellChatDB.mouse$interaction$ligand)]
write.csv(final.df,"path/to/output/1b_high_confidence_LR_pairs.at_least_one_core_deg.neighbor_based.top1_celltype.permu_test.csv",row.names = F)

# -------------------- 7. Step4 全局结果可视化 --------------------

library(parallel)
library(circlize)
library(igraph)
library(ggplot2)
library(paletteer)

final.df = read.csv("path/to/output/2_final.high_confidence_LR_pairs.one_methods_merged.csv",stringsAsFactors = F)
# Ke LoH
color.ls5 = c("#7CAE00","#ffdf00","#AD0000","#7F7F7F","#C09B00","#ff4500","#00BAE0","#A3A500","#000080","#90EE90",
              "#EA8331","#00B0F6","#FF0000","#704241","#FFFF00","#C77CFF","#daa520","#adff2f","#8D2282","#E76BF3",
              "#39B600","#9590FF","#FF6A98","#FA62DB","#4169E1","#FF9999","#B22222","#FF62BC","#FFD032","#A2FF00",
              "#006400","#c0c0c0","#004A00")
names(color.ls5) = 0:32

## 1 CC communication Graph ##
count.df = as.data.frame(table(final.df$Cluster_L,final.df$Cluster_R))
count.df = count.df[count.df$Freq>0,]
count.df = count.df[order(count.df$Var1,count.df$Var2),]
count.df$Var1 = as.character(count.df$Var1)
count.df$Var2 = as.character(count.df$Var2)
# table(count.df$Freq)
# 1 2 3 4 6
# 6 3 3 2 1

colnames(count.df) = c("cti","ctj","count")
levels.ls = c("0_oCD","1_iMFib","2_aLoH","3_Ukn","4_UEndo","5_oMFib","6_PTprog","7_cNP","8_NP","9_mMfib",
              "10_PT","11_Capil","12_Pod","13_oCFib","14_DT","15_UBT","16_CUkn","17_UFib","18_TC","19_LA",
              "20_dLoH","21_USMC","22_Uprog","23_RBC","24_Dprog","25_Podprog","26_Mesa","27_Corin","28_preDT","29_MD",
              "30_iCD","31_Utran","32_iCFib")
names(levels.ls) = 0:32
count.df$cti = levels.ls[count.df$cti]
count.df$ctj = levels.ls[count.df$ctj]
count.df$count = pmin(count.df$count,3)

oe.graph = graph.data.frame(count.df,directed = T)
pdf("path/to/output2_LR_pair_count.cell_type_network.pdf")
print(plot(oe.graph, edge.width = count.df$count,
           edge.color = colorRamp2(c(1, 3), c("grey","red"))(count.df$count),
           vertex.size = 5,vertex.label.cex=.7,vertex.label.dist=1))
print(plot(oe.graph, edge.width = count.df$count,
           edge.color = colorRamp2(c(1, 3), c("grey","red"))(count.df$count),
           vertex.size = 5,vertex.label.cex=.7,vertex.label.dist=1))
dev.off()

## 2 Chord Diagram ##
count.df = as.data.frame(table(final.df$Cluster_L,final.df$Cluster_R))
count.df = count.df[count.df$Freq>0,]
count.df = count.df[order(count.df$Var1,count.df$Var2),]
count.df$Var1 = as.character(count.df$Var1)
count.df$Var2 = as.character(count.df$Var2)
colnames(count.df) = c("cti","ctj","count")

pdf("path/to/output3_LR_pair_count.chordDiagram.pdf")
circos.clear()
chordDiagram(count.df,grid.col = color.ls5,transparency = 0.3,
             link.lwd = 1.5, # Line width
             link.lty = 1, # Line type
             link.border = ifelse(count.df$cti %in% c("13","32","5","9","1"),"black",NA), # Border color
             order = c("13","32","8","15","7","24",
                       "11","5","22","31"))
dev.off()

## 3 LR pair dot plot ##
plotDotLRs <- function(name.i,ct.ls1,ct.ls2){
  sub.df = final.df[(final.df$Cluster_L %in% ct.ls1) & (final.df$Cluster_R %in% ct.ls2),]
  sub.df = sub.df[order(factor(sub.df$Cluster_L,levels = ct.ls1),
                        factor(sub.df$Cluster_R,levels = ct.ls2)),]
  sub.df$ct_pair = paste0(sub.df$Cluster_L," -> ",sub.df$Cluster_R)
  sub.df$LR_pair = paste0(sub.df$Gene_L," - ",sub.df$Gene_R)
  sub.df$ct_pair = factor(sub.df$ct_pair,levels = unique(sub.df$ct_pair))
  sub.df$LR_pair = factor(sub.df$LR_pair,levels = unique(sub.df$LR_pair))
  sub.df$mergeFC = sub.df$permu_oe
  # sub.df$mergeFC[is.na(sub.df$mergeFC)] = sub.df$Wilco_oe[is.na(sub.df$mergeFC)]
  sub.df$logp = -1*log10(sub.df$permu_p + 0.0001)
  # sub.df$logp[is.na(sub.df$logp)] = (-1*log10(sub.df$Wilco_p + 0.0001))[is.na(sub.df$logp)]
  
  pdf(paste0("path/to/output4_LR_pair.dot_plot.",name.i,".pdf"),width = 2 + length(levels(sub.df$ct_pair))/2,height = 2 + length(levels(sub.df$LR_pair))/3)
  p1 <- ggplot(sub.df, aes(x = ct_pair, y = LR_pair, color = mergeFC, size = logp)) +
    geom_point(pch = 16, na.rm = TRUE) +
    scale_x_discrete(position = "top") + 
    # scale_color_gradientn(colors = paletteer_c("grDevices::Viridis", 30)) + 
    scale_radius(limits = c(1,4),range = c(3,8)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0,color = "#000000"), axis.text.y = element_text(color = "#000000"),
                       axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(p1)
  p1 <- ggplot(sub.df, aes(x = ct_pair, y = LR_pair, color = mergeFC, size = logp)) +
    geom_point(pch = 16, na.rm = TRUE) +
    scale_x_discrete(position = "top") + 
    scale_color_gradientn(colors = paletteer_c("grDevices::Viridis", 30),limits = c(0,max(sub.df$mergeFC))) + 
    scale_radius(limits = c(1,4),range = c(3,8)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0,color = "#000000"), axis.text.y = element_text(color = "#000000"),
                       axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(p1)
  dev.off()
}

plotDotLRs(name.i = "cortFib_progenitors",
           ct.ls1 = c("13","32","8","15","11"),
           ct.ls2 = c("13","32","8","15","11"))
plotDotLRs(name.i = "nephron",
           ct.ls1 = c("15","7","24","11"),
           ct.ls2 = c("15","7","24","11"))
plotDotLRs(name.i = "ureter",
           ct.ls1 = c("11","5","22","31"),
           ct.ls2 = c("11","5","22","31"))
plotDotLRs(name.i = "ureter_new",
           ct.ls1 = c("22","31"),
           ct.ls2 = c("22","31"))

## 4 manual revised LR pair dot plot ##
plotDotLRsNew <- function(name.i,must_include,may_include,focus=c("Vegfa","Vegfc"),only.focus=F){
  sub.df = Reduce(rbind,lapply(must_include, function(must_i){
    a = final.df[(final.df$Cluster_L %in% must_i) & (final.df$Cluster_R %in% c(must_include,may_include)),]
    b = final.df[(final.df$Cluster_L %in% c(must_include,may_include)) & (final.df$Cluster_R %in% must_i),]
    return(rbind(a,b))
  }))
  sub.df = sub.df[!duplicated(sub.df),]
  sub.df$ct_pair = paste0(sub.df$Cluster_L," -> ",sub.df$Cluster_R)
  sub.df$LR_pair = paste0(sub.df$Gene_L," - ",sub.df$Gene_R)
  sub.df$ct_pair = factor(sub.df$ct_pair,levels = unique(sub.df$ct_pair))
  
  LR.ls = Reduce(c,lapply(focus,grep,x=sub.df$LR_pair,value=T))
  if(only.focus){
    sub.df = sub.df[sub.df$LR_pair %in% LR.ls,]
    sub.df$LR_pair = factor(sub.df$LR_pair,levels = unique(LR.ls))
  }else{
    sub.df$LR_pair = factor(sub.df$LR_pair,levels = union(LR.ls,unique(sub.df$LR_pair)))
  }
  sub.df$mergeFC = sub.df$permu_oe
  sub.df$mergeFC[is.na(sub.df$mergeFC)] = sub.df$Wilco_oe[is.na(sub.df$mergeFC)]
  sub.df$logp = -1*log10(sub.df$permu_p + 0.0001)
  sub.df$logp[is.na(sub.df$logp)] = (-1*log10(sub.df$Wilco_p + 0.0001))[is.na(sub.df$logp)]
  
  pdf(paste0("path/to/output5_LR_pair.dot_plot.",name.i,".pdf"),width = 2 + length(levels(sub.df$ct_pair))/2,height = 2 + length(levels(sub.df$LR_pair))/3)
  p1 <- ggplot(sub.df, aes(x = ct_pair, y = LR_pair, color = mergeFC, size = logp)) +
    geom_point(pch = 16, na.rm = TRUE) +
    scale_x_discrete(position = "top") + 
    # scale_color_gradientn(colors = paletteer_c("grDevices::Viridis", 30)) + 
    scale_radius(limits = c(1,4),range = c(3,8)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0,color = "#000000"), axis.text.y = element_text(color = "#000000"),
                       axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(p1)
  p1 <- ggplot(sub.df, aes(x = ct_pair, y = LR_pair, color = mergeFC, size = logp)) +
    geom_point(pch = 16, na.rm = TRUE) +
    scale_x_discrete(position = "top") + 
    scale_color_gradientn(colors = paletteer_c("grDevices::Viridis", 30),limits = c(0,max(sub.df$mergeFC))) + 
    scale_radius(limits = c(1,4),range = c(3,8)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0,color = "#000000"), axis.text.y = element_text(color = "#000000"),
                       axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(p1)
  dev.off()
}

plotDotLRsNew(name.i = "Endo_system",
              must_include = c(11),
              may_include = c(32,24,7),
              focus=c("Itga","Itgb","Vegfa","Vegfc","Pdgfc","Notch"))
plotDotLRsNew(name.i = "Endo_system2",
              must_include = c(11),
              may_include = c(32,24,7),
              focus=c("Itga","Itgb","Vegfa","Vegfc","Pdgfc","Notch"),
              only.focus=T)

# 脚本执行完成
cat("脚本执行完毕。\n")
