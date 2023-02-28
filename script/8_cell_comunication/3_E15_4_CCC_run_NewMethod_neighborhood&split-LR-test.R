# Functions
get_neighbor_cells <- function(meta.df,dist.m,cti,ctj,permut=FALSE){
  if(permut){ # random cells within each cell type
    cellnum = dim(meta.df)[1]
    meta.df2 = meta.df[sample(1:cellnum,cellnum),]
    meta.df = meta.df2[order(meta.df2$Ann_v5),]
  }
  dist.df = data.frame("celli"=rep(meta.df$coord, times=dim(meta.df)[1]),
                       "cellj"=rep(meta.df$coord, each=dim(meta.df)[1]),
                       "cti"=rep(meta.df$Ann_v5, times=dim(meta.df)[1]),
                       "ctj"=rep(meta.df$Ann_v5, each=dim(meta.df)[1]),
                       "dist"=as.vector(dist.m),stringsAsFactors = F)
  dist.df = dist.df[dist.df$dist < 2 & dist.df$cti == cti & dist.df$ctj == ctj,]
  return(list(unique(dist.df$celli),unique(dist.df$cellj)))
}

get_neighborhood <- function(meta.df,dist.m,nei_cell,nbhood=FALSE,cti=NA,ctj=NA,permut=FALSE){
  if(!nbhood){
    nei_cell = c(nei_cell[[1]],nei_cell[[2]])
    dist.df = data.frame("celli"=rep(meta.df$coord, times=dim(meta.df)[1]),
                         "cellj"=rep(meta.df$coord, each=dim(meta.df)[1]),
                         "cti"=rep(meta.df$Ann_v5, times=dim(meta.df)[1]),
                         "ctj"=rep(meta.df$Ann_v5, each=dim(meta.df)[1]),
                         "dist"=as.vector(dist.m),stringsAsFactors = F)
    dist.df = dist.df[dist.df$dist < 2 & dist.df$cellj %in% nei_cell,]
    nei_cell= unique(dist.df$celli)
    return(nei_cell)
  }
  meta.df = meta.df[meta.df$coord %in% nei_cell,]
  if(permut){
    meta.df$Ann_v5 = sample(meta.df$Ann_v5,dim(meta.df)[1])
  }
  return(list(meta.df$coord[meta.df$Ann_v5==cti],meta.df$coord[meta.df$Ann_v5==ctj]))
}

get_LR_expr_score2 <- function(count.m,nei_cell,lr.df){
  LR_score1 = apply(count.m[lr.df$V1,nei_cell[[1]]], 1, mean)
  LR_score2 = apply(count.m[lr.df$V2,nei_cell[[2]]], 1, mean)
  names(LR_score1) = paste0(lr.df$V1,"_",lr.df$V2)
  names(LR_score2) = paste0(lr.df$V1,"_",lr.df$V2)
  return(list(LR_score1,LR_score2))
}

get_LR_expr_01 <- function(count.m,cti,ctj,lr.df){
  LR_01.1 = apply(count.m[lr.df$V1,meta.df$coord[meta.df$Ann_v5==cti]], 1, mean)
  LR_01.2 = apply(count.m[lr.df$V2,meta.df$coord[meta.df$Ann_v5==ctj]], 1, mean)
  LR_01 = pmin(LR_01.1,LR_01.2) # 1*(LR_01.1 > 0.1 & LR_01.2 > 0.1)
  names(LR_01) = paste0(lr.df$V1,"_",lr.df$V2)
  return(LR_01)
}

library(parallel)
load("Plot_E15_4/2_CCC_prepare.RData")

# 0 get distance between celli cellj
meta.df = meta.df[order(meta.df$Ann_v5),] # prepare work: shuffle cell coord while keep cell identity ; see also: get_neighbor_cells(permut=T)
distA.m = meta.df$nA %*% matrix(1,nrow = 1,ncol = dim(meta.df)[1]) - matrix(1,nrow = dim(meta.df)[1],ncol = 1) %*% t(meta.df$nA)
distB.m = meta.df$nB %*% matrix(1,nrow = 1,ncol = dim(meta.df)[1]) - matrix(1,nrow = dim(meta.df)[1],ncol = 1) %*% t(meta.df$nB)
dist.m = sqrt(distA.m**2 + distB.m**2)
dist.m = dist.m + diag(100,nrow = dim(meta.df)[1])
rm(list = c("distA.m","distB.m"))

# 1 Get neighborhood
neighbor_i = c("15","8") # c("8","15")
name_i = "UBT_NP" # "NP_UBT"
# for direct CM-UBT pairs
neighbor_i_cell = get_neighbor_cells(meta.df,dist.m,neighbor_i[1],neighbor_i[2])
meta.df$interact = "Others"
meta.df$interact[meta.df$coord %in% neighbor_i_cell[[1]]] = neighbor_i[1]
meta.df$interact[meta.df$coord %in% neighbor_i_cell[[2]]] = neighbor_i[2]
# for neighborhood of direct CM-UBT pairs
neighborhood_i = get_neighborhood(meta.df,dist.m,neighbor_i_cell)
meta.df$nbhood = "Others"
meta.df$nbhood[meta.df$coord %in% neighborhood_i] = "Local"
meta.df$nbhood[meta.df$coord %in% neighborhood_i & meta.df$Ann_v5 == neighbor_i[1]] = neighbor_i[1]
meta.df$nbhood[meta.df$coord %in% neighborhood_i & meta.df$Ann_v5 == neighbor_i[2]] = neighbor_i[2]
# for cell type label
meta.df$shape = "background"
meta.df$shape[meta.df$Ann_v5==neighbor_i[1]] = neighbor_i[1]
meta.df$shape[meta.df$Ann_v5==neighbor_i[2]] = neighbor_i[2]
library(ggplot2)
library("viridis")
pdf(paste0("Plot_E15_4/3_spatial_map_overview.",name_i,".pdf"),width = 7.6,height = 6.3)
color.ls = c("Others"="grey","Local"="green","L"="red","R"="blue")
names(color.ls)[3:4] = neighbor_i
shape.ls = c("background"=15,"L"=19,"R"=17)
names(shape.ls)[2:3] = neighbor_i
ggplot(meta.df,aes(x=nB,y=nA,color=interact)) + 
  geom_point(aes(shape=shape),size=1.5) + theme_minimal() + 
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
  expand_limits(y=c(0,96),x=c(0,96)) + 
  scale_colour_manual(values = color.ls) +
  scale_shape_manual(values = shape.ls)
ggplot(meta.df,aes(x=nB,y=nA,color=nbhood)) + 
  geom_point(aes(shape=shape),size=1.5) + theme_minimal() + 
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
  expand_limits(y=c(0,96),x=c(0,96)) + 
  scale_colour_manual(values = color.ls) +
  scale_shape_manual(values = shape.ls)
dev.off()

# 2 permutation test
neighbor_i_score = get_LR_expr_score2(correct_count.m,
                                      get_neighborhood(meta.df,dist.m,neighborhood_i,nbhood = T,
                                                       cti = neighbor_i[1],ctj = neighbor_i[2],permut = F),
                                      lr_pair.df)
neighbor_i_score.shuff = Reduce(function(i,j) return(list(cbind(i[[1]],j[[1]]),cbind(i[[2]],j[[2]]))),
                                mclapply(1:5000,function(i){
                                  cell_k = get_neighborhood(meta.df,dist.m,neighborhood_i,nbhood = T,
                                                            cti = neighbor_i[1],ctj = neighbor_i[2],permut = T)
                                  return(get_LR_expr_score2(correct_count.m,cell_k,lr_pair.df))
                                },mc.cores = 4))
neighbor_i_fc = (neighbor_i_score[[1]]/apply(neighbor_i_score.shuff[[1]],1,mean) +
                   neighbor_i_score[[2]]/apply(neighbor_i_score.shuff[[2]],1,mean))/2
neighbor_i_pv = cbind(apply((neighbor_i_score.shuff[[1]]-neighbor_i_score[[1]])>0,1,sum)/5000,
                      apply((neighbor_i_score.shuff[[2]]-neighbor_i_score[[2]])>0,1,sum)/5000)
neighbor_i_pv = apply(neighbor_i_pv, 1, mean) # max
neighbor_i_01 = get_LR_expr_01(correct_count.m,neighbor_i[1],neighbor_i[2],lr_pair.df)

result.df = data.frame("lr"=names(neighbor_i_01),
                       "fc"=neighbor_i_fc,
                       "pv"=neighbor_i_pv,
                       "exp"=neighbor_i_01,stringsAsFactors = F)
result.df[result.df$lr %in% c("Gdnf_Gfra1","Gdnf_Ret","Fgf10_Fgfr2"),]
result.df = result.df[result.df$pv<0.05 & result.df$exp>0.15,]
write.csv(result.df,file = paste0("Plot_E15_4/4_significant_LR_pairs.",name_i,".csv"),row.names=F)

# 3 2D plot
pdf(paste0("Plot_E15_4/5_spatial_map.",name_i,".pdf"),width = 7.6,height = 6.3)
for (i in 1:dim(result.df)[1]) {
  meta.df$Li = correct_count.m[strsplit(result.df$lr[i],split = "_")[[1]][1],meta.df$coord]
  meta.df$Ri = correct_count.m[strsplit(result.df$lr[i],split = "_")[[1]][2],meta.df$coord]
  meta.df$Li = pmin(meta.df$Li,sort(meta.df$Li,decreasing = T)[20])
  meta.df$Ri = pmin(meta.df$Ri,sort(meta.df$Ri,decreasing = T)[20])
  p1 = ggplot(meta.df,aes(x=nB,y=nA,color=Li)) + 
    geom_point(aes(shape=shape),size=1.5) + theme_minimal() + 
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",strsplit(result.df$lr[i],split = "_")[[1]][1])) + 
    expand_limits(y=c(0,96),x=c(0,96)) + 
    scale_color_viridis(option = "B",begin = 0.1,end = 0.9) +
    # scale_colour_gradient(low = "grey",high = "blue4",limits = c(0,sort(meta.df$Li,decreasing = T)[20])) +
    scale_shape_manual(values = shape.ls)
  print(p1)
  p1 = ggplot(meta.df,aes(x=nB,y=nA,color=Ri)) + 
    geom_point(aes(shape=shape),size=1.5) + theme_minimal() + 
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) + 
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",strsplit(result.df$lr[i],split = "_")[[1]][2])) + 
    expand_limits(y=c(0,96),x=c(0,96)) + 
    scale_color_viridis(option = "B",begin = 0.1,end = 0.9) +
    # scale_colour_gradient(low = "grey",high = "blue4",limits = c(0,sort(meta.df$Ri,decreasing = T)[20])) +
    scale_shape_manual(values = shape.ls)
  print(p1)
}
dev.off()


