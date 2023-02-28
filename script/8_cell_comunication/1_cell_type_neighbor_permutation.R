library(parallel)
library(ComplexHeatmap)
library(circlize)
library(igraph)

get_neighbor_count <- function(meta.df,dist.m,clus.ls,permut=FALSE){
  if(permut) meta.df$Ann_v5 = sample(meta.df$Ann_v5,dim(meta.df)[1])
  dist.df = data.frame("cti"=factor(rep(meta.df$Ann_v5, times=dim(meta.df)[1]),levels = clus.ls), # cautious 
                       "ctj"=factor(rep(meta.df$Ann_v5, each=dim(meta.df)[1]),levels = clus.ls), # !!!
                       "dist"=as.vector(dist.m))
  dist.df = dist.df[dist.df$dist < 2,]
  neighbor.df = as.data.frame(table(dist.df$cti,dist.df$ctj))
  colnames(neighbor.df) = c("cti","ctj","count")
  return(neighbor.df)
}

for (sample_i in c("E18_4")) { # c("E11_1","E11_2","E12_1","E12_2","E12_3","E12_4","E13_1","E13_2","E15_1","E15_2","E15_3","E15_4","E16_1","E18_1","E18_2","E18_3","E18_4")
  # 0 Filter data
  cut_off = ifelse(sample_i %in% c("E11_1","E11_2","E12_1","E12_2","E12_3","E12_4","E13_1","E13_2"),5,10)
  
  meta.df = read.csv("../74_Combine_E11_E12_E13_E15_E16_E18_Adult/Out/6_metadata.Annv5.csv",stringsAsFactors = F,row.names = 1)
  meta.df = meta.df[!(meta.df$Ann_v5 %in% c(3,16)),]
  meta.df = meta.df[meta.df$orig.ident == sample_i,]
  
  cluster.ls = table(meta.df$Ann_v5)
  cluster.ls = names(cluster.ls)[cluster.ls > cut_off]
  meta.df = meta.df[meta.df$Ann_v5 %in% cluster.ls,]
  
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
  system.time(random.m <- mclapply(1:times, function(i) get_neighbor_count(meta.df,dist.m,cluster.ls,permut=T)$count,mc.cores=5)) # cost much times!!!
  # user   system  elapsed
  # 2739.466  858.215  745.945
  random.m = Reduce(cbind,random.m)
  random.m = as.matrix(random.m)
  
  neighbor.df$expect = apply(random.m,1,mean)
  neighbor.df$pvalue = apply((random.m-neighbor.df$count)>0,1,sum)/times
  neighbor.df$oe = neighbor.df$count/neighbor.df$expect
  
  levels.ls = c("0_oCD","1_iMFib","2_aLoH","3_Ukn","4_UEndo","5_oMFib","6_PTprog","7_cNP","8_NP","9_mMfib",
                "10_PT","11_Capil","12_Pod","13_oCFib","14_DT","15_UBT","16_CUkn","17_UFib","18_TC","19_LA",
                "20_dLoH","21_USMC","22_Uprog","23_RBC","24_Dprog","25_Podprog","26_Mesa","27_Corin","28_preDT","29_MD",
                "30_iCD","31_Utran","32_iCFib")
  names(levels.ls) = 0:32
  levels(neighbor.df$cti) = levels.ls[cluster.ls]
  levels(neighbor.df$ctj) = levels.ls[cluster.ls]
  
  rm(list = c("dist.m","meta.df","random.m","times"))
  base::save.image(file = paste0("Out/",sample_i,"_0_cell_type_neighbor_permutation.RData"))
  write.csv(neighbor.df,file = paste0("Out/",sample_i,"_1_save_neighboring.csv"),row.names = F)
  
  ## 3a [optional] Confine to medullary cell types
  # cluster.ls = intersect(cluster.ls,c("5","9","1","2","27","20","0","30"))
  # neighbor.df = neighbor.df[(neighbor.df$cti %in% levels.ls[cluster.ls]) & (neighbor.df$ctj %in% levels.ls[cluster.ls]),]
  # neighbor.df$cti = as.character(neighbor.df$cti)
  # neighbor.df$ctj = as.character(neighbor.df$ctj)
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
    pdf(paste0("Plot/",sample_i,"_0_oe_dist_between_cell_type.",pcut_off,".pdf"),width = 2+length(cluster.ls)/2,height = 2+length(cluster.ls)/2)
    print(Heatmap(oe.m,cluster_rows=T,cluster_columns=T,rect_gp = gpar(col = "white", lwd = 1),
                  col = colorRamp2(c(0, 2, 5), c("blue", "white", "red")),
                  cell_fun = function(j, i, x, y, width, height, fill){
                    if(!is.na(oe.m[i, j])) if(oe.m[i, j] > 0) grid.text(sprintf("%.2f", oe.m[i, j]), x, y, gp = gpar(fontsize = 5))
                  }))
    dev.off()

    # 4 Network
    neighbor.df2 = neighbor.df[!duplicated(neighbor.df$expect) | !duplicated(neighbor.df$oe),]
    neighbor.df2 = neighbor.df2[neighbor.df2$cti != neighbor.df2$ctj,]
    neighbor.df2 = neighbor.df2[neighbor.df2$pvalue < pcut_off,]
    oe.graph = graph.data.frame(neighbor.df2,directed = F)
    pdf(paste0("Plot/",sample_i,"_1_neighboring_cell_type_network.",pcut_off,".pdf"))
    print(plot(oe.graph, edge.width = -log10(neighbor.df2$pvalue+0.00001),
               edge.color = colorRamp2(c(0, 4), c("grey","red"))(neighbor.df2$oe),
               vertex.size = 5,vertex.label.cex=.7,vertex.label.dist=1))
    print(plot(oe.graph, edge.width = -log10(neighbor.df2$pvalue+0.00001),
               edge.color = colorRamp2(c(0, 4), c("grey","red"))(neighbor.df2$oe),
               vertex.size = 5,vertex.label.cex=.7,vertex.label.dist=1))
    dev.off()
  }
}

