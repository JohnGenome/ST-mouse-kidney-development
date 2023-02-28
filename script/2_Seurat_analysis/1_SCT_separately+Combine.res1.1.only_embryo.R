library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(Matrix)

LoadPixelInTissue <- function(i){
  print(Sample.ls[i])
  count.m = read.table(data_path.ls[i])
  count.m = as.matrix(count.m)
  tmp = CreateSeuratObject(counts = t(count.m), project = Sample.ls[i], assay = "RNA")
  tmp$channel = channel.ls[i]
  tmp$stage = Stage.ls[i]
  tmp$coord = paste0(Sample.ls[i],"_",row.names(tmp@meta.data))
  tmp$nA = as.integer(str_split_fixed(tmp$coord,"_|x",4)[,3])
  tmp$nB = as.integer(str_split_fixed(tmp$coord,"_|x",4)[,4])
  
  tmp.subset = tmp[,tmp$nCount_RNA > Filter.m[i]]
  tmp.subset <- SCTransform(tmp.subset, assay = "RNA", verbose = FALSE,variable.features.n = 3000)
  return(tmp.subset)
}
# 1 Merge data
data_path.ls = c("...")
channel.ls = c("20um","20um",
               "20um","20um","20um","20um",
               "20um","20um",
               "20um","20um","20um","10um",
               "10um",
               "20um","20um","10um","10um")
Stage.ls = c("E11","E11",
             "E12","E12","E12","E12",
             "E13","E13",
             "E15","E15","E15","E15",
             "E16",
             "E18","E18","E18","E18")
Sample.ls = c("E11_1","E11_2",
              "E12_1","E12_2","E12_3","E12_4",
              "E13_1","E13_2",
              "E15_1","E15_2","E15_3","E15_4",
              "E16_1",
              "E18_1","E18_2","E18_3","E18_4")
Filter.m = c(1000,1000,
             2000,1000,2000,2000,
             1500,2000,
             1000,800,800,500,
             1000,
             1000,1000,800,800)

# Run 1 by 1
E11.1 = LoadPixelInTissue(1)
E11.2 = LoadPixelInTissue(2)
E12.1 = LoadPixelInTissue(3)
E12.2 = LoadPixelInTissue(4)
E12.3 = LoadPixelInTissue(5)
E12.4 = LoadPixelInTissue(6)
E13.1 = LoadPixelInTissue(7)
E13.2 = LoadPixelInTissue(8)
E15.1 = LoadPixelInTissue(9)
E15.2 = LoadPixelInTissue(10)
E15.3 = LoadPixelInTissue(11)
E15.4 = LoadPixelInTissue(12)
E16.1 = LoadPixelInTissue(13)
E18.1 = LoadPixelInTissue(14)
E18.2 = LoadPixelInTissue(15)
E18.3 = LoadPixelInTissue(16)
E18.4 = LoadPixelInTissue(17)

mKid.merge <- merge(E11.1, y = c(E11.2,E12.1,E12.2,E12.3,E12.4,E13.1,E13.2,E15.1,E15.2,E15.3,E15.4,E16.1,E18.1,E18.2,E18.3,E18.4), 
                    add.cell.ids = Sample.ls[1:17], project = "mKidney", merge.data = T)
rm(list = c("E11.1","E11.2","E12.1","E12.2","E12.3","E12.4","E13.1","E13.2","E15.1","E15.2","E15.3","E15.4",
            "E16.1","E18.1","E18.2","E18.3","E18.4"))

# 2 Clustering
# Seurat Issue 2814/5183
set.seed(2022)
VariableFeatures(mKid.merge[["SCT"]]) <- rownames(mKid.merge[["SCT"]]@scale.data) # 6480. In the merged object, the genes in the scale.data are the intersected genes
mKid.merge <- RunPCA(mKid.merge, assay = "SCT", verbose = T)
k1 = 30;k2 = 1.1
mKid.merge <- FindNeighbors(mKid.merge, reduction = "pca", dims = 1:k1)
mKid.merge <- FindClusters(mKid.merge, verbose = T, resolution = k2)
mKid.merge <- RunUMAP(mKid.merge, reduction = "pca", dims = 1:k1)

library(scales)
color.ls2 = c("#7CAE00","#ffdf00","#39B600","#7F7F7F","#C09B00","#ff4500","#00BAE0","#A3A500","#000080","#90EE90",
              "#EA8331","#00B0F6","#FF0000","#704241","#FFFF00","#C77CFF","#daa520","#adff2f","#C77CFF","#E76BF3",
              "#FF62BC","#9590FF","#FF6A98","#FA62DB")
names(color.ls2) = levels(mKid.merge$seurat_clusters)

pdf(paste0("Plot/PCA",k1,"_res",k2,"_1_cluster_UMAP2.pdf"))
DimPlot(mKid.merge, reduction = "umap", label = T,shuffle = T) + scale_color_manual(values = color.ls2)
dev.off()
pdf("Plot/1_cluster_UMAP_split.pdf",height = 4,width = 50)
DimPlot(mKid.merge, reduction = "umap",group.by = "Ann_v3_Fib", split.by = "orig.ident",label.size = 4) + scale_color_manual(values = color.ls)
dev.off()
pdf("Plot/1_cluster_UMAP_split_stage.pdf",height = 4,width = 20)
DimPlot(mKid.merge, reduction = "umap",group.by = "Ann_v3_Fib", split.by = "stage",label.size = 4) + scale_color_manual(values = color.ls)
dev.off()
pdf("Plot/1_cluster_UMAP_split_stage2.pdf",height = 4,width = 20)
DimPlot(mKid.merge, reduction = "umap",label = T, split.by = "stage",label.size = 2) + scale_color_manual(values = color.ls2)
dev.off()

# 3 Known markers
mKid.merge$mean_depth_sample = tapply(mKid.merge$nCount_SCT,mKid.merge$orig.ident,mean)[mKid.merge$orig.ident]
mKid.merge@assays[["SCT"]]@data = Matrix(t(t(as.matrix(mKid.merge@assays[["SCT"]]@counts))/mKid.merge$mean_depth_sample*10000),sparse=T)

pdf("Plot/2_marker_UMAP.CM.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Six2", "Cited1","Pax2","Crym"),slot = "data",order = T) # CM 2018 Scie + 2019 Dev
FeaturePlot(mKid.merge, features = c("Itga8"),slot = "data",order = T, pt.size = .8) # MM & CM 2018 Dev
VlnPlot(mKid.merge, features = c("Six2", "Cited1","Pax2","Crym","Itga8"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.G.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Nphs2","Nphs1","Podxl","Mafb"),slot = "data",order = T) # 2018 Scie + 2019 Dev
FeaturePlot(mKid.merge, features = c("Ptpro","Wt1","Podxl"),slot = "data",order = T) # 2018 Scie
VlnPlot(mKid.merge, features = c("Nphs2","Nphs1","Podxl","Mafb"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.PT.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Slc27a2","Lrp2","Slc5a2","Slc5a12"),slot = "data",order = T) # 2018 Scie PT*2 + S1*2
FeaturePlot(mKid.merge, features = c("Fxyd2","Atp11a","Slc13a3","Slc22a30"),slot = "data",order = T) # 2018 Scie S2*2 + S3*1 + 2021 NC PST
FeaturePlot(mKid.merge, features = c("Ace","Fbp1","Spp2","Slc3a1"),slot = "data",order = T) # 2019 Dev + 2018 Dev
FeaturePlot(mKid.merge, features = c("Slc34a1","Slc22a8","Slc17a3","Slc22a7"),slot = "data",order = T) # 2021 NC + 2018 Scie
FeaturePlot(mKid.merge, features = c("Slc16a9","Slc7a13","Slc13a3"),slot = "data",order = T) # 2018 Scie
VlnPlot(mKid.merge, features = c("Slc27a2","Lrp2","Slc5a2","Slc5a12","Fxyd2","Atp11a","Slc13a3","Slc22a30",
                                 "Gpx3","Slc3a1"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.LoH.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Slc12a1","Umod","Cldn16"),slot = "data",order = T) # 2018 Scie + 2018 Scie
FeaturePlot(mKid.merge, features = c("Fst","Aqp1","Slc14a2","Bst1"),slot = "data",order = T) # 2022-review descending LoH
VlnPlot(mKid.merge, features = c("Slc12a1","Umod","Cldn16","Fst","Aqp1","Slc14a2","Bst1"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.DT.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Slc12a3","Pvalb","Tmem52b","Slc8a1"),slot = "data",order = T) # 2018 Scie + 2018 Dev + 2018 Scie
VlnPlot(mKid.merge, features = c("Slc12a3","Pvalb","Tmem52b","Slc8a1"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.UB.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Ret","Slco4c1","Lhx1","Hnf1b"),slot = "data",order = T) # 2018 Dev UBT*2(2nd our defined) + 2018 Scie UB*2
VlnPlot(mKid.merge, features = c("Ret","Slco4c1","Lhx1","Hnf1b"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.CD.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Aqp2","Hsd11b2","Atp6v1g3","Atp6v0d2"),slot = "data",order = T) # 2018 Scie PC*2+IC*2
FeaturePlot(mKid.merge, features = c("Slc26a7","Slc26a4","Gata3"),slot = "data",order = T) # 2021 NC IC*2 + 2018 Dev CD
VlnPlot(mKid.merge, features = c("Aqp2","Hsd11b2","Atp6v1g3","Atp6v0d2"),slot = "data",stack = T)
VlnPlot(mKid.merge, features = c("Slc26a7","Slc26a4","Gata3"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.Pelvic.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Upk1b","Upk1a","Myh11","Myocd"),slot = "data",order = T) # 2018 Scie + 2 our defined
FeaturePlot(mKid.merge, features = c("Hba.a1","Hba.a2","Mrc1","F13a1"),slot = "data",order = T) # RBC + 2 our defined Macro
FeaturePlot(mKid.merge, features = c("Hba.a1","Hba.a2","Mrc1","F13a1"),slot = "data",order = T,max.cutoff = 10)
FeaturePlot(mKid.merge, features = c("C1qb","S100a8","Col1a1","Fabp4"),slot = "data",order = T) # Macro+Neutrophils+Strom+Endo
VlnPlot(mKid.merge, features = c("Hba.a1","Hba.a2"),slot = "data",y.max = 10, ncol = 1)
VlnPlot(mKid.merge, features = c("Mrc1","F13a1"),slot = "data",stack = T)
dev.off()
pdf("Plot/2_marker_UMAP.Endo.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Pecam1","Cdh5","Flt1","Egfl7"),slot = "data",order = T) # 2019 Dev + 2021 NC Endo
VlnPlot(mKid.merge, features = c("Pecam1","Cdh5","Flt1","Egfl7"),slot = "data",stack = T)
FeaturePlot(mKid.merge, features = c("Mgp","Fbln5","Kdr"),slot = "data",order = T)
VlnPlot(mKid.merge, features =  c("Mgp","Fbln5","Kdr"),slot = "data",stack = T,flip = T,idents = c(11,19))
VlnPlot(mKid.merge, features = c("Mgp","Fbln5","Edn1","Eln","Gja5", # large artery
                                 "Thbd","Calca","Aqp1","Sox17","Glul", # cap artery
                                 "Rgcc","Esm1","Igfbp3","Igfbp5","Kdr", # capillary
                                 "Igf1","Bgn","Nr2f2", # vein
                                 "Aplnr","Trp53l11", # angiogenic
                                 "Plat","Pi16","Lpl" # glomeruli
),slot = "data",stack = T,flip = T)
dev.off()
pdf("Plot/3_Vln_UMI_data.pdf",width = 10,height = 5)
VlnPlot(mKid.merge, features = "nCount_RNA")
VlnPlot(mKid.merge, features = "nCount_RNA",log = T)
VlnPlot(mKid.merge, features = "nCount_RNA",y.max = 10000)
dev.off()

# 3 Known markers 2022 Review
pdf("Plot/3_marker_UMAP.Pod.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Nphs1","Nphs2","Synpo","Cdkn1c"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Wt1","Foxc2","Mafb","Efnb2"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Nphs1","Nphs2","Synpo","Cdkn1c","Wt1","Foxc2","Mafb","Efnb2","Foxl1"),slot = "data",stack = T)
dev.off()
pdf("Plot/3_marker_UMAP.PT.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Slc34a1","Lrp2","Hxyd2","Hrsp12"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Acsm1","Acsm2","Cpt1a","Acox3"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc26a6","Slc9a3","Glud1","Pck1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Aqp8","Hnf4a","Ppara","Slc5a2"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc5a12","Adra1a","Slc6a19","Slc7a8"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc7a9","Atp11a","Slc13a3","Slc16a9"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc27a2","Slc7a13","Slc22a6","Slc1a1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Notch2","Lgr4","Havcr1","Krt20"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Hspa1a","Vcam1","Dcdc2a","Sema5a"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Slc34a1","Lrp2","Hxyd2","Hrsp12","Acsm1","Acsm2","Cpt1a","Acox3",
                                 "Slc26a6","Slc9a3","Glud1","Pck1","Aqp8","Hnf4a","Ppara","Slc5a2",
                                 "Slc5a12","Adra1a","Slc6a19","Slc7a8","Slc7a9","Atp11a","Slc13a3","Slc16a9",
                                 "Slc27a2","Slc7a13","Slc22a6","Slc1a1","Notch2","Lgr4","Havcr1","Krt20",
                                 "Hspa1a","Vcam1","Dcdc2a","Sema5a"),slot = "data",stack = T,flip = T) + NoLegend()
dev.off()
pdf("Plot/3_marker_UMAP.LoH.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Fst","Aqp1","Slc14a2","Bst1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Epha7","Cryab","Tshz2","Cald1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Lypd2","Mx2","Clcnka","Slc12a1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Umod","Tmem207","Foxq1","Cldn10"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Ptger3","Kcnj1","Enox1","Thsd4"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Mt2","Slc5a3","Nos1","Avpr1a"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Fst","Aqp1","Slc14a2","Bst1","Epha7","Cryab","Tshz2","Cald1","Lypd2","Mx2","Clcnka","Slc12a1","Umod",
                                 "Tmem207","Foxq1","Cldn10","Ptger3","Kcnj1","Enox1","Thsd4","Mt2","Slc5a3","Nos1","Avpr1a"),slot = "data",stack = T,flip = T) + NoLegend()
dev.off()
pdf("Plot/3_marker_UMAP.DT.data.pdf",width = 10,height = 10) # Almost the same markers between DCT1/2 CNT
FeaturePlot(mKid.merge, features = c("Pvalb","Slc12a3","Trpm7","Wnk1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Wnk4","Stk39","Calb1","Slc8a1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Egf","Trpm6","Cnnm2","Atp1a1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Atp1a2","Atp1a3","Atp1a4","Fxyd2"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc12a3","Trpm7","Wnk1","Wnk4"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Klhl3","Stk39","Calb1","Slc8a1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Egf","Trpm6","Cnnm2","Atp1a1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Atp1a2","Atp1a3","Atp1a4","Klk1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Trpv5","Trpm6","S100g","Atp2b1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Atp2b4","Scnn1b","Scnn1g","Kcne1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Fxyd2","Calb1","Slc8a1","Egf"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Klk1","Trpv5","Trpm6","S100g"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Atp2b1","Scnn1b","Scnn1g","Kcne1"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Pvalb","Slc12a3","Trpm7","Wnk1","Wnk4","Stk39","Calb1","Slc8a1","Egf","Trpm6","Cnnm2","Atp1a1",
                                 "Atp1a2","Atp1a3","Atp1a4","Fxyd2","Slc12a3","Trpm7","Wnk1","Wnk4","Klhl3","Stk39","Calb1","Slc8a1",
                                 "Egf","Trpm6","Cnnm2","Atp1a1","Atp1a2","Atp1a3","Atp1a4","Klk1","Trpv5","Trpm6","S100g","Atp2b1",
                                 "Atp2b4","Scnn1b","Scnn1g","Kcne1","Fxyd2","Calb1","Slc8a1","Egf","Klk1","Trpv5","Trpm6","S100g",
                                 "Atp2b1","Scnn1b","Scnn1g","Kcne1"),slot = "data",stack = T,flip = T) + NoLegend()
dev.off()
pdf("Plot/3_marker_UMAP.CD.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Scnn1b","Scnn1g","Aqp2","Avpr2"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Hsd11b2","Rhbg","Elf5","Fxyd4"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Aqp3","Apela","Kcne1","Npnt"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Kcnj10","Tcfcp2l1","Foxi","Atp6v1g3"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Atp6v0d2","Insr","Atp6v1b1","Atp4a"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc4a1","Aqp6","Kit","Adgrf5"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Mme","Slc26a4","Hmx2","Spink8"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Hsd11b2","Rhbg","Atp6v1g3","Atp6v0d2"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Insr","Atp6v1b1","Parm1","Sec23b"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Scnn1b","Scnn1g","Aqp2","Avpr2","Hsd11b2","Rhbg","Elf5","Fxyd4","Aqp3","Apela","Kcne1","Npnt",
                                 "Kcnj10","Tcfcp2l1","Foxi","Atp6v1g3","Atp6v0d2","Insr","Atp6v1b1","Atp4a","Slc4a1","Aqp6","Kit","Adgrf5",
                                 "Mme","Slc26a4","Hmx2","Spink8","Hsd11b2","Rhbg","Atp6v1g3","Atp6v0d2","Insr","Atp6v1b1","Parm1","Sec23b"),slot = "data",stack = T,flip = T) + NoLegend()
dev.off()
pdf("Plot/3_marker_UMAP.Immune.data.pdf",width = 12,height = 16)
FeaturePlot(mKid.merge, features = c("C1qa","C1qb","Itgam","Apoe","C1qc","Cd74","Ctss","Fcer1g","Aif1","Ms4a7","S100a8","S100a9"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Lyz2","Plac8","Ifitm3","Cebpb","Tyrobp","Lst1","Hp","Ifitm1","Hdc","Mcmpt8","Fcer1a","Csrp3"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Ms4a2","Cyp11a1","Cd200r3","Il6","Il4","Cd209a","Wfdc17","Mgl2","Ccl6","Ccl9","Alox5ap","Irf8"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Naaa","Plbd1","Cbfa2t3","Basp1","Rnase6","Wdfy3","Sept3","Ppm1m","Rab7b","Ly6d","Siglech","Cox6a2"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Sell","Ccr9","Runx2","Cd209d","Bcl11a","Lair1","Cd79a","Cd79b","Ms4a1","Ebf1","Cd22","Cd19"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Fcmr","Siglecg","Fcrl1","Cxcr6","Cd247","Nkg7","Lef1","Ms4a4b","Il7r","Ccr7","Klf2","Tcf7"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Dapl1","Satb1","Cd3d","Ccl5","Cd8b1","Cd8a","Hcst","Cd3g","Lck","Tnfrsf4","Capg","Ikzf2"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Izumo1r","Ifi27l2a","S100a4","Rgs1","Ltb","Tnfrsf18","Ly6c2","Gimap3","Tmsb10","Gimap4","Ctsw","Gzma"),slot = "data",order = T,ncol = 3, pt.size = .2)
FeaturePlot(mKid.merge, features = c("Cd7","Xcl1","Klrd1","Klrk1","Ncr1","Klre1","Il2rb"),slot = "data",order = T,ncol = 3, pt.size = .2)
VlnPlot(mKid.merge, features = c("C1qa","C1qb","Itgam","Apoe","C1qc","Cd74","Ctss","Fcer1g","Aif1","Ms4a7","S100a8","S100a9",
                                 "Lyz2","Plac8","Ifitm3","Cebpb","Tyrobp","Lst1","Hp","Ifitm1","Hdc","Mcmpt8","Fcer1a","Csrp3",
                                 "Ms4a2","Cyp11a1","Cd200r3","Il6","Il4","Cd209a","Wfdc17","Mgl2","Ccl6","Ccl9","Alox5ap","Irf8",
                                 "Naaa","Plbd1","Cbfa2t3","Basp1","Rnase6","Wdfy3","Sept3","Ppm1m","Rab7b","Ly6d","Siglech","Cox6a2",
                                 "Sell","Ccr9","Runx2","Cd209d","Bcl11a","Lair1","Cd79a","Cd79b","Ms4a1","Ebf1","Cd22","Cd19",
                                 "Fcmr","Siglecg","Fcrl1","Cxcr6","Cd247","Nkg7","Lef1","Ms4a4b","Il7r","Ccr7","Klf2","Tcf7",
                                 "Dapl1","Satb1","Cd3d","Ccl5","Cd8b1","Cd8a","Hcst","Cd3g","Lck","Tnfrsf4","Capg","Ikzf2",
                                 "Izumo1r","Ifi27l2a","S100a4","Rgs1","Ltb","Tnfrsf18","Ly6c2","Gimap3","Tmsb10","Gimap4","Ctsw","Gzma",
                                 "Cd7","Xcl1","Klrd1","Klrk1","Ncr1","Klre1","Il2rb"),slot= "data",stack = T,flip = T) + NoLegend()
dev.off()
pdf("Plot/3_marker_UMAP.SMC.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Serpine2","Fhl2","Des","Prkca"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Art3","Nt5e","Pdgfrb","Vim"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Tagln","Myh11","Pdgfrb","Tagln"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Myh11","Acta2","Gata3","Rergl"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Map3k7cl","Ren1","Akr1b7","Rgs5"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Serpine2","Fhl2","Des","Prkca","Art3","Nt5e","Pdgfrb","Vim","Tagln","Myh11","Pdgfrb","Tagln",
                                 "Myh11","Acta2","Gata3","Rergl","Map3k7cl","Ren1","Akr1b7","Rgs5"),slot = "data",stack = T,flip = T) + NoLegend()
dev.off()
pdf("Plot/3_marker_UMAP.Endo.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Nrp1","Cdh5","Eln","Plat"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Emcn","Tsapn7","Mapt","Kdr"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Smad6","Ehd3","Lpl","Flt1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Fbln2","Mgp","Trpv4","Bmx"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Kdr","Smad6","Ehd3","Lpl"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Flt1","Fbln2","Mgp","Trpv4"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Bmx","Sox17","Cxcl12","Gja5"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Edn1","Fbln5","Cldn5","Efnb2"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Klf4","Cryab","Gas6","Podxl"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Plvap","Bgn","Cd9","Nr2f2"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Fxyd2","Fxyd6","Igfbp7","Slc14a1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Aqp1","S100a4","S100a4"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Nrp1","Cdh5","Eln","Plat","Emcn","Tsapn7","Mapt","Kdr","Smad6","Ehd3","Lpl","Flt1",
                                 "Fbln2","Mgp","Trpv4","Bmx","Kdr","Smad6","Ehd3","Lpl","Flt1","Fbln2","Mgp","Trpv4",
                                 "Bmx","Sox17","Cxcl12","Gja5","Edn1","Fbln5","Cldn5","Efnb2","Klf4","Cryab","Gas6","Podxl",
                                 "Plvap","Bgn","Cd9","Nr2f2","Fxyd2","Fxyd6","Igfbp7","Slc14a1","Aqp1","S100a4"),slot = "data",stack = T,flip = T) + NoLegend()
dev.off()

# 4 Known markers 2022 Review
pdf("Plot/4_marker_UMAP.Develop.data.pdf",width = 10,height = 10)
FeaturePlot(mKid.merge, features = c("Six2","Cited1","Wnt4","Lhx1"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Pax8","Jag1","Wt1","Olfm3"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Mafb","Podxl","Cldn1","Hnf4a"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc3a1","Cdh1","Fgf8","Sox9"),slot = "data",order = T)
FeaturePlot(mKid.merge, features = c("Slc12a1","Irx1","Tfap2a","Gata3"),slot = "data",order = T)
VlnPlot(mKid.merge, features = c("Six2","Cited1","Wnt4","Lhx1","Pax8","Jag1","Wt1","Olfm3","Mafb",
                                 "Podxl","Cldn1","Hnf4a","Slc3a1","Cdh1","Fgf8","Sox9","Slc12a1", # "Pou3f3"
                                 "Irx1","Tfap2a","Gata3"),slot = "data",stack = T,flip = T) + NoLegend()
dev.off()

# 4 DE genes
mKid.merge <- NormalizeData(mKid.merge,assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)
de_markers <- FindAllMarkers(mKid.merge, assay = "RNA", only.pos = T,min.pct = 0.1,logfc.threshold = 0.5,min.diff.pct = 0.1)
de_markers %>% group_by(cluster) %>% slice_max(n=3,order_by = avg_log2FC) -> top3
de_markers = de_markers[order(de_markers$avg_log2FC,decreasing = T),]
de_markers = de_markers[order(de_markers$cluster),]
write.csv(de_markers,row.names = F,file = "Out/4_de_genes_FCorder.csv")
de_markers = read.csv("Out/4_de_genes.csv",stringsAsFactors=F)

