## Cell type 2D map ##
pdf("1_spatial_map.cluster.void.pdf")
mKid.merge@meta.data$focus = mKid.merge$Ann_v5
mKid.merge$focus[!(mKid.merge$focus %in% c(8,7,24,6,25,12,10,2,27,20,28,14,29))] = NA # optional
for (Samplei in Sample.ls) {
  p1 = ggplot(mKid.merge@meta.data[mKid.merge$orig.ident==Samplei,],aes(x=nB,y=nA,color=focus)) +
    geom_point(shape=15,size=1.9,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,96),x=c(0,96)) +
    scale_color_manual(values = color.ls5)
  print(p1)
}
dev.off()

## Gene*1 2D map ##
# v1 origin
pdf("2_spatial_map.gene_i.multi_stages.pdf")
mKid.merge$focus = mKid.merge@assays$SCT@data[gene_i,]
maxi = sort(mKid.merge$focus[mKid.merge$orig.ident=="E18_3"],decreasing = T)[30]
mini = 0
mKid.merge$focus = pmin(mKid.merge$focus,maxi)
mKid.merge$focus = pmax(mKid.merge$focus,mini)
p1 = ggplot(mKid.merge@meta.data[mKid.merge$orig.ident=="E18_3",],aes(x=nB,y=nA,color=focus)) +
  geom_point(shape=15,size=1.9,show.legend=T) + theme_minimal() +
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB","E18_3")) +
  expand_limits(y=c(0,96),x=c(0,96)) +
  scale_colour_gradient(low = "#F0F0F0",high = "blue4",limits = c(mini,maxi))
print(p1)
for (Samplei in Sample.ls) {
  p1 = ggplot(mKid.merge@meta.data[mKid.merge$orig.ident==Samplei,],aes(x=nB,y=nA,color=focus)) +
    geom_point(shape=15,size=1.9,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,96),x=c(0,96)) +
    scale_colour_gradient(low = "#F0F0F0",high = "blue4",limits = c(mini,maxi))
  print(p1)
}
dev.off()
# v2 layer boundary
pdf("2_spatial_map.gene_i.multi_stages.pdf")
mKid.merge$focus = mKid.merge@assays$SCT@data[gene_i,]
maxi = sort(mKid.merge$focus[mKid.merge$orig.ident=="E18_3"],decreasing = T)[30]
mini = 0
mKid.merge$focus = pmin(mKid.merge$focus,maxi)
mKid.merge$focus = pmax(mKid.merge$focus,mini)
p1 = ggplot(mKid.merge@meta.data[mKid.merge$orig.ident=="E18_3",],aes(x=nB,y=nA)) +
  geom_point(aes(color=focus),shape=15,size=1.9,show.legend=T) + theme_minimal() +
  geom_point(data=mKid.merge@meta.data[mKid.merge$orig.ident=="E18_3" & mKid.merge$depth>=0 & (mKid.merge$depth %% 1)==0,],shape=0,size=1.9,fill="#000000",show.legend=F) +
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB","E18_3")) +
  expand_limits(y=c(0,96),x=c(0,96)) +
  scale_colour_gradient(low = "#F0F0F0",high = "blue4",limits = c(mini,maxi))
print(p1)
for (Samplei in Sample.ls){
  p1 = ggplot(mKid.merge@meta.data[mKid.merge$orig.ident==Samplei,],aes(x=nB,y=nA)) +
    geom_point(aes(color=focus),shape=15,size=1.9,show.legend=F) + theme_void() +
    geom_point(data=mKid.merge@meta.data[mKid.merge$orig.ident==Samplei & mKid.merge$depth>=0 & (mKid.merge$depth %% 1)==0,],shape=0,size=1.9,fill="#000000",show.legend=F) +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,96),x=c(0,96)) +
    scale_colour_gradient(low = "#F0F0F0",high = "blue4",limits = c(mini,maxi))
  print(p1)
}
dev.off()

## Gene*2 2D map ##
mKid.merge$Grip1 = mKid.merge@assays[["SCT"]]@data["Grip1",]
mKid.merge$Grip1[!(mKid.merge$Ann_v5 %in% c(0,30))] = 0
mKid.merge$Slc14a2 = mKid.merge@assays[["SCT"]]@data["Slc14a2",]
mKid.merge$Slc14a2[!(mKid.merge$Ann_v5 %in% c(0,30))] = 0
mKid.merge$Merge2 = mKid.merge$Grip1
k=3/5
mKid.merge$Merge2[mKid.merge$Slc14a2 > k*mKid.merge$Grip1] = 
  -mKid.merge$Slc14a2[mKid.merge$Slc14a2 > k*mKid.merge$Grip1]
mKid.merge$Merge2[mKid.merge$Merge2 > 50] = 50
mKid.merge$Merge2[mKid.merge$Merge2 < -30] = -30
mKid.merge$Merge2[mKid.merge$Grip1 > 25 & mKid.merge$Slc14a2 > 15] = 70
# mKid.merge$Merge2[abs(mKid.merge$Merge2)<5] = 0

pdf("3_two_color_2D_map.Grip1_Slc14a2.max50_30.only_CD.E15_E16_E18.pdf")
ggplot(mKid.merge@meta.data[mKid.merge$orig.ident=="E18_3",],aes(x=nB,y=nA,color=Merge2)) +
  geom_point(shape=15,size=1.9) + theme_minimal() +
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB")) +
  expand_limits(y=c(0,96),x=c(0,96)) +
  scale_colour_gradientn(colours = c("blue4","#F0F0F0","red4","purple"),limits = c(-30,70),values = c(0,3/10,8/10,1))
for (Samplei in c("E15_4","E16_1","E18_3")) {
  p1=ggplot(mKid.merge@meta.data[mKid.merge$orig.ident==Samplei,],aes(x=nB,y=nA,color=Merge2)) +
    geom_point(shape=15,size=1.9,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,96),x=c(0,96)) +
    scale_colour_gradientn(colours = c("blue4","#F0F0F0","red4","purple"),limits = c(-30,70),values = c(0,3/10,8/10,1))
  print(p1)
}

## Ligation receptor 2D map ##
mKid.E15_4$shape = "background"
mKid.E15_4$shape[mKid.E15_4$dep_bin %in% 1:4] = "ligation"
mKid.E15_4$shape[mKid.E15_4$Ann_v5 %in% c(9)] = "receptor"
shape.ls = c("background"=15,"ligation"=19,"receptor"=17)
mKid.E15_4$Epha4[!(mKid.E15_4$Ann_v5 %in% c(13,32,5,9,1))] = 0
mKid.E15_4$Merge2 = mKid.E15_4$Efna5
k=2/3
mKid.E15_4$Merge2[mKid.E15_4$Epha4 > k*mKid.E15_4$Efna5] = 
  -mKid.E15_4$Epha4[mKid.E15_4$Epha4 > k*mKid.E15_4$Efna5]
mKid.E15_4$Merge2[mKid.E15_4$Merge2 > 40] = 40
mKid.E15_4$Merge2[mKid.E15_4$Merge2 < -30] = -30
mKid.E15_4$Merge2[mKid.E15_4$Efna5 > 20 & mKid.E15_4$Epha4 > 15] = 50

mKid.E15_4$Merge2[abs(mKid.E15_4$Merge2)<6] = 0
pdf("Plot/2_use_depth/E15_4/4_two_color_2D_map.Efna5_Epha4.max40_30.C9_receptor_3.min6.only_Fib_Epha4.pdf")
ggplot(mKid.E15_4@meta.data,aes(x=nB,y=nA,color=Merge2)) +
  geom_point(aes(shape=shape),size=1.9) + theme_minimal() +
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB")) +
  expand_limits(y=c(0,96),x=c(0,96)) +
  scale_colour_gradientn(colours = c("blue4","#F0F0F0","red4","purple"),limits = c(-30,50),values = c(0,3/8,7/8,1)) +
  scale_shape_manual(values = shape.ls)
ggplot(mKid.E15_4@meta.data,aes(x=nB,y=nA,color=Merge2)) +
  geom_point(aes(shape=shape),size=1.9,show.legend=F) + theme_void() +
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB")) +
  expand_limits(y=c(0,96),x=c(0,96)) +
  scale_colour_gradientn(colours = c("blue4","#F0F0F0","red4","purple"),limits = c(-30,50),values = c(0,3/8,7/8,1)) +
  scale_shape_manual(values = shape.ls)
dev.off()


mKid.merge$shape = "background"
mKid.merge$shape[mKid.merge$dep_bin %in% 1:4] = "ligation"
mKid.merge$shape[mKid.merge$Ann_v5 %in% c(9)] = "receptor"
shape.ls = c("background"=15,"ligation"=19,"receptor"=17)
mKid.merge$Efna5 = mKid.merge@assays[["SCT"]]@data["Efna5",]
mKid.merge$Epha4 = mKid.merge@assays[["SCT"]]@data["Epha4",]
mKid.merge$Epha4[!(mKid.merge$Ann_v5 %in% c(13,32,5,9,1))] = 0
mKid.merge$Merge2 = mKid.merge$Efna5
max1 = 25
min1 = 0
max2 = 18
min2 = 0
k=max2/max1
ranges = max1+max2-min1-min2+10
mKid.merge$Merge2[mKid.merge$Epha4 > k*mKid.merge$Efna5] = 
  -mKid.merge$Epha4[mKid.merge$Epha4 > k*mKid.merge$Efna5]
mKid.merge$Merge2[mKid.merge$Merge2 > max1] = max1
mKid.merge$Merge2[mKid.merge$Merge2 < -1*max2] = -1*max2
mKid.merge$Merge2[mKid.merge$Efna5 > max1/2 & mKid.merge$Epha4 > max2/2] = max1+10
mKid.merge$Merge2[mKid.merge$Merge2>0] = pmax(0,mKid.merge$Merge2[mKid.merge$Merge2>0]-min1)
mKid.merge$Merge2[mKid.merge$Merge2<0] = pmin(0,mKid.merge$Merge2[mKid.merge$Merge2<0]+min2)
pdf(paste0("Plot_Epha4_Efna5_pathway/4_two_color_2D_map.Efna5_Epha4.max",max1,"_",max2,".min",min1,"_",min2,".only_Fib_Epha4.pdf"))
ggplot(mKid.merge@meta.data[mKid.merge$orig.ident=="E15_4",],aes(x=nB,y=nA,color=Merge2)) +
  geom_point(aes(shape=shape),size=1.7) + theme_minimal() +
  scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB")) +
  expand_limits(y=c(0,96),x=c(0,96)) +
  scale_colour_gradientn(colours = c(paletteer_d("RColorBrewer::PiYG",direction = -1),"purple"),
                         limits = c(-max2+min2,max1-min1+10),values = c(0,(max2-min2)/ranges*1:5/5,(max1-min1)/ranges*1:5/5+(max2-min2)/ranges,1)) +
  scale_shape_manual(values = shape.ls)
for (samplei in c("E15_4","E16_1","E18_1","E18_3")) {
  p1 = ggplot(mKid.merge@meta.data[mKid.merge$orig.ident==samplei,],aes(x=nB,y=nA,color=Merge2)) +
    geom_point(aes(shape=shape),size=1.7,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB")) +
    expand_limits(y=c(0,96),x=c(0,96)) +
    scale_colour_gradientn(colours = c(paletteer_d("RColorBrewer::PiYG",direction = -1),"purple"),
                           limits = c(-max2+min2,max1-min1+10),values = c(0,(max2-min2)/ranges*1:5/5,(max1-min1)/ranges*1:5/5+(max2-min2)/ranges,1)) +
    scale_shape_manual(values = shape.ls)
  print(p1)
}
dev.off()
