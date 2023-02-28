## Decode coord ##
data.df = data.frame("orig.ident"=rep("E18_3",each=(96*106)),
                     "nA"=rep(rep(1:96,each=106),1),
                     "nB"=rep(rep(0:105,times=96),1),
                     "Ann_v5"=factor(33,levels = 0:33))
row.names(data.df) = paste0(data.df$orig.ident,"_",data.df$nA,"x",data.df$nB)
data.df[mKid.E18_3$coord,"Ann_v5"] = mKid.E18_3$Ann_v5
data.df$svgx = read.table("Plot/1_ai_edited_svg/0_spatial_map.E18_3.coord")$V1
data.df$svgy = read.table("Plot/1_ai_edited_svg/0_spatial_map.E18_3.coord")$V2

## Decode layer ##
L1.df = read.table("Plot/1_ai_edited_svg_E18/2_spatial_map.E18_3.L1.coord")
L2.df = read.table("Plot/1_ai_edited_svg_E18/2_spatial_map.E18_3.L2.coord")
L3.df = read.table("Plot/1_ai_edited_svg_E18/2_spatial_map.E18_3.L3.coord")
L4.df = read.table("Plot/1_ai_edited_svg_E18/2_spatial_map.E18_3.L4.coord")
L5.df = read.table("Plot/1_ai_edited_svg_E18/2_spatial_map.E18_3.L5.coord")
data.df$layer = factor(0,levels = 0:5)
data.df[paste0(data.df$svgx,data.df$svgy) %in% paste0(L1.df$V1,L1.df$V2),"layer"] = "1"
data.df[paste0(data.df$svgx,data.df$svgy) %in% paste0(L2.df$V1,L2.df$V2),"layer"] = "2"
data.df[paste0(data.df$svgx,data.df$svgy) %in% paste0(L3.df$V1,L3.df$V2),"layer"] = "3"
data.df[paste0(data.df$svgx,data.df$svgy) %in% paste0(L4.df$V1,L4.df$V2),"layer"] = "4"
data.df[paste0(data.df$svgx,data.df$svgy) %in% paste0(L5.df$V1,L5.df$V2),"layer"] = "5"
rm(list = c("L1.df","L2.df","L3.df","L4.df","L5.df"))

## Decode boundry ##
bound.df= read.table("Plot/1_ai_edited_svg_E18/3_spatial_map.E18_3.bound0.coord")
pelvis.df = read.table("Plot/1_ai_edited_svg_E18/3_spatial_map.E18_3.bound5.coord")
data.df$bound = factor(NA,levels = 0:5)
data.df[paste0(data.df$svgx,data.df$svgy) %in% paste0(bound.df$V1,bound.df$V2),"bound"] = "0"
data.df$bound[data.df$bound==0 & data.df$layer==2] = "1" 
data.df$bound[data.df$bound==0 & data.df$layer==3] = "2" 
data.df$bound[data.df$bound==0 & data.df$layer==4] = "3"
data.df$bound[data.df$bound==0 & data.df$layer==5] = "4"
data.df[paste0(data.df$svgx,data.df$svgy) %in% paste0(pelvis.df$V1,pelvis.df$V2),"bound"] = "5"
rm(list = c("bound.df","pelvis.df"))

# data.df[1:3,]
#          orig.ident nA nB Ann_v5  svgx svgy layer bound
# E18_3_1x0      E18_3  1  0     33 479.1 29.2     0  <NA>
# E18_3_1x1      E18_3  1  1     33 474.8 29.2     0  <NA>
# E18_3_1x2      E18_3  1  2     33 470.4 29.2     0  <NA>


