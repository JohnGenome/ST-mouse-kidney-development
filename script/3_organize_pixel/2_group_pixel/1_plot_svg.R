# extend nB from 1-96 to 0-105 for E16_1
data.df = data.frame("orig.ident"=rep(Sample.ls,each=(96*106)),
                     "nA"=rep(rep(1:96,each=106),17),
                     "nB"=rep(rep(0:105,times=96),17),
                     "Ann_v5"=factor(33,levels = 0:33))
row.names(data.df) = paste0(data.df$orig.ident,"_",data.df$nA,"x",data.df$nB)
data.df[mKid.merge$coord,"Ann_v5"] = mKid.merge$Ann_v5

color.ls4 = c(color.ls4,"33"="#f0ffff")
for (Samplei in Sample.ls) {
  svg(paste0("Plot/0_total_pixel_svg/0_spatial_map.",Samplei,".svg"))
  p1 = ggplot(data.df[data.df$orig.ident==Samplei,],aes(x=nB,y=nA,color=Ann_v5)) +
    geom_point(shape=15,size=1.5,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,96),x=c(0,105)) +
    scale_color_manual(values = color.ls4)
  print(p1)
  dev.off()
}



