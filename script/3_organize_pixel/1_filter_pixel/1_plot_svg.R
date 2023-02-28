for (Samplei in Sample.ls) {
  svg(paste0("1_spatial_map.",Samplei,".svg"))
  p1 = ggplot(mKid.merge@meta.data[mKid.merge$orig.ident==Samplei,],aes(x=nB,y=nA,color=Ann_v5)) +
    geom_point(shape=15,size=1.5,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,96),x=c(0,96)) +
    scale_color_manual(values = color.ls)
  print(p1)
  dev.off()
}
