final_cell.ls = c()
for (si in sample.ls) {
  cord.df = read.table(paste0("Filter_pixel_Annv5/origin/1_spatial_map.",si,".coord"))
  cord_crop.df = read.table(paste0("Filter_pixel_Annv5/outFinal/3_spatial_map.",si,".final.coord"))
  
  cord.df$cellid = row.names(mKid.merge@meta.data)[mKid.merge$orig.ident==si]
  cord.df$isfinal = ifelse(paste0(cord.df$V1,"_",cord.df$V2) %in% 
                             paste0(cord_crop.df$V1,"_",cord_crop.df$V2),1,0)
  final_cell.ls = c(final_cell.ls,cord.df$cellid[cord.df$isfinal==1])
}

mKid.merge$isfinal = 0
mKid.merge@meta.data[final_cell.ls,"isfinal"] = 1
