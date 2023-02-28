## Import depth ##
get_min_dist <- function(df1,df2){
  distA.m = df1$nA %*% matrix(1,nrow = 1,ncol = dim(df2)[1]) - matrix(1,nrow = dim(df1)[1],ncol = 1) %*% t(df2$nA)
  distB.m = df1$nB %*% matrix(1,nrow = 1,ncol = dim(df2)[1]) - matrix(1,nrow = dim(df1)[1],ncol = 1) %*% t(df2$nB)
  dist.m = sqrt(distA.m**2 + distB.m**2)
  return(apply(dist.m,1,min))
}

bd.df = data.df[!is.na(data.df$bound),]
bd0.df = bd.df[bd.df$bound==0,]
bd1.df = bd.df[bd.df$bound==1,]
bd2.df = bd.df[bd.df$bound==2,]
bd3.df = bd.df[bd.df$bound==3,]
bd4.df = bd.df[bd.df$bound==4,]
bd5.df = bd.df[bd.df$bound==5,]

data.df$dist0 = get_min_dist(data.df,bd0.df)
data.df$dist1 = get_min_dist(data.df,bd1.df)
data.df$dist2 = get_min_dist(data.df,bd2.df)
data.df$dist3 = get_min_dist(data.df,bd3.df)
data.df$dist4 = get_min_dist(data.df,bd4.df)
data.df$dist5 = get_min_dist(data.df,bd5.df)

data.df$layer_depth = NA
data.df$total_depth = NA
tmp = data.df[data.df$layer==1,]
data.df$layer_depth[data.df$layer==1] = tmp$dist0/(tmp$dist0+tmp$dist1)
tmp = data.df[data.df$layer==2,]
data.df$layer_depth[data.df$layer==2] = tmp$dist1/(tmp$dist1+tmp$dist2)
tmp = data.df[data.df$layer==3,]
data.df$layer_depth[data.df$layer==3] = tmp$dist2/(tmp$dist2+tmp$dist3)
tmp = data.df[data.df$layer==4,]
data.df$layer_depth[data.df$layer==4] = tmp$dist3/(tmp$dist3+tmp$dist4)
tmp = data.df[data.df$layer==5,]
data.df$layer_depth[data.df$layer==5] = tmp$dist4/(tmp$dist4+tmp$dist5)
data.df$layer = as.numeric(as.character(data.df$layer))
data.df$total_depth = data.df$layer_depth + data.df$layer - 1

# devide into 200 cell bins
mKid.merge$layer = 0
mKid.merge$depth = -1
mKid.merge$dep_bin = 0

cell.ls = intersect(mKid.merge$coord,row.names(data.df[data.df$layer!=0,]))
mKid.merge$layer[cell.ls] = data.df[cell.ls,"layer"]
mKid.merge$depth[cell.ls] = data.df[cell.ls,"total_depth"]

N = 0
for (layi in 1:5) {
  subi = mKid.merge@meta.data[mKid.merge$layer==layi & mKid.merge$orig.ident=="E18_3",]
  subi = subi[order(subi$depth),]
  Ni = round(dim(subi)[1]/200)
  subi$dep_bin = ceiling((1:dim(subi)[1])/200)
  subi$dep_bin = pmin(subi$dep_bin,Ni)
  subi$dep_bin = subi$dep_bin + N
  mKid.merge$dep_bin[subi$coord] = subi$dep_bin
  N = N + Ni
  print(Ni)
}
