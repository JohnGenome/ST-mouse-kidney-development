library(getopt)
command=matrix(c("order_path","1",1,"character",
                 "refe_path","2",1,"character",
                 "Log","l",1,"character",
                 "output_path","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$order_path) || is.null(args$refe_path) || is.null(args$Log)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

order_path = args$order_path
refe_path = args$refe_path
Log = args$Log
output_path = args$output_path

## 1
sink(Log,append=TRUE,split=FALSE)
cat("Using custom order\n")

order.df = read.table(order_path,stringsAsFactors = F,header = F)
b100 = read.table(refe_path,stringsAsFactors = F,header = F)$V1
nchannel = dim(order.df)[1]

comb.5050 = data.frame("b1"=rep(b100[order.df$V1],each=nchannel),
                       "b2"=rep(b100[order.df$V2],times=nchannel),
                       "i1"=rep(1:nchannel,each=nchannel),
                       "i2"=rep(1:nchannel,times=nchannel),stringsAsFactors = F)
comb.5050$comb = paste0(comb.5050$b2,comb.5050$b1)

write.table(comb.5050[,c("comb","i1","i2")],file = output_path,
            row.names = F,col.names = F,sep = "\t",quote = F)



