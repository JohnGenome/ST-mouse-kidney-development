library(getopt)
command=matrix(c("core","q",1,"integer",
                 "Log","l",1,"character",
                 "Path","d",1,"character",
                 "Input_name","i",1,"character",
                 "Barcode_name","b",1,"character",
                 "order_path","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$core) || is.null(args$Log) || is.null(args$Path)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

core_num = args$core
Log = args$Log
work_path = args$Path
inpu_name = args$Input_name
refe_path = args$Barcode_name
order_path = args$order_path

setwd(work_path)
library(parallel)
library(ComplexHeatmap)

sink(Log,append=TRUE,split=FALSE)
cat("\nStarting to evaluate cross-channel:",inpu_name,"\n")
seq.ls = read.table(inpu_name,stringsAsFactors = F)$V1
order.df = read.table(order_path,stringsAsFactors = F,header = F)
ref_barcode = read.table(refe_path,stringsAsFactors = F)$V1[order.df$V2]
# 1
get_barcode_num <- function(s1){
  sub1 = strsplit(s1,split = "ATCCACGTGCTTGAG")[[1]]
  return(length(sub1)-1)
}
check_consistent <- function(s1){
  sub1 = strsplit(s1,split = "ATCCACGTGCTTGAG")[[1]]
  barcode1 = unlist(lapply(sub1[-length(sub1)], function(x) substr(x, nchar(x)-8+1, nchar(x))))
  return(length(unique(barcode1)))
}

count.ls = unlist(lapply(seq.ls, get_barcode_num))
unique.ls = unlist(lapply(seq.ls, check_consistent))

cat("2.1 Total reads:",length(seq.ls),"\n")
print(table(unique.ls))
print(table(count.ls,unique.ls))


# 2
get_barcode_B <- function(s1,n=1){
  sub1 = strsplit(s1,split = "ATCCACGTGCTTGAG")[[1]]
  barcode1 = unlist(lapply(sub1[-length(sub1)], function(x) substr(x, nchar(x)-8+1, nchar(x))))
  return(rep(unique(barcode1),2)[n])
}

comb.1.ls = unlist(mclapply(seq.ls, get_barcode_B, mc.cores = core_num))
comb.2.ls = unlist(mclapply(seq.ls, get_barcode_B, n = 2, mc.cores = core_num))
comb.df = data.frame("B1"=comb.1.ls,"B2"=comb.2.ls,stringsAsFactors = F)

comb.df$B1 = factor(comb.df$B1,levels = ref_barcode)
comb.df$B2 = factor(comb.df$B2,levels = ref_barcode)
comb.df = comb.df[!(is.na(comb.df$B1) | is.na(comb.df$B2)),]
comb.m = as.data.frame.matrix(table(comb.df))
cat("2.2 Proper barcode reads:",dim(comb.df)[1],"\n")

comb.m = as.matrix(comb.m)
comb.m.filter = as.matrix(comb.m)
diag(comb.m.filter) = 0

pdf(paste0("2_multiple_barcode2_",gsub("_2ormore_linker2.read2.txt","",inpu_name),".pdf"),width = 10,height = 10)
Heatmap(comb.m,cluster_rows=F,cluster_columns=F)
Heatmap(comb.m.filter,cluster_rows=F,cluster_columns=F)
dev.off()




