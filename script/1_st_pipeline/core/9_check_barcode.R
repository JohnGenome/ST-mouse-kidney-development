library(getopt)
command=matrix(c("combine_barcode","b",1,"character",
                 "work_path","p",1,"character",
                 "file_name","n",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$combine_barcode) || is.null(args$work_path) || is.null(args$file_name)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

work_path = args$work_path
inpu_name = args$file_name
refe_path = args$combine_barcode

setwd(work_path)
library(ggplot2)
library(ComplexHeatmap)
# 1
ref_barcode = read.table(refe_path,stringsAsFactors = F)
our_barcode = read.table(inpu_name,stringsAsFactors = F)$V1
print("total barcode num")
length(our_barcode)

# 2
our_barcode = our_barcode[our_barcode %in% ref_barcode$V1]
print("barcode in ref list")
length(our_barcode)
our_barcode = factor(our_barcode,levels = ref_barcode$V1)
ref_barcode$count = table(our_barcode)

nchannel = max(ref_barcode$V2)
all_combi.m = matrix(0,nchannel,nchannel) # barcodeA*barcodeB
row.names(all_combi.m) = paste0("A",1:nchannel)
colnames(all_combi.m) = paste0("B",1:nchannel)
tmp <- lapply(1:dim(ref_barcode)[1], function(i) all_combi.m[ref_barcode$V2[i],ref_barcode$V3[i]] <<- ref_barcode$count[i])

# from B96 --> B1
all_combi.m = all_combi.m[,nchannel:1]
nwid = 10
if(nchannel>90){nwid = 15}
pdf("0_barcode_combination.pdf",width = nwid,height = nwid)
Heatmap(all_combi.m,cluster_rows=F,cluster_columns=F)
Heatmap(log2(all_combi.m+1),cluster_rows=F,cluster_columns=F)
dev.off()

library(circlize)
nlwd = 5
if(nchannel>90){nlwd = 2}
svg("1_barcode_combination.svg",width = 6,height = 6)
Heatmap(all_combi.m,cluster_rows=F,cluster_columns=F,rect_gp = gpar(col = "white", lwd = nlwd),
        col = colorRamp2(c(0, quantile(all_combi.m,0.9)/4, max(all_combi.m)), c("blue", "green", "red")),
        show_heatmap_legend = F,show_row_names=F,show_column_names=F)
dev.off()

