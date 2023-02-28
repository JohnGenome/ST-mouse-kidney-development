library(monocle3)
library(dplyr)
library(ggplot2)
library(scales)

expression_matrix = readRDS("Out/1_SCT.CD_system_expression_matrix.RDS")
cell_metadata = readRDS("Out/1_SCT.CD_system_cell_metadata.RDS")
gene_metadata = readRDS("Out/1_SCT.CD_system_gene_metadata.RDS")

## 1 CDS obeject ##
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, method = "PCA", num_dim = 50, norm_method = "size_only", scaling = TRUE)
cds <- reduce_dimension(cds, reduction_method = "UMAP") # c("UMAP", "tSNE", "PCA", "LSI", "Aligned")

## 2 Load Seurat embedding #
pca.df = as.matrix(read.csv("Out/3_PCA_embedding.csv",row.names = 1))
UMAP.df = as.matrix(read.csv("Out/3_UMAP_embedding.csv",row.names = 1))
cds@int_colData$reducedDims$PCA = pca.df
cds@int_colData$reducedDims$UMAP = UMAP.df

pdf("Plot/3_monocle3.CD_system.using_Seurat_embedding.pdf")
plot_cells(cds, label_groups_by_cluster=T,group_label_size = 4,  color_cells_by = "Ann_v5",cell_size = 0.65) + scale_color_manual(values = color.ls5)
plot_cells(cds, label_groups_by_cluster=T,group_label_size = 4,  color_cells_by = "orig.ident",cell_size = 0.65)
plot_cells(cds, label_groups_by_cluster=T,group_label_size = 4,  color_cells_by = "stage",cell_size = 0.65) +
  scale_color_manual(values = c("#FF0000","#CDCD00","#00CD00","#0000CD","#FFFF00","#c4aead"))
dev.off()

## 3 Partition ##
set.seed(2022)
cds <- cluster_cells(cds) # resolution = 10
pdf("Plot/4_monocle.CD_system.partition.pdf")
plot_cells(cds,group_label_size = 4, color_cells_by = "partition",cell_size = 0.65)
dev.off()

cds <- learn_graph(cds,learn_graph_control = list("euclidean_distance_ratio"=1,"geodesic_distance_ratio"=1/3,
                                                  "minimal_branch_len"=4,"orthogonal_proj_tip"=F,"prune_graph"=T,"rann.k"=25))
pdf("Plot/4_monocle.CD_system.trajectory.pdf")
plot_cells(cds,
           color_cells_by = "Ann_v5",
           label_groups_by_cluster=T,group_label_size = 4,cell_size = 0.65,
           label_leaves=FALSE,
           label_branch_points=FALSE) + scale_color_manual(values = color.ls5)
plot_cells(cds,
           color_cells_by = "Ann_v5",
           label_groups_by_cluster=T,group_label_size = 4,cell_size = 0.65,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=4) + scale_color_manual(values = color.ls5)
dev.off()

## 4 identify the root principal points ##
cell_ids <- which(colData(cds)[, "stage"] == "E11") # 120 cells
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex # return each cell's nearest vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name # get 161 vertex names
# root_pr_nodes <- root_pr_nodes[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
root_pr_nodes <- root_pr_nodes[c(2,93)]

cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
pdf("Plot/4_monocle.CD_system.trajectory3.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F, cell_size = 0.65,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F, cell_size = 0.65,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=1.5)
dev.off()

write.csv(as.matrix(cds@principal_graph_aux$UMAP$pseudotime),file = "Plot/4_monocle.CD_system.trajectory.pseudotime.csv")

