library(Seurat)

## export to monocle3 ##
# Matrix
expression_matrix = as.matrix(mKid.merge.CDs@assays$SCT@counts)
cell_metadata = mKid.merge.CDs@meta.data
gene_metadata = data.frame("gene_short_name"=row.names(mKid.merge.CDs@assays$SCT@counts),
                           row.names=row.names(mKid.merge.CDs@assays$SCT@counts))
saveRDS(expression_matrix,"Out/1_SCT.CD_system_expression_matrix.RDS")
saveRDS(cell_metadata,"Out/1_SCT.CD_system_cell_metadata.RDS")
saveRDS(gene_metadata,"Out/1_SCT.CD_system_gene_metadata.RDS")

# UMAP coordinate
write.csv(mKid.merge.CDs@reductions$pca@cell.embeddings,"Out/3_PCA_embedding.csv")
write.csv(mKid.merge.CDs@reductions$umap@cell.embeddings,"Out/3_UMAP_embedding.csv")
