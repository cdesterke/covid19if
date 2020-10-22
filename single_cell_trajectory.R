# R code to performed single cell trajectory

# input data import matrix and meta data export for Seurat object after cell selection
expression<-read.table("matrix.txt",h=T)
meta<-read.table("meta.txt,h=T)
library(monocle)

# Create monocle single cell object with matrix (expression) and meta data (meta)
mat<-as.matrix(expression)
colonnes<-as.data.frame(colnames(mat))
colonnes$group<-meta$orig.ident
colnames(colonnes)<-c("id","group")
row.names(colonnes)<-colonnes$id
pd <- new("AnnotatedDataFrame", data = colonnes)
genes<-as.data.frame(row.names(matrix))
colnames(genes)<-"gene_short_name"
row.names(genes)<-genes$gene_short_name
fd <- new("AnnotatedDataFrame", data = genes)
cds <- newCellDataSet(mat, phenoData = pd,featureData = fd)

# defined cell hierarchy on makers expression
cth <- newCellTypeHierarchy()
PSMB8_id <- row.names(subset(fData(cds), gene_short_name == "PSMB8"))
cth <- addCellType(cth, "PSMB8_HIGH", classify_func = function(x) { x[PSMB8_id,] >= 8})
cth <- addCellType(cth, "PSMB8_MEDIUM", classify_func = function(x) { x[PSMB8_id,] < 8 & x[PSMB8_id,] > 3} )
cth <- addCellType(cth, "PSMB8_LOW", classify_func = function(x) { x[PSMB8_id,] <= 3 })
cds <- classifyCells(cds, cth)
table(pData(cds)$CellType)

# object transformations
meta<-as.data.frame(colnames(table))
meta$group<-Idents(monocytes)
colnames(meta)<-c("identifier","group")
row.names(meta)<-meta$identifier
cds@expressionFamily = negbinomial.size()
cds@lowerDetectionLimit = 0
my_feat <- fData(cds)
my_feat$id<-my_feat$gene_short_name
head(my_feat)

# size factor and dispersion
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# dispersion table and variance 
disp_table <- dispersionTable(cds)
table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds, return_all = FALSE)
cds <- reduceDimension(cds, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = TRUE)

# clustering cells
cds <- clusterCells(cds)
table(pData(cds)$CellType)

# ggplot graph cell type
pie <- ggplot(pData(cds),
aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# cluster cells with hiearchy
cds <- clusterCells(cds,cth) 
plot_cell_clusters(cds, 1, 2, color = "Cluster") +  facet_wrap(~CellType)
plot_cell_clusters(cds, 1, 2, color = "CellType")
save(cds,file="monohiearchy.rda")

# differential expressed genes against cell type
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~CellType")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
my_ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
cds2 <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)

# dimensional reduction and gene selection
cds2 <- reduceDimension(cds, method = 'DDRTree')
gene_to_cluster <- row.names(diff_test_res)[order(diff_test_res$qval)][1:89] 
conca<-c("PSMB8",gene_to_cluster)
cds2 <- orderCells(cds2)

# build graphs on trajectory 
plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="FABP4")
plot_cell_trajectory(cds2, color_by = "CellType",)
plot_cell_trajectory(cds2, color_by = "group",)
plot_cell_trajectory(cds2, color_by = "group",markers="PSMB8",markers_linear = T,show_branch_points=T)
plot_cell_trajectory(cds2, color_by = "Pseudotime")
plot_cell_trajectory(cds2, color_by = "CellType") +  facet_wrap(~group)
my_pseudotime_cluster <- plot_pseudotime_heatmap(cds2[conca,],cores = 8,
show_rownames = TRUE,return_heatmap = TRUE,cluster_rows = TRUE)
plot_genes_in_pseudotime(cds2[c("PSMB8","CALR","TAP1"),],cell_size = 2, color_by = "group",ncol = 1)

# save data
write.table(pData(cds),file="phenotype.txt")
write.table(diff_test_res,file="difexpCD9trajectory.txt")
save(cds,file="monoclecds.rda")
save(cds2,file="monoclecdsfiltre.rda")

# BEAM 3 analysis
BEAM_res <- BEAM(cds2, branch_point = 3, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

sig_gene_names2 <- row.names(subset(BEAM_res, qval < 0.1))
my_ordering_genes <- row.names(BEAM_res)[order(BEAM_res$qval)][1:1000]
gene_to_cluster2 <- row.names(BEAM_res)[order(BEAM_res$qval)][1:75] 
conca2<-c("PSMB8",gene_to_cluster2)
plot_genes_branched_heatmap(cds2[conca2,],branch_point = 3,num_clusters = 2, cores = 1, use_gene_short_name = T, show_rownames = T)
genes <- row.names(subset(fData(cds2), gene_short_name %in% c("PSMB8", "MARCO", "CD68","TNFSF13","CALR")))
plot_genes_branched_pseudotime(cds2[genes,], branch_point = 2, color_by = "group", ncol = 1)






