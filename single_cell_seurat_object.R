# create merged Seurat object with h5 input files 

library(Seurat)
library(hdf5r)


C51_HD <- Read10X_h5('GSM4475048_C51_HD.h5')

C52_HD <- Read10X_h5('GSM4475049_C52_HD.h5')

C100_HD <- Read10X_h5('GSM4475050_C100_HD.h5')

C148_HD <- Read10X_h5('GSM4475051_C148_HD.h5')

C149_HD <- Read10X_h5('GSM4475052_C149_HD.h5')

C152_HD <- Read10X_h5('GSM4475053_C152_HD.h5')


C141_mildCOVID <- Read10X_h5('GSM4339769_C141_mildCOVID_BALF.h5')

C142_mildCOVID <- Read10X_h5('GSM4339770_C142_mildCOVID_BALF.h5')

C144_mildCOVID <- Read10X_h5('GSM4339772_C144_mildCOVID_BALF.h5')


C143_severeCOVID <- Read10X_h5('GSM4339771_C143_severeCOVID_BALF.h5')

C145_severeCOVID <- Read10X_h5('GSM4339773_C145_severeCOVID_BAFL.h5')

C146_severeCOVID <- Read10X_h5('GSM4339774_C146_severeCOVID_BAFL.h5')


control_51 <- CreateSeuratObject(counts = C51_HD, project = "control", min.cells = 3, min.features = 200)

control_52 <- CreateSeuratObject(counts = C52_HD, project = "control", min.cells = 3, min.features = 200)

control_100 <- CreateSeuratObject(counts = C100_HD, project = "control", min.cells = 3, min.features = 200)

control_148 <- CreateSeuratObject(counts = C148_HD, project = "control", min.cells = 3, min.features = 200)

control_149 <- CreateSeuratObject(counts = C149_HD, project = "control", min.cells = 3, min.features = 200)

control_152 <- CreateSeuratObject(counts = C152_HD, project = "control", min.cells = 3, min.features = 200)




mildcovid_141 <- CreateSeuratObject(counts = C141_mildCOVID, project = "mc", min.cells = 3, min.features = 200)

mildcovid_142 <- CreateSeuratObject(counts = C142_mildCOVID, project = "mc", min.cells = 3, min.features = 200)

mildcovid_144 <- CreateSeuratObject(counts = C144_mildCOVID, project = "mc", min.cells = 3, min.features = 200)


severecovid_143 <- CreateSeuratObject(counts = C143_severeCOVID, project = "sc", min.cells = 3, min.features = 200)

severecovid_145 <- CreateSeuratObject(counts = C145_severeCOVID, project = "sc", min.cells = 3, min.features = 200)

severecovid_146 <- CreateSeuratObject(counts = C146_severeCOVID, project = "sc", min.cells = 3, min.features = 200)



all <- merge(control_51, y = c(control_52, control_100, control_148, control_149, control_152 , mildcovid_141, mildcovid_142, mildcovid_144, severecovid_143, severecovid_145, severecovid_146 ), add.cell.ids = c("ct51", "ct52","ct100","ct148","ct149","ct152","mc141", "mc142", "mc144", "sc143","sc145","sc146"), project = "covid")


# normalization

list <- SplitObject(all, split.by = "orig.ident")
list <- lapply(X = list, FUN = function(x) {
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


anchors <- FindIntegrationAnchors(object.list = list, dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) <- "integrated"

# scale data and dimension reductions

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)

# plot outputs
## dotplot
jpeg(file = "DOTPLOTMARKERS.jpeg",  width = 700, height = 300)
DotPlot(combined, features = rev(c("TAP1","TAP2","TAPBPL","PSMB8","PSMB8-AS1","PSMB9","PSMB10","CALR","CANX")),
cols=c("aquamarine4","darkorchid3","darkorange"), dot.scale = 16, split.by = "orig.ident",assay = "RNA") +
RotatedAxis()
dev.off()

## feature plot
jpeg(file = "biplotCD163PSMB10.jpg",  width = 250, height = 250)
FeatureScatter(combined, feature1 = "rna_CD163", feature2 = "rna_PSMB10", slot="counts",cols=c("aquamarine4","darkorchid3","darkorange"))
dev.off()

## dimplot
jpeg(file = "dimplotumap.jpg",  width = 4, height = 4)
DimPlot(combined, reduction = "umap",cols=c("aquamarine4","darkorchid3","darkorange")) 
dev.off()

## violinplot
jpeg(file = "TAP1VLN.jpg", width = 300, height = 300)
VlnPlot(combined, features = c("TAP1"), slot = "data", log = TRUE,split.by= "orig.ident",pt.size=0.1,cols=c("aquamarine4","darkorchid3","darkorange"))
dev.off()

## biplot
jpeg(file = "biplotCD68&TAP1.jpg",  width = 250, height = 250)
FeatureScatter(combined, feature1 = "rna_CD68", feature2 = "rna_TAP1", slot="counts",cols=c("aquamarine4","darkorchid3","darkorange"))
dev.off()

# CD14+/CD68+ subseting in Seurat
m1mc<-WhichCells(all, idents = "mc", expression = CD14 > 5 & CD68 > 5 & PSMB8  >=3, slot ="counts")
m1sc<-WhichCells(all, idents = "sc", expression = CD14 > 5 & CD68 > 5 & PSMB8 < 3, slot ="counts")
length(m1sc)
mixmono<-c(m1mc,m1sc)
length(mixmono)
monocytes <- SubsetData(object = all, cells = mixmono)
save(monocytes,file="m1.rda")

expression <- GetAssayData(object = monocytes, slot = "counts")
table<-as.matrix(expression)
write.table(table,file="m1matrix.txt")
meta<-monocytes[[]]
write.table(meta,file="m1meta.txt")

### end of code






