library(tidyverse)
library(sctransform)
library(Seurat)

setwd("~/cellranger_out_temp/")
getwd()
sc_list <- list()
files <- c("./aCSF1R/outs/filtered_feature_bc_matrix/",
           "./D14/outs/filtered_feature_bc_matrix/",
           "./IgG/outs/filtered_feature_bc_matrix/",
           "./Naive/outs/filtered_feature_bc_matrix/")
sc_list <- lapply(files, Read10X)
names(sc_list) <- c("aCSF1R", "D14", "IgG", "Naive")


######## Read in data and SCTRANSFORM ############

sc_obj2 <- list()

for(i in 1:length(sc_list)) {
  
  sc_obj2[[i]] <- CreateSeuratObject(counts = sc_list[[i]], project = names(sc_list)[i], 
                                     min.cells = 3, min.features = 200)
  names(sc_obj2)[i] <- names(sc_list)[i]
  
}
class(sc_obj2[["Naive"]])

library(Seurat)
#all_dat <- merge(x = sc_obj2[["Naive"]], y = c(sc_obj2[["D14"]], sc_obj2[["IgG"]], sc_obj2[["aCSF1R"]]),
#                add.cell.ids = c("Naive", "D14", "IgG", "aCSF1R"), project = "schmid")
all_dat <- merge(x = sc_obj2[["Naive"]], y = c(sc_obj2[["IgG"]], sc_obj2[["aCSF1R"]]),
                 add.cell.ids = c("Naive", "IgG", "aCSF1R"), project = "schmid")

all_dat[["percent.mt"]] <- PercentageFeatureSet(all_dat, pattern = "^mt-")
#"nFeature_RNA", "nCount_RNA",

#VlnPlot(all_dat, features = c("nFeature_RNA"), ncol = 1) + ylim(1000,6500)

all_dat <- subset(all_dat, subset = nFeature_RNA < 6500 & nFeature_RNA > 1500 & percent.mt < 10)

#VlnPlot(all_dat, features = "nCount_RNA", ncol = 1) + ylim(500,50000)

all_dat <- SCTransform(all_dat)

all_dat <- RunPCA(all_dat) %>% 
  RunUMAP(dims = 1:50, n.neighbors = 30) %>%
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = 0.8)
ElbowPlot(all_dat, ndims = 50)

DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident")
DimPlot(all_dat, pt.size = 0.4, label = TRUE, label.size = 6)



############# REMOVE OUTLIER CLUSTERS ##############


#all_dat <- all_dat[,!all_dat$seurat_clusters %in% c(11:14)]
#all_dat <- all_dat[,!all_dat$seurat_clusters %in% c(16,17,19)]
all_dat <- all_dat[,!all_dat$seurat_clusters %in% c(13,14)]
#all_dat <- all_dat[,all_dat@reductions[["umap"]]@cell.embeddings[,2] < 10 & all_dat@reductions[["umap"]]@cell.embeddings[,1] > -10]

DimPlot(all_dat, group.by = "orig.ident")
DimPlot(all_dat, label = T)


############# RECLUSTER ###############

###This is now SCTransform all_dat without D14 from the start

all_dat <- RunPCA(all_dat) %>% 
  RunUMAP(dims = 1:50, n.neighbors = 30) %>%
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = 1)

DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident")
DimPlot(all_dat, pt.size = 0.4, label = TRUE, label.size = 6)

## Extract CAFs

#cafs <- all_dat[,all_dat$SCT_snn_res.1 %in% c(0,2,5,14,8)]

cafs <- all_dat[,all_dat$SCT_snn_res.1 %in% c(0,2,5,13)]
DimPlot(cafs, label = T)
not_cafs <- all_dat[,!all_dat$SCT_snn_res.1 %in% c(0,2,5,13)]
DimPlot(not_cafs)
not_cafs <- not_cafs[,not_cafs@reductions$umap@cell.embeddings[,2] < 2.5]


### For naive alone, data was filtered with more stringent cutoffs due to cells that had slightly odd clustering, potentially due to expresion issues. 
#The cutoffs for filtering for Naive alone were; subset = nFeature_RNA < 6500 & nFeature_RNA > 1500 & percent.mt < 10
### full details of the code can be found in Log_of_project_work


naive_alone <- not_cafs
naive_alone <- naive_alone[,naive_alone$orig.ident == "Naive"]
DimPlot(naive_alone, group.by = "orig.ident")


#n.neighbors was 30 res = 0.5 dims = 50
naive_alone <- RunPCA(naive_alone) %>% 
  FindNeighbors(dims = 1:10, k.param = 20) %>%
  RunUMAP(dims = 1:10, n.neighbors = 50, min.dist = 0.5) %>%
  FindClusters(resolution = 0.3)
ElbowPlot(naive_alone, ndims = 50) +ylim(0,35)
#FeaturePlot(naive_alone, features = "nFeature_SCT")

DimPlot(naive_alone, label = T, pt.size = 2)














VlnPlot(naive_alone, features = c("Gsn", "Clec3b", "Dpt", "Cd34", "Mfap4", "Entpd2", "Fbln2", "Col15a1")) # Fibroblast markers
ggsave("../project/ReRun_AfterMeetingWithMei/Fibroblast_Markers_1105.png", dpi = 1000, width = 10, height = 10, units = "in")
VlnPlot(naive_alone, features = c("Ecm1", "Vipr1", "Colec11", "Reln", "Lrat", "Pth1r", "Hgf", "Rgs5", "Ngfr")) # Hepatic stellate markers
ggsave("../project/ReRun_AfterMeetingWithMei/HepaticStellate_Markers_1105.png", dpi = 1000, width = 10, height = 10, units = "in")
VlnPlot(naive_alone, features = c("Tagin", "Acta2", "Tpm2", "Cnn1", "Myh11", "Pln")) # VSMC markers
ggsave("../project/ReRun_AfterMeetingWithMei/VascularSmoothMuscleCell_Markers_1105.png", dpi = 1000, width = 10, height = 10, units = "in")

naive_alone$CellType <- case_when(naive_alone$SCT_snn_res.0.3 %in% c(0,1) ~ "Fibroblast",
                                  naive_alone$SCT_snn_res.0.3 %in% c(5,2) ~ "HSC",
                                  naive_alone$SCT_snn_res.0.3 %in% c(3,4,6) ~ "VSMC")


Idents(naive_alone) <- naive_alone@meta.data$SCT_snn_res.0.3
DimPlot(naive_alone, pt.size = 0.4) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Naive_Clusters_UMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")

DimPlot(naive_alone, group.by = "CellType", pt.size = 0.4)  + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Naive_FibHSCVSMC_UMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")


VlnPlot(naive_alone, features = c("Gsn", "Clec3b", "Dpt", "Cd34", "Mfap4", "Entpd2", "Fbln2", "Col15a1"), group.by = "CellType")
ggsave("../project/ReRun_AfterMeetingWithMei/Fibroblast_Markers_CellType_1105.png", dpi = 1000, width = 10, height = 10, units = "in")
VlnPlot(naive_alone, features = c("Ecm1", "Vipr1", "Colec11", "Reln", "Lrat", "Pth1r", "Hgf", "Rgs5", "Ngfr"), group.by = "CellType")
ggsave("../project/ReRun_AfterMeetingWithMei/HepaticStellate_Markers_CellType_1105.png", dpi = 1000, width = 10, height = 10, units = "in")
VlnPlot(naive_alone, features = c("Tagin", "Acta2", "Tpm2", "Cnn1", "Myh11", "Pln"), group.by = "CellType")
ggsave("../project/ReRun_AfterMeetingWithMei/VascularSmoothMuscleCell_Markers_CellType_1105.png", dpi = 1000, width = 10, height = 10, units = "in")

DimPlot(naive_alone, group.by = "CellType", pt.size = 2)
ggsave("../project/ReRun_AfterMeetingWithMei/NaiveCellTypeIdentity_1305.png", dpi = 1000, width = 7, height = 7)

DimPlot(naive_alone, label = T, pt.size = 2)
ggsave("../project/ReRun_AfterMeetingWithMei/NaiveIndividClusters_1305.png", dpi = 1000, width = 7, height = 7)
naive_alone$Outlier <- if_else(naive_alone@reductions[["umap"]]@cell.embeddings[,2] < 1 & naive_alone$CellType == "Vascular_Smooth_Muscle", "VSCM_Outlier", "Not")

VlnPlot(naive_alone, features = c("Gdf2"))

Idents(naive_alone) <- naive_alone@meta.data$CellType

naive_alone_cluster_markers <- FindAllMarkers(naive_alone)
naive_alone_cluster_markers <- naive_alone_cluster_markers[abs(naive_alone_cluster_markers$avg_log2FC) > 0.5,]



naive_pubs_markers <- read_csv("./PublicationPlots/Naive_markers_for_heatmap_signature.csv", col_names = F)
colnames(naive_pubs_markers) <- c("Genes", "CellType")
naive_pubs_markers <- naive_pubs_markers[-1,]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

default_gg_cols <- gg_color_hue(4)


heatmaps_data_anno <- data.frame(Cell_Ident = Idents(naive_alone), 
                                 row.names = colnames(naive_alone))

annot_colors=list(Cell_Ident=c(Fibroblast=default_gg_cols[1],HSC=default_gg_cols[2], 
                        VSMC=default_gg_cols[3]),
                  Markers = c(Fibroblast=default_gg_cols[1],HSC=default_gg_cols[2], 
                                   VSMC=default_gg_cols[3]))

marker_anno <- data.frame(Markers = naive_pubs_markers$CellType, 
                          row.names = naive_pubs_markers$Genes)
marker_anno$Markers[marker_anno$Markers == "Fb"] <- "Fibroblast"

#Fix order of naive pub markers to Fibroblast/HSC/VSMC wheen you get back
naive_pubs_markers <- naive_pubs_markers[order(naive_pubs_markers$CellType),]

heatmap_naive <- naive_pubs_markers$Genes[naive_pubs_markers$Genes %in% naive_alone_cluster_markers$gene]
heatmaps_order_idx <- naive_alone@meta.data$CellType[order(naive_alone@meta.data$CellType, decreasing = F)] #This orders HSC, Fibroblast, VSMC
heatmaps_data <- naive_alone@assays$SCT@scale.data

#Remove Rgs5 as per Mei's request
heatmap_naive <- heatmap_naive[heatmap_naive != "Rgs5"]

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Naive_Mei_NoClust.pdf", width = 5, height = 4.5)
pheatmap(heatmaps_data[match(heatmap_naive, rownames(heatmaps_data)),heatmaps_order_idx], 
         main = "Heatmap 1 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = F, cluster_rows = F,
         show_colnames = F, clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,annotation_colors = annot_colors, angle_col = 45, fontsize = 6,
         annotation_row = marker_anno)
dev.off()




names(heatmap_naive)





write.csv(naive_alone_cluster_markers, "../project/ReRun_AfterMeetingWithMei/naive_alone_allmarkers_1305.csv")
naive_alone_fibro_markers <- FindMarkers(naive_alone, ident.1 = c(0,1), ident.2 = c(2,3,4,5,6))
write.csv(naive_alone_fibro_markers, "../project/ReRun_AfterMeetingWithMei/naive_alone_fibromarkers_1305.csv")
naive_alone_HSC_markers <- FindMarkers(naive_alone, ident.1 = c(5,2), ident.2 = c(0,3,4,1,6))
write.csv(naive_alone_HSC_markers, "../project/ReRun_AfterMeetingWithMei/naive_alone_HSCmarkers_1305.csv")
naive_alone_VSMC_markers <- FindMarkers(naive_alone, ident.1 = c(3,4,6), ident.2 = c(0,1,2,5))
write.csv(naive_alone_VSMC_markers, "../project/ReRun_AfterMeetingWithMei/naive_alone_VSMCmarkers_1305.csv")

#loop over umap and save them
#show violins by cell type vsmc hsc fibro
naive_genes <- c("Acta2", "Pdgfra", "Pdgfrb", "Ly6a")
for(i in 1:length(naive_genes)) {
  FeaturePlot(naive_alone, features = naive_genes[i])
  ggsave(paste0("../project/IgGaCSF1R_Final/Naive_", naive_genes[i], "_UMAP.png"), dpi = 1000, width = 7, height = 7, units = "in")
}


VlnPlot(naive_alone, features = naive_genes, group.by = "CellType")
ggsave("../project/Naive_Final_Script/Naive_genes_vlnplot.png", dpi = 1000, height = 10, width = 10, units = "in")
library(BiocManager)
library(remotes)
library(devtools)


install.packages("remotes")



remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
SeuratWrappers::ExportToCellbrowser(object = cafs, dir = "../project/IgGaCSF1R_Final/", 
                                    dataset.name = "IgG_aCSF1R_CAFS_0511", reductions = "umap")

ExportToCellbrowser

remotes::install_version("Seurat", version = "3.2.3")
library(Seurat)



naive_heatmap <- read_csv("./PublicationPlots/Naive_markers_heatmap.csv", col_names = F)
naive_heatmap <- naive_heatmap$X1

naive_heatmap <- naive_heatmap[naive_heatmap %in% naive_alone_cluster_markers$gene]


heatmap_data <- naive_alone@assays$SCT@scale.data
dim(heatmap_data)

heatmap_data <- heatmap_data[,order(naive_alone@meta.data$CellType)]

heatmaps_data_anno <- data.frame(as.factor(naive_alone@meta.data$CellType[order(naive_alone@meta.data$CellType)]), 
                                 row.names = colnames(heatmap_data))

colnames(heatmaps_data_anno) <- c("CellType")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

default_gg_cols <- gg_color_hue(4)

annot_colors=list(CellType=c(Fibroblast=default_gg_cols[1],HSC=default_gg_cols[2], 
                  VSMC=default_gg_cols[3]))

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_Naive_Mei_NoClust.pdf", width = 5, height = 4)
pheatmap(heatmap_data[rownames(heatmap_data) %in% naive_heatmap,], 
         main = "Heatmap 2 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = F, cluster_rows = T,
         show_colnames = F, clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,
         #annotation_colors = annot_colors, 
         angle_col = 45, fontsize_col = 6)
dev.off()


pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_Naive_Mei_Clust.pdf", width = 5, height = 4)
pheatmap(heatmap_data[rownames(heatmap_data) %in% naive_heatmap,], 
         main = "Heatmap 2 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = T, cluster_rows = T,
         show_colnames = F, clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,
         #annotation_colors = annot_colors, 
         angle_col = 45, fontsize_col = 6)
dev.off()


### Naive Average

heatmap1_avg <- data.frame(matrix(nrow = length(naive_heatmap), ncol = 3))
colnames(heatmap1_avg) <- c("Fibroblast", "VSMC", "HSC")

for(i in 1:nrow(heatmap1_avg)) {
  
  heatmap1_avg$Fibroblast[i] <- median(naive_alone@assays$SCT@scale.data[rownames(naive_alone@assays$SCT@scale.data) == naive_heatmap[i],naive_alone@meta.data$CellType == "Fibroblast"])
  heatmap1_avg$VSMC[i] <- median(naive_alone@assays$SCT@scale.data[rownames(naive_alone@assays$SCT@scale.data) == naive_heatmap[i],naive_alone@meta.data$CellType == "VSMC"])
  heatmap1_avg$HSC[i] <- median(naive_alone@assays$SCT@scale.data[rownames(naive_alone@assays$SCT@scale.data) == naive_heatmap[i],naive_alone@meta.data$CellType == "HSC"])
  rownames(heatmap1_avg)[i] <- naive_heatmap[i] 
  
}


pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_Naive_Mei_Clust_Average.pdf", width = 5, height = 4)
pheatmap(t(heatmap1_avg), 
         main = "Heatmap Naive Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "column", cluster_cols = T, cluster_rows = F,
         show_colnames = T, clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_colors = annot_colors, angle_col = 45, fontsize_col = 6, cellheight = 15)
dev.off()




Idents(naive_alone) <- naive_alone@meta.data$CellType
DimPlot(naive_alone)
de_naive_celltype <- FindAllMarkers(naive_alone)




p_annotop10 <- data.frame(CellType = naive_alone@meta.data$CellType[order(naive_alone@meta.data$CellType)], row.names = rownames(naive_alone@meta.data)[order(naive_alone@meta.data$CellType)])

rna_scale2 <- ScaleData(naive_alone, features = rownames(naive_alone))
rna_scale2 <- GetAssayData(rna_scale2, slot = "scale.data")
dim(rna_scale2)
rna_scale2 <- rna_scale2 %>% as.data.frame() %>% rownames_to_column(var = "gene")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

default_gg_cols <- gg_color_hue(3)



naive_de_top10_bycluster <- de_naive_celltype %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10genes_ordered <- data.frame(Genes = naive_de_top10_bycluster$gene)
top10genes_ordered <- top10genes_ordered %>% left_join(rna_scale2, by = c("Genes" = "gene"))
top10genes_ordered <- top10genes_ordered %>% column_to_rownames(var = "Genes")
top10genes_ordered <- top10genes_ordered[,order(naive_alone@meta.data$CellType)]

naive_cols <- default_gg_cols
names(naive_cols) <- c("Fibroblast", "HSC", "VSMC")
top10_colours <- list(CellType = naive_cols)

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Top10GenesLogFC_NaiveCellTypes_FixedColoursJuly.pdf", width = 5, height = 4)
pheatmap(top10genes_ordered, 
         scale = "row", cluster_rows = F, cluster_cols = F, 
         show_colnames = F, breaks = seq(-1.5,1.5,length.out = 101), 
         annotation = p_annotop10, clustering_method = "ward.D2", 
         fontsize_row = 8, annotation_colors = top10_colours, gaps_col = c(1398, (1398+475)), gaps_row = c(10,20,30))
dev.off()




naive_de_top20_bycluster <- de_naive_celltype %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20genes_ordered <- data.frame(Genes = naive_de_top20_bycluster$gene)
top20genes_ordered <- top20genes_ordered %>% left_join(rna_scale2, by = c("Genes" = "gene"))
top20genes_ordered <- top20genes_ordered %>% column_to_rownames(var = "Genes")
top20genes_ordered <- top20genes_ordered[,order(naive_alone@meta.data$CellType)]



pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Top20GenesLogFC_NaiveCellTypes_FixedColoursJuly.pdf", width = 5, height = 4)
pheatmap(top20genes_ordered, 
         scale = "row", cluster_rows = F, cluster_cols = F, 
         show_colnames = F, breaks = seq(-1.5,1.5,length.out = 101), 
         annotation = p_annotop10, clustering_method = "ward.D2", 
         fontsize_row = 8, annotation_colors = top10_colours, 
         gaps_col = c(1398, (1398+475)), gaps_row = c(20,40,60), fontsize = 6)
dev.off()

write.csv(top20genes_ordered, "~/cellranger_out_temp/PublicationPlots/Top20GenesLogFC_Naive_GeneList.csv")








