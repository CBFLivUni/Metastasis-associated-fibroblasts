library(tidyverse)
library(sctransform)
library(Seurat)

setwd("~/cellranger_out_temp/")

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

#all_dat <- merge(x = sc_obj2[["Naive"]], y = c(sc_obj2[["D14"]], sc_obj2[["IgG"]], sc_obj2[["aCSF1R"]]),
#                add.cell.ids = c("Naive", "D14", "IgG", "aCSF1R"), project = "schmid")
all_dat <- merge(x = sc_obj2[["Naive"]], y = c(sc_obj2[["IgG"]], sc_obj2[["aCSF1R"]]),
                 add.cell.ids = c("Naive", "IgG", "aCSF1R"), project = "schmid")

all_dat[["percent.mt"]] <- PercentageFeatureSet(all_dat, pattern = "^mt-")
#"nFeature_RNA", "nCount_RNA",

#VlnPlot(all_dat, features = c("nFeature_RNA"), ncol = 1) #+ ylim(1000,6500)

VlnPlot(all_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

all_dat <- subset(all_dat, subset = nFeature_RNA < 6500 & nFeature_RNA > 500 & percent.mt < 10)

VlnPlot(all_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#VlnPlot(all_dat, features = "nCount_RNA", ncol = 1) + ylim(500,50000)


all_dat <- SCTransform(all_dat)

all_dat <- RunPCA(all_dat) %>% 
  RunUMAP(dims = 1:50, n.neighbors = 30) %>%
  FindNeighbors(dims = 1:50) %>% 
  FindClusters(resolution = 0.8)
ElbowPlot(all_dat, ndims = 50)

DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident")
DimPlot(all_dat, pt.size = 0.4, label = TRUE, label.size = 6)

all_dat_umap_gene <- c("Acta2", "Pdgfra", "Pdgfrb", "Ly6a", "Col1a1", "Col1a2", "Lrat", "Myh11")
#for loop over umaps to save output

#for(i in 1:length(all_dat_umap_gene)) {
#FeaturePlot(all_dat, features = all_dat_umap_gene[i])
#  ggsave(paste0("../project/IgGaCSF1R_Final/Global_", all_dat_umap_gene[i], "_UMAP.png"), dpi = 1000, width = 7, height = 7, units = "in")
#}

############# REMOVE OUTLIER CLUSTERS ##############


#all_dat <- all_dat[,!all_dat$seurat_clusters %in% c(11:14)]
#all_dat <- all_dat[,!all_dat$seurat_clusters %in% c(16,17,19)]
all_dat <- all_dat[,!all_dat$seurat_clusters %in% c(13,15,16,17)]
all_dat <- all_dat[,all_dat@reductions[["umap"]]@cell.embeddings[,2] < 9 & all_dat@reductions[["umap"]]@cell.embeddings[,1] > -10]

DimPlot(all_dat, group.by = "orig.ident", pt.size = 1)
DimPlot(all_dat, label = T, label.size = 6)
FeaturePlot(all_dat, features = "Runx1", pt.size = .3)
DimPlot(all_dat[,all_dat$SCT_snn_res.1.2 %in% c(17)], pt.size = 2)
DimPlot(all_dat[,all_dat$SCT_snn_res.1.2 %in% c(6)], pt.size = 2, group.by = "orig.ident")
DimPlot(all_dat, pt.size = 2, label = T)


############# RECLUSTER ###############

###This is now SCTransform all_dat without D14 from the start

all_dat <- RunPCA(all_dat) %>% 
  RunUMAP(dims = 1:30, n.neighbors = 10) %>%
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 1.2)

DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident")
DimPlot(all_dat, pt.size = 0.4, label = TRUE, label.size = 6) + 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/GlobalUMAP_CAFSelection.pdf", dpi = 1000, width = 5.5, height = 4.5, units = "in")



DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident") + 
  scale_colour_manual(values = c("red3","red3", alpha("grey",0.3))) + 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/TumourBearingUMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

default_gg_cols <- gg_color_hue(4)


DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident") + 
  scale_colour_manual(values = c("red3","red3", alpha("grey",0.3))) + 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/TumourBearingUMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")

DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident") + 
  scale_colour_manual(values = c(alpha("lightgrey",0.3),alpha("lightgrey",0.3),"#332288")) + 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")

p1 <- p1$data
levels(p1$orig.ident)
ggsave("~/cellranger_out_temp/PublicationPlots/NonTumourBearingUMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")

all_dat$CellType_Label <- case_when(all_dat$SCT_snn_res.1.2 %in% c(0,2,5,8,17) ~ "MAFs",
                                    all_dat$SCT_snn_res.1.2 %in% c(3,14,18) ~ "SmoothMuscle",
                                    all_dat$SCT_snn_res.1.2 %in% c(11,1,6,15,12,9,7,4,13,10,16) ~ "Naive")

DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident") + 
  scale_colour_manual(values = c(alpha("lightgrey",0.3), "red3",alpha("lightgrey",0.3))) + 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/IgGRedUMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")

DimPlot(all_dat, pt.size = 0.4, label = F, group.by = "orig.ident") + 
  scale_colour_manual(values = c("seagreen4",alpha("lightgrey",0.3),alpha("lightgrey",0.3))) + 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/aCSF1RGreenUMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")


table(all_dat[,all_dat@meta.data$orig.ident == "aCSF1R"]$CellType_Label)
table(all_dat[,all_dat@meta.data$orig.ident == "IgG"]$CellType_Label) 
all_dat$MAFS <- case_when(all_dat$SCT_snn_res.1.2 %in% c(0,2,5,8,17) ~ "MAFs",
                                    all_dat$SCT_snn_res.1.2 %in% c(3,14,18) ~ "NotMAFs",
                                    all_dat$SCT_snn_res.1.2 %in% c(11,1,15,12,9,7,6,4,13,10,16) ~ "NotMAFs")

maf_notmaf_markers <- FindMarkers(all_dat, ident.1 = c(0,2,5,8,17), ident.2 = c(3,14,18,   11,1,6,15,12,9,7,4,13,10,16))
write.csv(maf_notmaf_markers, "./PublicationPlots/MAF_NOTMAF_Markers.csv")

FeaturePlot(all_dat, features = c("Ptprc"))

DimPlot(all_dat, pt.size = 1, label = F, group.by = "MAFS")

GlobalCAFS_vs_Everything <- FindMarkers(all_dat, ident.1 = c(0,2,5,8,17), ident.2 = c(3,14,18,   11,1,15,12,9,7,4,13,10,16))
GlobalSmoothMuscle_vs_Everything <- FindMarkers(all_dat, ident.1 = c(3,14,18), ident.2 = c(0,2,5,8,17,   11,1,15,12,9,7,4,13,10,16))
GlobalNaive_vs_Everything <- FindMarkers(all_dat, ident.1 = c(11,1,15,12,9,7,4,13,10,16), ident.2 = c(0,2,5,6,8,17,   11,1,15,12,9,7,4,13,10,16))

write.csv(GlobalCAFS_vs_Everything, "../project/IgGaCSF1R_Final/GlobalCAFS_vs_Everything.csv")
write.csv(GlobalSmoothMuscle_vs_Everything, "../project/IgGaCSF1R_Final/GlobalSmoothMuscle_vs_Everything.csv")
write.csv(GlobalNaive_vs_Everything, "../project/IgGaCSF1R_Final/GlobalNaive_vs_Everything.csv")


DimPlot(all_dat, pt.size = 1, label = F, group.by = "CellType_Label")

DimPlot(all_dat, pt.size = 1, label = F, group.by = "orig.ident")
ggsave("../project/IgGaCSF1R_Final/GlobalClustering.png", dpi = 1000, height = 10, width = 10, units = "in")
DimPlot(all_dat, pt.size = 1, label = TRUE, label.size = 6)
ggsave("../project/IgGaCSF1R_Final/GlobalClustering_w_labels.png", dpi = 1000, height = 10, width = 10, units = "in")
library(RColorBrewer)
brewcols <- brewer.pal(n = 11, name = "RdYlBu")
RdYlBu

FeaturePlot(all_dat, features = c("Ptprc"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Ptprc.pdf",  dpi = 1000, width = 5, height = 4, units = "in")
FeaturePlot(all_dat, features = c("Epcam"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Epcam.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = "Pecam1", cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Pecam1.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Col1a1"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Col1a1.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Col1a2"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Col1a2.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Col3a1"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Col3a1.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Acta2"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Acta2.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Spp1"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Spp1.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Hgf"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Hgf.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Pdgfra"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/PDGFRa.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Pdgfra"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/PDGFRa.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Cd34"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Cd34.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Ly6a"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Ly6a.pdf",  dpi = 1000, width = 5, height = 4, units = "in")

FeaturePlot(all_dat, features = c("Myh11"), cols = rev(brewcols)) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/Myh11.pdf",  dpi = 1000, width = 5, height = 4, units = "in")


cells_in_global <- table(all_dat$SCT_snn_res.1.2) %>% as.data.frame()
colnames(cells_in_global) <- c("Cluster", "Count")
write.csv(cells_in_global, "~/cellranger_out_temp/PublicationPlots/CellsInGlobal.csv")

Acta2
Spp1
Hgf
PDGFRa
CD34
LY6A
Myh11 




VlnPlot(all_dat, features = c("Col1a1", "Acta2"))

FeaturePlot(all_dat, pt.size = 1, features = c("Acta2"), cols = c("white", "darkred"))


caf_landscape <- c("Acta2", "Pdgfra", "Pdgfrb", "Ly6a", "Postn", "Tnc", "Mif", "Spp1", "Col1a1", "Col1a2", "Has2", "Hgf", "Lrat", "Myh11", "Gsn", "Cxcl12", "Gas6", "Fbln2")


for(i in 1:length(caf_landscape)) {
print(FeaturePlot(all_dat, pt.size = 1, features = caf_landscape[i], cols = c("white", "darkred")))
}

## Extract CAFs

#cafs <- all_dat[,all_dat$SCT_snn_res.1 %in% c(0,2,5,14,8)]

cafs <- all_dat[,all_dat$SCT_snn_res.1.2 %in% c(0,2,5,8,17)]
DimPlot(cafs, label = T)
DimPlot(cafs, pt.size = 1, label = F, group.by = "orig.ident")

#Check with Mei incase he wants to include 6
cafs <- all_dat[,all_dat$SCT_snn_res.1.2 %in% c(0,2,5,6,8,17)]
DimPlot(cafs, label = T)
DimPlot(cafs, pt.size = 1, label = F, group.by = "orig.ident")







#### CAFS & IGG ALONE

### For this all_dat needs to be filtered with less stringent cutoffs to retain as much aCSF1R & IgG as possible

cafs <- cafs[,cafs$orig.ident != "Naive"]

cafs <- RunPCA(cafs) %>% 
  RunUMAP(dims = 1:20, n.neighbors = 50) %>%
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.3)

DimPlot(cafs, pt.size = 0.8, label = F) + xlab("UMAP 1") + ylab("UMAP 2")
table(cafs@meta.data$orig.ident, cafs@meta.data$SCT_snn_res.0.3)
  #scale_colour_manual(values = c("red3","seagreen4","#332288","purple3"))
ggsave("~/cellranger_out_temp/PublicationPlots/AllCAFSUMAP_originalcolours.pdf", dpi = 1000, width = 5, height = 4, units = "in")

DimPlot(cafs, pt.size = 0.8, label = F, group.by = "orig.ident") + 
  scale_colour_manual(values = c(alpha("grey",0.3), "#332288")) + 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/IgGCAFSRedUMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")


DimPlot(cafs, pt.size = 0.8, label = F, group.by = "orig.ident") + 
  scale_colour_manual(values = c("red3",alpha("grey",0.2)))+ 
  ggtitle("") + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/aCSF1RCAFSUMAP.pdf", dpi = 1000, width = 5, height = 4, units = "in")

Idents(cafs) <- cafs@meta.data$SCT_snn_res.0.3
DimPlot(cafs[,cafs@meta.data$orig.ident == "IgG"], pt.size = 0.8) + xlim(-6.5,4) + ylim(-7,4.2) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/IgG_only_UMAPCAFS.pdf", dpi = 1000, width = 5, height = 4, units = "in")

DimPlot(cafs[,cafs@meta.data$orig.ident == "aCSF1R"], pt.size = 0.8) + xlim(-6.5,4) + ylim(-7,4.2) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("~/cellranger_out_temp/PublicationPlots/aCSF1R_only_UMAPCAFS.pdf", dpi = 1000, width = 5, height = 4, units = "in")



meis_heatmap <- c("Lpl", "Sparcl1", "Cepbd", "Itga8", "Lrat", "Hgf", "Socs3", "Cthrc1", "Spp1", 
                  "Acta2", "Tnc", "Col1a1", "Col1a2", "Vcan", "Postn", "Has2", 
                  "Gsn", "C3", "Gas6", "Ly6a", "Pi16", "Cxcl12", "Ube2c", "Hmgb2", "Mki67", "Cdk1")

library(pheatmap)

p_anno <- data.frame(periMAF = if_else(cafs@meta.data$SCT_snn_res.0.3 == 0, 1, 0),
                     myMAF = if_else(cafs@meta.data$SCT_snn_res.0.3 == 1, 1, 0),
                     iMAF = if_else(cafs@meta.data$SCT_snn_res.0.3 == 2, 1, 0),
                     cycMAF = if_else(cafs@meta.data$SCT_snn_res.0.3 == 3, 1, 0), row.names = colnames(cafs@assays$SCT@scale.data))




marker_anno <- data.frame(Markers = rep("", length(rownames(meis_heatmap_df))), row.names = rownames(meis_heatmap_df))
marker_anno[rownames(marker_anno) %in% c("Sparcl1", "Lpl", "Lrat", "Hgf", "Itga8", "Socs3"),] <- c("periMAF")
marker_anno[rownames(marker_anno) %in% c("Pi16", "Ly6a", "C3", "Gsn", "Cxcl12", "Gas6"),] <- c("iMAF")
marker_anno[rownames(marker_anno) %in% c("Mki67", "Hmgb2", "Ube2c", "Cdk1"),] <- c("cycMAF")
marker_anno[rownames(marker_anno) %in% c("Col1a2", "Col1a1", "Spp1", "Has2", "Vcan", "Tnc", "Cthrc1", "Postn", "Acta2"),] <- c("myMAF")

cafs$MAF_Type <- case_when(cafs$SCT_snn_res.0.3 == 0 ~ "periMAF",
                           cafs$SCT_snn_res.0.3 == 1 ~ "myMAF",
                           cafs$SCT_snn_res.0.3 == 2 ~ "iMAF",
                           cafs$SCT_snn_res.0.3 == 3 ~ "cycMAF")

p_anno <- data.frame(MAF = as.factor(cafs$MAF_Type), row.names = colnames(cafs@assays$SCT@scale.data))

annot_colors=list(MAF=c(periMAF=default_gg_cols[1],myMAF=default_gg_cols[2], 
                        iMAF=default_gg_cols[3],cycMAF=default_gg_cols[4]),
                  Markers = c(periMAF=default_gg_cols[1],myMAF=default_gg_cols[2], 
                              iMAF=default_gg_cols[3],cycMAF=default_gg_cols[4]))

pheatmap(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) %in% meis_heatmap,order(cafs@meta.data$MAF_Type)], 
         main = "MAF Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = F,
         show_colnames = F, clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean", annotation = p_anno, 
         annotation_row = marker_anno,cutree_cols = 4, cutree_rows = 4,annotation_colors = annot_colors)






meis_heatmap_df <- cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) %in% meis_heatmap,]
avg_scaled_data <- data.frame(matrix(nrow = nrow(meis_heatmap_df), ncol = 4))
colnames(avg_scaled_data) <- c("periMAF", "myMAF", "iMAF", "cycMAF")
mean_i <- c()
sd_i <- c()
for(i in 1:nrow(meis_heatmap_df)) {
rownames(avg_scaled_data)[i] <- rownames(meis_heatmap_df)[i]  

mean_i <- c(mean_i, mean(meis_heatmap_df[i,]))
names(mean_i) <- rownames(meis_heatmap_df)[i]

sd_i <- c(sd_i, sd(meis_heatmap_df[i,]))
names(sd_i) <- rownames(meis_heatmap_df)[i]
  
# avg_scaled_data$periMAF[i] <- (mean(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "periMAF"]) - mean_i[i]) / sd_i[i]
# avg_scaled_data$myMAF[i] <- (mean(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "myMAF"]) - mean_i[i]) / sd_i[i]
# avg_scaled_data$iMAF[i] <- (mean(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "iMAF"]) - mean_i[i]) / sd_i[i]
# avg_scaled_data$cycMAF[i] <- (mean(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "cycMAF"]) - mean_i[i]) / sd_i[i]

avg_scaled_data$periMAF[i] <- median(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "periMAF"])
avg_scaled_data$myMAF[i] <- median(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "myMAF"])
avg_scaled_data$iMAF[i] <- median(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "iMAF"]) 
avg_scaled_data$cycMAF[i] <- median(meis_heatmap_df[i,cafs@meta.data$MAF_Type == "cycMAF"]) 

}



avg_scaled_data <- avg_scaled_data[order(marker_anno$Markers),]
avg_scaled_data <- avg_scaled_data %>% rownames_to_column(var = "Gene")
marker_anno_join <- marker_anno %>% rownames_to_column(var = "Gene")

avg_scaled_data <- avg_scaled_data %>% left_join(marker_anno_join)
avg_scaled_data <- avg_scaled_data %>% arrange(desc(Markers))
marker_anno_join <- marker_anno_join[match(avg_scaled_data$Gene, marker_anno_join$Gene),]
rownames(marker_anno_join) <- marker_anno_join$Gene
marker_anno_join$Gene <- NULL

avg_scaled_data <- avg_scaled_data %>% column_to_rownames(var = "Gene")

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/HeatmapOrdered.pdf", width = 5, height = 4)
pheatmap(t(avg_scaled_data[,1:4]), 
         main = "MAF Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "column", cluster_cols = F, cluster_rows = F,
         show_colnames = T, clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = marker_anno_join,annotation_colors = annot_colors, angle_col = 45, fontsize_col = 6, cellheight = 15)
dev.off()


pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/HeatmapClustered.pdf", width = 5, height = 4)
pheatmap(t(avg_scaled_data[,1:4]), 
         main = "MAF Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "column", cluster_cols = T, cluster_rows = F,
         show_colnames = T, clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = marker_anno_join,annotation_colors = annot_colors, angle_col = 45, fontsize_col = 6, cellheight = 15)
dev.off()



cafs_de_nb <- FindAllMarkers(cafs, test.use = "negbinom")
cafs_de <- FindAllMarkers(cafs)
cafs_de_auc <- FindAllMarkers(cafs, test.use = "roc")
write.csv(cafs_de_auc, "./PublicationPlots/CAFS_perclusterAUC.csv")

cafs_de_auc <- cafs_de_auc %>% arrange(cluster, desc(myAUC))
experimental_celltype <- colnames(cafs)
experimental_celltype <- gsub("_.*", "", experimental_celltype)
unique(experimental_celltype)
Idents(cafs) <- experimental_celltype
cafs_de_acsfvsigg <- FindMarkers(cafs[,cafs@meta.data$SCT_snn_res.0.3 %in% c("0","1")], ident.1 = "aCSF1R", ident.2 = "IgG", test.use = "roc")

(DimPlot(cafs[,cafs@meta.data$SCT_snn_res.0.3 %in% c("0","1")])+ xlab("UMAP 1") + ylab("UMAP 2")) + 
(FeaturePlot(cafs[,cafs@meta.data$SCT_snn_res.0.3 %in% c("0","1")], features = "Uba52", cols = rev(brewcols), split.by = "ident") + xlab("UMAP 1") + ylab("UMAP 2")) +
  plot_layout(widths=c(1,2))




cafs_de_auc_filt <- cafs_de_auc[cafs_de_auc$myAUC > 0.7,]
cafs_de_auc_filt <- split.data.frame(cafs_de_auc_filt, f = cafs_de_auc_filt$cluster)

FeaturePlot(cafs[,cafs@meta.data$SCT_snn_res.0.3 %in% c("0","1")], features = c("Nme2"))

Idents(cafs) <- cafs@meta.data$SCT_snn_res.0.3
cell_type_cafs <- gsub("_.*", "", colnames(cafs))
Idents(cafs) <- cell_type_cafs

cafs_celltypemarkers <- FindAllMarkers(cafs[,cafs@meta.data$SCT_snn_res.0.3 %in% c("0")], test.use = "roc")
cafs_celltypemarkers <- cafs_celltypemarkers %>% arrange(cluster, desc(myAUC))
DimPlot(cafs)

heatmap1 <- read_csv("./PublicationPlots/Heatmap 1_Mei_signature.csv", col_names = T)
heatmap1_annot <- heatmap1
heatmap1 <- heatmap1$`Heatmap 1 Selected markers_Tuveson`
cafs_de_filt <- cafs_de[abs(cafs_de$avg_log2FC) > 0.5,]
heatmap1 <- heatmap1[heatmap1 %in% cafs_de_filt$gene]
heatmaps_data <- cafs@assays$SCT@scale.data
heatmaps_data <- ScaleData(cafs, features = rownames(cafs))

heatmaps_order_idx <- order(cafs@meta.data$MAF_Type, decreasing = T) #This orders periMAF, iMAF, myMAF, cycMAF
cafs@meta.data$MAF_Type[heatmaps_order_idx]

#heatmaps_data <- heatmaps_data[,order(cafs@meta.data$MAF_Type)]
default_gg_cols

### Heatmap 1
cafs@meta.data$MAF_Type
heatmaps_data_anno <- data.frame(as.factor(cafs@meta.data$MAF_Type[heatmaps_order_idx]), 
                                 row.names = colnames(heatmaps_data[,heatmaps_order_idx]))
table(heatmaps_data_anno$as.factor.cafs.meta.data.MAF_Type.heatmaps_order_idx..)
heatmap1_annot <- data.frame(Markers1 = heatmap1_annot$`Signature ID`, row.names =heatmap1_annot$`Heatmap 1 Selected markers_Tuveson`)

annot_colors$Markers1 <- annot_colors$MAF[1:3]
names(annot_colors$Markers1) <- c("myCAF", "iCAF", "apCAF")
colnames(heatmaps_data_anno) <- "MAF"
class(heatmaps_data_anno$MAF)
heatmaps_data_anno$MAF <- factor(heatmaps_data_anno$MAF, levels = c("periMAF", "myMAF", "iMAF", "cycMAF"))

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_1_Mei_NoClust.pdf", width = 5, height = 4.5)
pheatmap(heatmaps_data[match(heatmap1, rownames(heatmaps_data)),heatmaps_order_idx], 
         main = "Heatmap 1 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = F, cluster_rows = F,
         show_colnames = F, clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,annotation_colors = annot_colors, angle_col = 45, fontsize = 6,
         annotation_row = heatmap1_annot)
dev.off()

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_1_Mei_Clust.pdf", width = 5, height = 4)
pheatmap(heatmaps_data[rownames(heatmaps_data) %in% heatmap1,], 
         main = "Heatmap 1 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = T, cluster_rows = T,
         show_colnames = F, clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,annotation_colors = annot_colors, angle_col = 45, fontsize_col = 6)
dev.off()

### Heatmap 1 Avg

heatmap1_avg <- data.frame(matrix(nrow = length(heatmap1), ncol = 4))
colnames(heatmap1_avg) <- c("periMAF", "myMAF", "iMAF", "cycMAF")

for(i in 1:nrow(heatmap1_avg)) {
  
  heatmap1_avg$periMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap1[i],cafs@meta.data$MAF_Type == "periMAF"])
  heatmap1_avg$myMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap1[i],cafs@meta.data$MAF_Type == "myMAF"])
  heatmap1_avg$iMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap1[i],cafs@meta.data$MAF_Type == "iMAF"])
  heatmap1_avg$cycMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap1[i],cafs@meta.data$MAF_Type == "cycMAF"])
  
  rownames(heatmap1_avg)[i] <- heatmap1[i] 
  
}
match(colnames(cafs@assays$SCT@scale.data), rownames(cafs@meta.data))


pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_1_Mei_NoClust_Average.pdf", width = 5, height = 4)
pheatmap(t(heatmap1_avg), 
         main = "Heatmap 1 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "column", cluster_cols = F, cluster_rows = F,
         show_colnames = T, clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_colors = annot_colors, angle_col = 45, fontsize = 5, cellheight = 15,
         annotation_col = heatmap1_annot)
dev.off()


### Heatmap 2


heatmap2 <- read_csv("./PublicationPlots/Heatmap 2_Mei_signature.csv", col_names = T)
heatmap2_annot <- heatmap2
heatmap2 <- heatmap2$`Heatmap 2 Selected markers_Schwabbe`
heatmap2 <- heatmap2[heatmap2 %in% cafs_de_filt$gene]
heatmap2_annot <- data.frame(Markers2 = heatmap2_annot$`Signature ID`, row.names =heatmap2_annot$`Heatmap 2 Selected markers_Schwabbe`)
annot_colors$Markers2 <- annot_colors$MAF[1:3]
names(annot_colors$Markers2) <- c("myCAF", "iCAF", "mesCAF")
unique(heatmap2_annot$Markers2)
pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_2_Mei_NoClust.pdf", width = 5, height = 4)
pheatmap(heatmaps_data[rownames(heatmaps_data) %in% heatmap2,], 
         main = "Heatmap 2 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = F, cluster_rows = F,
         show_colnames = F, clustering_method = "ward.D2", 
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,annotation_colors = annot_colors, angle_col = 45, fontsize_col = 6,
         annotation_row = heatmap2_annot)
dev.off()

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_2_Mei_Clust.pdf", width = 5, height = 4.5)
pheatmap(heatmaps_data[match(heatmap2, rownames(heatmaps_data), nomatch = 0),heatmaps_order_idx], 
         main = "Heatmap 2 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = F, cluster_rows = F,
         show_colnames = F, clustering_method = "ward.D2", show_rownames = F,
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,annotation_colors = annot_colors, annotation_row = heatmap2_annot, angle_col = 45, fontsize = 6)
dev.off()

png(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_2_Mei_NoClust.png", width = 5, height = 4.5, units = "in", res = 1000)
pheatmap(heatmaps_data[match(heatmap2, rownames(heatmaps_data), nomatch = 0),heatmaps_order_idx], 
         main = "Heatmap 2 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "row", cluster_cols = F, cluster_rows = F,
         show_colnames = F, clustering_method = "ward.D2", show_rownames = F,
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_col = heatmaps_data_anno,annotation_colors = annot_colors, annotation_row = heatmap2_annot, angle_col = 45, fontsize = 6)
dev.off()







### Heatmap 2 avg

heatmap2_avg <- data.frame(matrix(nrow = length(heatmap2), ncol = 4))
colnames(heatmap2_avg) <- c("periMAF", "myMAF", "iMAF", "cycMAF")

for(i in 1:nrow(heatmap2_avg)) {
  
  heatmap2_avg$periMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap2[i],cafs@meta.data$MAF_Type == "periMAF"], na.rm = T)
  heatmap2_avg$myMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap2[i],cafs@meta.data$MAF_Type == "myMAF"],na.rm = T)
  heatmap2_avg$iMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap2[i],cafs@meta.data$MAF_Type == "iMAF"],na.rm =T)
  heatmap2_avg$cycMAF[i] <- median(cafs@assays$SCT@scale.data[rownames(cafs@assays$SCT@scale.data) == heatmap2[i],cafs@meta.data$MAF_Type == "cycMAF"],na.rm=T)
  
  rownames(heatmap2_avg)[i] <- heatmap2[i] 
  
}
pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_2_Mei_Clust_Average.pdf", width = 5, height = 4)
pheatmap(t(heatmap2_avg[!is.na(rowSums(heatmap2_avg)),]), 
         main = "Heatmap 2 Markers", 
         breaks = seq(-1.5,1.5, length.out = 101), scale = "column", cluster_cols = F, cluster_rows = F,
         show_colnames = F, clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",
         annotation_colors = annot_colors, angle_col = 45, fontsize = 5, cellheight = 15,
         annotation_col = heatmap2_annot)
dev.off()










# optional colours - , color=colorRampPalette(c("darkgreen", "white", "darkred"))(100)
sum(cafs@meta.data$orig.ident == "aCSF1R")

DimPlot(cafs, label = T, pt.size = 2)
ggsave("../project/IgGaCSF1R_Final/IgG_aCSF1R_Cafs_Clusters.png", dpi = 1000, height = 10, width = 10, units = "in")
DimPlot(cafs, group.by = "orig.ident", pt.size = 2)
ggsave("../project/IgGaCSF1R_Final/IgG_aCSF1R_Cafs_OrigIdent.png", dpi = 1000, height = 10, width = 10, units = "in")
table(cafs$SCT_snn_res.0.3, cafs$orig.ident)



DimPlot(cafs, label = T)
DimPlot(cafs, group.by = "orig.ident")

cafs@meta.data %>% group_by(orig.ident, SCT_snn_res.0.3) %>% summarise(n = n())

cafs_de <- FindAllMarkers(cafs)
cafs_de_rnaassay <- FindAllMarkers(cafs, assay = "RNA", test.use = "DESeq2")
cafs_de_cluster1markers <- FindMarkers(cafs, ident.1 = 1, ident.2 = 0)
write.csv(cafs_de, "../project/IgGaCSF1R_Final/aCSF1R_IgG_cafs_all_markers_by_cluster.csv")

cafs@meta.data$celltype_cluster <- paste0(cafs@meta.data$orig.ident, "_", cafs@meta.data$SCT_snn_res.0.3)
cafs$celltype_cluster <- factor(cafs$celltype_cluster, levels = c("IgG_0","aCSF1R_0", "IgG_1", "aCSF1R_1", "IgG_2", "aCSF1R_2", "IgG_3", "aCSF1R_3"))

VlnPlot(cafs, features = c("Hsp90b1", "Sparcl1", "Hmgb1", "Rps21", "Erh", "Cbx3"), group.by = "celltype_cluster", ncol = 3)
VlnPlot(cafs, features = c("Hsp90b1", "Sparcl1", "Hmgb1", "Rps21", "Erh", "Cbx3"))



unique(cafs@meta.data$celltype_cluster)






#plots for Mei's talk

caf_landscape <- c("Acta2", "Pdgfra", "Pdgfrb", "Ly6a", "Postn", "Tnc", "Mif", "Spp1", "Col1a1", "Col1a2", "Has2", "Hgf", "Lrat", "Myh11", "Gsn", "Cxcl12", "Gas6", "Fbln2")


for(i in 1:length(caf_landscape)) {
  FeaturePlot(cafs, features = caf_landscape[i])
  ggsave(paste0("../project/IgGaCSF1R_Final/CAFS_", caf_landscape[i], "_UMAP.png"), dpi = 1000, width = 7, height = 7, units = "in")
}



VlnPlot(cafs, features = caf_landscape[1:6])
ggsave("../project/IgGaCSF1R_Final/caf_landscape_markers_1.png", dpi = 1000, height = 10, width = 10, units = "in")
VlnPlot(cafs, features = caf_landscape[7:12])
ggsave("../project/IgGaCSF1R_Final/caf_landscape_markers_2.png", dpi = 1000, height = 10, width = 10, units = "in")
VlnPlot(cafs, features = caf_landscape[13:18])
ggsave("../project/IgGaCSF1R_Final/caf_landscape_markers_3.png", dpi = 1000, height = 10, width = 10, units = "in")











cafs_0_igg_vs_acsf <- FindMarkers(cafs, group.by = "celltype_cluster", ident.1 = "aCSF1R_0", ident.2 = "IgG_0")
cafs_0_igg_vs_acsf_downsample <- FindMarkers(cafs, group.by = "celltype_cluster", ident.1 = "aCSF1R_0", ident.2 = "IgG_0", max.cells.per.ident = 260)
cafs_0_igg_vs_acsf <- cafs_0_igg_vs_acsf[cafs_0_igg_vs_acsf$p_val_adj < 0.05,]
cafs_0_igg_vs_acsf_downsample <- cafs_0_igg_vs_acsf_downsample[cafs_0_igg_vs_acsf_downsample$p_val_adj < 0.05,]

cafs@meta.data %>% group_by(cafs$SCT_snn_res.0.3) %>% summarise(n = n())



cafs_1_igg_vs_acsf <- FindMarkers(cafs, group.by = "celltype_cluster", ident.1 = "aCSF1R_1", ident.2 = "IgG_1", max.cells.per.ident = 150)
cafs_1_igg_vs_acsf <- cafs_1_igg_vs_acsf[cafs_1_igg_vs_acsf$p_val_adj < 0.05,]
write.csv(cafs_1_igg_vs_acsf, "../project/IgGaCSF1R_Final/CAFS_Cluster1_aCSF1RvsIgG.csv")

VlnPlot(cafs, features = c("Fos", "Jun", "Socs3"), ncol = 3, split.by = "orig.ident")
VlnPlot(cafs, features = c("Uba52"), ncol = 1, split.by = "orig.ident")
write.csv(cafs_0_igg_vs_acsf, "../project/IgGaCSF1R_Final/CAFS_Cluster0_aCSF1RvsIgG.csv")


install.packages("enrichR")
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- ("GO_Molecular_Function_2015")
enriched <- enrichr(rownames(cafs_1_igg_vs_acsf), dbs)


write.csv(enriched[1:50,], "../project/IgGaCSF1R_Final/CAFS_Cluster0_BiologicalProcesses.csv")


write.csv(cafs_de, "../project/CAFS_noD14.csv")


caf_markers <- c("Lrat", "Spp1", "Socs3")
RidgePlot(cafs, features = caf_markers, ncol = 3)
VlnPlot(cafs, features = "Acta2", ncol = 1, assay = "SCT")
FeaturePlot(cafs, features = caf_markers)

jakstat_genes <- c("Cdkn1a", "Cebpb", "Fos", "Jun", "Rap2a", "Socs3", "Stat3")
RidgePlot(cafs, features = jakstat_genes, ncol = 3)

VlnPlot(cafs, features = "Col8a1", ncol = 1, split.by = "orig.ident")
cafs_de_list <- cafs_de %>% split.data.frame(f = cafs_de$cluster)

cafs_de_0vs1 <- FindMarkers(cafs, ident.1 = 0, ident.2 = 1)



for(i in 1:length(cafs_de_list)) {
  
  write.csv(cafs_de_list[[i]], paste0("../project/IgGaCSF1R_Final/CAFS_IgGaCSF1R_Cluster_", i-1, ".csv"))
  print(i-1)
}


iCAFmyCAF_markers <- read_csv("../project/IgGaCSF1R_Final/myCAFiCAFmesCAF_genes.csv")
myCAF_filt <- iCAFmyCAF_markers$myCAF[!is.na(iCAFmyCAF_markers$myCAF)] 
myCAF_filt <- myCAF_filt[myCAF_filt %in% cafs_de$gene[abs(cafs_de$avg_log2FC) > 1]]

iCAF_filt <- iCAFmyCAF_markers$iCAF[!is.na(iCAFmyCAF_markers$iCAF)] 
iCAF_filt <- iCAF_filt[iCAF_filt %in% cafs_de$gene[abs(cafs_de$avg_log2FC) > 1]]

mesCAF_filt <- iCAFmyCAF_markers$mesCAF[!is.na(iCAFmyCAF_markers$mesCAF)] 
mesCAF_filt <- mesCAF_filt[mesCAF_filt %in% cafs_de$gene[abs(cafs_de$avg_log2FC) > 0.5]]

DoHeatmap(cafs, features = myCAF_filt, group.bar = T, assay = "SCT", slot = "scale.data")
DoHeatmap(cafs, features = iCAF_filt, group.bar = T, assay = "SCT", slot = "scale.data")
DoHeatmap(cafs, features = mesCAF_filt, group.bar = T, assay = "SCT", slot = "scale.data")
scale.dat <- cafs@assays$SCT@scale.data

library(pheatmap)
cafs_meta_ordered <- cafs@meta.data
cafs_meta_ordered <- cafs_meta_ordered %>% arrange(SCT_snn_res.0.3)
pheatmap_annot <- data.frame(Cluster = cafs_meta_ordered$SCT_snn_res.0.3,
                             aCSF1R = if_else(cafs_meta_ordered$orig.ident == "IgG", 0,1),
                             IgG = if_else(cafs_meta_ordered$orig.ident == "IgG", 1, 0),
                             row.names = rownames(cafs_meta_ordered))

my_cols <- colorRampPalette(colors = c("darkgreen", "white", "red"))(101)

cafs@meta.data$celltype_encode <- if_else(cafs$orig.ident == "IgG", 1, 0)

heatmap_breaks <- c(sum(cafs_meta_ordered$SCT_snn_res.0.3 == 0), 
                    sum(cafs_meta_ordered$SCT_snn_res.0.3 == 1),
                    sum(cafs_meta_ordered$SCT_snn_res.0.3 == 2),
                    sum(cafs_meta_ordered$SCT_snn_res.0.3 == 3))

heatmap_breaks <- cumsum(heatmap_breaks)

other_markers_from_mei <- read_csv("../project/IgGaCSF1R_Final/Schmid v Tuveson.csv")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

default_gg_cols <- gg_color_hue(4)


newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotdf$category))))
mycolors <- newCols(length(unique(annotdf$category)))
names(mycolors) <- unique(annotdf$category)
mycolors <- list(category = mycolors)

names(default_gg_cols) <- factor(c(0,1,2,3), levels = c("0", "1", "2", "3"))
default_gg_cols <- list(Cluster = default_gg_cols)

other_markers_from_mei$Associated.Gene.Name -> other_markers_from_mei
rna_scale <- cafs@assays$SCT@scale.data

rna_scale <- rna_scale[,match(rownames(cafs_meta_ordered), colnames(rna_scale))] 
pheatmap(rna_scale[rownames(rna_scale) %in% myCAF_filt,], main = "myCAF markers abs(avglogFC) > 1",annotation_colors = default_gg_cols, breaks = seq(-2,2, length.out = 101), gaps_col = heatmap_breaks[1:3], color = my_cols, annotation = pheatmap_annot, scale = "row", cluster_cols = FALSE, show_colnames = F)
pheatmap(rna_scale[rownames(rna_scale) %in% iCAF_filt,], main = "iCAF markers abs(avglogFC) > 1",annotation_colors = default_gg_cols, breaks = seq(-2,2, length.out = 101), gaps_col = heatmap_breaks[1:3], color = my_cols, annotation = pheatmap_annot, scale = "row", cluster_cols = FALSE, show_colnames = F)
pheatmap(rna_scale[rownames(rna_scale) %in% mesCAF_filt,], main = "mesCAF markers abs(avglogFC) > 0.5",annotation_colors = default_gg_cols, breaks = seq(-2,2, length.out = 101), gaps_col = heatmap_breaks[1:3], color = my_cols, annotation = pheatmap_annot, scale = "row", cluster_cols = FALSE, show_colnames = F)


pheatmap(rna_scale[rownames(rna_scale) %in% other_markers_from_mei,],annotation_colors = default_gg_cols, main = "Tuveson Markers", breaks = seq(-2,2, length.out = 101), gaps_col = heatmap_breaks[1:3], color = my_cols, annotation = pheatmap_annot, scale = "row", cluster_cols = FALSE, show_colnames = F)




igg_alone_cafs <- cafs[,cafs$orig.ident == "IgG"]
DimPlot(igg_alone_cafs)
igg_alone_cafs <-RunPCA(igg_alone_cafs) %>% 
  RunUMAP(dims = 1:20, n.neighbors = 50) %>%
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.3)

DimPlot(igg_alone_cafs, label = T, pt.size = 2)
ggsave("../project/IgGaCSF1R_Final/IgG_alone_clusters.png", dpi = 1000, height = 10, width = 10, units = "in")
DimPlot(igg_alone_cafs, label = F, pt.size = 2, group.by = "orig.ident")

igg_alone_markers <- FindAllMarkers(igg_alone_cafs)
write.csv(igg_alone_markers, "../project/IgGaCSF1R_Final/IgGAlone_markers.csv")
DimPlot(all_dat)


jakstat_genes[jakstat_genes %in% rownames(cafs_0_igg_vs_acsf)]

show_props <-data.frame(CellType = cafs$orig.ident,
                        Cluster = cafs$SCT_snn_res.0.3)


show_props <- show_props %>% group_by(CellType, Cluster) %>% summarise(proportion = n())

show_props$proportion2 <- if_else(show_props$CellType == "aCSF1R", ( show_props$proportion / sum(show_props$proportion[show_props$CellType == "aCSF1R"]) )*100, (show_props$proportion / sum(show_props$proportion[show_props$CellType == "IgG"]) ) * 100)

my_gg_theme <- theme(axis.title.y = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(size = 0.25, linetype = 'dashed',
                                      colour = "grey"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'dashed',
                                      colour = "grey"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5))
show_props$CellType <- factor(show_props$CellType, levels = c("IgG", "aCSF1R"))
ggplot(show_props, aes(x = Cluster, y = proportion2, fill = CellType)) +
  geom_bar(alpha = 0.8, position = position_dodge(), stat = "identity", col = "black") +
  my_gg_theme + ylab("% of Cells") 




ggplot(cafs@meta.data[cafs$orig.ident == "aCSF1R",], aes(x = SCT_snn_res.0.3, fill = SCT_snn_res.0.3)) +
  geom_density() + my_gg_theme + xlab("Cluster") + ylab("Density")

cafs$MAF_Type <- case_when(cafs$SCT_snn_res.0.3 == 0 ~ "periMAFs",
                           cafs$SCT_snn_res.0.3 == 1 ~ "myMAFs",
                           cafs$SCT_snn_res.0.3 == 2 ~ "iMAF",
                           cafs$SCT_snn_res.0.3 == 3 ~ "cycMAF")

DimPlot(cafs, group.by = "MAF_Type")

FeaturePlot(cafs, features = "Spp1", pt.size = 4) + xlab("UMAP 1") + ylab("UMAP 2") + scale_color_gradient2(low = c(alpha("darkgrey",0.3)),
                                                                                                            high = c(rev(brewcols)[10]))
                                                                                                                                                                                                                    high = c(rev(brewcols)[10]))


cols = rev(brewcols)[5:10]

FeaturePlot(cafs, features = "Spp1", pt.size = 4) + xlab("UMAP 1") + ylab("UMAP 2") + scale_color_gradient2(low = c(alpha("darkgrey",0.3)),
                                                                                                            high = c(rev(brewcols)[10]))
                                                                                                                    
scale_colour_manual(values = c(alpha("darkgrey",0.3), "red3",alpha("lightgrey",0.3)))





p_annotop10 <- data.frame(MAF = cafs@meta.data$MAF_Type[order(desc(cafs@meta.data$MAF_Type))], row.names = rownames(cafs@meta.data)[order(desc(cafs@meta.data$MAF_Type))])

rna_scale2 <- ScaleData(cafs, features = rownames(cafs))
rna_scale2 <- GetAssayData(rna_scale2, slot = "scale.data")
dim(rna_scale2)
rna_scale2 <- rna_scale2 %>% as.data.frame() %>% rownames_to_column(var = "gene")

top10genes_ordered <- data.frame(Genes = cafs_de_top10_bycluster$gene)
top10genes_ordered <- top10genes_ordered %>% left_join(rna_scale2, by = c("Genes" = "gene"))
top10genes_ordered_genes <- top10genes_ordered[,1]
top10genes_ordered <- top10genes_ordered[,-1]

Idents(cafs) <- cafs@meta.data$MAF_Type
DoHeatmap(cafs, group.by = "MAF_Type", features = cafs_de_top10_bycluster$gene)

#rna_scale2 <- rna_scale2[rownames(rna_scale2) %in% cafs_de_top10_bycluster$gene ,]

dim(rna_scale2)
rna_scale2 <- rna_scale2[,order(desc(cafs@meta.data$MAF_Type))]



match(rownames(cafs_de_top10_bycluster), rownames(rna_scale2))

cafs_de_top10_bycluster <- cafs_de %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10genes_ordered <- data.frame(Genes = cafs_de_top10_bycluster$gene)
top10genes_ordered <- top10genes_ordered %>% left_join(rna_scale2, by = c("Genes" = "gene"))
top10genes_ordered <- top10genes_ordered %>% column_to_rownames(var = "Genes")
top10genes_ordered <- top10genes_ordered[,order(desc(cafs@meta.data$MAF_Type))]

maf_cols <- default_gg_cols
names(maf_cols) <- c("periMAF", "myMAF", "iMAF", "cycMAF")
top10_colours <- list(MAF = maf_cols)

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Top10GenesLogFC_MAFTypes_FixedColoursJuly.pdf", width = 5, height = 4)
pheatmap(top10genes_ordered, 
         scale = "row", cluster_rows = F, cluster_cols = F, 
         show_colnames = F, breaks = seq(-1.5,1.5,length.out = 101), 
         annotation = p_annotop10, clustering_method = "ward.D2", 
         fontsize_row = 8, annotation_colors = top10_colours, gaps_col = c(877, (877+701), (877+701+184)), gaps_row = c(10,20,30))
dev.off()




cafs_de_top20_bycluster <- cafs_de %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top20genes_ordered <- data.frame(Genes = cafs_de_top20_bycluster$gene)
top20genes_ordered <- top20genes_ordered %>% left_join(rna_scale2, by = c("Genes" = "gene"))
top20genes_ordered <- top20genes_ordered %>% column_to_rownames(var = "Genes")
top20genes_ordered <- top20genes_ordered[,order(desc(cafs@meta.data$MAF_Type))]

maf_cols <- default_gg_cols
names(maf_cols) <- c("periMAF", "myMAF", "iMAF", "cycMAF")
top10_colours <- list(MAF = maf_cols)

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Top20GenesLogFC_MAFTypes_FixedColoursJuly.pdf", width = 5, height = 4)
pheatmap(top20genes_ordered, 
         scale = "row", cluster_rows = F, cluster_cols = F, 
         show_colnames = F, breaks = seq(-1.5,1.5,length.out = 101), 
         annotation = p_annotop10, clustering_method = "ward.D2", 
         fontsize_row = 8, annotation_colors = top10_colours, 
         gaps_col = c(877, (877+701), (877+701+184)), gaps_row = c(20,40,60), fontsize = 6)
dev.off()

write.csv(top20genes_ordered, "~/cellranger_out_temp/PublicationPlots/Top20GenesLogFC_GeneList.csv")









rna_top10_peri <- rowMeans(rna_scale2[,p_annotop10$MAF == "periMAF"])
rna_top10_myMaf <- rowMeans(rna_scale2[,p_annotop10$MAF == "myMAF"])
rna_top10_iMaf <- rowMeans(rna_scale2[,p_annotop10$MAF == "iMAF"])
rna_top10_cycMaf <- rowMeans(rna_scale2[,p_annotop10$MAF == "cycMAF"])

head(rna_top10_cycMaf)
head(rna_top10_iMaf)

top_10_avg <- data.frame(periMAF = rna_top10_peri,
                         myMAF = rna_top10_myMaf,
                         iMAF = rna_top10_iMaf,
                         cycMAF = rna_top10_cycMaf, row.names = names(rna_top10_peri))

pdf(file = "/home/rstudio/cellranger_out_temp/PublicationPlots/Heatmap_Top10Avg.pdf", width = 5, height = 4)
pheatmap(t(top_10_avg), scale = "column", clustering_method = "ward.D2",
         show_colnames = T, cluster_rows = F,
         clustering_distance_cols = "euclidean",clustering_distance_rows = "euclidean",
         annotation_colors = annot_colors, angle_col = 45, fontsize = 5, cellheight = 15)
dev.off()







