library(tidyverse)
library(Seurat)

load('data/00_process_output.Rdata')
# get clusters from seuratCC object
DimPlot(seuratCC, reduction = "umap", group.by = 'seurat_clusters', label = TRUE, label.box = TRUE) +
  scale_color_manual(values = pals::alphabet2() %>% unname()) +
  scale_fill_manual(values = pals::alphabet2() %>% unname() )
# 0, 1, 2, 5, 6
# pulls barcodes
bcs <- row.names( seuratCC@meta.data %>% filter(seurat_clusters %in% c('0','1','2','5','6')) )

# rm(seurat)
# rm(seuratCC)
# rm(seurat.markers)
# rm(seurat.markersCC)

# reload
seurat.data  <- Read10X(data.dir = '~/data/sc_mouse_OFC_e11/aggr/outs/count/filtered_feature_bc_matrix')
seurat <- CreateSeuratObject(counts = seurat.data, project = "e11_mouse_eye_cup", min.cells = 3, min.features = 200)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

seurat@meta.data$replicate <- str_extract(colnames(seurat), '\\d')
seurat@meta.data$orig.ident <- paste0(seurat@meta.data$orig.ident, '__', seurat@meta.data$replicate)
Idents(seurat) <- seurat@meta.data$orig.ident


seurat <- subset(seurat, subset = nFeature_RNA > 200  & percent.mt < 10)



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


m.s.genes <-s.genes %>% str_to_title()
m.g2m.genes <- g2m.genes %>% str_to_title()
seurat <- CellCycleScoring(seurat, s.features = m.s.genes, g2m.features = m.g2m.genes)


seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score
options(future.globals.maxSize = 8000 * 1024^2)

seuratSUBSET <- seurat[,bcs]
genes_of_interest <- c('Ntn1','Vax1','Sox1','Fzd8','Shtn1','Slitrk1','Dlk1','Pmel','Mitf','Pax2','Pax6','Sox2','Vsx2','Vax2','Aldh1a3','Bmp4','Tbx5','Aldh1a1','Igf2','Jag1')

plan(strategy = "multicore", workers = 10)
seuratCCSUB <- ScaleData(seuratSUBSET, features = rownames(seuratSUBSET))
seuratCCSUB <- FindVariableFeatures(seuratCCSUB, selection.method = "vst", nfeatures = 1000)
seuratCCSUB <- RunPCA(seuratCCSUB, features = c(genes_of_interest,
                                                VariableFeatures(object = seuratCCSUB) %>%
                                                  grep('^Hist|^Rpl|^Rps|^mt-|^Hbb',.,value =TRUE, invert = TRUE)) %>%
                        unique(),
                      npcs = 50)
DimPlot(seuratCCSUB, reduction = "pca")


seuratCCSUB <- FindNeighbors(seuratCCSUB, dims = 1:20)
seuratCCSUB <- FindClusters(seuratCCSUB, resolution = 1.2)
seuratCCSUB <- RunUMAP(seuratCCSUB, dims = 1:20, min.dist = 0.3)



# CT assignment from 02_figures.Rmd

ct_cluster_CC <- rbind(c(0, "RPC"),
                       c(1, 'RPC'),
                       c(2, 'RPE Precursor'),
                       c(3, 'Fibroblast'),
                       c(4, 'Fibroblast'),
                       c(5, 'Optic Fissure'),
                       c(6, 'Proliferating RPC'),
                       c(7, 'Epithelial'),
                       c(8, 'Red Blood Cell'),
                       c(9, 'Corneal Progenitor'),
                       c(10, 'Blood Vessel'),
                       c(11, "PBMC"),
                       c(12, "NK/T"),
                       c(13, "Lens"
                       )) %>% data.frame()
colnames(ct_cluster_CC) <- c('seurat_clusters','CellType')
seuratCC@meta.data$CellType <- seuratCC@meta.data %>% dplyr::select(seurat_clusters) %>% left_join(ct_cluster_CC) %>% pull(CellType)
bc_ct <- seuratCC@meta.data %>% as_tibble(rownames = 'Barcode') %>%  dplyr::select(Barcode, CellType)
seuratCCSUB@meta.data$CellType <- seuratCCSUB@meta.data %>% as_tibble(rownames = 'Barcode') %>% left_join(bc_ct, by = 'Barcode') %>% pull(CellType)


seurat.markersCCSUB <- list()
for (i in seuratCCSUB@meta.data$seurat_clusters %>% unique() %>% sort()){
  seurat.markersCCSUB[[i]] <- FindMarkers(seuratCCSUB, ident.1 = i, only.pos = TRUE, logfc.threshold = 0.2) %>% as_tibble(rownames = 'Gene')
}
seuratCCSUB <- NormalizeData(seuratCCSUB, normalization.method = "LogNormalize", scale.factor = 10000)
save(seuratCCSUB, seurat.markersCCSUB, file = 'data/03_nr_rpe_of_filtering.freeze_2023.11.01.Rdata')
#'
#' #'
#' dot_size <- 0.4
#' cowplot::plot_grid(plotlist = list(DimPlot(seuratCCSUB, reduction = "umap", label = TRUE, label.box = TRUE) +
#'                                      scale_color_manual(values = pals::alphabet2() %>% unname()) +
#'                                      scale_fill_manual(values = pals::alphabet2() %>% unname()),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap",  features = 'CellType' ) +
#'                                      scale_color_manual(values = pals::alphabet() %>% unname()) +
#'                                      scale_fill_manual(values = pals::alphabet() %>% unname()),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Ntn1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Vax1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Sox1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Fzd8'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Shtn1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Slitrk1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Smoc1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Dlk1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Pmel'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Mitf'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Pax2'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Pax6'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Sox2'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Vsx2'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Vax2'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Aldh1a3'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Aldh1a1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Bmp4'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Tbx5'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Foxd1'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Foxn4'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Foxn4'), pt.size = dot_size),
#'                                    FeaturePlot(seuratCCSUB, reduction = "umap", features =
#'                                                  c('Tbx2'), pt.size = dot_size)
#' ),
#'                    ncol = 4)
#' #'
#' #'
#' #' genes <- c('Ntn1','Vax1','Vax2','Pax2','Pax6',
#' #'            'Pitx2',
#' #'            'Bmp4','Yap1', 'Zfp503', 'Pmel', 'Abca4', 'Tenm3', 'Aldh1a3',
#' #'            'Dkk1','Tbx5',
#' #'            #'Hmga2', 'Aldh1a1', 'Rps26', 'Trpm1', 'Serf2',
#' #'            'Rbp1','Slitrk1', 'Slit2')
#' #'
# genes <- c('Ntn1','Pmel','Pax2','Pax6','Vax1','Sox2','Bmp4','Vsx2', "Vax2","Tbx5", 'Shtn1','Dlk1','Fezf1','Sox1','Vim','Zic1')
# FeaturePlot(seuratCCSUB, reduction = "umap", features =
#               genes, pt.size = 0.2)
#
#
# library(ComplexHeatmap)
# counts <- Seurat::FetchData(seuratCCSUB, vars = c('seurat_clusters', 'CellType', genes)) %>% arrange(CellType,seurat_clusters)
# ct_color <- pals::glasbey(n = length(unique(counts$CellType)))
# names(ct_color) <- counts$CellType %>% unique() %>% sort()
#
# clus_color <- pals::alphabet(n = length(unique(counts$seurat_clusters)))
# names(clus_color) <- counts$seurat_clusters %>% unique() %>% sort()
#
# col_anno <- HeatmapAnnotation(df = data.frame(seurat_clusters = counts$seurat_clusters %>% factor(),
#                                               CellType = counts$CellType %>% as.factor()),
#                               col = list(seurat_clusters = clus_color,CellType = ct_color))
#
# hm_counts <- counts[,c(3:ncol(counts))] %>% as.matrix() %>% t()
# colnames(hm_counts) <- NULL
#
# col_fun <- circlize::colorRamp2(c(0,2), c("#440154FF","#FDE725FF"))
# Heatmap((hm_counts), top_annotation = col_anno, cluster_columns = FALSE, col = col_fun, name = ' ' , use_raster = FALSE)

