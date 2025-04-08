options(stringsAsFactors = F)
suppressMessages(library("optparse"))
suppressMessages(library("CellTrek"))
suppressMessages(library("dplyr"))
suppressMessages(library("Seurat"))
suppressMessages(library("viridis"))
suppressMessages(library("ConsensusClusterPlus"))
suppressMessages(library("tidyverse"))
option_list <- list(
    make_option(c("-a","--scRNA_celltype_rds"),type="character", default=NULL, help="input scRNA_celltype.rds"),
    make_option(c("-b","--Spatial_rds"),type="character", default=NULL, help="input Spatial.rds"),
    make_option(c("-c","--ST_pdf"),type="character", default=NULL, help="output Spatial annotation file"),
    make_option(c("-d","--sc_st_coembedding"),type="character", default=NULL, help="output sc_st_coembedding file"),
    make_option(c("-e","--sc_st_SpatialDimPlot"),type="character", default=NULL, help="output sc_st_SpatialDimPlot file"),
    make_option(c("-f","--celltype"),type="character", default=NULL, help="select celltype"),
    make_option(c("-g","--heatmap"),type="character", default=NULL, help="output heatmap"),
    make_option(c("-m","--celltype_featureplot_module"),type="character", default=NULL, help="output celltype_featureplot_module file"),
    make_option(c("-i","--celltype_spatialplot_module"),type="character", default=NULL, help="output celltype_spatialplot_module file"),
    make_option(c("-j","--gene"),type="character", default=NULL, help="select gene"),
    make_option(c("-k","--gene_spatial"),type="character", default=NULL, help="output gene spatial infomation"),
    make_option(c("-l","--gene_spatial_module"),type="character", default=NULL, help="output gene_spatial_module infomation")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$scRNA_celltype_rds)|| is.null(opt$Spatial_rds)|| is.null(opt$ST_pdf)|| is.null(opt$sc_st_coembedding)|| is.null(opt$sc_st_SpatialDimPlot)|| 
    is.null(opt$celltype)|| is.null(opt$heatmap)||is.null(opt$celltype_featureplot_module)||is.null(opt$celltype_spatialplot_module)||
    is.null(opt$gene)|| is.null(opt$gene_spatial)||is.null(opt$gene_spatial_module)){
  print_help(opt_parser);
  stop("Please provide -a scRNA_celltype_rds, -b Spatial_rds ,-c ST_pdf,-d sc_st_coembedding,-e sc_st_SpatialDimPlot ,-f celltype,
    -g heatmap,-m celltype_featureplot_module ,-i celltype_spatialplot_module,-j gene,-k gene_spatial,-l gene_spatial_module", call.=FALSE);
}

scRNA_rds=opt$scRNA_celltype_rds
SpatialRds=opt$Spatial_rds
STpdf=opt$ST_pdf
Coembedding=opt$sc_st_coembedding
SpatialDimPlot=opt$sc_st_SpatialDimPlot
Celltype=opt$celltype
Heatmap=opt$heatmap
CelltypeFeatureplot=opt$celltype_featureplot_module
CelltypeSpatialplot=opt$celltype_spatialplot_module
Gene=opt$gene
GeneSpatial=opt$gene_spatial
GeneSpatialModule=opt$gene_spatial_module
##############################################################################################################################
st <- readRDS(SpatialRds)
sc <- readRDS(scRNA_rds)
#cells.sub1 <- subset(sc@meta.data,tissue == c("Brain Myeloid Cortex"))
#cells.sub2 <- subset(sc@meta.data,tissue == c("Brain Non-Myeloid Cortex"))
#cells.sub = rbind(cells.sub1,cells.sub2)
#sc <- subset(sc,cells = row.names(cells.sub))

## Rename the cells/spots with syntactically valid names
st <- RenameCells(st, new.names=make.names(Cells(st)))
sc <- RenameCells(sc, new.names=make.names(Cells(sc)))
Idents(sc) <- "celltype"
## Visualize the ST data
SpatialDimPlot(st)
ggsave(STpdf)
## Visualize the scRNA-seq data
#DimPlot(sc, label = T, label.size = 4.5)
#ggsave("result/sc.pdf")

## Cell charting using CellTrek
brain_traint <- CellTrek::traint(st_data=st, sc_data=sc, sc_assay='RNA', cell_names='celltype')
DimPlot(brain_traint, group.by = "type")
ggsave(Coembedding)
#After coembedding, we can chart single cells to their spatial locations. Here, we use the non-linear interpolation (intp = T, intp_lin=F) approach to augment the ST spots.
brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', sc_data=sc, sc_assay = 'RNA', 
                                   reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                   dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek
#After cell charting, we can interactively visualize the CellTrek result using celltrek_vis

brain_celltrek$celltype <- factor(brain_celltrek$celltype, levels=sort(unique(brain_celltrek$celltype)))
Idents(brain_celltrek) <- "celltype"
SpatialDimPlot(brain_celltrek)
ggsave(SpatialDimPlot)
# table(brain_celltrek$celltype)

#       Astrocytes Endothelial cells         Microglia           Neurons 
#             1068               599              2198                91 
#         NK cells  Oligodendrocytes     Stromal cells 
#               52              1492               225

## Cell colocalization analysis
#glut_cell <- c('Microglia', 'Oligodendrocytes', 'Astrocytes', 'Neurons')
#names(glut_cell) <- make.names(glut_cell)
#brain_celltrek_glut <- subset(brain_celltrek, subset=celltype %in% glut_cell)
#brain_celltrek_glut$celltype <- factor(brain_celltrek_glut$celltype, levels=glut_cell)

#Then we can use scoloc module to perform colocalization analysis.
#brain_sgraph_KL <- CellTrek::scoloc(brain_celltrek_glut, col_cell='celltype', use_method='KL', eps=1e-50)
## We extract the minimum spanning tree (MST) result from the graph
#brain_sgraph_KL_mst_cons <- brain_sgraph_KL$mst_cons
#rownames(brain_sgraph_KL_mst_cons) <- colnames(brain_sgraph_KL_mst_cons) <- glut_cell[colnames(brain_sgraph_KL_mst_cons)]
## We then extract the metadata (including cell types and their frequencies)
#brain_cell_class <- brain_celltrek@meta.data %>% dplyr::select(id=celltype) %>% unique
#brain_celltrek_count <- data.frame(freq = table(brain_celltrek$celltype))
#brain_cell_class_new <- merge(brain_cell_class, brain_celltrek_count, by.x ="id", by.y = "freq.Var1")
#CellTrek::scoloc_vis(brain_sgraph_KL_mst_cons, meta_data=brain_cell_class)

## Spatial-weighted gene co-expression analysis within the cell type of interest
#基于 CellTrek 结果，我们可以使用 SCoexp 模块进一步研究感兴趣的细胞类型内的共表达模式。 
#在这里，我们将使用共识聚类（CC）方法以Microglia为例。 首先从图表结果中提取Microglia
brain_celltrek_cellsub <- subset(brain_celltrek, subset=celltype== Celltype)
brain_celltrek_cellsub@assays$RNA@scale.data <- matrix(NA, 1, 1)
#brain_celltrek_cellsub$cluster <- gsub('L5 IT VISp ', '', brain_celltrek_cellsub$cluster)
#DimPlot(brain_celltrek_l5, group.by = 'cluster')
brain_celltrek_cellsub<- FindVariableFeatures(brain_celltrek_cellsub)
vst_df <- brain_celltrek_cellsub@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(brain_celltrek_cellsub[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^mt-', rownames(brain_celltrek_cellsub), value=T)
rp_gene <- grep('^Rpl|^Rps', rownames(brain_celltrek_cellsub), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]
#
brain_celltrek_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=brain_celltrek_cellsub, assay='RNA', approach='cc', sigm=140, avg_cor_min=.4, zero_cutoff=3, min_gen=40, max_gen=400)

brain_celltrek_k =c()
for(i in 1:length(names(brain_celltrek_scoexp_res_cc$gs))){
    data_frame=data.frame(gene=c(brain_celltrek_scoexp_res_cc$gs[[i]]), G=names(brain_celltrek_scoexp_res_cc$gs)[i])
    brain_celltrek_k=rbind(brain_celltrek_k,data_frame)
}
brain_celltrek_k=brain_celltrek_k %>%magrittr::set_rownames(.$gene) %>% dplyr::select(-1)

pheatmap::pheatmap(brain_celltrek_scoexp_res_cc$wcor[rownames(brain_celltrek_k), rownames(brain_celltrek_k)], 
                   clustering_method='ward.D2', annotation_row=brain_celltrek_k, show_rownames=F, show_colnames=F, 
                   treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                   color=viridis(10), main='spatial co-expression',filename=Heatmap)
brain_celltrek_cellsub<- AddModuleScore(brain_celltrek_cellsub, features=brain_celltrek_scoexp_res_cc$gs, name='CC_', nbin=10, ctrl=50, seed=42)
## First we look into the coexpression module based on the sCCcRNA-seq embeddingCC
FeaturePlot(brain_celltrek_cellsub, grep('CC_', colnames(brain_celltrek_cellsub@meta.data), value=T))
ggsave(CelltypeFeatureplot)
SpatialFeaturePlot(brain_celltrek_cellsub, grep('CC_', colnames(brain_celltrek_cellsub@meta.data), value=T))
ggsave(CelltypeSpatialplot)
SpatialFeaturePlot(brain_celltrek_cellsub,features=Gene)
ggsave(GeneSpatial)
SpatialFeaturePlot(brain_celltrek_cellsub, ,features=Gene,grep('CC_', colnames(brain_celltrek_cellsub@meta.data), value=T))
ggsave(GeneSpatialModule)

