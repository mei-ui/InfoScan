suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))
option_list <- list(
                make_option(c("-a","--finall_fpkm"),type="character", default=NULL, help="input finally_fpkm file"),
                make_option(c("-b","--group_chose"),type="character", default=NULL, help="input id file"),
		            make_option(c("-c","--scRNA_rds"),type="character", default=NULL, help="input scRMA.rds file"),
                make_option(c("-e","--variable_features"),type="character", default=NULL, help="output VariableFeatures.pdf file"),
                make_option(c("-f","--pca"),type="character", default=NULL, help="output pca.pdf file"),
                make_option(c("-g","--cell_cluster"),type="character", default=NULL, help="output cell_cluster.csv file"),
                make_option(c("-l","--tsne"),type="character", default=NULL, help="output tSNE.pdf file"),
                make_option(c("-i","--umap"),type="character", default=NULL, help="output UMAP.pdf file"),
		            make_option(c("-m","--lncRNA_heatmap"),type="character", default=NULL, help="output top10 lncRNA marker gene file"),
		            make_option(c("-n","--tsne_tissue"),type="character", default=NULL, help="output tSNE_tissue.pdf file"),
		            make_option(c("-o","--umap_tissue"),type="character", default=NULL, help="output umap_tissue.pdf file"),
                make_option(c("-j","--resolution_chose"),type="character", default=NULL, help="input resolution_chose"),
                make_option(c("-k","--pca_num"),type="character", default=NULL, help="input pca number"),
                make_option(c("-t","--outputdir"),type="character", default=NULL, help="input outputdir")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$finall_fpkm)|| is.null(opt$group_chose)|| is.null(opt$variable_features)|| is.null(opt$pca)|| is.null(opt$cell_cluster)|| is.null(opt$tsne)|| 
    is.null(opt$umap)|| is.null(opt$scRNA_rds)||is.null(opt$lncRNA_heatmap)||is.null(opt$tsne_tissue)||is.null(opt$umap_tissue)||is.null(opt$resolution_chose)
    ||is.null(opt$pca_num)||is.null(opt$outputdir)){
  print_help(opt_parser);
  stop("Please provide -a finall_fpkm, -b group_chose ,-e VariableFeatures ,-f pca,-g cell_cluster and -l tSNE,-i UMAP,-m lncRNA_marker,-c scRNA_rds,-n tsne_tissue,
      -o umap_tissue,-j resolution_chose,-k pca_num", call.=FALSE);
}
Finallfpkm=opt$finall_fpkm
SCRNA=opt$scRNA_rds
VariableFeatures=opt$variable_features
PCA=opt$pca
CELLcluster=opt$cell_cluster
TSNE=opt$tsne
UMAP=opt$umap
LncRNAheatmap=opt$lncRNA_heatmap
TSNE_tissue=opt$tsne_tissue
UMAP_tissue=opt$umap_tissue
group=opt$group_chose
Resolution=opt$resolution_chose
PcaNum=opt$pca_num
OutDir=opt$outputdir
##加载数据
finally_fpkm <- read.table(Finallfpkm,header=T)
name<-as.data.frame(table(finally_fpkm$Geneid))
name<-filter(name,Freq>1)
if(nrow(name)!=0){finally_fpkm_1<-finally_fpkm[-which(finally_fpkm$Geneid%in%c(as.character(name$Var1))),]}else {finally_fpkm_1<-finally_fpkm}
protein_coding<-finally_fpkm_1[which(finally_fpkm_1$type%in%c('protein_coding')),]
annotated_lncRNA<-finally_fpkm_1[which(finally_fpkm_1$type%in%c('annotated_lncRNA')),]
unannotated_lncRNA<-finally_fpkm_1[which(finally_fpkm_1$type%in%c('unannotated_lncRNA')),]
unannotated_lncRNA$Geneid<-stri_replace_all_regex(unannotated_lncRNA$Geneid,'_','-')
unannotated_protein<-finally_fpkm[which(finally_fpkm$type%in%c('unannotated_coding_transcript')),]
unannotated_protein$Geneid<-stri_replace_all_regex(unannotated_protein$Geneid,'_','-')

protein_matrix<-rbind(as.data.frame(protein_coding),as.data.frame(unannotated_protein))
protein_id<-as.character(protein_matrix$Geneid)
lncRNA_matrix = rbind(annotated_lncRNA,unannotated_lncRNA)
lncRNA_id<-as.character(lncRNA_matrix$Geneid)
scRNA<-readRDS(SCRNA)
lncRNA<-scRNA[which(rownames(scRNA)%in%c(lncRNA_id)),]
protein<-scRNA[which(rownames(scRNA)%in%c(protein_id)),]
##06寻找高变基因
lncRNA <- NormalizeData(lncRNA, normalization.method = "LogNormalize", scale.factor = 10000)
lncRNA <- FindVariableFeatures(lncRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(lncRNA), 10)
plot1 <- VariableFeaturePlot(lncRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom")
ggsave(VariableFeatures, plot = plot, width = 8, height = 6)

##07标准化
scale.genes <-  rownames(lncRNA)
lncRNA<- ScaleData(lncRNA, features = scale.genes)
# 结果存储在scRNA[["RNA"]]@scale.data中

##08PCA降维
lncRNA<- RunPCA(lncRNA,features=VariableFeatures(lncRNA))
plot1 <- DimPlot(lncRNA, reduction = "pca")
plot2 <- ElbowPlot(lncRNA, ndims=20, reduction="pca")
plotc <- plot1+plot2
ggsave(PCA, plot = plotc, width = 8, height = 4)

##09细胞聚类
PcaNum=as.numeric(PcaNum)
PcaNum
pc.num=1:PcaNum
Resolution=as.numeric(Resolution)
Resolution
lncRNA <- FindNeighbors(lncRNA, dims = pc.num)
lncRNA <- FindClusters(lncRNA, resolution=Resolution)
table(lncRNA@meta.data$seurat_clusters)
metadata <- lncRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,CELLcluster,row.names = F)

##10降维之TSNE
set.seed(888)
lncRNA = RunTSNE(lncRNA, dims = pc.num,perplexity=10,check_duplicates = FALSE)
embed_tsne <- Embeddings(lncRNA, 'tsne')
plot1 = DimPlot(lncRNA, reduction = "tsne",label = F)
ggsave(TSNE, plot = plot1, width = 8, height = 7)
p1 = DimPlot(lncRNA, reduction = "tsne",label = F,group.by=group)
ggsave(TSNE_tissue, plot = p1, width = 10, height = 7)
##UMAP
lncRNA <- RunUMAP(lncRNA, dims = pc.num,perplexity=10)
embed_umap <- Embeddings(lncRNA, 'umap')
plot2 = DimPlot(lncRNA, reduction = "umap",label = F)
ggsave(UMAP, plot = plot2, width = 8, height = 7)
p2 = DimPlot(lncRNA, reduction = "umap",label = F,group.by=group)
ggsave(UMAP_tissue, plot = p2, width = 10, height = 7)

diff.wilcox = FindAllMarkers(lncRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_genes = CaseMatch(search = as.vector(top10$gene), match = rownames(lncRNA))
plot1 = DoHeatmap(lncRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave(LncRNAheatmap, plot=plot1, width=15, height=18)
out_dir=paste0(OutDir,"/snakemake/result/11_data_analysis/cluster/lncRNA/lncRNA.rds")
saveRDS(lncRNA,file=out_dir)