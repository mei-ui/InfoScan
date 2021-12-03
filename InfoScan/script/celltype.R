suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))
suppressMessages(library('SingleR'))
suppressMessages(library('pheatmap'))

option_list <- list(
    make_option(c("-a","--scRNA_rds"),type="character", default=NULL, help="input scRNA.rds"),
    make_option(c("-b","--ref_all"),type="character", default=NULL, help="input singleR annotation file"),
    make_option(c("-c","--ref_imm"),type="character", default=NULL, help="input singleR annotation file"),
		make_option(c("-d","--cell_type"),type="character", default=NULL, help="output celltype annotation file"),
		make_option(c("-e","--tSNE_celltype"),type="character", default=NULL, help="output tSNE-celltype.pdf file"),
		make_option(c("-f","--UMAP_celltype"),type="character", default=NULL, help="output UMAP-celltype.pdf file"),
    make_option(c("-g","--conservation_score"),type="character", default=NULL, help="input gene conservation score file"),
		make_option(c("-m","--scRNA_celltype_rds"),type="character", default=NULL, help="output scRNA_celltype.rds file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$scRNA_rds)|| is.null(opt$ref_all)|| is.null(opt$ref_imm)|| is.null(opt$cell_type)|| is.null(opt$tSNE_celltype)|| 
    is.null(opt$UMAP_celltype)|| is.null(opt$scRNA_celltype_rds)||is.null(opt$conservation_score)){
  print_help(opt_parser);
  stop("Please provide -ascRNA_rds, -b refall ,-c refe_imm,-d celltype,-e tSNE_celltype ,-f UMAP_celltype,-g conservation_score,-m scRNA_celltype_rds", call.=FALSE);
}

Refall=opt$ref_all
Refimm=opt$ref_imm
Cell_Type=opt$cell_type
TSNE=opt$tSNE_celltype
UMAP=opt$UMAP_celltype
SCrds=opt$scRNA_rds
scRNA_celltype=opt$scRNA_celltype_rds
Conservation = opt$conservation_score
#数据处理
scRNA<-readRDS(SCrds)
ref.se=get(load(Refall))
ref.se_1=get(load(Refimm))
testdata <- GetAssayData(scRNA, slot="data")  
#testdata <- GetAssayData(scRNA,slot="scale.data")
clusters <- scRNA@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = list(ref.se,ref.se_1), labels = list(ref.se$label.main,ref.se_1$label.main), 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = 1)
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,Cell_Type,row.names = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

#scRNAsub@meta.data[which(scRNAsub@meta.data$seurat_clusters == 17),'celltype'] <- "Epithelial cells"

p1 = DimPlot(scRNA, group.by="celltype", label=F, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype", label=F, label.size=5, reduction='umap')
p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave(TSNE, p1, width=7 ,height=6)
ggsave(UMAP, p2, width=7 ,height=6)
saveRDS(scRNA,file=scRNA_celltype)

Conservation_score<-read.table(Conservation,header=T)
Conservation_score$gene_name<-stri_replace_all_regex(Conservation_score$gene_name,'_','-')
for( i in unique(scRNA@meta.data$celltype)){
  markers_df <- FindMarkers(object = scRNA, ident.1 = i, min.pct = 0.25,group.by = 'celltype')
  write.csv(markers_df, paste0('snakemake/result/11_data_analysis/celltype/',i,'_marker_','.csv'), row.names = T)
  print(x = head(markers_df))
  markers_genes =  rownames(head(x = markers_df, n = 15))
  VlnPlot(object = scRNA, features =markers_genes,log =T ,pt.size = 0.1, ncol = 5,group.by = 'celltype')
  ggsave(filename=paste0('snakemake/result/11_data_analysis/celltype/',i,'_VlnPlot_','.pdf'),width=15, height=10)
  FeaturePlot(object = scRNA, features=markers_genes ,ncol = 5,reduction = "tsne")
  ggsave(filename=paste0('snakemake/result/11_data_analysis/celltype/',i,'_FeaturePlot_','.pdf'),width=15, height=8)

  markers_df$Geneid = rownames(markers_df)
  unannotated_marker<-na.omit(str_match(markers_df$Geneid,pattern='TCONS-[0-9]*'))
  unannotated_marker<-as.data.frame(unannotated_marker)
  #print(x = head(unannotated_marker))
  unannotated_marker<-inner_join(markers_df,unannotated_marker,by=c("Geneid"="V1"))
  print(x = head(unannotated_marker))
  #unannotated_marker<-filter(unannotated_marker,p_val<="0.05")
  unannotated_marker<-inner_join(unannotated_marker,Conservation_score,by=c("Geneid"="gene_name"))
  print(x = head(unannotated_marker))
  write.csv(unannotated_marker, paste0('snakemake/result/11_data_analysis/celltype/',i,'_unannotated_marker_','.csv'), row.names = T)
  unannotated_gene = unannotated_marker$Geneid
  if(length(unannotated_gene)<=5){
  p1<-VlnPlot(scRNA,features = unannotated_gene,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0('snakemake/result/11_data_analysis/celltype/',i,'_unannotated_marker_V_','.pdf'), p1, width=16,height=4)
  p2 <- FeaturePlot(scRNA,features = unannotated_gene,
                  reduction = "tsne",label = T,ncol =5)
  ggsave(paste0('snakemake/result/11_data_analysis/celltype/',i,'_unannotated_marker_F_','.pdf'), p2, width=18 ,height=4)}
  if(length(unannotated_gene)<=10&length(unannotated_gene)>5){
  p1<-VlnPlot(scRNA,features = unannotated_gene,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0('snakemake/result/11_data_analysis/celltype/',i,'_unannotated_marker_V_','.pdf'), p1, width=16,height=8)
  p2 <- FeaturePlot(scRNA,features = unannotated_gene,
                  reduction = "tsne",label = T,ncol =5)
  ggsave(paste0('snakemake/result/11_data_analysis/celltype/',i,'_unannotated_marker_F_','.pdf'), p2, width=18 ,height=8)}
  if(length(unannotated_gene)<=15&length(unannotated_gene)>10){
  p1<-VlnPlot(scRNA,features = unannotated_gene,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0('snakemake/result/11_data_analysis/celltype/',i,'_unannotated_marker_V_','.pdf'), p1, width=16,height=12)
  p2 <- FeaturePlot(scRNA,features = unannotated_gene,
                  reduction = "tsne",label = T,ncol =5)
  ggsave(paste0('snakemake/result/11_data_analysis/celltype/',i,'_unannotated_marker_F_','.pdf'), p2, width=18 ,height=12)}


}