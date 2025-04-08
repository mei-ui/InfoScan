suppressMessages(library('reticulate'))
python_env=Sys.which("celltypist")
python_env
python_env=as.character(python_env)
python_env=gsub("celltypist$","python",python_env)
use_python(python_env,required = T)
suppressMessages(library("optparse"))
suppressMessages(library("Seurat"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))

option_list <- list(
    make_option(c("-a","--scRNA_rds"),type="character", default=NULL, help="input scRNA.rds"),
    make_option(c("-b","--model"),type="character", default=NULL, help="input celltypist annotation file"),
		make_option(c("-d","--cell_type"),type="character", default=NULL, help="output celltype annotation file"),
		make_option(c("-e","--tSNE_celltype"),type="character", default=NULL, help="output tSNE-celltype.pdf file"),
		make_option(c("-f","--UMAP_celltype"),type="character", default=NULL, help="output UMAP-celltype.pdf file"),
    make_option(c("-g","--conservation_score"),type="character", default=NULL, help="input gene conservation score file"),
		make_option(c("-m","--scRNA_celltype_rds"),type="character", default=NULL, help="output scRNA_celltype.rds file"),
    make_option(c("-t","--outputdir"),type="character", default=NULL, help="input outputdir")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$scRNA_rds)|| is.null(opt$model)|| is.null(opt$cell_type)|| is.null(opt$tSNE_celltype)|| 
    is.null(opt$UMAP_celltype)|| is.null(opt$scRNA_celltype_rds)||is.null(opt$conservation_score)||is.null(opt$outputdir)){
  print_help(opt_parser);
  stop("Please provide -a scRNA_rds, -b model ,-d celltype,-e tSNE_celltype ,-f UMAP_celltype,-g conservation_score,-m scRNA_celltype_rds", call.=FALSE);
}

Model=opt$model
Cell_Type=opt$cell_type
TSNE=opt$tSNE_celltype
UMAP=opt$UMAP_celltype
SCrds=opt$scRNA_rds
scRNA_celltype=opt$scRNA_celltype_rds
Conservation = opt$conservation_score
OutDir=opt$outputdir
#载入python 模块
scanpy = import("scanpy")
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")
#加载参考数据集
#model_type = list.files("snakemake/script/Celltypist")
#model_type1 = paste("snakemake/script/Celltypist/",model_type,sep="")
#names(model_type) = str_split(string = model_type,pattern = "\\.", simplify = T)[,1]

#model_list = lapply(model_type, function(x){
#  celltypist$models$Model$load(model = x)})

model=celltypist$models$Model$load(model = Model)
#数据处理
scRNA<-readRDS(SCrds)
####seurat to celltypist
adata = scanpy$AnnData(X = numpy$array(as.matrix(t(as.matrix(scRNA[['RNA']]@counts)))),
                       obs = pandas$DataFrame(scRNA@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(scRNA[['RNA']]@counts),
                                                         row.names = rownames(scRNA[['RNA']]@counts)))
)
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
#细胞亚群预测和可视化
### 1. Immune_All_High
Model
predictions = celltypist$annotate(adata, model = model, majority_voting = T)
## 把这些信息加入到seurat对象中去
scRNA  = AddMetaData(scRNA, predictions$predicted_labels$majority_voting, col.name ="celltype") 


celltype_list = predictions$predicted_labels

write.csv(celltype_list,Cell_Type,row.names = T)

#scRNAsub@meta.data[which(scRNAsub@meta.data$seurat_clusters == 17),'celltype'] <- "Epithelial cells"
#scRNA$celltype=scRNA$CellType
p1 = DimPlot(scRNA, group.by="celltype", label=F, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype", label=F, label.size=5, reduction='umap')
#p3 = plotc <- p1+p2+ plot_layout(guides = 'collect')
ggsave(TSNE, p1)
ggsave(UMAP, p2)
saveRDS(scRNA,file=scRNA_celltype)
Model
#Conservation_score<-read.table(Conservation,header=T)
#Conservation_score$gene_name<-stri_replace_all_regex(Conservation_score$gene_name,'_','-')
if(length(unique(scRNA@meta.data$celltype))>1){
  for( i in unique(scRNA@meta.data$celltype)){
  markers_df <- FindMarkers(object = scRNA, ident.1 = i, min.pct = 0.25,group.by = 'celltype',only.pos=T,min.cells.group=1)
  write.csv(markers_df, paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_marker','.csv'), row.names = T)
  print(x = head(markers_df))
  markers_gene = markers_df  %>% subset(p_val<0.05)
  markers_gene = markers_gene %>% top_n(n = 15, wt = avg_log2FC)
  markers_gene = markers_gene[order(-markers_gene$avg_log2FC),]
  #markers_genes =  rownames(head(x = markers_df, n = 15))
  markers_genes = rownames(markers_gene)
  VlnPlot(object = scRNA, features =markers_genes,log =T ,pt.size = 0.1, ncol = 5,group.by = 'celltype')
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_VlnPlot','.pdf'),width=15, height=10,limitsize = FALSE)
  FeaturePlot(object = scRNA, features=markers_genes ,ncol = 5,reduction = "tsne")
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_FeaturePlot','.pdf'),width=15, height=8,limitsize = FALSE)
  #unannotated marker
  unannotated_markers_gene = markers_df  %>% subset(p_val<0.05)
  a<-unannotated_markers_gene[str_match(rownames(unannotated_markers_gene),pattern='TCONS-[0-9]*'),]
  a=na.omit(a)
  a=a[order(-a$avg_log2FC),]
  a = rownames(a)
  if(length(a)<=10&length(a)>0){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=8,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=8,limitsize = FALSE)
  }
  if(length(a)>10&length(a)<=20){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=12,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=12,limitsize = FALSE)
  }
  if(length(a)>20&length(a)<=30){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=16,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=16,limitsize = FALSE)
  }
  if(length(a)>30&length(a)<=40){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=20,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=20,limitsize = FALSE)
  }
  if(length(a)>40&length(a)<=50){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=24,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=24,limitsize = FALSE)
  }
  if(length(a)>50&length(a)<=60){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=28,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=28,limitsize = FALSE)
  }
  if(length(a)>60&length(a)<=70){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=32,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=32,limitsize = FALSE)
  }
  if(length(a)>70&length(a)<=80){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=36,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=36,limitsize = FALSE)
  }
  if(length(a)>80&length(a)<=90){
  p1<-VlnPlot(scRNA,features = a,
            pt.size = 0.1,group.by = "celltype",ncol = 5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_VlnPlot','.pdf'),p1,width=16,height=40,limitsize = FALSE)
  p2 <- FeaturePlot(scRNA,features = a,reduction = "tsne",label = T,ncol =5)
  ggsave(filename=paste0(OutDir,'/snakemake/result/11_data_analysis/celltype/',i,'_unannotated_FeaturePlot','.pdf'),p2,width=18 ,height=40,limitsize = FALSE)
  }
  }
}