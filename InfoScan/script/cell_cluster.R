suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))
suppressMessages(library('sctransform'))

option_list <- list(
		make_option(c("-a","--finall_fpkm"),type="character", default=NULL, help="input finally fpkm expression matrix"),
		make_option(c("-b","--all_id"),type="character", default=NULL, help="input id corresponding table"),
		make_option(c("-c","--before_qc"),type="character", default=NULL, help="output before_qc.pdf"),
		make_option(c("-k","--pear_plot"),type="character", default=NULL, help="output pearplot_before_qc.pdf"),
		make_option(c("-d","--after_qc"),type="character", default=NULL, help="output after_qc.pdf"),
		make_option(c("-e","--variable_features"),type="character", default=NULL, help="output VariableFeatures.pdf"),
		make_option(c("-f","--pca"),type="character", default=NULL, help="output pca.pdf"),
		make_option(c("-g","--cell_cluster"),type="character", default=NULL, help="output cell_cluster.csv"),
		make_option(c("-l","--tsne"),type="character", default=NULL, help="output tSNE.pdf"),
		make_option(c("-i","--umap"),type="character", default=NULL, help="output UMAP.pdf"),
		make_option(c("-m","--scRNA_rds"),type="character", default=NULL, help="output UMAP.pdf"),
		make_option(c("-n","--min_gene"),type="character", default=NULL, help="input QC thresold-minGene"),
		make_option(c("-p","--max_gene"),type="character", default=NULL, help="input QC thresold-maxGene"),
		make_option(c("-q","--pct_mt"),type="character", default=NULL, help="input QC thresold-pctMT"),
		make_option(c("-j","--meta_data"),type="character", default=NULL, help="input metadata file"),
		make_option(c("-t","--pca_num"),type="character", default=NULL, help="input pca number"),
		make_option(c("-s","--pct_ercc"),type="character", default=NULL, help="input QC thresold-pctERCC"),
		make_option(c("-r","--resolution_chose"),type="character", default=NULL, help="input resolution_chose")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$finall_fpkm)||is.null(opt$min_gene) || is.null(opt$all_id)|| is.null(opt$max_gene)|| is.null(opt$before_qc)|| 
		is.null(opt$pear_plot)|| is.null(opt$after_qc)|| is.null(opt$variable_features)|| is.null(opt$pca)|| is.null(opt$cell_cluster)|| 
		is.null(opt$tsne)|| is.null(opt$umap)|| is.null(opt$scRNA_rds)||is.null(opt$pct_mt)||is.null(opt$meta_data)||is.null(opt$pct_ercc)||is.null(opt$resolution_chose)){
  print_help(opt_parser);
  stop("Please provide -a finall fpkm, -b all_id ,-c before_qc,-d after_qc,-e VariableFeatures ,-f pca,-g cell_cluster and -l tSNE,
 				-i UMAP,-j meta_data ,-k pearplot,-m scRNA_rds,-n min_gene,-p max_gene,-q pct_mt,-t pca_num,-s pct_ercc,-r resolution_chose", call.=FALSE);
}

Finallfpkm=opt$finall_fpkm
ALLid=opt$all_id
MetaData=opt$meta_data
BeforeQC=opt$before_qc
Pearplot=opt$pear_plot
AfterQC=opt$after_qc
VariableFeatures=opt$variable_features
PCA=opt$pca
CELLcluster=opt$cell_cluster
TSNE=opt$tsne
UMAP=opt$umap
SCrds=opt$scRNA_rds
MinGene=opt$min_gene
MaxGene=opt$max_gene
pctmt=opt$pct_mt
PcaNum=opt$pca_num
pctERCC=opt$pct_ercc
Resolution=opt$resolution_chose
#数据处理
finally_fpkm<-read.table(Finallfpkm,header=T)
all_id<-read.table(ALLid,header=T)
name<-as.data.frame(table(finally_fpkm$Geneid))
name<-filter(name,Freq>1)
if(nrow(name)!=0){finally_fpkm_1<-filter(finally_fpkm,Geneid!=name$Var1)}else {finally_fpkm_1<-finally_fpkm}
finally_count<-finally_fpkm_1[,7:(ncol(finally_fpkm_1)-1)]
rownames(finally_count)<-finally_fpkm_1$Geneid
dim(finally_count)
dat_1<-t(finally_count)
dat_1<-as.data.frame(dat_1)
dat_1$sample<-rownames(dat_1)
dat_1<-as.data.frame(dat_1$sample)
colnames(dat_1)<-"id"
##01创建seurat对象
scRNA<-CreateSeuratObject(finally_count,projiect='scRNA')
table(scRNA@meta.data$orig.ident)
# add metadata to Seurat object 这里的metadata需要sample为列名的一列细胞名称
metadata = read.table(MetaData,header=T,sep='\t')
metadata<-unique(metadata)
rownames(metadata)<-metadata$sample
metadata<-metadata %>% select(sample,everything())
meta_name<-matrix(0,1,ncol(metadata))
for (i in 1:ncol(metadata)){
	tmp=paste("group",i,sep='_')
	meta_name[1,i]=tmp
}
meta_name<-as.character(meta_name)
colnames(metadata)<-meta_name
meatdata<-inner_join(metadata,dat_1,by=c("group_1"="id"))
scRNA <- AddMetaData(object = scRNA, metadata = metadata)
#scRNA = AddMetaData(object = scRNA, percent.ercc, col.name = "percent.ercc")
##02计算质控指标
#计算细胞中线粒体基因比例
#rownames(scRNA)
is.mt <- grep("^mt-", rownames(scRNA))
#rownames(scRNA)[is.mt]
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
scRNA[["percent.ercc"]] <- PercentageFeatureSet(scRNA, pattern = "^ERCC-")
#计算ERCC的比例
#scRNA[["percent.ercc"]] <- PercentageFeatureSet(scRNA, pattern = "^ERCC-")
#head(scRNA@meta.data)
p1 <- VlnPlot(scRNA,
		features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ercc'), 
		pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
		ncol =4) +
		theme(axis.title.x=element_blank(),
					axis.text.x=element_blank(), 
					axis.ticks.x=element_blank())
ggsave(BeforeQC, plot = p1, width = 12, height = 6)
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3<- FeatureScatter(scRNA,feature1 = "nFeature_RNA",feature2 = "percent.mt")
plot4<- FeatureScatter(scRNA,feature1 = "nFeature_RNA",feature2 = "percent.ercc")
pearplot <- CombinePlots(plots = list(plot1, plot2,plot3,plot4), nrow=2, legend="none")
ggsave(Pearplot,plot=pearplot, width = 12, height = 8)
##03设置质控标准
#print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=5000", "pctMT=20"))
pctMT=as.numeric(pctmt)
minGene=as.numeric(MinGene)
maxGene=as.numeric(MaxGene)
pctERCC=as.numeric(pctERCC)
pctMT
minGene
maxGene
pctERCC
scRNA
##04数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.ercc < pctERCC)
scRNA
violin <-VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ercc'), 
                  pt.size = 0.1, 
                  ncol = 4) + 
   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave(AfterQC, plot = violin, width = 12, height = 6) 
##05Normalization
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
#scRNA <- SCTransform(scRNA, vars.to.regress = "percent.mt",verbose = FALSE)
##06寻找高变基因
#dir.create("cluster")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
ggsave(VariableFeatures, plot = plot, width = 8, height = 6) 

##07标准化
scale.genes <-  rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes, vars.to.regress = c("nCount_RNA","percent.ercc"))
# 结果存储在scRNA[["RNA"]]@scale.data中

##08PCA降维
scRNA <- RunPCA(scRNA,features=VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
ggsave(PCA, plot = plotc, width = 8, height = 4)

##09细胞聚类
PcaNum=as.numeric(PcaNum)
PcaNum
pc.num=1:PcaNum
Resolution=as.numeric(Resolution)
Resolution
scRNA <- FindNeighbors(scRNA, dims = pc.num) 
scRNA <- FindClusters(scRNA, resolution=Resolution)
table(scRNA@meta.data$seurat_clusters)
metadata_1 <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata_1), cluster_ID=metadata_1$seurat_clusters)
write.csv(cell_cluster,CELLcluster,row.names = F)

##10降维之TSNE
set.seed(888)
scRNA = RunTSNE(scRNA, dims = pc.num)
embed_tsne <- Embeddings(scRNA, 'tsne')
#write.csv(embed_tsne,'cluster/embed_tsne.csv')
plot1 = DimPlot(scRNA, reduction = "tsne",label = T) 
ggsave(TSNE, plot = plot1, width = 8, height = 7)

##UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
#write.csv(embed_umap,'cluster/embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap",label = T) 
ggsave(UMAP, plot = plot2, width = 8, height = 7)

saveRDS(scRNA,file=SCrds)
