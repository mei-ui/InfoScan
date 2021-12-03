suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))

option_list <- list(
		make_option(c("-a","--scRNA_rds"),type="character", default=NULL, help="input scRNA_rds file "),
 		make_option(c("-b","--diff_genes_wilcox"),type="character", default=NULL, help="output diff_genes_wilcox.csv"),
		make_option(c("-c","--top10_diff_genes_wilcox"),type="character", default=NULL, help="output top10_diff_genes_wilcox.csv"),
 		make_option(c("-d","--unannotated_lncRNA_vloplot"),type="character", default=NULL, help="output unannotated_lncRNA_vloplot.pdf file"),
		make_option(c("-e","--unannotated_lncRNA_featureplot"),type="character", default=NULL, help="output unannotated_lncRNA_featureplot.pdf file"),
		make_option(c("-f","--unannotated_lncRNA"),type="character", default=NULL, help="input unannotated_lncRNA.csv file"),
		make_option(c("-g","--unannotated_lncRNA_marker"),type="character", default=NULL, help="output unannotated_lncRNA_marker.csv file"),
		make_option(c("-i","--top10_marker"),type="character", default=NULL,help="output top10 marker gene heatmap file"),
		make_option(c("-k","--unannotated_protein"),type="character", default=NULL, help="input unannotated_ptotein.csv file")
) 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$scRNA_rds)|| is.null(opt$diff_genes_wilcox) || is.null(opt$top10_diff_genes_wilcox)|| is.null(opt$unannotated_lncRNA_vloplot)||
		is.null(opt$unannotated_lncRNA_featureplot) || is.null(opt$unannotated_lncRNA)|| is.null(opt$unannotated_lncRNA_marker)||
		is.null(opt$top10_marker)||is.null(opt$unannotated_protein)){
  print_help(opt_parser);
  stop("Please provide -a scRNA_rds, -b diff_genes_wilcox ,-c top10_diff_genes_wilcox,-d unannotated_lncRNA_vloplot，
  	-e unannotated_lncRNA_featureplot ,-f unannotated_lncRNA,-g unannotated_lncRNA_marker ,-i top10_marker ", call.=FALSE);
}
SCRNA=opt$scRNA_rds
DIFFGENE=opt$diff_genes_wilcox
TOP10=opt$top10_diff_genes_wilcox
VLOPLOT=opt$unannotated_lncRNA_vloplot
Featureplot=opt$unannotated_lncRNA_featureplot
UNANNOTATED=opt$unannotated_lncRNA
UNANNOTATED_PROTEIN=opt$unannotated_protein
LNCRNA=opt$unannotated_lncRNA_marker
TOP10MARKER=opt$top10_marker

##11 Cluster 差异基因
scRNA<-readRDS(SCRNA)
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_genes = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA)) 
plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave(TOP10MARKER, plot=plot1, width=15, height=18) 
write.csv(all.markers, DIFFGENE, row.names = F)
write.csv(top10, TOP10, row.names = F)
scRNA
a<-str_match(all.markers$gene,pattern='TCONS-[0-9]*')
a<-as.data.frame(a)
a<-na.omit(a)
lncRNA_marker<-all.markers[which(all.markers$gene %in% c(as.character(a$V1))),]
unannotated_lncRNA<-read.csv(UNANNOTATED)
unannotated_protein<-read.csv(UNANNOTATED_PROTEIN)
unannotated_lncRNA<-rbind(unannotated_lncRNA,unannotated_protein)
unannotated_lncRNA$Geneid<-stri_replace_all_regex(unannotated_lncRNA$Geneid,'_','-')
unannotated_lncRNA_eset<-unannotated_lncRNA[,7:(ncol(unannotated_lncRNA)-1)]
unannotated_lncRNA_eset$FPKM<-apply(unannotated_lncRNA_eset,1,sum)
unannotated_lncRNA$FPKM<-unannotated_lncRNA_eset$FPKM
lncRNA<-inner_join(unannotated_lncRNA,lncRNA_marker,by=c('Geneid'='gene'))
lncRNA<-lncRNA[,which(colnames(lncRNA)%in%c('Geneid','Chr','Start','End','Strand','Length','type','p_val','avg_logFC','pct.1','pct.2','p_val_adj','cluster','FPKM'))]
write.csv(lncRNA,LNCRNA,row.names = F)
scRNA
#####unannotated lncRNA
lncRNA_marker_gene<-lncRNA_marker %>% select(gene, everything()) %>% subset(p_val<0.05)
lncRNA_marker_gene_1 <- c(lncRNA_marker_gene$gene)
lncRNA_marker_gene_1<-unique(lncRNA_marker_gene_1)
#vlnplot展示
if(length(lncRNA_marker_gene_1)<=10){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=8)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=8)

}
if(10<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=20){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=12)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=12)

}
if(20<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=30){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=16)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=16)

}
if(30<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=40){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=20)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=20)

}
if(40<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=50){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=24)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=24)

}
if(50<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=60){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=28)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=28)

}
if(60<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=70){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=32)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=32)

}
if(70<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=80){
	p1<-VlnPlot(scRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=36)
	p2 <- FeaturePlot(scRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=36)

}
