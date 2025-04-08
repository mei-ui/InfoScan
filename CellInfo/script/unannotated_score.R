suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("Seurat"))
suppressMessages(library("patchwork"))

option_list <- list(
		make_option(c("-a","--scRNA"),type="character", default=NULL, help="input scRNA RDS"),
		make_option(c("-i","--unannotated_marker"),type="character", default=NULL, help="input unannotated_marker.csv file"),
		make_option(c("-I","--unannotated_score"),type="character", default=NULL,help="input unannotated conservation score file"),
		make_option(c("-k","--unannotated_conservation"),type="character", default=NULL, help="output unannotated_marker_score.csv file"),
		make_option(c("-f","--fpkm_filter"),type="character", default=NULL, help="input fpkm filter threshold"),
		make_option(c("-d","--unannotated_lncRNA_vloplot"),type="character", default=NULL, help="output unannotated_lncRNA_vloplot.pdf file"),
		make_option(c("-e","--unannotated_lncRNA_featureplot"),type="character", default=NULL, help="output unannotated_lncRNA_featureplot.pdf file"),
		make_option(c("-c","--conservation_filter"),type="character", default=NULL, help="input conservation filter threshold")
) 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$unannotated_marker)|| is.null(opt$unannotated_score) || is.null(opt$unannotated_conservation)||
	is.null(opt$fpkm_filter)||is.null(opt$conservation_filter)||is.null(opt$scRNA)|| is.null(opt$unannotated_lncRNA_vloplot)
	||is.null(opt$unannotated_lncRNA_featureplot)){
  print_help(opt_parser);
  stop("Please provide -a scRNA.rds,-i unannotated_marker, -I unannotated_score ,-k unannotated_conservation,-f fpkm_filter,-c conservation_filter
  		-d unannotated_lncRNA_vloplot,-e unannotated_lncRNA_featureplot", call.=FALSE);
}
sc_RNA=opt$scRNA
marker=opt$unannotated_marker
score=opt$unannotated_score
conservation=opt$unannotated_conservation
FPKMfilter=opt$fpkm_filter
CONSERVATION=opt$conservation_filter
Featureplot=opt$unannotated_lncRNA_featureplot
VLOPLOT=opt$unannotated_lncRNA_vloplot

SCRNA<-readRDS(sc_RNA)
table1<-read.csv(marker)
table2<-read.table(score,header=T)
table2$transcript_id<-stri_replace_all_regex(table2$transcript_id,'_','-')
table2<-table2[,-which(colnames(table2)%in%c("type"))]
table3<-inner_join(table1,table2,by=c("Geneid"="transcript_id"))
CONSERVATION
table3<-filter(table3,conserved_score>CONSERVATION)
FPKMfilter
table3<-filter(table3,TPM>FPKMfilter)
write.csv(table3,conservation,row.names = F)

lncRNA_marker_gene_1 <- c(table3$Geneid)
lncRNA_marker_gene_1<-unique(lncRNA_marker_gene_1)
#vlnplot展示
if(length(lncRNA_marker_gene_1)<=10){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=8)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=8)

}
if(10<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=20){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=12)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=12)

}
if(20<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=30){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=16)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=16)

}
if(30<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=40){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=20)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=20)

}
if(40<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=50){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=24)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=24)

}
if(50<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=60){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=28)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=28)

}
if(60<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=70){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=32)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=32)

}
if(70<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=80){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=36)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=36)

}
if(80<length(lncRNA_marker_gene_1)&length(lncRNA_marker_gene_1)<=90){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=36)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=40)

}
if(90<length(lncRNA_marker_gene_1)){
	p1<-VlnPlot(SCRNA,features = lncRNA_marker_gene_1,
            pt.size = 0.1,group.by = "seurat_clusters",ncol = 5)
	ggsave(VLOPLOT, p1, width=16,height=100,limitsize = FALSE)
	p2 <- FeaturePlot(SCRNA,features = lncRNA_marker_gene_1,
                  reduction = "tsne",label = T,ncol =5)
	ggsave(Featureplot, p2, width=18 ,height=100,limitsize = FALSE)

}