suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggsci"))
suppressMessages(library('cummeRbund'))
suppressMessages(library('pheatmap'))

option_list <- list(
        make_option(c("-a","--finall_fpkm"),type="character", default=NULL, help="input finally_fpkm file"),
        make_option(c("-i","--protein_coding_heatmap"),type="character", default=NULL, help="output protein_coding_heatmap.gtf file"),
        make_option(c("-I","--annotated_lncRNA_heatmap"),type="character", default=NULL, help="annotated_lncRNA_heatmap.pdf file"),
        make_option(c("-c","--celltype_rds"),type="character", default=NULL, help="input scRNA celltype rds file"),
	    make_option(c("-o","--unannotated_lncRNA_heatmap"),type="character", default=NULL, help="unannotated_lncRNA_heatmap.pdf file"),
	    make_option(c("-O","--unannotated_protein_heatmap"),type="character", default=NULL, help="unannotated_protein_heatmap.pdf file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$finall_fpkm)|| is.null(opt$protein_coding_heatmap)||is.null(opt$celltype_rds)||is.null(opt$annotated_lncRNA_heatmap)||
    is.null(opt$unannotated_lncRNA_heatmap)||is.null(opt$unannotated_protein_heatmap)){
  print_help(opt_parser);
  stop("Please provide -a finall_fpkm ,-i protein_coding_heatmap and -I annotated_lncRNA_heatmap,-c celltype_rds,-o unannotated_lncRNA_heatmap ,
        -O unannotated_protein_heatmap", call.=FALSE);
}

finallfpkm=opt$finall_fpkm
PCH=opt$protein_coding_heatmap
ALH=opt$annotated_lncRNA_heatmap
CR=opt$celltype_rds
ULH=opt$unannotated_lncRNA_heatmap
UPH=opt$unannotated_protein_heatmap

##加载数据

##load data
finally_fpkm <- read.table(finallfpkm,header=T)
name<-as.data.frame(table(finally_fpkm$Geneid))
name<-filter(name,Freq>1)
finally_fpkm_1<-filter(finally_fpkm,Geneid!=name$Var1)
scRNA<-readRDS(CR)
celltype<-as.data.frame(scRNA$celltype)
colnames(celltype)<-'celltype'
celltype$sample<-rownames(celltype)
fpkm<-finally_fpkm_1[,7:(ncol(finally_fpkm_1)-1)]
rownames(fpkm)<-finally_fpkm_1$Geneid
fpkm_1<-as.data.frame(t(fpkm))
fpkm_1$sample<-rownames(fpkm_1)
fpkm_1<-inner_join(fpkm_1,celltype,by=c('sample'='sample'))
fpkm_group = fpkm_1 %>% group_by(celltype) 
n<-nrow(table(fpkm_1$celltype))
x<-group_split(fpkm_group)
data<-x[[1]]
data<-as.data.frame(data)                               
data_eset<-data[,1:(ncol(data)-2)]
data_1<-apply(data_eset,2,mean)
name<-as.character(unique(data$celltype))
data_1<-as.data.frame(data_1)
colnames(data_1)<-name
for (i in 2:as.numeric(n)) {
            data_2<-x[[i]]
            data_2<-as.data.frame(data_2)
            data_eset_2<-data_2[,1:(ncol(data_2)-2)]
            data_3<-apply(data_eset_2,2,mean)
            name_2<-as.character(unique(data_2$celltype))
            data_3<-as.data.frame(data_3)
            colnames(data_3)<-name_2
            data_1<-cbind(data_1,data_3)
}
fpkm_2<-data_1
fpkm_2$gene_id<-row.names(fpkm_2)
id<-finally_fpkm_1[,which(colnames(finally_fpkm_1)%in%c("Geneid","type"))]
finally_fpkm_3<-inner_join(id,fpkm_2,by=c('Geneid'='gene_id'))
#protein_coding
protein_coding<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('protein_coding')),]
protein_coding<-as.data.frame(protein_coding)
protein_coding_eset<-protein_coding[,-which(colnames(protein_coding)%in%c("Geneid","type"))]
rownames(protein_coding_eset)<-protein_coding$Geneid
protein_coding_eset<- protein_coding_eset[apply(protein_coding_eset, 1, function(x) sd(x)!=0),]
pheatmap(protein_coding_eset, 
         show_rownames = F, 
         show_colnames = TRUE, 
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = F,
         border_color = NA,
         main = "protein-coding-heatmap",
         filename=PCH)
#annotated lncRNA 
annotated_lncRNA<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('annotated_lncRNA')),]
annotated_lncRNA<-as.data.frame(annotated_lncRNA)
annotated_lncRNA_eset<-annotated_lncRNA[,-which(colnames(annotated_lncRNA)%in%c("Geneid","type"))]
rownames(annotated_lncRNA_eset)<-annotated_lncRNA$Geneid
annotated_lncRNA_eset<- annotated_lncRNA_eset[apply(annotated_lncRNA_eset, 1, function(x) sd(x)!=0),]
pheatmap(annotated_lncRNA_eset, 
         show_rownames = F, 
         show_colnames = TRUE, 
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = F,
         border_color = NA,
         main = "annotated-lncRNA-heatmap",
         filename=ALH)
#unannotated lncRNA
unannotated_lncRNA<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('unannotated_lncRNA')),]
unannotated_lncRNA<-as.data.frame(unannotated_lncRNA)
unannotated_lncRNA_eset<-unannotated_lncRNA[,-which(colnames(unannotated_lncRNA)%in%c("Geneid","type"))]
rownames(unannotated_lncRNA_eset)<-unannotated_lncRNA$Geneid
unannotated_lncRNA_eset<- unannotated_lncRNA_eset[apply(unannotated_lncRNA_eset, 1, function(x) sd(x)!=0),]
pheatmap(unannotated_lncRNA_eset, 
         show_rownames = F, 
         show_colnames = TRUE, 
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = F,
         border_color = NA,
         main = "unannotated-lncRNA-heatmap",
         filename=ULH)

#unannotated protein coding 
unannotated_protein<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('unannotated_coding_transcript')),]
unannotated_protein<-as.data.frame(unannotated_protein)
unannotated_protein_eset<-unannotated_protein[,-which(colnames(unannotated_protein)%in%c("Geneid","type"))]
rownames(unannotated_protein_eset)<-unannotated_protein$Geneid
unannotated_protein_eset<- unannotated_protein_eset[apply(unannotated_protein_eset, 1, function(x) sd(x)!=0),]
pheatmap(unannotated_protein_eset, 
         show_rownames = F, 
         show_colnames = TRUE, 
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = F,
         border_color = NA,
         main = "unannotated-coding-transcript-heatmap",
         filename=UPH)