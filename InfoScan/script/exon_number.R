suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggsci"))

option_list <- list(
		make_option(c("-a","--fpkm_matrix"),type="character", default=NULL, help="input finally fpkm"),
		make_option(c("-c","--all_gtf"),type="character", default=NULL, help="input all.gtf"),
		make_option(c("-d","--protein_gtf"),type="character", default=NULL, help="input protein_coding.gtf"),
		make_option(c("-e","--annotated_lncRNA"),type="character", default=NULL, help="input annotated_lncRNA.gtf"),
		make_option(c("-f","--unannotated_lncRNA"),type="character", default=NULL, help="input unannotated_lncRNA.gtf"),
		make_option(c("-g","--unprotein"),type="character", default=NULL, help="input unannotated_protein.gtf"),
    make_option(c("-q","--exon_filter"),type="character", default=NULL, help="input exon number filter"),
		make_option(c("-o","--exon_number"),type="character", default=NULL, help="input exon_number.pdf"),
    make_option(c("-O","--finally_fpkm"),type="character", default=NULL, help="output exon number filter result")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$fpkm_matrix)|| is.null(opt$all_gtf) || is.null(opt$exon_number)||is.null(opt$protein_gtf)||
    is.null(opt$annotated_lncRNA)|| is.null(opt$unannotated_lncRNA)|| is.null(opt$unprotein)||is.null(opt$exon_filter)||is.null(opt$finally_fpkm)){
  print_help(opt_parser);
  stop("Please provide -a fpkm_matrix,-c all_gtf and -d protein_gtf,-e annotated_lncRNA,-f unannotated_lncRNA ,-g unprotein,-q exon_filter, -o exon_number,-O finally_fpkm ", call.=FALSE);
}

fpkmMatrix=opt$fpkm_matrix
allgtf=opt$all_gtf
filterExon=opt$exon_filter
exonnumber=opt$exon_number
proteingtf=opt$protein_gtf
annotatedgtf=opt$annotated_lncRNA
unannotatedlncgtf=opt$unannotated_lncRNA
unannotatedproteingtf=opt$unprotein
FPKM=opt$finally_fpkm

#Read-in data
finally_fpkm <- read.table(fpkmMatrix,header=T)
finally_count_1<-finally_fpkm 

gtf_all<-rtracklayer::import(allgtf)
gtf_all<-as.data.frame(gtf_all)


protein_coding<-finally_count_1[which(finally_count_1$type%in%c('protein_coding')),]
protein<-as.data.frame(protein_coding$Geneid)
colnames(protein)<-'gene'
protein_gtf<-inner_join(gtf_all,protein,by=c('gene_name'='gene'))
protein_gtf<-protein_gtf[,-which(colnames(protein_gtf)%in%c('ref_gene_id'))]
protein_gtf_1<-filter(protein_gtf,type=="exon")
protein_exon_number<-as.data.frame(table(protein_gtf_1$transcript_id))
protein_exon<-filter(protein_exon_number,Freq> 0)
protein_exon_id<-as.data.frame(protein_exon$Var1)
colnames(protein_exon_id)<-'id'
protein_gtf_2<-inner_join(protein_gtf,protein_exon_id,by=c('transcript_id'='id'))
export(protein_gtf_2,proteingtf,'gtf')
protein_exon$group<-'protein_coding'

annotated_lncRNA<-finally_count_1[which(finally_count_1$type%in%c('annotated_lncRNA')),]
annotated<-as.data.frame(annotated_lncRNA$Geneid)
colnames(annotated)<-'gene'
annotated_gtf<-inner_join(gtf_all,annotated,by=c('gene_name'='gene'))
annotated_gtf<-annotated_gtf[,-which(colnames(annotated_gtf)%in%c('ref_gene_id'))]
annotated_gtf_1<-filter(annotated_gtf,type=="exon")
annotated_exon_number<-as.data.frame(table(annotated_gtf_1$transcript_id))
annotated_exon<-filter(annotated_exon_number,Freq> 0)
annotated_exon_id<-as.data.frame(annotated_exon$Var1)
colnames(annotated_exon_id)<-'id'
annotated_gtf_2<-inner_join(annotated_gtf,annotated_exon_id,by=c('transcript_id'='id'))
annotated_exon$group<-'annotated_lncRNA'
export(annotated_gtf_2,annotatedgtf,'gtf')

unannotated_lncRNA<-finally_count_1[which(finally_count_1$type%in%c('unannotated_lncRNA')),]
unannotated<-as.data.frame(unannotated_lncRNA$Geneid)
colnames(unannotated)<-'gene'
unannotated_lncRNA_gtf<-inner_join(gtf_all,unannotated,by=c('transcript_id'='gene'))
unannotated_lncRNA_gtf_1<-filter(unannotated_lncRNA_gtf,type=="exon")
unannotated_lncRNA_exon_number<-as.data.frame(table(unannotated_lncRNA_gtf_1$transcript_id))
unannotated_lncRNA_exon<-filter(unannotated_lncRNA_exon_number,Freq> filterExon)
unannotated_lncRNA_exon_id<-as.data.frame(unannotated_lncRNA_exon$Var1)
colnames(unannotated_lncRNA_exon_id)<-'id'
unannotated_lncRNA_gtf_2<-inner_join(unannotated_lncRNA_gtf,unannotated_lncRNA_exon_id,by=c('transcript_id'='id'))
unannotated_lncRNA_exon$group<-'unannotated_lncRNA'
export(unannotated_lncRNA_gtf_2,unannotatedlncgtf,'gtf')

unannotated_protein<-finally_count_1[which(finally_count_1$type%in%c('unannotated_coding_transcript')),]
unannotated_protein<-as.data.frame(unannotated_protein$Geneid)
colnames(unannotated_protein)<-'gene'
unannotated_protein_gtf<-inner_join(gtf_all,unannotated_protein,by=c('transcript_id'='gene'))
unannotated_protein_gtf_1<-filter(unannotated_protein_gtf,type=="exon")
unannotated_protein_exon_number<-as.data.frame(table(unannotated_protein_gtf_1$transcript_id))
unannotated_protein_exon<-filter(unannotated_protein_exon_number,Freq> 0)
unannotated_protein_exon_id<-as.data.frame(unannotated_protein_exon$Var1)
colnames(unannotated_protein_exon_id)<-'id'
unannotated_protein_gtf_2<-inner_join(unannotated_protein_gtf,unannotated_protein_exon_id,by=c('transcript_id'='id'))
unannotated_protein_exon$group<-'unannotated_protein'
export(unannotated_protein_gtf_2,unannotatedproteingtf,'gtf')
exon_number<-rbind(protein_exon,annotated_exon,unannotated_lncRNA_exon,unannotated_protein_exon)
dup<-subset(exon_number,exon_number$Freq<50)

annotated_gtf<-rbind(annotated_gtf_2,protein_gtf_2)
annotated_id<-as.data.frame(unique(annotated_gtf$gene_name))
colnames(annotated_id)<-'id'
annotated_fpkm<-inner_join(finally_fpkm,annotated_id,by=c("Geneid"="id"))
unannotated_id =rbind(unannotated_lncRNA_exon_id,unannotated_protein_exon_id)
unannotated_fpkm<-inner_join(finally_fpkm,unannotated_id,by=c("Geneid"="id"))
finally_fpkm_1<-rbind(annotated_fpkm,unannotated_fpkm)
ercc.index = grep(pattern = "^ERCC-", x = finally_fpkm$Geneid, value = FALSE)
ercc_fpkm<-finally_fpkm[ercc.index,]
finally_fpkm_2<-rbind(finally_fpkm_1,ercc_fpkm)
write.table(finally_fpkm_2,FPKM,sep='\t',quote=F,row.names=F)

ggplot(dup,mapping=aes(x=Freq,fill = group))+
  geom_bar()+
  labs(title = "exon number")+
  labs(fill="gene-type")+#修改图例的名字
  scale_fill_aaas()+
  scale_x_continuous(breaks = seq(1, 50, by = 3))+
  ylab('Total Transcript Number')+#修改y轴名称
  xlab('exon number')+#修改x轴名称
  theme_bw()
ggsave(exonnumber,width = 8,height = 4)
