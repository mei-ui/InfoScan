suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("Seurat"))
suppressMessages(library("ggsci"))
suppressMessages(library("rtracklayer"))
option_list <- list(
      make_option(c("-a","--conserved_score"),type="character", default=NULL, help="input conservation Score txt file"),
      make_option(c("-b","--all_id"),type="character", default=NULL, help="input id Correspondence Table"),
      make_option(c("-c","--conservation_density"),type="character", default=NULL, help="output conservation_density.pdf"),
      make_option(c("-d","--conservation_ecdf"),type="character", default=NULL, help="output conservation_ecdf.pdf"),
      make_option(c("-e","--all_gtf"),type="character", default=NULL, help="input all_gene.gtf"),
			make_option(c("-f","--UnannotatedScore"),type="character", default=NULL, help="output unannotated gene conservedscore file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$conserved_score)|| is.null(opt$all_id) || is.null(opt$conservation_density) || is.null(opt$conservation_ecdf)||
    is.null(opt$all_gtf)||is.null(opt$UnannotatedScore)){
  print_help(opt_parser);
  stop("Please provide -a conserved_score, -b all_id ,-c conservation_density,-d conservation_ecdf -e all_gtf -f UnannotatedScore", call.=FALSE);
}
conservedscore=opt$conserved_score
allid=opt$all_id
DENSITY=opt$conservation_density
ECDF=opt$conservation_ecdf
GTF=opt$all_gtf
unannotated_score=opt$UnannotatedScore
#
conservedscore<-read.table(conservedscore)
all_id<-read.table(allid,header=T)
gtf_all<-rtracklayer::import(GTF)
gtf_df<-as.data.frame(gtf_all)
a<-str_match(gtf_df$transcript_id,pattern='TCONS_[0-9]*')
a<-as.data.frame(a)
a<-na.omit(a)
gtf_df_1<-inner_join(gtf_df,a,by=c("transcript_id"="V1"))
gtf_df_1$gene_name<-gtf_df_1$transcript_id
gtf_df_1<-unique(gtf_df_1)

id<-gtf_df[,which(colnames(gtf_df)%in%c("transcript_id","gene_name"))]
id<-na.omit(id)
id<-unique(id)

id_1<-gtf_df_1[,which(colnames(gtf_df_1)%in%c("transcript_id","gene_name"))]
id_1<-na.omit(id_1)
id_1<-unique(id_1)

id_2<-rbind(id,id_1)
conservedscore <-inner_join(conservedscore,id_2,by=c("V1"="transcript_id"))

conservedscore_1<-inner_join(conservedscore,all_id,by=c('gene_name'='id'))
conservedscore_1<-unique(conservedscore_1)
#protein_coding
protein_coding<-conservedscore_1[which(conservedscore_1$type%in%c('protein_coding')),]
protein_coding_1<-protein_coding[,which(colnames(protein_coding)%in%c('V1','V6','type'))]
colnames(protein_coding_1)<-c('transcript_id','conserved_score','type')
#annotated_lncRNA
annotated_lncRNA<-conservedscore_1[which(conservedscore_1$type%in%c('annotated_lncRNA')),]
annotated_lncRNA_1<-annotated_lncRNA[,which(colnames(annotated_lncRNA)%in%c('V1','V6','type'))]
colnames(annotated_lncRNA_1)<-c('transcript_id','conserved_score','type')
#unannotated_lncRNA
unannotated_lncRNA<-conservedscore_1[which(conservedscore_1$type%in%c('unannotated_lncRNA')),]
unannotated_lncRNA_1<-unannotated_lncRNA[,which(colnames(unannotated_lncRNA)%in%c('V1','V6','type'))]
colnames(unannotated_lncRNA_1)<-c('transcript_id','conserved_score','type')
#unannotated_protein
unannotated_protein<-conservedscore_1[which(conservedscore_1$type%in%c('unannotated_coding_transcript')),]
unannotated_protein_1<-unannotated_protein[,which(colnames(unannotated_protein)%in%c('V1','V6','type'))]
colnames(unannotated_protein_1)<-c('transcript_id','conserved_score','type')
df_1<-rbind(unannotated_lncRNA,unannotated_protein)
df_1<-df_1[,which(colnames(df_1)%in%c("gene_name","V1","V6","type"))]
colnames(df_1)<-c("transcript_id","conserved_score","gene_name","type")
write.table(df_1,unannotated_score,quote=F,col.names=T,row.names=F)
df<-rbind(protein_coding_1,annotated_lncRNA_1,unannotated_lncRNA_1,unannotated_protein_1)
ggplot(df,aes(conserved_score,color=type)) + 
  xlab("Conserved Score") + 
  geom_density(alpha=1) + 
	labs(title = "gene conservation")+
  theme_bw()+scale_color_aaas()
ggsave(DENSITY,width = 6, height = 4)

ggplot(df,aes(conserved_score,color=type)) + 
  xlab("Conserved Score") + 
  stat_ecdf(alpha=1) + 
	labs(title = "gene conservation")+
  theme_bw()+scale_color_aaas()
ggsave(ECDF,width = 6, height = 4)

