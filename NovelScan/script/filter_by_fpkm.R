suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggsci"))
option_list <- list(
			make_option(c("-i","--annotated_count"),type="character", default=NULL, help="input annotated_count.txt"),
			make_option(c("-I","--unannotated_count"),type="character", default=NULL, help="input unannotated_count.txt"),
      make_option(c("-q","--fpkm_threshold"),type="character", default=NULL, help="input fpkm threshold"),
      make_option(c("-b","--all_id"),type="character", default=NULL, help="input all id txt file"),
      make_option(c("-o","--fpkm_matrix"),type="character", default=NULL, help="output fpkm expression file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$annotated_count)|| is.null(opt$unannotated_count) || is.null(opt$fpkm_matrix)||is.null(opt$all_id) 
  ||is.null(opt$fpkm_threshold)){
  print_help(opt_parser);
  stop("Please provide -o fpkm_matrix, -b all_id ,-i annotated_count,-I unannotated_count ,-q fpkm_threshold", call.=FALSE);
}

annoCount=opt$annotated_count
unannoCount=opt$unannotated_count
filter_q=opt$fpkm_threshold
allid=opt$all_id
FPKM=opt$fpkm_matrix
#filter_before=opt$before_filter
#filter_after=opt$after_filter

#读入数据
#anno_count
anno_count<-read.table(annoCount,header=T,quote='\t',skip=1)
sampleName<-colnames(anno_count[,7:ncol(anno_count)])
sampleName<-str_remove_all(sampleName,pattern="\\.bam")
sampleName<-str_remove_all(sampleName,pattern="^.*3_HISAT2_aligned.")
names(anno_count)[7:ncol(anno_count)] <- sampleName
#unanno_count
unanno_count<-read.table(unannoCount,header=T,quote='\t',skip=1)
sampleName<-colnames(unanno_count[,7:ncol(unanno_count)])
sampleName<-str_remove_all(sampleName,pattern="\\.bam")
sampleName<-str_remove_all(sampleName,pattern="^.*3_HISAT2_aligned.")
names(unanno_count)[7:ncol(unanno_count)] <- sampleName
unanno_count<-as.data.frame(unanno_count)
unanno_count$Start<-as.character(unanno_count$Start)
unanno_count$End<-as.character(unanno_count$End)
#all id
all_id<-read.table(allid,header = F)
colnames(all_id) <- c("id","type")

finally_count<- rbind(anno_count,unanno_count)
finally_count<-unique(finally_count)
#数据预处理
finally_eset<-finally_count[,7:ncol(finally_count)]
finally_eset<-as.data.frame(finally_eset)
rownames(finally_eset)<-finally_count$Geneid
dat=finally_eset
dat_1<-dat
dat_1$gene<-row.names(dat_1)
#read count 转 FPKM
finally_id<-finally_count[,1:6]
dat_1<-inner_join(dat_1,finally_id,by=c('gene'='Geneid'))
finally_count_1<-inner_join(dat_1,all_id,by=c('gene'='id'))
kb<-finally_count$Length / 1000
rpk <- dat / kb
fpkm <- t(t(rpk)/colSums(dat) * 10^6)
#head(fpkm)
fpkm<-as.data.frame(fpkm) 
fpkm<-fpkm[apply(fpkm,1, function(x) sum(x>=1) > filter_q),] 
fpkm$gene_id<-row.names(fpkm)
finally_fpkm<-inner_join(finally_id,fpkm,by=c('Geneid'='gene_id'))
finally_fpkm<-inner_join(finally_fpkm,all_id,by=c('Geneid'='id'))
write.table(finally_fpkm,FPKM,sep='\t',quote=F,row.names=F)