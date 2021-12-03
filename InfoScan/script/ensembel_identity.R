
suppressMessages(library("rtracklayer"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("tidyverse"))
suppressMessages(library("optparse"))
option_list <- list(
    make_option(c("-i","--gtf_file"),type="character", default=NULL,help="输入gtf文件"),
    make_option(c("-o","--outfile_prefix"),type="character", default=NULL, help="输出含有id对应关系的txt文件"),
    make_option(c("-O","--unannotated_gene_id"),type="character", default=NULL, help="unannotated_gene_id")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gtf_file) || is.null(opt$outfile_prefix) || is.null(opt$unannotated_gene_id) ){
  print_help(opt_parser);
  stop("Please provide -i gtf_file,-O unannotated_gene_id and -o outfile_prefix option", call.=FALSE);
}
gtffile=opt$gtf_file
outfile=opt$outfile_prefix
unannotatedid=opt$unannotated_gene_id
gtf1 <- rtracklayer::import(gtffile)
gtf_df <- as.data.frame(gtf1)
lncRNA_1<-gtf_df
id<-lncRNA_1[,which(colnames(lncRNA_1)%in%c("gene_id","transcript_id","gene_name","oId"))]
id<-unique(id)
id_1<-na.omit(id)
unannotated<-id[-which(c(id$transcript_id) %in% c(id_1$transcript_id)),]
unknow_gene_id<-unannotated$transcript_id
write.table(id,outfile,quote=F,row.names=F)
write.table(unknow_gene_id,unannotatedid,quote=F,col.names=F,row.names=F)
