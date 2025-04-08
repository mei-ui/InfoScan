suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))
option_list <- list(
                make_option(c("-a","--finall_fpkm"),type="character", default=NULL, help="input finally fpkm expression matrix"),
		            make_option(c("-c","--lncRNA_id"),type="character", default=NULL, help="input lncRNA id"),
                make_option(c("-i","--co_expression"),type="character", default=NULL, help="output fpkm expression result"),
                make_option(c("-I","--row_number"),type="character", default=NULL, help="output  the line of unannotated gene in the file") 
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$finall_fpkm)|| is.null(opt$lncRNA_id)||is.null(opt$co_expression)|| is.null(opt$row_number)){
  print_help(opt_parser);
  stop("Please provide -a finall_fpkm ,-c unannotated_marker ,-i co_expression ,-I row_number", call.=FALSE);
}
finallfpkm=opt$finall_fpkm
Unannotated_marker=opt$lncRNA_id
COexpression=opt$co_expression
Rownumber=opt$row_number
#un_ConservedScore=opt$unannotated_conservedscore
#conservedscore=opt$conserved_score
##加载数据
finally_fpkm<-read.table(finallfpkm,header=T)
name<-as.data.frame(table(finally_fpkm$Geneid))
name<-filter(name,Freq>1)
if(nrow(name)!=0){finally_fpkm_1<-finally_fpkm[-which(finally_fpkm$Geneid%in%c(as.character(name$Var1))),]}else {finally_fpkm_1<-finally_fpkm}
finally_count<-finally_fpkm_1[,7:(ncol(finally_fpkm_1)-1)]
rownames(finally_count)<-finally_fpkm_1$Geneid
fpkm=finally_count
#计算ERCC比例，并去掉外源基因ERCC
ercc.index = grep(pattern = "^ERCC-", x = rownames(x = fpkm), value = FALSE) # 获取 index
# 删除 fpkm 矩阵中的 ERCC
fpkm = fpkm[-ercc.index,]
dim(fpkm)
rownames(fpkm)<-stri_replace_all_regex(rownames(fpkm),'_','-')
fpkm$geneID<-row.names(fpkm)
fpkm$geneName<-row.names(fpkm)
fpkm<-fpkm %>% select(geneID,geneName,everything())
write.table(fpkm,COexpression,col.names=F,quote=F,row.names=F)

unannotated_marker<-as.character(Unannotated_marker)
#move_to_first <- function(df, n) df[c(n,setdiff(seq_len(nrow(df)), n)), ]

#for (i in 1:length(unannotated_marker_id)){
    #j=1
#	  row.number<-which(fpkm$geneID==unannotated_marker_id[i])
#		fpkm_1<-move_to_first(fpkm,row.number)
#}
row.number<- which(fpkm$geneID%in%c(unannotated_marker))
row.number<-as.data.frame(row.number)
row.number<-row.number-1
write.table(row.number,Rownumber,col.names=F,quote=F,row.names=F)

