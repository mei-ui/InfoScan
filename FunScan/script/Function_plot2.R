suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("optparse"))
suppressMessages(library("ggsci"))

option_list <- list(
    make_option(c("-i","--id"),type="character", default=NULL, help="input gene id"),
    make_option(c("-o","--output"),type="character", default=NULL, help="output result"),
    make_option(c("-d","--outputdir"),type="character", default=NULL, help="input output dir")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$id)|| is.null(opt$output)||is.null(opt$outputdir)){
  print_help(opt_parser);
  stop("Please provide -i gene id, -o output file", call.=FALSE);
}

gene_id=opt$id
Out=opt$output
OutDir=opt$outputdir
gocc=paste(OutDir,"/snakemake/result/11_data_analysis/coexpression/",gene_id,"/GO_CC.txt",sep="")
data1<-read.table(gocc,header=T)
if (nrow(data1)!=0){
  tmp<-arrange(data1,data1$`log10.qval.`)
  tmp$number <- factor(rev(1:nrow(tmp)))
  tmp$termName=tolower(tmp$termName)
  tmp$`log10.qval.`=round(tmp$`log10.qval.`,2)
  #tmp$termName=str_to_title(tmp$termName)
  #tmp<-arrange(tmp,desc(tmp$commonGeneNum))
  tmp$type='GO_CC'
  if(nrow(tmp)<=10){tmp=tmp}else{tmp<-tmp[1:10,]}
  go_cc=tmp
}

gobp=paste(OutDir,"/snakemake/result/11_data_analysis/coexpression/",gene_id,"/GO_BP.txt",sep="")
data2<-read.table(gobp,header=T)
if (nrow(data2)!=0){
  tmp<-arrange(data2,data2$`log10.qval.`)
  tmp$number <- factor(rev(1:nrow(tmp)))
  tmp$termName=tolower(tmp$termName)
  tmp$`log10.qval.`=round(tmp$`log10.qval.`,2)
  #tmp$termName=str_to_title(tmp$termName)
  #tmp<-arrange(tmp,desc(tmp$commonGeneNum))
  tmp$type='GO_BP'
  if(nrow(tmp)<=10){tmp=tmp}else{tmp<-tmp[1:10,]}
  go_bp=tmp
}

gomf=paste(OutDir,"/snakemake/result/11_data_analysis/coexpression/",gene_id,"/GO_MF.txt",sep="")
data3<-read.table(gomf,header=T)
if (nrow(data3)!=0){
  tmp<-arrange(data3,data3$`log10.qval.`)
  tmp$number <- factor(rev(1:nrow(tmp)))
  tmp$termName=tolower(tmp$termName)
  tmp$`log10.qval.`=round(tmp$`log10.qval.`,2)
  #tmp$termName=str_to_title(tmp$termName)
  #tmp<-arrange(tmp,desc(tmp$commonGeneNum))
  tmp$type='GO_MF'
  if(nrow(tmp)<=10){tmp=tmp}else{tmp<-tmp[1:10,]}
  go_mf=tmp
}

kegg=paste(OutDir,"/snakemake/result/11_data_analysis/coexpression/",gene_id,"/KEGG.txt",sep="")
data4<-read.table(kegg,header=T)
if (nrow(data4)!=0){
  tmp<-arrange(data4,data4$`log10.qval.`)
  tmp$number <- factor(rev(1:nrow(tmp)))
  tmp$termName=tolower(tmp$termName)
  tmp$`log10.qval.`=round(tmp$`log10.qval.`,2)
  #tmp$termName=str_to_title(tmp$termName)
  #tmp<-arrange(tmp,desc(tmp$commonGeneNum))
  tmp$type='KEGG'
  if(nrow(tmp)<=10){tmp=tmp}else{tmp<-tmp[1:10,]}
  kegg=tmp
}

if (nrow(data1)!=0 && nrow(data2)!=0 && nrow(data3)!=0 && nrow(data4)!=0) {data_plot=rbind(go_cc,go_bp,go_mf,kegg)}
if (nrow(data1)!=0 && nrow(data2)!=0 && nrow(data3)!=0 && nrow(data4)==0) {data_plot=rbind(go_cc,go_bp,go_mf)}
if (nrow(data1)==0 && nrow(data2)!=0 && nrow(data3)!=0 && nrow(data4)!=0) {data_plot=rbind(go_bp,go_mf,kegg)}
if (nrow(data1)==0 && nrow(data2)==0 && nrow(data3)!=0 && nrow(data4)!=0) {data_plot=rbind(go_mf,kegg)}
if (nrow(data1)==0 && nrow(data2)==0 && nrow(data3)==0 && nrow(data4)!=0) {data_plot=kegg}
if (nrow(data1)!=0 && nrow(data2)==0 && nrow(data3)==0 && nrow(data4)!=0) {data_plot=rbind(go_cc,kegg)}
if (nrow(data1)!=0 && nrow(data2)!=0 && nrow(data3)==0 && nrow(data4)!=0) {data_plot=rbind(go_cc,go_bp,kegg)}
if (nrow(data1)!=0 && nrow(data2)==0 && nrow(data3)!=0 && nrow(data4)!=0) {data_plot=rbind(go_cc,go_mf,kegg)}
if (nrow(data1)!=0 && nrow(data2)!=0 && nrow(data3)==0 && nrow(data4)==0) {data_plot=rbind(go_cc,go_bp)}
if (nrow(data1)!=0 && nrow(data2)==0 && nrow(data3)!=0 && nrow(data4)==0) {data_plot=rbind(go_cc,go_mf)}
if (nrow(data1)==0 && nrow(data2)!=0 && nrow(data3)!=0 && nrow(data4)==0) {data_plot=rbind(go_bp,go_mf)}
if (nrow(data1)==0 && nrow(data2)!=0 && nrow(data3)==0 && nrow(data4)!=0) {data_plot=rbind(go_bp,kegg)}
if (nrow(data1)!=0 && nrow(data2)==0 && nrow(data3)==0 && nrow(data4)==0) {data_plot=go_cc}
if (nrow(data1)==0 && nrow(data2)!=0 && nrow(data3)==0 && nrow(data4)==0) {data_plot=go_bp}
if (nrow(data1)==0 && nrow(data2)==0 && nrow(data3)!=0 && nrow(data4)==0) {data_plot=go_mf}
if (nrow(data1)==0 && nrow(data2)==0 && nrow(data3)==0 && nrow(data4)==0) {data_plot=go_mf}

data_plot<-data_plot[order(data_plot$type,data_plot$`log10.qval.`,decreasing =F),]
data_plot$termName<-factor(data_plot$termName,levels = data_plot$termName)
tab_class<-as.data.frame(table(data_plot$type))
col_times<-tab_class$Freq
col_key<-pal_futurama()(nrow(tab_class))

ggplot(data_plot,aes(x=-log10.qval.,y=termName,-log10.qval.))+
  geom_col(aes(fill=type),alpha=0.7,width = 0.8)+
  scale_fill_manual(values = col_key)+
  geom_text(aes(x=-log10.qval.,y=termName,label=-log10.qval.))+
  theme_test()+
  theme(legend.position = 'right')+
  theme(axis.text.y = element_text(face = 'plain',size=14,colour =rep(col_key,times=col_times)))+
  guides(fill = guide_legend(reverse=TRUE))+
  labs(x='-log10qval',y=NULL,fill='Pathway')

#outFile=paste("snakemake/result/11_data_analysis/coexpression/",gene_id,"/",Out,sep="")
ggsave(Out,width = 15,height = 12)