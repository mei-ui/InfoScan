suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggrepel"))
suppressMessages(library("optparse"))
suppressMessages(library("ggsci"))
option_list <- list(
    make_option(c("-i","--id"),type="character", default=NULL, help="input gene set enrichment result"),
    make_option(c("-o","--output"),type="character", default=NULL, help="input singleR annotation file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$id)|| is.null(opt$output)){
  print_help(opt_parser);
  stop("Please provide -i gene id, -o output file", call.=FALSE);
}

gene_id=opt$id
Out=opt$output

data1<-read.table(gene_id,header=T)
if (nrow(data1)!=0){
  tmp<-arrange(data1,data1$`log10.qval.`)
  tmp$number <- factor(rev(1:nrow(tmp)))
  tmp$termName=tolower(tmp$termName)
  tmp$`log10.qval.`=round(tmp$`log10.qval.`,2)
  #tmp$termName=str_to_title(tmp$termName)
  #tmp<-arrange(tmp,desc(tmp$commonGeneNum))
  tmp$type='pathway'
  if(nrow(tmp)<=20){tmp=tmp}else{tmp<-tmp[1:20,]}
}

data_plot=tmp

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
ggsave(Out,width = 15,height = 8)