suppressMessages(library("openxlsx"))
suppressMessages(library("tidyverse"))
suppressMessages(library("optparse"))
option_list <- list(
                make_option(c("-a","--coexpression_path"),type="character", default=NULL, help="输入coexpression gene path")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$coexpression_path)) {
  print_help(opt_parser);
  stop("Please provide -a coexxpression_path", call.=FALSE);
}
copath=opt$coexpression_path



#setwd('/public/home/meisq/02data/10Mus_brain/data_2/3_11_M/snakemake/11_data_analysis/coexpression')
#path<- c("/public/home/meisq/02data/10Mus_brain/data_2/3_11_M/snakemake/11_data_analysis/coexpression/")
path<-c(copath)
file_names<- list.files(path)

for (i in 1:length(file_names)) {
  
  name<-gsub(".txt","",file_names[i])
  data<-read.table(paste0(path,file_names[i]),header = T, stringsAsFactors = FALSE)
	querygene<-unique(data$queryGene)
	coexpsymbol<-data$coexpSymbol
	coexgene<-c(querygene,coexpsymbol)
	coexgene<-t(coexgene)
	data.out<-paste(path,querygene,'.txt',sep="")
	write.table(coexgene,data.out,quote=F,col.names=F,row.names=F,sep='\t')
}
	
