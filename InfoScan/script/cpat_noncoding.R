suppressMessages(library("optparse"))
suppressMessages(library("tidyverse"))

option_list <- list(
        make_option(c("-i","--cpatoutput_file"),type="character", default=NULL,help="input cpat result"),
		make_option(c("-I","--cpatoutput_nonORF_file"),type="character", default=NULL,help="input cpat result"),
        make_option(c("-o","--outfile_prefix"),type="character", default=NULL, help="output coding transcript"),
        make_option(c("-p","--coding_probability_cutoff"),type="character",default=NULL, help="input coding probability (CP) cutoff")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$cpatoutput_file)|| is.null(opt$outfile_prefix) || is.null(opt$cpatoutput_nonORF_file)|| is.null(opt$coding_probability_cutoff)){
  print_help(opt_parser);
  stop("Please provide -i cpat_file,-I nonORF_file and -o outfile_prefix,-p coding_probability_cutoff option", call.=FALSE);
}
cpatoutput=opt$cpatoutput_file
nonORF=opt$cpatoutput_nonORF_file
outfile=opt$outfile_prefix
CPcutoff=opt$coding_probability_cutoff

cpat<-read_tsv(cpatoutput)
cpat<-as.data.frame(cpat)
cpat$ID<-str_remove_all(cpat$ID,pattern='_[A-Z]*_[0-9]')
noncoding = cpat %>% group_by(ID) %>% top_n(n = 1, wt = Coding_prob)
noncoding<-as.data.frame(noncoding)
noncoding<-noncoding[noncoding$Coding_prob<CPcutoff,]
noncoding_1<-as.data.frame(noncoding$ID)
colnames(noncoding_1)<-'id'
no_ORF<-read.table(nonORF)
colnames(no_ORF)<-'id'
cpat<-rbind(noncoding_1,no_ORF)
write.table(cpat,outfile,quote=F,col.names=F,row.names=F)
