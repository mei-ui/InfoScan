suppressMessages(library("optparse"))
suppressMessages(library("tidyverse"))
option_list <- list(
    make_option(c("-i","--cpatoutput_file"),type="character", default=NULL,help="input cpat result "),
    make_option(c("-o","--outfile_prefix"),type="character", default=NULL, help="output coding transcript"),
    make_option(c("-p","--coding_probability_cutoff"),type="character",default=NULL, help="input coding probability (CP) cutoff")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$cpatoutput_file)|| is.null(opt$outfile_prefix)|| is.null(opt$coding_probability_cutoff) ){
  print_help(opt_parser);
  stop("Please provide -i cpatoutput_file and -o outfile_prefix ,-p coding_probability_cutoff option", call.=FALSE);
}
cpatoutput=opt$cpatoutput_file
outfile=opt$outfile_prefix
CPcutoff=opt$coding_probability_cutoff

cpat<-read_tsv(cpatoutput)
cpat<-as.data.frame(cpat)
cpat$ID<-str_remove_all(cpat$ID,pattern='_[A-Z]*_[0-9]')
coding = cpat %>% group_by(ID) %>% top_n(n = 1, wt = Coding_prob)
coding<-as.data.frame(coding)
#print(CPcutoff)
#CPcutoff <- as.integer(CPcutoff)
#print(CPcutoff)
coding<-coding[coding$Coding_prob>CPcutoff,]
coding_1<-as.data.frame(coding$ID)
colnames(coding_1)<-'id'
write.table(coding_1,outfile,quote=F,col.names=F,row.names=F)
