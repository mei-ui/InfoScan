#/public/home/meisq/miniconda3/pkgs/r-base-3.6.3-hd23ff56_6/lib/R/bin/Rscript
suppressMessages(library("rtracklayer"))
suppressMessages(library("SummarizedExperiment"))
suppressMessages(library("tidyverse"))
suppressMessages(library("optparse"))
option_list <- list(
    make_option(c("-i","--gtf_file"),type="character", default=NULL,help="输入gtf文件"),
    make_option(c("-I","--txt_file"),type="character", default=NULL, help="输入txt文件"),
    make_option(c("-o","--outfile_prefix"),type="character", default=NULL, help="输出含有txt文件的gtf文件")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gtf_file)|| is.null(opt$txt_file) || is.null(opt$outfile_prefix) ){
  print_help(opt_parser);
  stop("Please provide -i gtf_file, -I txt_file and -o outfile_prefix option", call.=FALSE);
}

gtffile=opt$gtf_file
txtfile=opt$txt_file
outfile=opt$outfile_prefix

gtf1 <- rtracklayer::import(gtffile)
gtf_df <- as.data.frame(gtf1)
test<-read.table(txtfile)
lncRNA_1<-inner_join(gtf_df,test,by=c("transcript_id"='V1'))
export(lncRNA_1,outfile,'gtf')
