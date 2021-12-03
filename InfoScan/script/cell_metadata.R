suppressMessages(library("openxlsx"))
suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))
option_list <- list(
                make_option(c("-f","--pca"),type="character", default=NULL, help="output pca.pdf file"),
                make_option(c("-l","--tsne1"),type="character", default=NULL, help="output tSNE1.pdf file"),
                make_option(c("-i","--umap1"),type="character", default=NULL, help="output UMAP1.pdf file"),
                make_option(c("-j","--group_chose"),type="character", default=NULL, help="input chose group"),
                make_option(c("-m","--scRNA_rds"),type="character", default=NULL, help="input seurat file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$pca)|| is.null(opt$tsne1)|| is.null(opt$umap1)|| is.null(opt$scRNA_rds)|| is.null(opt$group_chose)){
  print_help(opt_parser);
  stop("Please provide -f pca and -l tSNE1,-i UMAP1,-m scRNA_rds,-j group_chose", call.=FALSE);
}
PCA=opt$pca
TSNE1=opt$tsne1
UMAP1=opt$umap1
SCrds=opt$scRNA_rds
group=opt$group_chose
scRNA<-readRDS(SCrds)
##PCA
plot1 <- DimPlot(scRNA, reduction = "pca",group.by=group)
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca")
plotc <- plot1+plot2
ggsave(PCA, plot = plotc, width = 8, height = 4)
#TSNE
plot1 = DimPlot(scRNA, reduction = "tsne",label = F,group.by=group)
ggsave(TSNE1, plot = plot1, width = 10, height = 7)
#UMAP
plot2 = DimPlot(scRNA, reduction = "umap",label =F,group.by=group)
ggsave(UMAP1, plot = plot2, width = 10, height = 7)

