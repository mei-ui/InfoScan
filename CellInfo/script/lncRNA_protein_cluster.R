suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))
option_list <- list(
                make_option(c("-a","--scRNA_rds"),type="character", default=NULL, help="input scRNA rds file"),
                make_option(c("-b","--lncRNA_id"),type="character", default=NULL, help="input lncRNA cluster id file"),
		            make_option(c("-c","--protein_id"),type="character", default=NULL, help="input protein cluster id file"),
                make_option(c("-e","--TSNE_lncRNA"),type="character", default=NULL, help="output TSNE_lncRNA.pdf file"),
                make_option(c("-f","--UMAP_lncRNA"),type="character", default=NULL, help="output UMAP_lncRNA.pdf file"),
		            make_option(c("-n","--TSNE_protein"),type="character", default=NULL, help="output TSNE_protein.pdf file"),
		            make_option(c("-o","--UMAP_protein"),type="character", default=NULL, help="output UMAP_protein.pdf file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$scRNA_rds)|| is.null(opt$lncRNA_id)|| is.null(opt$protein_id)|| is.null(opt$TSNE_lncRNA)|| is.null(opt$UMAP_lncRNA)|| 
	is.null(opt$TSNE_protein)|| is.null(opt$UMAP_protein)){
  print_help(opt_parser);
  stop("Please provide -a scRNA_rds, -b lncRNA_id ,-c protein_id,-e TSNE_lncRNA ,-f UMAP_lncRNA,-n TSNE_protein,
      -o UMAP_protein", call.=FALSE);
}

SC=opt$scRNA_rds
lncRNA=opt$lncRNA_id
protein=opt$protein_id
tsne_lncRNA=opt$TSNE_lncRNA
umap_lncRNA=opt$UMAP_lncRNA
tsne_protein=opt$TSNE_protein
umap_protein=opt$UMAP_protein

scRNA = readRDS(SC)
lncRNA_cluster =  read.csv(lncRNA)
colnames(lncRNA_cluster) = c("cell_ID","lncRNA")
rownames(lncRNA_cluster) = lncRNA_cluster$cell_ID
scRNA = AddMetaData(object = scRNA, metadata = lncRNA_cluster)

protein_cluster =  read.csv(protein)
colnames(protein_cluster) = c("cell_ID","protein")
rownames(protein_cluster) = protein_cluster$cell_ID
scRNA = AddMetaData(object = scRNA, metadata = protein_cluster)

plot1 = DimPlot(scRNA, reduction = "tsne",label=T, label.size=2,label.box = TRUE,group.by="lncRNA",pt.size = 1)
ggsave(tsne_lncRNA, plot = plot1, width = 8, height = 7)
plot1 = DimPlot(scRNA, reduction = "umap",label=T, label.size=2,label.box = TRUE,group.by="lncRNA",pt.size = 1)
ggsave(umap_lncRNA, plot = plot1, width = 8, height = 7)

plot1 = DimPlot(scRNA, reduction = "tsne",label=T, label.size=2,label.box = TRUE,group.by="lncRNA",pt.size = 1)
ggsave(tsne_protein, plot = plot1, width = 8, height = 7)
plot1 = DimPlot(scRNA, reduction = "umap",label=T, label.size=2,label.box = TRUE,group.by="lncRNA",pt.size = 1)
ggsave(umap_protein, plot = plot1, width = 8, height = 7)