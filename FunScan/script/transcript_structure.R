suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("patchwork"))
suppressMessages(library("Seurat"))
suppressMessages(library("dplyr"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("ggsci"))
suppressMessages(library("RColorBrewer"))
source("snakemake/script/ggtranscript-master/R/geom_range.R")
source("snakemake/script/ggtranscript-master/R/add_exon_number.R")
source("snakemake/script/ggtranscript-master/R/add_utr.R")
source("snakemake/script/ggtranscript-master/R/data.R")
source("snakemake/script/ggtranscript-master/R/geom_half_range.R")
source("snakemake/script/ggtranscript-master/R/geom_intron.R")
source("snakemake/script/ggtranscript-master/R/geom_junction.R")
source("snakemake/script/ggtranscript-master/R/geom_junction_label_repel.R")
source("snakemake/script/ggtranscript-master/R/ggtranscript-package.R")
source("snakemake/script/ggtranscript-master/R/globals.R")
source("snakemake/script/ggtranscript-master/R/shorten_gaps.R")
source("snakemake/script/ggtranscript-master/R/to_diff.R")
source("snakemake/script/ggtranscript-master/R/to_intron.R")
source("snakemake/script/ggtranscript-master/R/utils.R")
#setwd("/public/home/meisq/02data/10Mus_brain/data_2/3_9_M/1020output")
option_list <- list(
    	make_option(c("-a","--lncRNA_class_unannotated"),type="character", default=NULL, help="input lncRNA_class_unannotated file"),
		make_option(c("-b","--mRNA_class_unannotated"),type="character", default=NULL, help="input mRNA_class_unannotated file"),
    	make_option(c("-c","--lncRNA_gtf_unannotated"),type="character", default=NULL, help="input lncRNA_gtf_unannotated file"),	
		make_option(c("-d","--mRNA_gtf_unannotated"),type="character", default=NULL, help="input mRNA_gtf_unannotated file"),
		make_option(c("-e","--All_gtf"),type="character", default=NULL, help="input All_gtf file"),
		make_option(c("-f","--transcript_id"),type="character", default=NULL, help="input transcript_id file"),
        make_option(c("-g","--scRNA"),type="character", default=NULL, help="input scRNA RDS file"),
        make_option(c("-o","--output"),type="character", default=NULL, help="output")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$lncRNA_class_unannotated)|| is.null(opt$mRNA_class_unannotated)||is.null(opt$lncRNA_gtf_unannotated)||is.null(opt$mRNA_gtf_unannotated)||
    is.null(opt$All_gtf)||is.null(opt$transcript_id)||is.null(opt$scRNA)||is.null(opt$output)){
  print_help(opt_parser);
  stop("Please provide -a lncRNA_class_unannotated ,-b mRNA_class_unannotated and -c lncRNA_gtf_unannotated,-d mRNA_gtf_unannotated,
        -e All_gtf,-f transcript_id,-g scRNA,-o output", call.=FALSE);
}
lncRNA_class=opt$lncRNA_class_unannotated
mRNA_class=opt$mRNA_class_unannotated
lncRNA_gtf=opt$lncRNA_gtf_unannotated
mRNA_gtf=opt$mRNA_gtf_unannotated
all_gtf=opt$All_gtf
SCRNA=opt$scRNA
Transcript_Id=opt$transcript_id
Output=opt$output

unannotated_lncRNA_class = read.table(lncRNA_class,header=T)
unannotated_lncRNA_class$type2="unannotated_lncRNA"
unannotated_mRNA_class = read.table(mRNA_class,header=T)
unannotated_mRNA_class$type2="unannotated_coding_transcripts"
class_code=rbind(unannotated_lncRNA_class,unannotated_mRNA_class)

unannotated_lncRNA_gtf= rtracklayer::import(lncRNA_gtf)
unannotated_lncRNA_gtf<-as.data.frame(unannotated_lncRNA_gtf)
unannotated_mRNA_gtf= rtracklayer::import(mRNA_gtf)
unannotated_mRNA_gtf<-as.data.frame(unannotated_mRNA_gtf)
unannotated_gtf = rbind(unannotated_lncRNA_gtf,unannotated_mRNA_gtf)
unannotated_gtf_1 = inner_join(unannotated_gtf,class_code,by=c("transcript_id"="cuff_gene_id"))
#head(unannotated_gtf_1)
gtf_all = rtracklayer::import(all_gtf)
gtf_all =as.data.frame(gtf_all)

if("transcript_name" %in% colnames(gtf_all)){
    if("gene_name" %in% colnames(gtf_all)){
      sln <- c('seqnames','start','end','width','strand','type','gene_id','gene_name','transcript_id','transcript_name')
    }else{
      sln <- c('seqnames','start','end','width','strand','type','gene_id','transcript_id','transcript_name')
    }
  }else{
    if("gene_name" %in% colnames(gtf_all)){
      sln <- c('seqnames','start','end','width','strand','type','gene_id','gene_name','transcript_id')
    }else{
      sln <- c('seqnames','start','end','width','strand','type','gene_id','transcript_id')
    }
  }
# select data
unannotated_gtf_1$transcript_id<-stri_replace_all_regex(unannotated_gtf_1$transcript_id,'_','-')
Transcript_Id<-stri_replace_all_regex(Transcript_Id,'_','-')
unannotated_transcript_gene = filter(unannotated_gtf_1,transcript_id==Transcript_Id)
gene=as.character(unique(unannotated_transcript_gene$ref_gene_id))
# if(gene=="-"){gene=NULL}
if (!is.null(gene)) {
  if (gene == "-") {
    gene <- NULL
  }
}

head(gene)
# use gene_id stands for gene_name
if(is.null(gene)==FALSE){   
    if("gene_name" %in% colnames(gtf_all)){
      # filter gene by gene name
      myGene <- gtf_all%>%
        dplyr::filter(gene_name %in% .env$gene) %>%
        dplyr::filter(type != 'gene') %>%
        dplyr::select(sln)
    }else{
      # filter gene by gene name
      myGene <- gtf_all%>%
        dplyr::select(sln) %>%
        dplyr::mutate(gene_name = gene_id) %>%
        dplyr::filter(gene_name %in% .env$gene) %>%
        dplyr::filter(type != 'gene')
    }
    unannotated_transcript_gene=filter(unannotated_gtf_1,ref_gene_id==gene)
    max = c()
    for(i in 1:nrow(unannotated_transcript_gene)){
    	if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="j"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="novel isoform transcript"
    		max = rbind(data,max)}
    	if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="i"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="intron transcript"
    		max = rbind(data,max)}
        if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="u"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="intergenic transcript"
    		max = rbind(data,max)}
        if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="o"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="novel isoform transcript"
    		max = rbind(data,max)}
        if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="x"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="novel isoform transcript"
    		max = rbind(data,max)}
    }
    unannotated_transcript_gene=max
    unannotated_transcript_gene$gene_name=gene
    unannotated_transcript_gene$gene_id=unique(myGene$gene_id)
    unannotated_transcript_gene$transcript_name=unannotated_transcript_gene$transcript_id

    unannotated_transcript_gene = unannotated_transcript_gene[,c('seqnames','start','end','width','strand','type','gene_id',
                                                                'gene_name','transcript_id','transcript_name','transcript_type')]
    unannotated_transcript_gene$strand=unique(myGene$strand)
    myGene$transcript_type="annotated isoform"
    myGene = rbind(myGene,unannotated_transcript_gene)
    myData <- myGene
    # 筛选外显子
    exons <- myData %>% dplyr::filter(type == "exon")
    exons_rescaled <- shorten_gaps(
      exons = exons,
      introns = to_intron(exons, "transcript_name"),
      group_var = "transcript_name"
    )
    # 绘图
    # let's split these for plotting
    exons_rescaled_exons <- exons_rescaled %>% dplyr::filter(type == "exon")
    exons_rescaled_introns <- exons_rescaled %>% dplyr::filter(type == "intron")

    p=ggplot(exons_rescaled_exons,aes(xstart = start,xend = end,y = transcript_name)) +
        geom_range(aes(fill = transcript_type)) +
        geom_intron(data = exons_rescaled_introns,aes(strand = strand),arrow.min.intron.length = 300)+
        scale_fill_manual(values = brewer.pal(length(unique(exons_rescaled_exons$transcript_type)), "Set3"))

    #transcript及其宿主基因的表达情况
    scRNA=readRDS(SCRNA)
    lnc = unique(unannotated_transcript_gene$transcript_id)
    lnc = c(gene,lnc)
    if(length(lnc)==3){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=3)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=3)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,1.5))
    ggsave(Output,p3,width = 12,height = 12)}
    if(length(lnc)==2){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=2)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=2)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,2))
    ggsave(Output,p3,width = 9,height = 9)}
    if(length(lnc)==4){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=2)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=2)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1.5,2))
    ggsave(Output,p3,width = 9,height = 15)}
    if(length(lnc)>4&length(lnc)<=6){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=3)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=3)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,2))
    ggsave(Output,p3,height = 20,width = 12)}
    if(length(lnc)>6&length(lnc)<=10){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=5)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=5)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,2))
    ggsave(Output,p3,height = 24,width = 18)}
    if(length(lnc)>10&length(lnc)<=15){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=5)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=5)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,2))
    ggsave(Output,p3,height = 28,width = 18)}
    if(length(lnc)>15&length(lnc)<=20){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=5)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=5)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,2))
    ggsave(Output,p3,height = 36,width = 18)}
    if(length(lnc)>20&length(lnc)<=25){
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype',ncol=5)
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T,ncol=5)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,2))
    ggsave(Output,p3,height = 44,width = 18)}
}
if(is.null(gene)){
    max = c()
    for(i in 1:nrow(unannotated_transcript_gene)){
    	if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="j"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="novel isoform transcript"
    		max = rbind(data,max)}
    	if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="i"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="intron transcript"
    		max = rbind(data,max)}
        if(as.character(unannotated_transcript_gene[i,][,"class_code"])=="u"){
    		data = unannotated_transcript_gene[i,]
    		data$transcript_type="intergenic transcript"
    		max = rbind(data,max)}
    }
    unannotated_transcript_gene=max
    unannotated_transcript_gene$gene_name=unannotated_transcript_gene$transcript_id
    unannotated_transcript_gene$gene_id=unannotated_transcript_gene$transcript_id
    unannotated_transcript_gene$transcript_name=unannotated_transcript_gene$transcript_id

    unannotated_transcript_gene = unannotated_transcript_gene[,c('seqnames','start','end','width','strand','type','gene_id',
                                                            'gene_name','transcript_id','transcript_name','transcript_type')]
    myGene = unannotated_transcript_gene
    myData <- myGene
    # 筛选外显子
    exons <- myData %>% dplyr::filter(type == "exon")
    exons_rescaled <- shorten_gaps(
      exons = exons,
      introns = to_intron(exons, "transcript_name"),
      group_var = "transcript_name"
    )
    # 绘图
    # let's split these for plotting
    exons_rescaled_exons <- exons_rescaled %>% dplyr::filter(type == "exon")
    exons_rescaled_introns <- exons_rescaled %>% dplyr::filter(type == "intron")

    p=ggplot(exons_rescaled_exons,aes(xstart = start,xend = end,y = transcript_name)) +
        geom_range(aes(fill = transcript_type)) +
        geom_intron(data = exons_rescaled_introns,aes(strand = strand),arrow.min.intron.length = 300)+
        scale_fill_manual(values = brewer.pal(length(unique(exons_rescaled_exons$transcript_type)), "Set3"))

    #transcript及其宿主基因的表达情况
    scRNA=readRDS(SCRNA)
    lnc = unique(unannotated_transcript_gene$transcript_id)
    p1=VlnPlot(object = scRNA, features =lnc,log =T ,pt.size = 0.1,group.by = 'celltype')
    p2=FeaturePlot(scRNA,features = lnc,reduction = "tsne",label = T)
    p3=p+p1+p2+plot_layout(ncol = 1, heights = c(2,1,2))
    ggsave(Output,p3,width = 9,height = 12)
}