suppressMessages(library("stringr"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggsci"))
suppressMessages(library('cummeRbund'))

option_list <- list(
    	make_option(c("-a","--finally_fpkm"),type="character", default=NULL, help="input finally_fpkm file"),
			make_option(c("-i","--cell_specificity"),type="character", default=NULL, help="output cell_specificity.gtf file"),
    	make_option(c("-I","--tissue_specificity"),type="character", default=NULL, help="tissue_specificity.pdf file"),	
			make_option(c("-c","--celltype_rds"),type="character", default=NULL, help="input scRNA celltype rds file"),
			make_option(c("-o","--cell_ecdf"),type="character", default=NULL, help="output cell_specificity_ecdf.gtf file"),
			make_option(c("-O","--tissue_ecdf"),type="character", default=NULL, help="output tissue_specificity_ecdf.gtf file"),
      make_option(c("-d","--group_chose"),type="character", default=NULL, help="input chose group")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$finally_fpkm)|| is.null(opt$cell_specificity)||is.null(opt$tissue_specificity)||is.null(opt$celltype_rds)||
    is.null(opt$cell_ecdf)||is.null(opt$tissue_ecdf)||is.null(opt$group_chose)){
  print_help(opt_parser);
  stop("Please provide -a finally_fpkm ,-i cell_specificity and -I tissue_specificity,-c celltype_rds,-o cell_specificity_ecdf,
        -O tissue_specificity_ecdf,-d group_chsoe", call.=FALSE);
}
finallfpkm=opt$finally_fpkm
CS=opt$cell_specificity
TS=opt$tissue_specificity
CR=opt$celltype_rds
CSE=opt$cell_ecdf
TSE=opt$tissue_ecdf
groupChose=opt$group_chose
#计算jsd函数
tissuespecificity<-function(fpkms,logMode=T,pseudocount=1,relative=FALSE,...){
  if(logMode){
    fpkms<-log10(fpkms+pseudocount)
  }
  fpkms<-t(makeprobs(t(fpkms)))
  d<-diag(ncol(fpkms))
  res<-apply(d,MARGIN=1,function(q){
    JSdistFromP(fpkms,q)
  })
  colnames(res)<-paste(colnames(fpkms),"_spec",sep="")
  
  if(relative){
    res<-res/max(res)
  }
  1-res
}
##load data
finally_fpkm <- read.table(finallfpkm,header=T)
name<-as.data.frame(table(finally_fpkm$Geneid))
name<-filter(name,Freq>1)
if(nrow(name)!=0){finally_fpkm_1<-finally_fpkm[-which(finally_fpkm$Geneid%in%c(as.character(name$Var1))),]}else {finally_fpkm_1<-finally_fpkm}

############celltype specificity###########
#load data
scRNA<-readRDS(CR)
if(length(unique(scRNA@meta.data$celltype))>1){
  celltype<-as.data.frame(scRNA$celltype)
  colnames(celltype)<-'celltype'
  celltype$sample<-rownames(celltype)
  fpkm<-finally_fpkm_1[,7:(ncol(finally_fpkm_1)-1)]
  rownames(fpkm)<-finally_fpkm_1$Geneid
  fpkm_1<-as.data.frame(t(fpkm))
  fpkm_1$sample<-rownames(fpkm_1)
  fpkm_1<-inner_join(fpkm_1,celltype,by=c('sample'='sample'))
  fpkm_group = fpkm_1 %>% group_by(celltype) 
  x<-group_split(fpkm_group)
  n=length(x)
  data<-x[[1]]
  data<-as.data.frame(data)								
  data_eset<-data[,1:(ncol(data)-2)]
  data_1<-apply(data_eset,2,mean)
  name<-as.character(unique(data$celltype))
  data_1<-as.data.frame(data_1)
  colnames(data_1)<-name
  for (i in 2:as.numeric(n)) {
  			data_2<-x[[i]]
  			data_2<-as.data.frame(data_2)
  			data_eset_2<-data_2[,1:(ncol(data_2)-2)]
  			data_3<-apply(data_eset_2,2,mean)
  			name_2<-as.character(unique(data_2$celltype))
  			data_3<-as.data.frame(data_3)
  			colnames(data_3)<-name_2
  			data_1<-cbind(data_1,data_3)
  }
  fpkm_2<-data_1
  fpkm_2$gene_id<-row.names(fpkm_2)
  id<-finally_fpkm_1[,which(colnames(finally_fpkm_1)%in%c("Geneid","type"))]
  finally_fpkm_3<-inner_join(id,fpkm_2,by=c('Geneid'='gene_id'))
  #protein_coding
  protein_coding<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('protein_coding')),]
  protein_coding<-as.data.frame(protein_coding)
  protein_coding_eset<-protein_coding[,-which(colnames(protein_coding)%in%c("Geneid","type"))]
  rownames(protein_coding_eset)<-protein_coding$Geneid
  protein_coding_js.distance<-tissuespecificity(protein_coding_eset)
  protein_coding_js.distance<-as.data.frame(apply(protein_coding_js.distance,1,max))
  colnames(protein_coding_js.distance)<-'value'
  protein_coding_js.distance$group<-'protein coding'
  #annotated lncRNA
  annotated_lncRNA<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('annotated_lncRNA')),]
  annotated_lncRNA<-as.data.frame(annotated_lncRNA)
  annotated_lncRNA_eset<-annotated_lncRNA[,-which(colnames(annotated_lncRNA)%in%c("Geneid","type"))]
  rownames(annotated_lncRNA_eset)<-annotated_lncRNA$Geneid
  annotated_lncRNA_js.distance<-tissuespecificity(annotated_lncRNA_eset)
  annotated_lncRNA_js.distance<-as.data.frame(apply(annotated_lncRNA_js.distance,1,max))
  colnames(annotated_lncRNA_js.distance)<-'value'
  annotated_lncRNA_js.distance$group<-'annotated lncRNA'
  #unannotated lncRNA
  unannotated_lncRNA<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('unannotated_lncRNA')),]
  unannotated_lncRNA<-as.data.frame(unannotated_lncRNA)
  unannotated_lncRNA_eset<-unannotated_lncRNA[,-which(colnames(unannotated_lncRNA)%in%c("Geneid","type"))]
  rownames(unannotated_lncRNA_eset)<-unannotated_lncRNA$Geneid
  unannotated_lncRNA_js.distance<-tissuespecificity(unannotated_lncRNA_eset)
  unannotated_lncRNA_js.distance<-as.data.frame(apply(unannotated_lncRNA_js.distance,1,max))
  colnames(unannotated_lncRNA_js.distance)<-'value'
  unannotated_lncRNA_js.distance$group<-'unannotated lncRNA'
  #unannotated coding transcript
  unannotated_protein<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('unannotated_coding_transcript')),]
  unannotated_protein<-as.data.frame(unannotated_protein)
  unannotated_protein_eset<-unannotated_protein[,-which(colnames(unannotated_protein)%in%c("Geneid","type"))]
  rownames(unannotated_protein_eset)<-unannotated_protein$Geneid
  unannotated_protein_js.distance<-tissuespecificity(unannotated_protein_eset)
  unannotated_protein_js.distance<-as.data.frame(apply(unannotated_protein_js.distance,1,max))
  colnames(unannotated_protein_js.distance)<-'value'
  unannotated_protein_js.distance$group<-'unannotated coding transcript'
  df<-rbind(protein_coding_js.distance,annotated_lncRNA_js.distance,unannotated_lncRNA_js.distance,unannotated_protein_js.distance)
  
  ggplot(df,aes(value,color=group)) + 
    xlab("JS Specificity Score") + 
    geom_density(alpha=1) + 
  	labs(title = "celltype specificity")+
    theme_bw()+scale_color_aaas()
  ggsave(CS,width = 6, height = 5)
  
  ggplot(df,aes(value,color=group)) +
    xlab("JS Specificity Score") +
    stat_ecdf(alpha=1) +
    labs(title = "celltype specificity")+  theme_bw()+scale_color_aaas()
  ggsave(CSE,width = 6, height = 5)
  ggsave("snakemake/result/celltype_specificity_ecdf.png",width = 6, height = 4)
}
##################### tissue specificity ######################
#load data
tissue = eval(parse(text = paste0("scRNA","$",groupChose)))
tissue<-as.data.frame(tissue)
colnames(tissue)<-"tissue"
tissue$sample<-rownames(tissue)
fpkm<-finally_fpkm_1[,7:(ncol(finally_fpkm_1)-1)]
rownames(fpkm)<-finally_fpkm_1$Geneid
fpkm_1<-as.data.frame(t(fpkm))
fpkm_1$sample<-rownames(fpkm_1)
fpkm_1<-inner_join(fpkm_1,tissue,by=c('sample'='sample'))
fpkm_group = fpkm_1 %>% group_by(tissue) 
#n<-nrow(table(fpkm_1$tissue))
x<-group_split(fpkm_group)
n = length(x)
data<-x[[1]]
data<-as.data.frame(data)
data_eset<-data[,1:(ncol(data)-2)]
data_1<-apply(data_eset,2,mean)
name<-as.character(unique(data$tissue))
data_1<-as.data.frame(data_1)
colnames(data_1)<-name
for (i in 2:as.numeric(n)) {
			data_2<-x[[i]]
			data_2<-as.data.frame(data_2)
			data_eset_2<-data_2[,1:(ncol(data_2)-2)]
			data_3<-apply(data_eset_2,2,mean)
			name_2<-as.character(unique(data_2$tissue))
			data_3<-as.data.frame(data_3)
			colnames(data_3)<-name_2
			data_1<-cbind(data_1,data_3)
}
fpkm_2<-data_1
fpkm_2$gene_id<-row.names(fpkm_2)
id<-finally_fpkm_1[,which(colnames(finally_fpkm_1)%in%c("Geneid","type"))]
finally_fpkm_3<-inner_join(id,fpkm_2,by=c('Geneid'='gene_id'))
#protein_coding
protein_coding<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('protein_coding')),]
protein_coding<-as.data.frame(protein_coding)
protein_coding_eset<-protein_coding[,-which(colnames(protein_coding)%in%c("Geneid","type"))]
rownames(protein_coding_eset)<-protein_coding$Geneid
protein_coding_js.distance<-tissuespecificity(protein_coding_eset)
protein_coding_js.distance<-as.data.frame(apply(protein_coding_js.distance,1,max))
colnames(protein_coding_js.distance)<-'value'
protein_coding_js.distance$group<-'protein coding'
#annotated lncRNA
annotated_lncRNA<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('annotated_lncRNA')),]
annotated_lncRNA<-as.data.frame(annotated_lncRNA)
annotated_lncRNA_eset<-annotated_lncRNA[,-which(colnames(annotated_lncRNA)%in%c("Geneid","type"))]
rownames(annotated_lncRNA_eset)<-annotated_lncRNA$Geneid
annotated_lncRNA_js.distance<-tissuespecificity(annotated_lncRNA_eset)
annotated_lncRNA_js.distance<-as.data.frame(apply(annotated_lncRNA_js.distance,1,max))
colnames(annotated_lncRNA_js.distance)<-'value'
annotated_lncRNA_js.distance$group<-'annotated lncRNA'
#unannotated lncRNA
unannotated_lncRNA<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('unannotated_lncRNA')),]
unannotated_lncRNA<-as.data.frame(unannotated_lncRNA)
unannotated_lncRNA_eset<-unannotated_lncRNA[,-which(colnames(unannotated_lncRNA)%in%c("Geneid","type"))]
rownames(unannotated_lncRNA_eset)<-unannotated_lncRNA$Geneid
unannotated_lncRNA_js.distance<-tissuespecificity(unannotated_lncRNA_eset)
unannotated_lncRNA_js.distance<-as.data.frame(apply(unannotated_lncRNA_js.distance,1,max))
colnames(unannotated_lncRNA_js.distance)<-'value'
unannotated_lncRNA_js.distance$group<-'unannotated lncRNA'
#unannotated protein coding 
unannotated_protein<-finally_fpkm_3[which(finally_fpkm_3$type%in%c('unannotated_coding_transcript')),]
unannotated_protein<-as.data.frame(unannotated_protein)
unannotated_protein_eset<-unannotated_protein[,-which(colnames(unannotated_protein)%in%c("Geneid","type"))]
rownames(unannotated_protein_eset)<-unannotated_protein$Geneid
unannotated_protein_js.distance<-tissuespecificity(unannotated_protein_eset)
unannotated_protein_js.distance<-as.data.frame(apply(unannotated_protein_js.distance,1,max))
colnames(unannotated_protein_js.distance)<-'value'
unannotated_protein_js.distance$group<-'unannotated_coding_transcript'
df<-rbind(protein_coding_js.distance,annotated_lncRNA_js.distance,unannotated_lncRNA_js.distance,unannotated_protein_js.distance)

p1=ggplot(df,aes(value,color=group)) + 
  xlab("JS Specificity Score") + 
  geom_density(alpha=1) + 
	labs(title = "tissue specificity")+
  theme_bw()+scale_color_aaas()
ggsave(TS,p1,width = 6, height = 5)

p2=ggplot(df,aes(value,color=group)) +
  xlab("JS Specificity Score") +
  stat_ecdf(alpha=1) +
  labs(title = "tissue specificity")+  theme_bw()+scale_color_aaas()
ggsave(TSE,p2,width = 6, height = 5)
ggsave("snakemake/result/tissue_specificity_ecdf.png",p2,width = 6, height = 4)
if(length(unique(scRNA@meta.data$celltype))==1){
  ggsave(CS,p1,width = 6, height = 5)
  ggsave(CSE,p2,width = 6, height = 5)
}