suppressMessages(library("stringr"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggsci"))
suppressMessages(library('Seurat'))
suppressMessages(library('pheatmap'))
suppressMessages(library('ComplexHeatmap'))

option_list <- list(
        make_option(c("-a","--finall_fpkm"),type="character", default=NULL, help="input finally_fpkm file"),
        make_option(c("-b","--group_chose"),type="character", default=NULL, help="input id file"),
        make_option(c("-c","--celltype_rds"),type="character", default=NULL, help="input scRNA celltype rds file"),
	    make_option(c("-o","--unannotated_lncRNA_heatmap"),type="character", default=NULL, help="unannotated_lncRNA_heatmap.pdf file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$finall_fpkm)|| is.null(opt$group_chose) ||is.null(opt$celltype_rds)||
    is.null(opt$unannotated_lncRNA_heatmap)){
  print_help(opt_parser);
  stop("Please provide -a finall_fpkm, -b group_chose ,-i protein_coding_heatmap and -I annotated_lncRNA_heatmap,-c celltype_rds,-o unannotated_lncRNA_heatmap ,
        -O unannotated_protein_heatmap", call.=FALSE);
}

finallfpkm=opt$finall_fpkm
CR=opt$celltype_rds
ULH=opt$unannotated_lncRNA_heatmap
groupChose=opt$group_chose
##加载数据

##load data
finally_fpkm <- read.table(finallfpkm,header=T)
name<-as.data.frame(table(finally_fpkm$Geneid))
name<-filter(name,Freq>1)
if(nrow(name)!=0){finally_fpkm_1<-finally_fpkm[-which(finally_fpkm$Geneid%in%c(as.character(name$Var1))),]}else {finally_fpkm_1<-finally_fpkm}
scRNA<-readRDS(CR)
#数据预处理
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
n=length(x)
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
finally_fpkm_3=finally_fpkm_3[which(finally_fpkm_3$type%in%c('protein_coding','annotated_lncRNA','unannotated_lncRNA','unannotated_coding_transcript')),]
finally_fpkm_eset = finally_fpkm_3[,-which(colnames(finally_fpkm_3)%in%c("Geneid","type"))]
group=finally_fpkm_3[,which(colnames(finally_fpkm_3)%in%c("Geneid","type"))]
rownames(finally_fpkm_eset)<-finally_fpkm_3$Geneid
finally_fpkm_eset<- finally_fpkm_eset[apply(finally_fpkm_eset, 1, function(x) sd(x)!=0),]
Group=as.data.frame(rownames(finally_fpkm_eset))
colnames(Group)='id'
Group=inner_join(Group,group,by=c("id"="Geneid"))
Group1=Group[,"type",drop=F]
rownames(Group1)=Group$id
finally_fpkm_eset=log2(finally_fpkm_eset+1)#标准化处理
Group1$type = factor(Group1$type,levels = c('protein_coding','annotated_lncRNA','unannotated_lncRNA','unannotated_coding_transcript'))

col_key1<-pal_futurama()(ncol(finally_fpkm_eset))
col_anno=as.data.frame(colnames(finally_fpkm_eset))
colnames(col_anno)=c("tissue")
top_annotation = HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col_key1),
                   labels = colnames(finally_fpkm_eset),
                   labels_gp = gpar(col = "black", fontsize = 8)),annotation_height=unit(0.4, "cm"))
col_key2=pal_futurama()(length(table(Group1$type)))
row_annotation = rowAnnotation(cluster = anno_block(gp = gpar(fill = col_key2)),width=unit(0.4, "cm"))
#pdf("snakemake/result/test.pdf")
pdf(ULH)
p1=Heatmap(finally_fpkm_eset,#表达矩阵
#     col = colorRampPalette(c("navy","white","firebrick3"))(100),#颜色定义
    show_row_names = F,#不展示行名
    row_split = Group1,#用group信息将热土分开，以group聚类
    column_title = "Tissue Heatmap",#不显示列标题
    show_column_names = F,#不显示列名
    row_title_rot = 0,
    column_split=col_anno,
    show_parent_dend_line = FALSE,#隐藏虚线
    top_annotation = top_annotation,#顶部分组信息
    left_annotation = row_annotation,
    gap = unit(2, "mm"),
    show_column_dend = FALSE,#不显示列的树状图
    row_title_gp = gpar(col = col_key2),
    border = TRUE,#添加边界
    heatmap_legend_param = list(title = "Expression"),
    cluster_rows = TRUE)
draw(p1)
dev.off()