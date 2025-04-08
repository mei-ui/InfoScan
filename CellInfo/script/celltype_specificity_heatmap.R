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
        make_option(c("-o","--heatmap_plot"),type="character", default=NULL, help="output heatmap.gtf file"),
        make_option(c("-c","--celltype_rds"),type="character", default=NULL, help="input scRNA celltype rds file")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$finall_fpkm)|| is.null(opt$heatmap_plot)||is.null(opt$celltype_rds)){
  print_help(opt_parser);
  stop("Please provide -a finall_fpkm ,-i protein_coding_heatmap and,-c celltype_rds,-o unannotated_lncRNA_heatmap ,
        -O unannotated_protein_heatmap", call.=FALSE);
}

finallfpkm=opt$finall_fpkm
HP=opt$heatmap_plot
CR=opt$celltype_rds

##加载数据

##load data
finally_fpkm <- read.table(finallfpkm,header=T)
name<-as.data.frame(table(finally_fpkm$Geneid))
name<-dplyr::filter(name,Freq>1)
if(nrow(name)!=0){finally_fpkm_1<-finally_fpkm[-which(finally_fpkm$Geneid%in%c(as.character(name$Var1))),]}else {finally_fpkm_1<-finally_fpkm}
if(nrow(finally_fpkm_1)>60000){finally_fpkm_1 <- finally_fpkm_1[sample(nrow(finally_fpkm_1),60000),]}
scRNA<-readRDS(CR)
if(length(unique(scRNA@meta.data$celltype))>1){
    celltype<-as.data.frame(scRNA$celltype)
    colnames(celltype)<-'celltype'
    celltype$sample<-rownames(celltype)
    fpkm<-finally_fpkm_1[,7:(ncol(finally_fpkm_1)-1)]
    id = finally_fpkm_1[,"Geneid",drop=F]
    rownames(fpkm)<-finally_fpkm_1$Geneid
    fpkm_1<-as.data.frame(t(fpkm))
    fpkm_1$sample<-rownames(fpkm_1)
    fpkm_1<-inner_join(fpkm_1,celltype,by=c('sample'='sample'))
    fpkm_group = fpkm_1 %>% group_by(celltype) 
    #n<-nrow(table(fpkm_1$celltype))
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
    Group1$type = factor(Group1$type,levels = c('protein_coding','annotated_lncRNA','unannotated_lncRNA','unannotated_coding_transcript')) # nolint # nolint
    
    col_key1<-pal_futurama()(ncol(finally_fpkm_eset))
    col_anno=as.data.frame(colnames(finally_fpkm_eset))
    colnames(col_anno)=c("celltype")
    top_annotation = HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = col_key1),
                       labels = colnames(finally_fpkm_eset),
                       labels_gp = gpar(col = "black", fontsize = 8)),annotation_height=unit(0.4, "cm"))

    col_key2=pal_futurama()(length(table(Group1$type)))
    row_annotation = rowAnnotation(cluster = anno_block(gp = gpar(fill = col_key2)),width=unit(0.4, "cm"))
    #pdf("snakemake/result/test.pdf")
    pdf(HP)
    p1=Heatmap(finally_fpkm_eset,#表达矩阵
  #     col = colorRampPalette(c("navy","white","firebrick3"))(100),#颜色定义
        show_row_names = F,#不展示行名
        row_split = Group1,#用group信息将热土分开，以group聚类
        column_title = "Celltype Heatmap",#不显示列标题
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
}
if(length(unique(scRNA@meta.data$celltype))==1){
    p1 = DimPlot(scRNA, group.by="celltype", label=F, label.size=5, reduction='tsne')
    ggsave(HP, p1, width=7 ,height=6)
}
