suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggsci))

# 定义命令行参数
option_list <- list(
    make_option(c("-a", "--fpkm_matrix"), type="character", help="Input FPKM matrix"),
    make_option(c("-b", "--all_id"), type="character", help="Input all ID file"),
    make_option(c("-c", "--protein_coding"), type="character", help="Output protein_coding.csv"),
    make_option(c("-d", "--annotated_lncRNA"), type="character", help="Output annotated_lncRNA.csv"),
    make_option(c("-e", "--unannotated_lncRNA"), type="character", help="Output unannotated_lncRNA.csv"),
    make_option(c("-f", "--unannotated_protein"), type="character", help="Output unannotated_protein.csv"),
    make_option(c("-g", "--fpkm_density"), type="character", help="Output FPKM density plot"),
    make_option(c("-j", "--fpkm_boxplot"), type="character", help="Output FPKM boxplot"),
    make_option(c("-p", "--fpkm_boxplot_png"), type="character", help="Output FPKM boxplot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查参数是否完整
if (any(sapply(opt, is.null))) {
    print_help(opt_parser)
    stop("Missing required arguments.", call.=FALSE)
}

# 读取数据
all_id <- read.table(opt$all_id, header=FALSE, col.names=c("id", "type"))
fpkm_matrix <- read.table(opt$fpkm_matrix, header=TRUE)

# 处理不同类型基因的数据并保存
process_and_save <- function(data, gene_type, output_file) {
    filtered_data <- data %>% filter(type == gene_type)
    write.csv(filtered_data, output_file, row.names=FALSE)
    log_fpkm <- log2(filtered_data[, 7:(ncol(filtered_data) - 1)] + 1)
    log_fpkm$fpkm <- rowSums(log_fpkm)
    log_fpkm$gene_id <- rownames(log_fpkm)
    log_fpkm %>% select(gene_id, fpkm) %>% mutate(group = gene_type)
}

fpkm_density_list <- list(
    process_and_save(fpkm_matrix, "protein_coding", opt$protein_coding),
    process_and_save(fpkm_matrix, "annotated_lncRNA", opt$annotated_lncRNA),
    process_and_save(fpkm_matrix, "unannotated_lncRNA", opt$unannotated_lncRNA),
    process_and_save(fpkm_matrix, "unannotated_coding_transcript", opt$unannotated_protein)
)

fpkm_density <- bind_rows(fpkm_density_list)

# 绘制 FPKM 密度图
fpkm_density_plot <- ggplot(fpkm_density, aes(fpkm, color=group)) +
    xlab("log2(FPKM+1)") +
    stat_ecdf(alpha=1) +
    labs(title="Gene Expression") +
    theme_bw() +
    scale_color_aaas()
ggsave(opt$fpkm_density, fpkm_density_plot, width=7, height=5)

# 计算四分位数并绘制箱线图
ylim1 <- boxplot.stats(fpkm_density$fpkm)$stats[c(1,5)]

fpkm_boxplot <- ggplot(fpkm_density, aes(x=reorder(group, -fpkm, median), y=fpkm, fill=group)) +
    geom_boxplot(alpha=0.7, width=0.3, position=position_dodge(0.2), outlier.shape=NA) +
    scale_y_continuous(name="log2(FPKM+1)") +
    coord_cartesian(ylim = ylim1 * 1.05) +
    scale_x_discrete(name="Group") +
    ggtitle("Boxplot of Expression") +
    theme_bw() +
    scale_fill_aaas() +
    theme(plot.title = element_text(size=14, face="bold"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(size=6, face="bold"),
          legend.position = "none") +
    stat_boxplot(geom='errorbar', width=0.15, aes(color=group)) +
    stat_summary(fun=mean, geom='point', shape=23, size=2, fill='white')

ggsave(opt$fpkm_boxplot, fpkm_boxplot, width=7, height=5)
ggsave(opt$fpkm_boxplot_png, fpkm_boxplot, width=7, height=5)