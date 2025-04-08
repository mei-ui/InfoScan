suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggsci))

# Define options
option_list <- list(
  make_option(c("-a", "--fpkm_matrix"), type="character", help="Input FPKM matrix"),
  make_option(c("-b", "--all_id"), type="character", help="Input all ID file"),
  make_option(c("-o", "--length_distribution"), type="character", help="Output length_distribution.pdf"),
  make_option(c("-c", "--before_filter"), type="character", help="Output Before FPKM filter"), 
  make_option(c("-d", "--after_filter"), type="character", help="Output After FPKM filter"),
  make_option(c("-e", "--length_distribution_png"), type="character", help="Output length_distribution.png"),
  make_option(c("-f", "--before_filter_png"), type="character", help="Output Before FPKM filter"), 
  make_option(c("-g", "--after_filter_png"), type="character", help="Output After FPKM filter")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (any(sapply(opt, is.null))) {
  print_help(opt_parser)
  stop("Please provide all required inputs.", call.=FALSE)
}

# Load data
all_id <- read.table(opt$all_id, header=FALSE, col.names=c("id", "type"))
fpkm_matrix <- read.table(opt$fpkm_matrix, header=TRUE)

# Filter by gene types and extract relevant columns
gene_types1 <- c("protein_coding", "annotated_lncRNA", "unannotated_lncRNA", "unannotated_coding_transcript")
table(fpkm_matrix$type)
length_data <- fpkm_matrix %>%
  filter(type %in% gene_types1) %>%
  select(Geneid, Length, type)

# Plot Length Distribution
ggplot(length_data, aes(Length, color=type)) + 
  xlab("Length") + 
  geom_density(alpha=1) + 
  scale_color_aaas() +
  scale_x_continuous(limits = c(0, 15000)) +
  theme_bw()
ggsave(opt$length_distribution, width=7, height=5)
ggsave(opt$length_distribution_png, width=7, height=5)

# Function to generate pie chart
generate_pie_chart <- function(data, title, output_file) {
  summary_data <- data %>%
    count(type, name="Freq") %>%
    arrange(desc(Freq)) %>%
    mutate(Label = paste0(type, " (", Freq, " - ", round(Freq / sum(Freq) * 100, 2), "%)"))
  
  p <- ggplot(summary_data, aes(x="", y=Freq, fill=type)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y") +
    labs(x="", y="", title=title) +
    theme_void() +
    theme(legend.title=element_blank(), legend.position="right") +
    scale_fill_discrete(breaks=summary_data$type, labels=summary_data$Label)
  
  ggsave(output_file, p, width=7, height=5)
}

# Generate Pie Charts (Before and After FPKM filtering)
gene_types2 <- c("protein_coding", "annotated_lncRNA", "pseudo_gene", "other_rna", 
                "unannotated_lncRNA", "unannotated_coding_transcript")
generate_pie_chart(all_id %>% filter(type %in% gene_types2), "Before FPKM and Exon Number Filter", opt$before_filter)
generate_pie_chart(fpkm_matrix %>% filter(type %in% gene_types2), "After FPKM and Exon Number Filter", opt$after_filter)

generate_pie_chart(all_id %>% filter(type %in% gene_types2), "Before FPKM and Exon Number Filter", opt$before_filter_png)
generate_pie_chart(fpkm_matrix %>% filter(type %in% gene_types2), "After FPKM and Exon Number Filter", opt$after_filter_png)
