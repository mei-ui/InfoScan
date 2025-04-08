suppressMessages(library("stringr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("stringi"))
suppressMessages(library("optparse"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggsci"))

option_list <- list(
  make_option(c("-a", "--fpkm_matrix"), type="character", default=NULL, help="input finally fpkm"),
  make_option(c("-c", "--all_gtf"), type="character", default=NULL, help="input all.gtf"),
  make_option(c("-d", "--protein_gtf"), type="character", default=NULL, help="output protein_coding.gtf"),
  make_option(c("-e", "--annotated_lncRNA"), type="character", default=NULL, help="output annotated_lncRNA.gtf"),
  make_option(c("-f", "--unannotated_lncRNA"), type="character", default=NULL, help="output unannotated_lncRNA.gtf"),
  make_option(c("-g", "--unprotein"), type="character", default=NULL, help="output unannotated_protein.gtf"),
  make_option(c("-q", "--exon_filter"), type="integer", default=1, help="exon number filter"),
  make_option(c("-o", "--exon_number"), type="character", default=NULL, help="output exon_number.pdf"),
  make_option(c("-p", "--exon_number_png"), type="character", default=NULL, help="output exon_number.png"),
  make_option(c("-O", "--finally_fpkm"), type="character", default=NULL, help="output filtered fpkm result")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check for missing arguments
required_args <- c("fpkm_matrix", "all_gtf", "exon_number", "protein_gtf", 
                   "annotated_lncRNA", "unannotated_lncRNA", "unprotein", "exon_filter", "finally_fpkm")
if (any(sapply(required_args, function(arg) is.null(opt[[arg]])))) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

# Read input data
finally_fpkm <- read.table(opt$fpkm_matrix, header=TRUE, sep="\t")
all_gtf <- read.table(opt$all_gtf, header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
colnames(all_gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Extract gene names and transcript IDs
gtf_all <- all_gtf %>%
  filter(feature == "exon") %>%
  mutate(gene_name = str_extract(attribute, "gene_name [^;]+") %>% str_remove("gene_name ") %>% str_trim(),
         gene_name = ifelse(gene_name == "nan" | is.na(gene_name), 
                            str_extract(attribute, "transcript_id [^;]+") %>% str_remove("transcript_id ") %>% str_trim(), 
                            gene_name),
         transcript_id = str_extract(attribute, "transcript_id [^;]+") %>% str_remove("transcript_id ") %>% str_trim())
head(gtf_all)
process_gtf <- function(fpkm_data, type_label, gtf_data, exon_filter=0) {
  selected_genes <- fpkm_data %>% filter(type == type_label) %>% pull(Geneid)
  gtf_filtered <- gtf_data %>% filter(gene_name %in% selected_genes)
  exon_counts <- gtf_filtered %>% count(transcript_id, name="Freq") %>% filter(Freq > exon_filter)
  gtf_final <- gtf_filtered %>% filter(transcript_id %in% exon_counts$transcript_id)
  #write.table(gtf_final, file=output_gtf, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  exon_counts %>% mutate(group = type_label)
}
# Process each gene type
protein_exon <- process_gtf(finally_fpkm, "protein_coding", gtf_all)
annotated_exon <- process_gtf(finally_fpkm, "annotated_lncRNA", gtf_all)
unannotated_lncRNA_exon <- process_gtf(finally_fpkm, "unannotated_lncRNA", gtf_all, opt$exon_filter)
head(unannotated_lncRNA_exon)
unannotated_protein_exon <- process_gtf(finally_fpkm, "unannotated_coding_transcript", gtf_all)

gtf_all2 <- all_gtf %>%
  mutate(gene_name = str_extract(attribute, "gene_name [^;]+") %>% str_remove("gene_name ") %>% str_trim(),
         gene_name = ifelse(gene_name == "nan" | is.na(gene_name), 
                            str_extract(attribute, "transcript_id [^;]+") %>% str_remove("transcript_id ") %>% str_trim(), 
                            gene_name),
         transcript_id = str_extract(attribute, "transcript_id [^;]+") %>% str_remove("transcript_id ") %>% str_trim())
head(gtf_all2)
process_gtf2 <- function(fpkm_data, type_label, gtf_data, output_gtf,unannotated_lncRNA_id = NULL) {
  # 规范属性字段生成函数
  build_attribute <- function(gene_id = "na", transcript_id, gene_name) {
    elements <- c(
      if(!is.na(gene_id)) sprintf('gene_id "%s"', gene_id),
      if(!is.na(transcript_id)) sprintf('transcript_id "%s"', transcript_id),
      if(!is.na(gene_name) && gene_name != "nan") sprintf('gene_name "%s"', gene_name)
    )
    paste(elements[!is.na(elements)], collapse = "; ") %>% paste0(";")
  }

  selected_genes <- fpkm_data %>% filter(type == type_label) %>% pull(Geneid)
  gtf_filtered <- gtf_data %>% filter(gene_name %in% selected_genes)
  if (is.null(unannotated_lncRNA_id)) {
    gtf_final = gtf_filtered[,c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute","transcript_id","gene_name")]
    } else{
      gtf_final <- gtf_filtered %>% filter(transcript_id %in% unannotated_lncRNA_id)
      gtf_final = gtf_final[,c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute","transcript_id","gene_name")]
    }
  # 重建规范化的GTF结构
   gtf_final_1 <- gtf_final %>%
    mutate(
      attribute = pmap_chr(list(gene_id = "na",transcript_id, gene_name), build_attribute),
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    select(seqname, source, feature, start, end, score, strand, frame, attribute) %>%
    arrange(seqname, start)  # 按基因组坐标排序

  write.table(gtf_final_1, file=output_gtf, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  return(gtf_final_1)
}
protein_gtf_1 <- process_gtf2(finally_fpkm, "protein_coding", gtf_all2, opt$protein_gtf)
head(protein_gtf_1)
annotated_gtf_1 <- process_gtf2(finally_fpkm, "annotated_lncRNA", gtf_all2, opt$annotated_lncRNA)
unannotated_lncRNA_gtf_1 <- process_gtf2(finally_fpkm, "unannotated_lncRNA", gtf_all2, opt$unannotated_lncRNA,unique(unannotated_lncRNA_exon$transcript_id))
head(unannotated_lncRNA_gtf_1)
unannotated_protein_gtf_1<- process_gtf2(finally_fpkm, "unannotated_coding_transcript", gtf_all2, opt$unprotein)
# Combine exon numbers
exon_number <- bind_rows(protein_exon, annotated_exon, unannotated_lncRNA_exon, unannotated_protein_exon)
dup <- exon_number %>% filter(Freq < 50)

# Merge FPKM results
table(finally_fpkm$type)
annotated_fpkm <- finally_fpkm %>% filter(type %in% c("protein_coding","annotated_lncRNA","pseudo_gene","other_rna"))
unannotated_fpkm <- finally_fpkm %>% filter(Geneid %in% unique(c(unannotated_lncRNA_exon$transcript_id, unannotated_protein_exon$transcript_id)))
finally_fpkm_result <- bind_rows(annotated_fpkm, unannotated_fpkm)

# Add ERCC genes
ercc_fpkm <- finally_fpkm %>% filter(str_detect(Geneid, "^ERCC-"))
finally_fpkm_final <- bind_rows(finally_fpkm_result, ercc_fpkm)
table(finally_fpkm_final$type)
write.table(finally_fpkm_final, opt$finally_fpkm, sep='\t', quote=FALSE, row.names=FALSE)

# Generate exon count plot
ggplot(dup, aes(x=Freq, fill=group)) +
  geom_bar() +
  labs(title="Exon Number", fill="Gene Type") +
  scale_fill_aaas() +
  scale_x_continuous(breaks=seq(1, 50, by=3)) +
  ylab("Total Transcript Number") +
  xlab("Exon Number") +
  theme_bw()
ggsave(opt$exon_number, width=8, height=4)
ggsave(opt$exon_number_png, width=8, height=4)
