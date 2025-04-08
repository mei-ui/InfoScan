suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(ggsci))

option_list <- list(
  make_option(c("-a", "--conserved_score"), type="character", help="input conservation Score txt file"),
  make_option(c("-b", "--all_id"), type="character", help="input id Correspondence Table"),
  make_option(c("-c", "--conservation_density"), type="character", help="output conservation_density.pdf"),
  make_option(c("-d", "--conservation_ecdf"), type="character", help="output conservation_ecdf.pdf"),
  make_option(c("-p", "--conservation_ecdf_png"), type="character", help="output conservation_ecdf.png"),
  make_option(c("-e", "--all_gtf"), type="character", help="input all_gene.gtf")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(sapply(opt, is.null))) {
  print_help(opt_parser)
  stop("Please provide all required inputs", call.=FALSE)
}

# Load data
conservedscore <- read.table(opt$conserved_score)
all_id <- read.table(opt$all_id, header=FALSE, col.names=c("id", "type"))
all_gtf <- read.table(opt$all_gtf, header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE, 
                      col.names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

# Extract gene names and transcript IDs
gtf_all <- all_gtf %>%
  filter(feature == "exon") %>%
  mutate(
    gene_name = str_extract(attribute, "gene_name [^;]+") %>% str_remove("gene_name ") %>% str_trim(),
    gene_name = ifelse(
      is.na(gene_name) | gene_name == "nan" | gene_name == "", 
      str_extract(attribute, "transcript_id [^;]+") %>% str_remove("transcript_id ") %>% str_trim(), 
      gene_name
    ),
    transcript_id = str_extract(attribute, "transcript_id [^;]+") %>% str_remove("transcript_id ") %>% str_trim()
  )

# Extract valid transcript IDs
gtf_filtered <- gtf_all %>%
  filter(!is.na(transcript_id) & str_detect(transcript_id, "MSTRG\\.[0-9]+")) %>%
  select(transcript_id, gene_name) %>%
  distinct()

# Merge IDs
id_combined <- bind_rows(gtf_filtered, gtf_all %>% select(transcript_id, gene_name) %>% distinct())
conservedscore <- inner_join(conservedscore, id_combined, by=c("V1"="transcript_id"))
conservedscore <- inner_join(conservedscore, all_id, by=c("gene_name"="id")) %>% distinct()

# Categorize and format data
extract_data <- function(df, type_filter) {
  df %>% filter(type %in% type_filter) %>% select(V1, V6, type) %>%
    rename(transcript_id = V1, conserved_score = V6)
}

protein_coding <- extract_data(conservedscore, "protein_coding")
annotated_lncRNA <- extract_data(conservedscore, "annotated_lncRNA")
unannotated_lncRNA <- extract_data(conservedscore, "unannotated_lncRNA")
unannotated_protein <- extract_data(conservedscore, "unannotated_coding_transcript")

df <- bind_rows(protein_coding, annotated_lncRNA, unannotated_lncRNA, unannotated_protein)

# Generate plots
ggplot(df, aes(conserved_score, color=type)) + 
  xlab("Conserved Score") + 
  geom_density(alpha=1) + 
  labs(title="Gene Conservation") +
  theme_bw() + scale_color_aaas()
ggsave(opt$conservation_density, width=6, height=4)

ggplot(df, aes(conserved_score, color=type)) + 
  xlab("Conserved Score") + 
  stat_ecdf(alpha=1) + 
  labs(title="Gene Conservation") +
  theme_bw() + scale_color_aaas()
ggsave(opt$conservation_ecdf, width=6, height=4)
ggsave(opt$conservation_ecdf_png, width=6, height=4)
