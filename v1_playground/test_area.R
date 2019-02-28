
setwd("C:/OConnell Lab/R analyses/Shiny Apps/Cell_identity_predictor/CIPR_deploy/v1_playground")


reference_log <- readRDS("data/immgen_v1_reference_log.rds")

ref_annotation <- readRDS("data/imm_annot.rds") 


dat <- read.csv("data/Trimmed_cluster_signatures.csv",
                check.names = F,
                strip.white = T,
                stringsAsFactors = F)




sel_clst <- dat %>%
  filter(!!rlang::sym(cluster_column) == "1") %>%
  select_(.dots = c(gene_column, logFC_column))


# Merge SCseq cluster log FC value with immgen log FC for shared genes
merged <- merge(sel_clst, reference_log, by.x = gene_column, by.y = ref_gene_column)





reference <- read.csv("../../misc/ImmGenData20180420_18_48_38.csv", check.names=FALSE, strip.white = TRUE, stringsAsFactors = F)


reference <- read.csv("../../misc/Immgen_v2_020419.csv", check.names=FALSE, strip.white = TRUE, stringsAsFactors = F)


reference <- read.csv("../../misc/immgen_combined.csv", check.names=FALSE, strip.white = TRUE, stringsAsFactors = F, row)

reference <- readRDS("../../misc/immgen_combined.rds")
