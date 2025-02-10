library(ggplot2)
library(dplyr)
library(ggsci)

setwd("/shared/projects/project_dputhier/043_irla_elife_2022/progs/downstream_analyses")


# Get info about gene biotypes 
# in the whole genome. Only select
# gene_biotypes to which at least 500 genes are associated 
# (thus the most common ones)

biotypes <- read.table("../../output/gene_info/biotypes.tsv", sep="\t", he=T)

biotypes <- biotypes[, c("gene_name", "gene_biotype")]
biotypes$Context <- "Whole\nGenome"

biotypes <- biotypes %>%
  group_by(gene_biotype) %>%
  filter(n() >= 500)

# Get Aire regulated

aire_reg <- read.table("../../output/comparison/RipmOVAvsOTII_2/RipmOVAvsOTII_2_DESeq2_EdgeR_pval_and_norm_count_log2.txt", sep="\t", row=1, head=TRUE) 
aire_reg <- aire_reg[-grep("XLOC", rownames(aire_reg)), ]
aire_reg <- aire_reg[rownames(aire_reg) %in% biotypes$gene_name, ]
aire_reg <- aire_reg[order(aire_reg$DESeq2_padj),]

# Get Aire induced 

aire_pos <- aire_reg[aire_reg$DESeq2_log2FoldChange < 0, ]
aire_pos <- aire_pos[aire_pos$DESeq2_padj < 0.05 & !is.na(aire_pos$DESeq2_padj),]
aire_pos <- rownames(aire_pos)

hit <- match(aire_pos, biotypes$gene_name)
df_pos <- data.frame(gene_name=aire_pos,  gene_biotype=biotypes$gene_biotype[hit], Context="RIPmOVA\ninduced")

# Get Aire repressed

aire_neg <- aire_reg[aire_reg$DESeq2_log2FoldChange > 0, ]
aire_neg <- aire_neg[aire_neg$DESeq2_padj < 0.05 & !is.na(aire_neg$DESeq2_padj),]
aire_neg <- rownames(aire_neg)

hit <- match(aire_neg, biotypes$gene_name)
df_neg <- data.frame(gene_name=aire_neg,  gene_biotype=biotypes$gene_biotype[hit], Context="RIPmOVA\nrepressed")

# Bind info about all genes and those 
# that are R regulated

all <- rbind(biotypes, df_neg, df_pos)


# Plot a diagram 

ggplot(data=all, mapping=aes(x=Context, fill=gene_biotype))  +
  geom_bar(position="fill") + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values=ggsci::pal_bmj()(length(table(biotypes$gene_biotype)))) +
  labs(x="Gene set", legend="Gene\nBiotype")

