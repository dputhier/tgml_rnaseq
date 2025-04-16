# Penser à indiquer le sha du dernier commit
# de la branche main/master
# devtools::install_github("dputhier/scigenex")

library(scigenex)
library(ggplot2)

d <-  read.table("/shared/projects/project_dputhier/043_irla_elife_2022/output/comparison/ALL_2/ALL_2_DESeq2_norm_count_table_log2.txt", sep="\t", head=TRUE, row=1)
colnames(d) <- paste(gsub("\\..*", "", colnames(d)), "_", 1:ncol(d), sep="")

# Delete unknown genes
d <- d[-grep("XLOC", rownames(d)), ]

# Try to limit the expression matrix
# to well quantified genes
d <- d[rowSums(d) > 40, ]

# a <- 1
# b <- 10
# fct <- 2
# plot(d[, a], d[, b], pch=".")
# lfc <- d[,b] - d[,a]
# text(d[abs(lfc) > fct, a], 
#      d[abs(lfc) > fct, b], 
#      label=rownames(d)[abs(lfc) > fct], 
#      cex=0.5,
#      pos=sample(1:4, sum(abs(lfc) > fct), rep=T))
# pairs(d, pch=".")

# Select genes that are well correlated with others
# and probably part of a cluster

gs <- select_genes(d, k = 90, noise_level = 0.75, 
                   distance_method = "pearson", 
                   no_anti_cor = TRUE)

# Construct a graph and partition it with MCL.
gc <- gene_clustering(gs, 
                      method="reciprocal_neighborhood", 
                      inflation = 1.6, 
                      threads = 4)

# Filter the gene module by size.
fgc <- scigenex::filter_cluster_size(gc, min_cluster_size = 15)

# Filter the gene module by std dev 
# (trying to restrict to those having most dispersion/variance)
fgc <- scigenex::filter_cluster_sd(fgc, min_sd = 1.5)
clust_size(fgc)

# Compute gene that are the most representative 
# of each cluster
fgc <-  top_genes(fgc, top = 75)

# Create a vector of classes
class <- gsub("_.*", "", colnames(d))
class <- as.numeric(as.factor(class))
class <- setNames(class, colnames(d))

# Plotting the most representative genes
plot_ggheatmap(fgc, use_top_genes = TRUE, ident = class)

# Plotting the mean expression profile
# inside each gene module.
sort(class)
plot_profiles(fgc, ident=class)

# Same in linear scale.
fgc_lin <- fgc
fgc_lin@dbf_output$center <- 2^fgc_lin@dbf_output$center
plot_profiles(fgc_lin, ident=class)

# Search genes
which_clust(fgc, "Aire")
grep_clust(fgc, "Fox")

plot_ggheatmap(fgc[1,], 
               use_top_genes = TRUE, 
               ident = class, hide_gene_name = F)  + 
  theme(axis.text = element_text(size=5))

cat(paste0(sort(fgc[11,]@gene_clusters[[1]]), collapse = "\n"))


# DNA Binding
go_0003677 <- top_by_go(fgc[1:3,], 
                        go_id = "GO:0003677", 
                        species="mmusculus",
                        gene_id="mgi_symbol")

go_0007166 <- top_by_go(fgc, 
                        go_id = "GO:0003677", 
                        species="mmusculus",
                        gene_id="mgi_symbol")

go_0007166 <- top_by_go(fgc, "GO:0007166", 
                        species="mmusculus",
                        gene_id="mgi_symbol")


go_0005125 <- top_by_go(fgc, "GO:0005125", 
                        species="mmusculus",
                        gene_id="mgi_symbol")


plot_ggheatmap(fgc, 
               use_top_genes = TRUE, 
               ident = class, hide_gene_name = F)  + 
  theme(axis.text = element_text(size=5))
