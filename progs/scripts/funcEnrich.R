#!/usr/bin/env Rscript 

##------------------------------------------------------------------------------
## I Malki
## Created:  2025-01-10 10:00:00
## Last modification:  2025-01-13 16:29:30  
##------------------------------------------------------------------------------
# funcEnrich.R  -i gene_matrix.diff


doc <- '
Usage: funcEnrich.R -i <input> 

Options:
  -i <INFILE>    the input file (gene expression matrix from Deseq2)
  -h --help      Show this help message.
'

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------


dir.create(Sys.getenv("R_LIBS_USER"), showWarnings=FALSE, recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path


install.packages( "docopt", repos = "http://cran.us.r-project.org")
library("docopt")

# Installez Bioconductor si nécessaire
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  suppressMessages(install.packages("BiocManager", quiet = TRUE))
}

# Chargez les packages nécessaires en supprimant les messages
suppressMessages({
  suppressWarnings({
    # Installer ggplot2 si nécessaire
    if (!require("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2", quiet = TRUE)
    }
    
    # Installer clusterProfiler si nécessaire
    if (!require("clusterProfiler", quietly = TRUE)) {
      BiocManager::install("clusterProfiler", quiet = TRUE)
    }
    
    # Installer dplyr si nécessaire
    if (!require("dplyr", quietly = TRUE)) {
      BiocManager::install("dplyr", quiet = TRUE)
    }

    # Installer DOSE si nécessaire
    if (!require("DOSE", quietly = TRUE)) {
      BiocManager::install("DOSE", quiet = TRUE)
    }
    
    # Installer org.Mm.eg.db si nécessaire
    if (!require("org.Mm.eg.db", quietly = TRUE)) {
      BiocManager::install("org.Mm.eg.db", quiet = TRUE)
    }
    
    # Installer enrichplot si nécessaire
    if (!require("enrichplot", quietly = TRUE)) {
      BiocManager::install("enrichplot", quiet = TRUE)
    }
  })
})


## -----------------------------------------------------------------------------
my_opts <- docopt(doc)

geneList <- my_opts$"-i"
## -----------------------------------------------------------------------------
# Récupérer le dossier parent de l'input
parent_dir <- dirname(geneList)

# Créer un nouveau dossier dans le dossier parent
output_dir <- file.path(parent_dir, "FuncEnrich")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

cat('starting analysis')
#Recuperer en forme de liste les genes significatif , l'ensemblID et la valeur de FC
dat<- read.table(geneList)
dat <- na.omit(dat)

up<- rownames(dat[dat$DESeq2_log2FoldChange > 0 & dat$DESeq2_padj <0.05 ,])
up <- sort(up)

down<- rownames(dat[dat$DESeq2_log2FoldChange < 0 & dat$DESeq2_padj <0.05 ,])
down <- sort (down)


write.csv(up , file = file.path(output_dir,paste0(basename(parent_dir),"_UP_REGULATED.txt")) , row.names= FALSE , quote =FALSE )
write.csv(down , file = file.path(output_dir,paste0(basename(parent_dir),"_DOWN_REGULATED.txt")) , row.names= FALSE , quote =FALSE)

# Enrichissement fonctionnel

#Down regulated
edo_down <- enrichGO(gene = down,
              OrgDb = org.Mm.eg.db,
              keyType = "SYMBOL",  
              ont = "BP",
              pvalueCutoff = 0.5,
              qvalueCutoff = 0.5)

edo_up <- enrichGO(gene = up,
              OrgDb = org.Mm.eg.db,
              keyType = "SYMBOL",  
              ont = "BP",
              pvalueCutoff = 0.5,
              qvalueCutoff = 0.5)

#Barplot

barplot_down <- barplot(edo_down,showCategory = 10 )
barplot_up <- DOSE::barplot(edo_up,showCategory = 10 )
#Dotplot
dotplot_down <- dotplot(edo_down,showCategory=10)
dotplot_up <- dotplot(edo_up,showCategory=10)


#upsetplot
upset_down <- upsetplot(edo_down, showCategory=10)
upset_up <- upsetplot(edo_up, showCategory=10)




#Gene concept network

#p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList_significant , cex_label_gene=0.5, cex_label_category=0.5)
#p3 <- cnetplot(edox, foldChange=geneList_significant, circular = TRUE, colorEdge = TRUE, cex_label_category=0.5 , cex_label_gene=0.5)




# Heatmap functionnal classification
#h1 <- heatplot(edox, showCategory=10)
#h2 <- heatplot(edox,foldChange = geneList_significant, showCategory = 10)
#cowplot::plot_grid(h1, h2, ncol=1, labels=LETTERS[1:2])
plots <- list(barplot_down=barplot_down,barplot_up=barplot_up, dotplot_down=dotplot_down,dotplot_up=dotplot_up ,upsetplot_down=upset_down,upsetplot_up=upset_up)

for (plot_name in names(plots)) {
  output_file <- file.path(output_dir, paste0(basename(parent_dir),"_",plot_name, ".png"))
  
  # Ouvrir le périphérique graphique PNG avec une résolution élevée
  png(output_file, width = 2493, height = 1413 , res=300)
  
  # Ajuster les paramètres graphiques
  par(mar = c(5, 5, 4, 2) + 0.1) # Marges : bas, gauche, haut, droite
  
  # Afficher le graphique
  print(plots[[plot_name]])
  
  # Fermer le périphérique graphique et enregistrer l'image
  dev.off()
  
  # Afficher un message indiquant que le graphique a été enregistré
  cat("Graphique enregistré dans :", output_file, "\n")
}


write.csv(edo_up@result , file = file.path(output_dir,paste0(basename(parent_dir),"_GENE_ONTOLOGY_UPREGULATED_RAW.txt")) , row.names = FALSE , quote=FALSE)


write.csv(edo_down@result , file = file.path(output_dir,paste0(basename(parent_dir),"_GENE_ONTOLOGY_DOWNREGULATED_RAW.txt")) , row.names = FALSE , quote = FALSE)