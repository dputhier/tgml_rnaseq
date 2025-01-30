#!/usr/bin/env Rscript 

##------------------------------------------------------------------------------
## I Malki
## Created:  2025-01-10 10:00:00
## Last modification:  2025-01-13 16:29:30  
##------------------------------------------------------------------------------
# funcEnrich.R  -i gene_matrix.diff


doc <- '
Usage: funcEnrich.R -i <input> -l <list>

Options:
  -i <INFILE>    The input file (gene expression matrix from Cuffdiff/Deseq2)
  -l <list>      List of 2 columns with term and gene
  -h --help      Show this help message.
'

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------


dir.create(Sys.getenv("R_LIBS_USER"), showWarnings=FALSE, recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path


install.packages("docopt", repos = "http://cran.us.r-project.org" , quiet=TRUE)
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
    # Installer tidyr si nécessaire
    if (!require("tidyr", quietly = TRUE)) {
      install.packages("tidyr", quiet = TRUE)
    }
	  
    # Installer dplyr si nécessaire
    if (!require("dplyr", quietly = TRUE)) {
      install.packages("dplyr", quiet = TRUE)
    }
	  
    # Installer readxl si nécessaire
    if (!require("readxl", quietly = TRUE)) {
      install.packages("readxl", quiet = TRUE)
    }    
    # Installer clusterProfiler si nécessaire
    if (!require("clusterProfiler", quietly = TRUE)) {
      BiocManager::install("clusterProfiler", quiet = TRUE)
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
path_list <- my_opts$"-l"
## -----------------------------------------------------------------------------
# Récupérer le dossier parent de l'input
parent_dir <- dirname(geneList)

# Créer un nouveau dossier dans le dossier parent
output_dir <- file.path(parent_dir, "FuncEnrichWithList")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


cat('starting analysis')

#Recuperer la liste des genes et pathway associé
if (grepl("\\.xlsx$", path_list)) {  # Vérifie si le fichier est un .xlsx
  TRA <- readxl::read_xlsx(path_list)
} else {  # Sinon, on suppose un fichier texte
  TRA <- read.table(path_list, header = TRUE , sep=",")
}

TRA<-TRA %>%
  pivot_longer(cols = everything(),  values_to = "Gene" , names_to = "Pathway",) %>%
  filter(Gene != "")

#Recuperer en forme de liste les genes significatif , l'ensemblID et la valeur de FC
geneL <- read.table(geneList , header = TRUE)
gene_significant <- geneL[geneL$significant == "yes", ]
geneList_significant <- setNames(gene_significant$log2.fold_change., gene_significant$gene)
#Recuperer seulement avec un |FC|>2
de <- names(geneList_significant)[abs(geneList_significant) > 2]

# Enrichissement fonctionnel
edo <- enricher(gene = de,
              pAdjustMethod = "BH",
              pvalueCutoff = 1, TERM2GENE=TRA,
              qvalueCutoff = 1 )

#Barplot

barplot <- barplot(edo,showCategory = 40)+theme_minimal()

#Dotplot
dotplot <- dotplot(edo ,showCategory=40)+theme_minimal()


#upsetplot
upset <- upsetplot(edo)+theme_minimal()

plots <- list(barplot=barplot, dotplot=dotplot ,upsetplot=upset )

for (plot_name in names(plots)) {
  output_file <- file.path(output_dir, paste0(plot_name, ".png"))
  
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