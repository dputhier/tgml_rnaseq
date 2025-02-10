#!/usr/bin/env Rscript 

##------------------------------------------------------------------------------
## I Malki
## Created:  2025-01-17 14:36:00
## Last modification:  2025-01-17 14:36:30  
##------------------------------------------------------------------------------
# patchwork.R -i plots


# Définir la documentation pour les options de commande
doc <- '
Usage: funcEnrich.R -i <input>... [options]

Options:
  -i <input>    Input files (one or more files).
  -h --help     Show this help message and exit.
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
	# Installer patchwork si nécessaire
    if (!require("patchwork", quietly = TRUE)) {
      install.packages("patchwork", quiet = TRUE)
    }
 })
})


# Parser les arguments avec docopt
args <- docopt(doc)

# Afficher les arguments pour vérification
print(args)

# Obtenir les fichiers d'entrée
input_files <- args$`<input>`

# Vérifier qu'au moins un fichier est fourni
if (length(input_files) == 0) {
  stop("Error: No input files provided. Use -i <input> to specify input files.")
}


generate_plot <- function(file) {
  # Exemple : lecture des données
  plot = readRDS(file)
}



plots= lapply(input_file , generate_plot)
combined_plot <- Reduce(`+`, plots)


parent_dir <- dirname(input_file)

# Créer un nouveau dossier dans le dossier parent
output_dir <- file.path(parent_dir, "Plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

output_file <- file.path(output_dir , "Plots.png")
ggsave(output_file, plot = combined_plot, width = 10, height = 8, dpi = 300)
cat("Plot saved to:", output_file, "\n")

