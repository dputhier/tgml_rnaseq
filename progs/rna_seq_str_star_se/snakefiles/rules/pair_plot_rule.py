rule diagnostic_plot:
        input: "output/quantification_known_genes/gene_counts_known_and_novel_mini.txt"
        output: "output/diagnostic_plot/diagnostic.pdf"
        threads: 1
        params: mem="4G"
        shell: """ 
            module load r/4.4.1
            echo "suppressWarnings(suppressMessages(load.bioc(geneplotter))); \
            library(geneplotter) ;\
            data <- read.table('{input}',  sep='\t', header=T, row.names=1) ;\
            data <- data[rowSums(data) > 0, ] ; \
            data <- log2(data + 1); \
            pdf('{output}'); \
                dev.null <- apply(data, 2, hist, border='white', col='blue'); \
                boxplot(data, color='blue', pch=16); \
                myDisplayFunction <- function(x,y){{ smoothScatter(x, y , pch='.', add=TRUE) }}; \
                pairs(data, upper.panel=myDisplayFunction,  lower.panel=NULL); \
            dev.off(); " | R --slave
      
        """