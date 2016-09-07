rule diagnostic_plot:
        input: "output/quantification_known_genes/gene_counts_known_and_novel_mini.txt"
        output: "output/diagnostic_plot/diagnostic.pdf"
        threads: 1
        run: R("""
            load.bioc <- function(x) {{ 
                x <- as.character(substitute(x)) 
                if(isTRUE(x %in% .packages(all.available=TRUE))) {{ 
                    eval(parse(text=paste("require(", x, ")", sep=""))) 
                }} else {{ 
                    eval(parse(text="source('http://bioconductor.org/biocLite.R'"))
                    eval(parse(text=paste("biocLite('", x, "')", sep="")))
                }} 
            }} 
            
            suppressWarnings(suppressMessages(load.bioc(geneplotter)))

            library(geneplotter)

            data <- read.table("{input}", 
                                sep="\t", 
                                header=T, 
                                row.names=1)
            data <- data[rowSums(data) > 0, ]
            data <- log2(data + 1)
            
            pdf("{output}")
                dev.null <- apply(data, 2, hist, border="white", col="blue")
                boxplot(data, color="blue", pch=16)
                myDisplayFunction <- function(x,y){{ smoothScatter(x, y , pch=".", add=TRUE) }}
                pairs(data, upper.panel=myDisplayFunction,  lower.panel=NULL)
            dev.off()
      
        """)