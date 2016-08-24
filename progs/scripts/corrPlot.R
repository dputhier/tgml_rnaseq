#!/usr/bin/Rscript 

##------------------------------------------------------------------------------
## D puthier
## Created: Fri Jun 26 00:16:06 BRT 2015
## Last modification: Wed Aug 13 17:52:46 CEST 2014
##------------------------------------------------------------------------------
# corrPlot.R  -i my_matrix.txt

"usage: coorPlot.R -i <INFILE>  [-o <OUTPREFIX> <EXPANSION> -w <WIDTH> -y <HEIGHT> -l <LABEL>	-t <LOG>] 
		
		options:
		-i <INFILE>		the input file (genes as row, samples as column).
		-o <OUTPREFIX>		A prefix for output (e.g. /tmp/prefix) [default: ./].
		-w <WIDTH>      png width [default: 6].
		-y <HEIGHT>    	png height [default: 6].
		-x <EXPANSION>		Character expansion [default: 1]
		-l <LABEL>		Sample labels (A,B,C,...).
		-t <LOG>    Perform log2 transformation (and add 1 as pseudo-count) [default: 0]
		" -> doc



## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------
# http://r.789695.n4.nabble.com/Install-package-automatically-if-not-there-td2267532.html

load.fun <- function(x) { 
	x <- as.character(substitute(x)) 
	if(isTRUE(x %in% .packages(all.available=TRUE))) { 
		suppressWarnings(eval(parse(text=paste("require(", x, ")", sep=""))))
	} else { 
		eval(parse(text=paste("install.packages('", x, "', repos = 'http://cran.us.r-project.org')", sep=""))) 
		suppressWarnings(eval(parse(text=paste("require(", x, ")", sep=""))))
	} 
} 

load.bioc <- function(x) { 
	x <- as.character(substitute(x)) 
	if(isTRUE(x %in% .packages(all.available=TRUE))) { 
		eval(parse(text=paste("require(", x, ")", sep=""))) 
	} else { 
		eval(parse(text="source('http://bioconductor.org/biocLite.R'"))
		eval(parse(text=paste("biocLite('", x, "')", sep="")))
		suppressWarnings(eval(parse(text=paste("require(", x, ")", sep=""))))
	} 
} 

suppressMessages(load.bioc("geneplotter"))

suppressMessages(load.fun("docopt"))


## -----------------------------------------------------------------------------
## Arg parsing
## -----------------------------------------------------------------------------



my_opts <- docopt(doc)


if(length(commandArgs()) < 5 ){
	print(commandArgs())
	cat("Use -h for more informations\n")
	q(status=1)
}


## -----------------------------------------------------------------------------
## Loading libs
## -----------------------------------------------------------------------------
suppressMessages(load.fun(corrplot))

## -----------------------------------------------------------------------------
## Parsing args
## -----------------------------------------------------------------------------


myfile <- my_opts$"-i"
expansion <- as.double(my_opts$"-x")
outprefix <- my_opts$"-o"
pngheigth <- as.numeric(my_opts$"-y")
pngwidth <- as.numeric(my_opts$"-w")
label <- my_opts$"-l"


## -----------------------------------------------------------------------------
## Loading files
## -----------------------------------------------------------------------------


d <- read.table(myfile, sep="\t", 
		header=TRUE, quote="", 
		comment.char="", row=1)
d <- d[apply(d,1,sum) > 0, ]

if(is.null(label)){
	label <- colnames(d)	
}

if(my_opts$"-t" == "1"){
	cat("Performing log2 transformation")
	d <- log2(d + 1)
}



## -----------------------------------------------------------------------------
## Parameters
## -----------------------------------------------------------------------------

colnames(d) <- label

corr <- cor(d, method="pearson")

png(file.path(outprefix,"CoorPlot_circle.png"), 
		units = 'in', res = 300, width=pngwidth, height=pngheigth)
corrplot(corr, method = "circle", order = "hclust", hclust.method ="average", type="lower", tl.cex=expansion,
		tl.col="black", tl.srt = 45)
devnull <- dev.off()


png(file.path(outprefix,"CoorPlot_square.png"), 
		units = 'in', res = 300, width=pngwidth, height=pngheigth)
corrplot(corr, method = "square", order = "hclust", hclust.method ="average", type="lower", tl.cex=expansion,
		tl.col="black", tl.srt = 45)
devnull <- dev.off()

png(file.path(outprefix,"CoorPlot_ellipe.png"), 
		units = 'in', res = 300, width=pngwidth, height=pngheigth)
corrplot(corr, method = "ellipse", order = "hclust", hclust.method ="average", type="lower", tl.cex=expansion,
		tl.col="black", tl.srt = 45)
devnull <- dev.off()

png(file.path(outprefix,"CoorPlot_square.png"), 
		units = 'in', res = 300, width=pngwidth, height=pngheigth)
corrplot(corr, method = "number", order = "hclust", hclust.method ="average", type="lower", tl.cex=expansion,
		tl.col="black", tl.srt = 45)
devnull <- dev.off()

png(file.path(outprefix,"CoorPlot_pie.png"), 
		units = 'in', res = 300, width=pngwidth, height=pngheigth)
corrplot(corr, method = "pie", order = "hclust", hclust.method ="average", type="lower", tl.cex=expansion,
		tl.col="black", tl.srt = 45)
devnull <- dev.off()

wb <- c("white","black")

png(file.path(outprefix,"CoorPlot_piewb.png"), 
		units = 'in', res = 300, width=pngwidth, height=pngheigth)
corrplot(corr, method = "pie", col = wb, bg="gold2", order = "hclust", hclust.method ="average", type="lower", 
		tl.cex=expansion, tl.col="black", tl.srt = 45)
devnull <- dev.off()

png(file.path(outprefix,"Pairs_plot.png"), 
		units = 'in', res = 300, width=pngwidth, height=pngheigth)
plotFun <- function(x,y){ dns <- densCols(x,y); points(x,y, col=dns, pch=".") }
pairs(d, upper.panel=plotFun, lower.panel = NULL)
devnull <- dev.off()
