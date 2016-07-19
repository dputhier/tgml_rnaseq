#!/usr/bin/Rscript 

##------------------------------------------------------------------------------
## D puthier
## Created: Tue Jul  8 15:09:18 CEST 2014
## Last modification: Wed Aug 13 17:52:46 CEST 2014
##------------------------------------------------------------------------------
#todo: http://stackoverflow.com/questions/20260434/test-significance-of-clusters-on-a-pca-plot
# Perform PCA with various R libraries (ADE4, ggplot2,...). Several output are proposed.
# Example usage: pca_plot.R -i matrix.txt -f pheno.txt


"usage: pca_plot.R -i <INFILE> -f <FILEPHENO> [-o <OUTPREFIX> -z -d <DELETE> -x <EXPANSION> -l <LABEL>  -t -q -w <WIDTH> -y <HEIGHT> -P <PALETTE>] 

options:
 -i <INFILE>		the input file (genes as row, samples as column).
 -f <FILEPHENO>		A file with two or three columns (no header): label, phenotype, color (html).
 -o <OUTPREFIX>		A prefix for output (e.g. /tmp/prefix) [default: ./].
 -d <DELETE>		Delete this string from sample names when plotting (e.g .txt) [default:]
 -x <EXPANSION>		Character expansion [default: 1]
 -t             	Indicate that the FILEPHENO is transposed. The row 1,2 and 3 corresponds to label, phenotype, color.
 -v       		Run in verbose mode.
 -q     		Perform TMM normalization.
 -z      	  	Compute log2(data+1).
 -w <WIDTH>    	png width [default: 6].
 -y <HEIGHT>    	png height [default: 6].
 -P <PALETTE>    	A palette: rob (red, orange, blue), rwb YlOrBr, jet, or rainbow. (need -f whitout third column) [default: rob].
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
	} 
} 

## -----------------------------------------------------------------------------
## Arg parsing
## -----------------------------------------------------------------------------

suppressMessages(load.fun("docopt"))

my_opts <- docopt(doc)


## -----------------------------------------------------------------------------
## Loading libs
## -----------------------------------------------------------------------------
suppressMessages(library(affy))
suppressMessages(load.fun(ade4))
suppressMessages(load.fun(rgl))
suppressMessages(load.fun(ggplot2))
suppressMessages(load.fun(proto))
suppressMessages(load.fun(grid))
suppressMessages(load.fun(ellipse))
suppressMessages(load.fun(plyr))


## -----------------------------------------------------------------------------
## Loading files
## -----------------------------------------------------------------------------


myfile <- my_opts$"-i"
pngwidth <- as.numeric(my_opts$"-w")
pngheigth <- as.numeric(my_opts$"-y")
expansion <- as.numeric(my_opts$"-x")
mypalette <- my_opts$"-P"
outprefix <- my_opts$"-o"
fileinfo <- my_opts$"-f"


d <- read.table(myfile, sep="\t", 
		header=TRUE, quote="", 
		comment.char="", row=1)

## -----------------------------------------------------------------------------
## Selecting expressed genes
## -----------------------------------------------------------------------------

d <- d[rowSums(d) > 0,]

## -----------------------------------------------------------------------------
## Log2 transformation
## -----------------------------------------------------------------------------


if(my_opts$"-z"){
	cat(" --> Performing log transformation.\n")
	d <- log2(d+1)
}

if(is.null(mypalette))
	mypalette="rob"

## -----------------------------------------------------------------------------
## Normalization
## -----------------------------------------------------------------------------


### Let's implement a function for normalization (~TMM)

estimSf <- function (cts){
	
	# Compute the geometric mean
	geomMean <- function(x) exp(1/length(x)*sum(log(x)))
	
	# Compute the geometric mean over the line
	gm.mean  <-  apply(cts, 1, geomMean)
	
	# Zero values are set to NA (avoid subsequentcdsdivision by 0)
	gm.mean[gm.mean == 0] <- NA
	
	# Divide each line by its corresponding geometric mean
	# sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
	# MARGIN: 1 or 2 (line or columns)
	# STATS: a vector of length nrow(x) or ncol(x), depending on MARGIN
	# FUN: the function to be applied
	cts <- sweep(cts, 1, gm.mean, FUN="/")
	
	# Compute the median over the columns
	med <- apply(cts, 2, median, na.rm=TRUE)
	head(cts)
	# Return the scaling factor
	return(med)
	
}

if(my_opts$"-q"){
	
	cat(" --> Normalizing.\n")
	png(file.path(outprefix, "PCA_Boxplot_before.png"))
	boxplot(log2(d+1), pch=".")
	devnull <-  dev.off()
	
	png(file.path(outprefix, "PCA_Density_before.png"))
	plotDensity(log2(d+1))
	devnull <-  dev.off()
	
	sf <- estimSf(d)
	cat(" --> Scaling factors:", sf, "\n")
	d <- round(sweep(d, 2, sf, "/"),0)
	
	
	png(file.path(outprefix, "PCA_Boxplot_after.png"))
	boxplot(log2(d+1), pch=".")
	devnull <- dev.off()	
	
	png(file.path(outprefix, "PCA_Density_after.png"))
	plotDensity(log2(d+1))
	devnull <-  dev.off()
}


## -----------------------------------------------------------------------------
## Parameters
## -----------------------------------------------------------------------------


if(my_opts$"-t" == FALSE | is.null(my_opts$"-t") == FALSE){
	
	
	info <- read.table(fileinfo, sep="\t")
	label   <- as.character(info[,1])
	mypheno <- as.character(info[,2])
	
	res <- try(mycolor <- as.character(info[,3]) ,silent = TRUE)
	if(class(res) == "try-error"){
		mycolor <- NULL
	}

}else{
	info <- as.matrix(read.table(fileinfo, sep="\t", comment.char=""))
	label   <- info[1,]
	mypheno <- as.character(info[2,])
	res <- try(mycolor <- as.character(info[3,]) ,silent = TRUE)
	if(class(res) == "try-error"){
		mycolor <- NULL
	}
}


nb <- length(levels(as.factor(mypheno)))

if(is.null(mycolor)){

	mycolor <- as.factor(mypheno)
	
	if(mypalette=="rob"){
		Lab.palette <- colorRampPalette(c("red", "orange", "blue"), space="Lab")
		levels(mycolor) <- Lab.palette(nb)
			
	}else if(mypalette=="rwb"){
		Lab.palette <- colorRampPalette(c("red", "white", "blue"), space="Lab")
		levels(mycolor) <- Lab.palette(nb)
	}else if(mypalette=="YlOrBr"){
		YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
		Lab.palette <- colorRampPalette(YlOrBr, space="Lab")
		levels(mycolor) <- Lab.palette(nb)
		
	}else if(mypalette=="jet"){
		jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
						"#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		Lab.palette <- colorRampPalette(jet, space="Lab")
		levels(mycolor) <- Lab.palette(nb)
		
	}else if(mypalette=="rainbow"){
		levels(mycolor) <- rainbow(length(levels(mycolor)))
	}else{
		levels(mycolor) <- rainbow(length(levels(mycolor)))	
	}

	mycolor <- as.character(mycolor)
}


names(mypheno) <- mycolor
names(mycolor) <- mypheno

## -----------------------------------------------------------------------------
## Redefine ADE4 functions
## -----------------------------------------------------------------------------

scatterutil.eti2 <- function (x, y, label, clabel, boxes = TRUE, coul = rep(1, length(x)), 
		horizontal = TRUE, bg = "white") 
{
	if (length(label) == 0) 
		return(invisible())
	if (is.null(label)) 
		return(invisible())
	if (any(label == "")) 
		return(invisible())
	cex0 <- par("cex") * clabel
	for (i in 1:(length(x))) {
		cha <- as.character(label[i])
		cha <- paste(" ", cha, " ", sep = "")
		x1 <- x[i]
		y1 <- y[i]
		xh <- strwidth(cha, cex = cex0)
		yh <- strheight(cha, cex = cex0) * 5/3
		if (!horizontal) {
			tmp <- scatterutil.convrot90(xh, yh)
			xh <- tmp[1]
			yh <- tmp[2]
		}
		if (boxes) {
			rect(x1 - xh/2, y1 - yh/2, x1 + xh/2, y1 + yh/2, 
					col = bg[i], border = "white")
		}
		if (horizontal) {
			text(x1, y1, cha, cex = cex0, col = "white")
		}
		else {
			text(x1, y1, cha, cex = cex0, col = "white", srt = 90)
		}
	}
}


fac2disj<- function(fac, drop = FALSE) {
	## Returns the disjunctive table corrseponding to a factor
	n <- length(fac)
	fac <- as.factor(fac)
	if(drop)
		fac <- factor(fac)
	x <- matrix(0, n, nlevels(fac))
	x[(1:n) + n * (unclass(fac) - 1)] <- 1
	dimnames(x) <- list(names(fac), as.character(levels(fac)))
	return(data.frame(x, check.names = FALSE))
}
#-------------------------------------------------------------------------------
s.class2 <- function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
		cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
		clabel = 1, cpoint = 1, pch = 20, col = rep(1, length(levels(fac))), 
		xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, origin = c(0, 
				0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
		cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE, bg="white") 
{
	head(dfxy)
	opar <- par(mar = par("mar"))
	par(mar = c(0.1, 0.1, 0.1, 0.1))
	on.exit(par(opar))
	dfxy <- data.frame(dfxy)
	if (!is.data.frame(dfxy)) 
		stop("Non convenient selection for dfxy")
	if (any(is.na(dfxy))) 
		stop("NA non implemented")
	if (!is.factor(fac)) 
		stop("factor expected for fac")
	dfdistri <- fac2disj(fac) * wt
	coul <- col
	w1 <- unlist(lapply(dfdistri, sum))
	dfdistri <- t(t(dfdistri)/w1)
	coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
	cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
	if (nrow(dfxy) != nrow(dfdistri)) 
		stop(paste("Non equal row numbers", nrow(dfxy), nrow(dfdistri)))
	coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
			xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
			cgrid = cgrid, include.origin = include.origin, origin = origin, 
			sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
			contour = contour, area = area, add.plot = add.plot)
	if (cpoint > 0) 
		for (i in 1:ncol(dfdistri)) {
			pch <- rep(pch, length = nrow(dfxy))
			points(coo$x[dfdistri[, i] > 0], coo$y[dfdistri[, 
									i] > 0], pch = pch[dfdistri[, i] > 0], cex = par("cex") * 
							cpoint, col = coul[i])
		}
	if (cstar > 0) 
		for (i in 1:ncol(dfdistri)) {
			scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar, 
					coul[i])
		}
	if (cellipse > 0) 
		for (i in 1:ncol(dfdistri)) {
			scatterutil.ellipse(coo$x, coo$y, dfdistri[, i], 
					cellipse = cellipse, axesell = axesell, coul[i])
		}
	if (clabel > 0) 
		scatterutil.eti2(coox, cooy, label, clabel, coul = "white", bg=bg)
	box()
	invisible(match.call())
}

## -----------------------------------------------------------------------------
## 2D ACP with ADE4
## -----------------------------------------------------------------------------


cat(" --> Performing PCA with ADE4.\n")

acp <- dudi.pca(t(d),  scannf=FALSE, nf=2)

if(ncol(acp$li) >= 2){
	# With Sample names
	png(file.path(outprefix, "PCA_dim_genes_2D_samples.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	
	head(acp$li)
	if(ncol(acp$li) > 2){
		s.class2(acp$li,
				as.factor(1:length(label)),
				cellipse=1.5,
				cstar=0,
				label=label,
				clabel=expansion,
				col="white",
				bg=mycolor)
	}

	devnull <- dev.off()
	
	
	centroids <- aggregate(cbind(x,y)~group,
						data.frame( x=acp$li[,1], 
									y=acp$li[,2], 
									group=mypheno),
						mean)
	centroids <- data.frame(centroids, mycolor=mycolor[as.character(centroids$group)])


# Ellipse Star and class name
	png(file.path(outprefix, "PCA_dim_genes_2D_classes.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	
	
		s.class(centroids[ ,c("x","y")],
				as.factor(centroids$group),
				cellipse=1,
				cstar=1,
				label=as.character(centroids$group),
				clabel=expansion,
				col=as.character(centroids$mycolor),
				cpoint=expansion)
	
	devnull <- dev.off()
	
	# Ellipse Star and class name and sample name
	png(file.path(outprefix, "PCA_dim_genes_2D_overlay.png"), width=pngwidth, height=pngheigth)



		s.class2(acp$li,
			as.factor(1:length(label)),
			cellipse=0,
			cstar=0,
			label=label,
			clabel=expansion,
			col="white",
			bg=mycolor)
	
	s.class(centroids[ ,c("x","y")],
			as.factor(centroids$group),
			cellipse=1,
			cstar=1,
			label=as.character(centroids$group),
			clabel=expansion,
			col=as.character(centroids$mycolor),
			cpoint=expansion,
			add.plot=TRUE)
	
	devnull <- dev.off()



	# Ebouli Val propre
	png(file.path(outprefix,"PCA_Eigen_val.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
		barplot(cumsum(acp$eig)/sum(acp$eig))
	devnull <- dev.off()


	## -----------------------------------------------------------------------------
	## Defining function for PCA with ggplot2
	## -----------------------------------------------------------------------------
	# from: http://stackoverflow.com/questions/7660893/boxed-geom-text-with-ggplot2/9814794#9814794
	
	
	btextGrob <- function (label,x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
			just = "centre", hjust = NULL, vjust = NULL, rot = 0, check.overlap = FALSE, 
			default.units = "npc", name = NULL, gp = gpar(), vp = NULL, expand_w, expand_h, box_gp = gpar()) {
		if (!is.unit(x)) 
			x <- unit(x, default.units)
		if (!is.unit(y)) 
			y <- unit(y, default.units)
		grob(label = label, x = x, y = y, just = just, hjust = hjust, 
				vjust = vjust, rot = rot, check.overlap = check.overlap, 
				name = name, gp = gp, vp = vp, cl = "text")
		tg <- textGrob(label = label, x = x, y = y, just = just, hjust = hjust, 
				vjust = vjust, rot = rot, check.overlap = check.overlap)
		w <- unit(rep(1, length(label)), "strwidth", as.list(label))
		h <- unit(rep(1, length(label)), "strheight", as.list(label))
		rg <- rectGrob(x=x, y=y, width=expand_w*w, height=expand_h*h,
				gp=box_gp)
		
		gTree(children=gList(rg, tg), vp=vp, gp=gp, name=name)
	}
	
	GeomTextbox <- proto(ggplot2:::GeomText, {
				objname <- "textbox"
				
				draw <- function(., data, scales, coordinates, ..., parse = FALSE, na.rm = FALSE,
						expand_w = 1.2, expand_h = 2, bgcol = "grey50", bgfill = "white", bgalpha = 1) {
					data <- remove_missing(data, na.rm, 
							c("x", "y", "label"), name = "geom_textbox")
					lab <- data$label
					if (parse) {
						lab <- parse(text = lab)
					}
					
					with(coord_transform(coordinates, data, scales),
							btextGrob(lab, x, y, default.units="native", 
									hjust=hjust, vjust=vjust, rot=angle, 
									gp = gpar(col = alpha(colour, alpha), fontsize = size * .pt,
											fontfamily = family, fontface = fontface, lineheight = lineheight),
									box_gp = gpar(fill = bgfill, alpha = bgalpha, col = bgcol),
									expand_w = expand_w, expand_h = expand_h)
					)
				}
				
			})
	
	geom_textbox <- function (mapping = NULL, data = NULL, stat = "identity", position = "identity", 
			parse = FALSE,  ...) { 
		GeomTextbox$new(mapping = mapping, data = data, stat = stat,position = position, 
				parse = parse, ...)
	}
	
	
	suppressMessages(load.fun(proto))
	suppressMessages(load.fun(ggplot2))
	
	StatEllipse <- proto(ggplot2:::Stat,
			{
				required_aes <- c("x", "y")
				default_geom <- function(.) GeomPath
				objname <- "ellipse"
				
				calculate_groups <- function(., data, scales, ...){
					.super$calculate_groups(., data, scales,...)
				}
				calculate <- function(., data, scales, level = 0.75, segments = 51,...){
					dfn <- 2
					dfd <- length(data$x) - 1
					if (dfd < 3){
						ellipse <- rbind(c(NA,NA))	
					} else {
						require(MASS)
						v <- cov.trob(cbind(data$x, data$y))
						shape <- v$cov
						center <- v$center
						radius <- sqrt(dfn * qf(level, dfn, dfd))
						angles <- (0:segments) * 2 * pi/segments
						unit.circle <- cbind(cos(angles), sin(angles))
						
						ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
					}
					
					ellipse <- as.data.frame(ellipse)
					colnames(ellipse) <- c("x","y")
					return(ellipse)
				}
			}
	)
	
	stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
		StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
	}

	## -----------------------------------------------------------------------------
	## PCA with ggplot2
	## -----------------------------------------------------------------------------
	
	
	colnames(acp$li) <- c("x","y")
	
	scores <- acp$li



	df <- data.frame(scores, row.names=label, group=as.factor(mypheno), colour=as.factor(mycolor))
	
	
	# calculating the ellipses by df$group
	# create an empty dataframe
	df_ell <- data.frame()
	# for each level in df$groups 
	for(g in unique(df$group)){
		# create 100 points per variable around the mean of each group
		df_ell <- rbind(df_ell, cbind(as.data.frame(with(df[df$group==g,], ellipse(cor(x, y),
												scale=c(sd(x),sd(y)),
												centre=c(mean(x),mean(y))
										)
								)
						),
						group=g))
	}
	
	
	tmp <- as.character(df_ell$group) 
	ucolor <- unique(as.character(mycolor))
	upheno <- unique(as.character(mypheno))
	for(i in 1:length(ucolor))
		tmp[tmp == upheno[i]] <- ucolor[i]
	df_ell$colour <- tmp
	
	cat(" --> Performing PCA with ggplot2.\n")
	cat(" -->  --> With Sample Names.\n")
	
	png(file.path(outprefix,"PCA_ggplot_sampleName.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)

		p <- ggplot(df, aes(x = x, y = y)) + geom_blank()
		p <- p + theme_bw()
		
		find_hull <- function(df) df[chull(df$x, df$y), ]
		for(i in ucolor){
			df_ell_sub <- df_ell[as.character(df_ell$colour)==i,]
			df_sub <- df[as.character(df$colour)==i,]
			cur_col <- as.character(unique(df_sub$colour))
			hulls <- ddply(df_sub, "colour", find_hull)
			p <- p + geom_polygon(data = hulls, aes(x=x, y=y), colour=NA, fill=cur_col, alpha=.5)
			p <- p + geom_polygon(data=df_ell_sub, aes(x=x, y=y), colour=cur_col, fill=cur_col, alpha=.1)
			p <- p + geom_textbox(data = df_sub, aes(x = x, y = y),label = rownames(df_sub), size=expansion, expand = 1, bgcol = cur_col, bgfill = cur_col, bgalpha = 0.9)
		}
		
		p <- p + geom_text(label=rownames(df),colour = "white", size=expansion)
		print(p)

	devnull <-  dev.off()


	cat(" -->  --> With Class Names.\n")
	
	png(file.path(outprefix,"PCA_ggplot_ClassName.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)

		p <- ggplot(df, aes(x = x, y = y)) + geom_blank()
		p <- p + theme_bw()
		
		find_hull <- function(df) df[chull(df$x, df$y), ]
		for(i in ucolor){
			df_ell_sub <- df_ell[as.character(df_ell$colour)==i,]
			df_sub <- df[as.character(df$colour)==i,]	
			cur_col <- as.character(unique(df_sub$colour))
			hulls <- ddply(df_sub, "colour", find_hull)
			p <- p + geom_polygon(data = hulls, aes(x=x, y=y), colour=NA, fill=cur_col, alpha=.5)
			p <- p + geom_polygon(data=df_ell_sub, aes(x=x, y=y), colour=cur_col, fill=cur_col, alpha=.1)
			
		}
		
		class_name <- centroids$group
		
		
		gg <- merge(df,aggregate(cbind(mean.x=x,mean.y=y)~group,df,mean),by="group")
	
		for(i in ucolor){
			gg_sub <- gg[as.character(gg$colour)==i,]
			cur_col <- as.character(unique(gg_sub$colour))
			p <- p + geom_point(data=gg_sub, aes(x,y), colour=cur_col)
			p <- p + geom_point(data=gg_sub, aes(x=mean.x,y=mean.y),size=0.1, colour=cur_col)
			p <- p + geom_segment(data=gg_sub, aes(x=mean.x, y=mean.y, xend=x, yend=y), colour=cur_col)
		}
	
	
		p <- p + geom_textbox(data = centroids, aes(x = x, y = y), 
				label = as.character(centroids$group), size = expansion, 
				expand = 1, bgcol = as.character(centroids$mycolor), 
				bgfill = as.character(centroids$mycolor), bgalpha = 0.9)
		
		p <- p + geom_text(data = centroids, aes(x = x, y = y), label = as.character(centroids$group), colour = "white", size=expansion)
	
		print(p)

	devnull <- dev.off()


	## -----------------------------------------------------------------------------
	## MDS 
	## -----------------------------------------------------------------------------
	
	cat(" --> Performing MDS with ggplot2.\n")
	
	suppressMessages(load.fun(MASS))
	
	dissimilarity <- 1 - cor(d)
	distance <- as.dist(dissimilarity)
	
	fit <- isoMDS(distance, k=2) # k is the number of dim
	
	x <- fit$points[,1]
	y <- fit$points[,2]
	
	
	df <- data.frame(x, y, row.names=label, 
			group=as.character(mypheno), 
			colour=as.character(mycolor))
	
	
	centroids <- aggregate(cbind(x,y)~group,
			data.frame( x=df$x, 
					y=df$y, 
					group=df$group),
			mean)
	centroids <- data.frame(centroids, mycolor=mycolor[as.character(centroids$group)])
		
	
	# calculating the ellipses by df$group
	# create an empty dataframe
	df_ell <- data.frame()
	# for each level in df$groups 
	for(g in unique(df$group)){
		# create 100 points per variable around the mean of each group
		df_ell <- rbind(df_ell, cbind(as.data.frame(with(df[df$group==g,], ellipse(cor(x, y),
												scale=c(sd(x),sd(y)),
												centre=c(mean(x),mean(y))
										)
								)
						),
						group=g))
	}
	
	
	png(file.path(outprefix,"MDS_ggplot_sampleName.png"), width=pngwidth, height=pngheigth)

	p <- ggplot(df, aes(x = x, y = y)) + geom_blank()
	p <- p + theme_bw()
	
	find_hull <- function(df) df[chull(df$x, df$y), ]
	for(i in ucolor){
		df_ell_sub <- df_ell[as.character(df_ell$colour)==i,]
		df_sub <- df[as.character(df$colour)==i,]
		cur_col <- as.character(unique(df_sub$colour))
		hulls <- ddply(df_sub, "colour", find_hull)
		p <- p + geom_polygon(data = hulls, aes(x=x, y=y), colour=NA, fill=cur_col, alpha=.5)
		p <- p + geom_polygon(data=df_ell_sub, aes(x=x, y=y), colour=cur_col, fill=cur_col, alpha=.1)
		p <- p + geom_textbox(data = df_sub, aes(x = x, y = y), 
				label = rownames(df_sub), size = expansion, expand = 1, bgcol = cur_col, bgfill = cur_col, bgalpha = 0.9)
	}
	
	p <- p + geom_text(label=rownames(df),colour = "white", size = expansion)
	print(p)	
	devnull <-  dev.off()



	png(file.path(outprefix,"MDS_ggplot_ClassName.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)

	p <- ggplot(df, aes(x = x, y = y)) + geom_blank()
	p <- p + theme_bw()
	
	find_hull <- function(df) df[chull(df$x, df$y), ]
	for(i in ucolor){
		df_ell_sub <- df_ell[as.character(df_ell$colour)==i,]
		df_sub <- df[as.character(df$colour)==i,]	
		cur_col <- as.character(unique(df_sub$colour))
		hulls <- ddply(df_sub, "colour", find_hull)
		p <- p + geom_polygon(data = hulls, aes(x=x, y=y), colour=NA, fill=cur_col, alpha=.5)
		p <- p + geom_polygon(data=df_ell_sub, aes(x=x, y=y), colour=cur_col, fill=cur_col, alpha=.1)
		
	}
	
	class_name <- centroids$group
	
	gg <- merge(df,aggregate(cbind(mean.x=x,mean.y=y)~group,df,mean),by="group")
	
	for(i in ucolor){
		gg_sub <- gg[as.character(gg$colour)==i,]
		cur_col <- as.character(unique(gg_sub$colour))
		p <- p + geom_point(data=gg_sub, aes(x,y), colour=cur_col)
		p <- p + geom_point(data=gg_sub, aes(x=mean.x,y=mean.y),size=0.1, colour=cur_col)
		p <- p + geom_segment(data=gg_sub, aes(x=mean.x, y=mean.y, xend=x, yend=y), colour=cur_col)
	}
	
	
	p <- p + geom_textbox(data = centroids, aes(x = x, y = y), 
			label = as.character(centroids$group), size = expansion, 
			expand = 1, bgcol = as.character(centroids$mycolor), 
			bgfill = as.character(centroids$mycolor), bgalpha = 0.9)
	
	p <- p + geom_text(data = centroids, aes(x = x, y = y), label = as.character(centroids$group), colour = "white", size=expansion)
	
	print(p)

	
	devnull <- dev.off()
}else{
	
	png(file.path(outprefix, "PCA_dim_genes_2D_samples.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()
	png(file.path(outprefix, "PCA_dim_genes_2D_classes.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()	
	png(file.path(outprefix, "PCA_dim_genes_2D_overlay.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()	
	png(file.path(outprefix,"PCA_Eigen_val.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()	
	png(file.path(outprefix,"PCA_ggplot_sampleName.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()	
	png(file.path(outprefix,"PCA_ggplot_ClassName.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()	
	png(file.path(outprefix,"MDS_ggplot_sampleName.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()
	png(file.path(outprefix,"MDS_ggplot_ClassName.png"), 
			units = 'in', res = 300, width=pngwidth, height=pngheigth)
	plot.new()
	text(0.5,0.5, "not enough sample")
	dev.off()	
}

#tsne marche très bien sur les gènes.
# ici on l'applique aux échantillons
#tsne_iris = tsne(t(d),  perplexity=50)