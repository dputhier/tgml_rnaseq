library(ggplot2)
library(ggthemes)
library(ggrepel)
library(knitr)

# Create a markdown links ![name](link) from a path
# pos  define the place 'name' when the path is splitted using
# "/'
vector2_md_link <- function(x, chunknb=3, insert=T){
	v <- c()
	for(i in x){
		nm <- unlist(strsplit(i, "/"))[chunknb]
		if(insert){
			v <- append(v, paste("![", nm ,"](", i, ")", sep=""))
		}else{
			v <- append(v, paste("[", nm ,"](", i, ")", sep=""))
		}
	}
	return(v)
}

vector2_html_img <- function(x, pos=5, width=100, nm=NULL){
	v <- c()
	for(i in x){
		if(is.null(nm))
			nm <- unlist(strsplit(i, "/"))[pos]
		
		v <- append(v, paste("<img src='", 
							i, "' width=",
							width, " alt='",
							nm, "'>", sep=""))
	}
	return(v)
}

find_img_and_dotable <- function(glob=NULL, width=300, ncol=3, pos=3, title=NULL){
	table_list <- list()
	nb_table <- 1
	img <- Sys.glob(glob)
	img_html <- vector2_html_img(img, pos=pos, width=width)
	
	while(length(img_html) >= ncol){{
			tmp_html <- img_html[1:ncol]; 
			img_html <- img_html[-c(1:ncol)]
			tmp_mat <- matrix(tmp_html, ncol=ncol, byrow=T)
			colnames(tmp_mat) <- rep(title, ncol)
			table_list [[nb_table]] <- tmp_mat
			nb_table <- nb_table + 1
		}}
	
	tmp_mat <- matrix(img_html, ncol=length(img_html), byrow=T)
	colnames(tmp_mat) <- rep(title, ncol(tmp_mat))
	table_list [[nb_table]] <- tmp_mat
	
	for(i in 1:length(table_list)){
			print(kable(table_list[[i]]))
		}
}

maplot <- function(R,G, title="bla", M_pts_lim=2, M_txt_lim=3, gene_names=NULL, cex=5){

	M <- G - R
	A <- G + R
	M_pts_lim <- abs(M) > M_pts_lim
	M_txt_lim <- abs(M) > M_txt_lim

	d <- data.frame(M=M, A=A, 
					Mptslim=factor(M_pts_lim), 
					M_txt_lim=factor(M_txt_lim),
					gene_names=gene_names)
	plot <- ggplot(data=d, aes(A, M, color=Mptslim))
	plot <- plot + geom_point(size=0.8, alpha=0.5,na.rm=T) 
	plot <- plot + labs(title=title)

	if(any(M_txt_lim)){
		plot <- plot + geom_text_repel(data=d[M_txt_lim,], 
										aes(label=gene_names),
										color = 'gray25',
										size=cex)
	}
	plot <- plot + theme_minimal()
	plot <- plot + theme_bw()
	plot <- plot + geom_hline(aes(yintercept = 0), color="gray")
	plot <- plot + geom_hline(aes(yintercept = 1), color="gray")
	plot <- plot + geom_hline(aes(yintercept = -1), color="gray")
	plot <- plot + stat_smooth(se = FALSE, color="violet")
	plot <- plot + theme(legend.position="none")
	return(plot)
	
}

maplot_pval <- function(R,G, title="bla", pval=NULL, thresh=0.05, gene_names=NULL, cex=5){
	
	M <- G - R
	A <- G + R

	
	pval <- ifelse(pval < thresh, TRUE,FALSE)
	pval.save <- pval

	d <- data.frame(M=M, A=A, 
			pval=factor(pval, levels=c("TRUE", "FALSE")), 
			gene_names=gene_names)
	plot <- ggplot(data=d, aes(A, M, color=pval))
	plot <- plot + geom_point(size=0.8, alpha=0.5,na.rm=T) 
	plot <- plot + labs(title=title)
	

	plot <- plot + geom_text_repel(data=d[pval,], 
				                   aes(label=gene_names[pval.save]),
								   color = 'gray25',
								   size=cex)

	plot <- plot + theme_minimal()
	plot <- plot + theme_bw()
	plot <- plot + geom_hline(aes(yintercept = 0), color="gray")
	plot <- plot + geom_hline(aes(yintercept = 1), color="gray")
	plot <- plot + geom_hline(aes(yintercept = -1), color="gray")
	plot <- plot + stat_smooth(se = FALSE, color="violet")
	plot <- plot + theme(legend.position="none")
	return(plot)
	
}

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
