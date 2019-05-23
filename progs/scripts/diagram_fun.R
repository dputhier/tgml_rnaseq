library("RColorBrewer")


#' barplot from a dataframe (single column)
#'
#' @param d a dataframe
#' @param column the column number
#' @return a plot
#' @examples
#' 


df_barplot <- function(d, 
					column=1, 
					col=NULL, 
					beside=FALSE){

	# d is a dataframe
	x <- as.matrix(d[,column])
	
	# color
	if(is.null(col)){
		col <- brewer.pal(length(x), "Dark2")
	}

	barplot(x, 
			col="white", 
			border="white", 
			beside=beside,
			names.arg=rownames(d), 
			las=1)
	grid(col="darkgrey")
	barplot(x, col=col, border=col, 
			beside=beside,
			add=TRUE, 
			las=1)

}



