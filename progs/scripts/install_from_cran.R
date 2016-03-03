load.fun <- function(x) {
		x <- as.character(substitute(x)) 
		if(isTRUE(x %in% .packages(all.available=TRUE))) {
				eval(parse(text=paste("require(", x, ")", sep=""))) 
			} else { 
				eval(parse(text=paste("install.packages('", x, "', 
												repos = 'http://cran.us.r-project.org')", sep=""))) 
				eval(parse(text=paste("require(", x, ")", sep=""))) 
			} 
	} 