
# Create a markdown links ![name](link) from a path
# pos  define the place 'name' when the path is splitted using
# "/'
vector2_md_link <- function(x, chunknb=3){
	v <- c()
	for(i in x){
		nm <- unlist(strsplit(i, "/"))[chunknb]
		v <- append(v, paste("![", nm ,"](", i, ")", sep=""))
	}
	return(v)
}

vector2_html_img <- function(x, pos=5, width=100){
	v <- c()
	for(i in x){
		nm <- unlist(strsplit(i, "/"))[pos]
		v <- append(v, paste("<img src='", 
							i, "' width=",
							width, " alt='",
							nm, "'>", sep=""))
	}
	return(v)
}