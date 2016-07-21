rule report:
    """
    Generate a report.
    """
    input:  code="output/code/Snakefile.py"

    params: wdir=config["workingdir"], \
            user=config["user"], \
            title = config["exp_title"],
            smp = config["samples"]

    output: html="output/report/report.html",
            rmd="output/report/report.html"
    run:
        import os
        import itertools
        import datetime
        today = datetime.date.today()

        ALL = ""
    
        #-----------------------------------------------------------------------
        # Create a header
        #-----------------------------------------------------------------------

        HEAD = """
---
title: {title}
author: {user}
date: '`r format(Sys.time(), "%d %B, %Y")`'
output: html_document
---

<style>
h3 {{
    background-color: gray !important;
    padding: 10px;
    color: white !important;
}}

td,tr,table,th {{  /* Table  */
   font-size: 12px;
}}

div.section h2 {{
    font-size: 12px !important;
    font-weight: bold;
}}

</style>

```{{r, echo=FALSE}}
# IMPORTANT NOTE: the code will be interpreted 
# relative to output/report

# Set wd
setwd('{path}')

# Source function
source("../../progs/scripts/rfunc.R")

# Libraries
library(knitr)
```
        """.format(title=params.title,
                   user=params.user,
                   path=os.path.join(params.wdir, 'output', 'report'))
    
        SAMPLES = params.smp.split()

        ALL += HEAD


        #-----------------------------------------------------------------------
        # fastq quality control
        #-----------------------------------------------------------------------


        img_type = ["adapter_content", 
                    "duplication_levels", 
                    "kmer_profiles", 
                    "per_base_n_content", 
                    "per_base_quality", 
                    "per_base_sequence_content",
                    "per_sequence_gc_content", 
                    "per_sequence_quality", 
                    "per_tile_quality",
                    "sequence_length_distribution"]
    
        QUALITY = """


### Quality control ({read}, {smp})

Quality control of sample: **{smp}**. Statistics were performed using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
More information on diagram interpretation can be obtained [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

```{{r, echo=FALSE, results='asis'}}
table_list <- list()
nb_table <- 1
img <- Sys.glob("../../output/fastqc_raw/{smp}/*{read}.fq*/Images/*.png")
img_html <- vector2_html_img(img, pos=3, width=300)

while(length(img_html) >= 3){{
    tmp_html <- img_html[1:3]; 
    img_html <- img_html[-c(1,2,3)]
    tmp_mat <- matrix(tmp_html, ncol=3, byrow=T)
    colnames(tmp_mat) <- gsub(".png.*", "", unlist(lapply(strsplit(tmp_html, "/"), "[", 8)))
    table_list [[nb_table]] <- tmp_mat
    nb_table <- nb_table + 1
}}

tmp_mat <- matrix(img_html, ncol=length(img_html), byrow=T)
colnames(tmp_mat) <- gsub(".png.*", "", unlist(lapply(strsplit(img_html, "/"), "[", 8)))
table_list [[nb_table]] <- tmp_mat

for(i in 1:length(table_list)){{
    print(kable(table_list[[i]]))
}}



```

--------------------------------------------------------------------------------

        """

        for i in list(itertools.product(["R1","R2"], SAMPLES)):
            ALL +=QUALITY.format(read=i[0], 
                                 smp=i[1])
    

        #-----------------------------------------------------------------------
        # Mapping stats
        #-----------------------------------------------------------------------

        MAPPING_STATS = """


### Mapping statistics control ({smp})

Mapping statistics for sample **{smp}**. Performed with R.

```{{r, echo=FALSE, results='asis'}}
img <- Sys.glob("../../output/mapping_stats/{smp}.stats.png")
img_html <- vector2_html_img(img, pos=5, width=300)

cat(img_html)

d <- read.table("../../output/mapping_stats/{smp}.stats.txt", head=T,row=1)
rownames(d) <- c("Total reads", "Trimmed reads", "mapped reads")
kable(d, align="l", digits = 2, padding=0)

```

        """

        for i in SAMPLES:
            ALL +=MAPPING_STATS.format(smp=i)  

        #-----------------------------------------------------------------------
        # Various correlation spots
        #-----------------------------------------------------------------------

        CORRELATION = """


### Various correlation plots

```{r, echo=FALSE, results='asis'}

table_list <- list()
nb_table <- 1
img <- Sys.glob("../../output/corr_plot/*png")
img_html <- vector2_html_img(img, pos=3, width=300)

while(length(img_html) >= 3){{
    tmp_html <- img_html[1:3]; 
    img_html <- img_html[-c(1,2,3)]
    tmp_mat <- matrix(tmp_html, ncol=3, byrow=T)
    colnames(tmp_mat) <- rep("Correlation diagram", 3)
    table_list [[nb_table]] <- tmp_mat
    nb_table <- nb_table + 1
}}

tmp_mat <- matrix(img_html, ncol=length(img_html), byrow=T)
colnames(tmp_mat) <- rep("Correlation diagram", length(img_html))
table_list [[nb_table]] <- tmp_mat

for(i in 1:length(table_list)){{
    print(kable(table_list[[i]]))
}}


img <- Sys.glob("../../output/corr_plot/Pairs_plot.png")
img_html <- vector2_html_img(img, pos=3, width=900)
cat(img_html)



```

        """


        ALL +=CORRELATION
   
        #-----------------------------------------------------------------------
        # Write
        #-----------------------------------------------------------------------

        fh = open("output/report/report.Rmd", "w")
        fh.write(ALL)
        fh.close()
        fh2 = open("report.Rmd.bla", "w")
        fh2.write(ALL)
        fh2.close()


        #-----------------------------------------------------------------------
        # Process Rmd to html
        #-----------------------------------------------------------------------
        
        CMD = """ cd output/report; Rscript -e "rmarkdown::render('report.Rmd')" """
        
        os.system(CMD)
        