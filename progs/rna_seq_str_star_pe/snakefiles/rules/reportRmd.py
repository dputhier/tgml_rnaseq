import os

rule report:
    """
    Generate a report.
    """

    input:  code="output/code/Snakefile.py"
    params: wdir=config["workingdir"], \
            user=config["user"], \
            title = config["exp_title"],
            smp = config["samples"],
            comp = config["comparison"],
            house = config["housekeeping"],
            class_code = config["cuffmerge"]["selected_class_code"],
            classes = config["classes"],
            threshold_maplot = config["report"]["threshold_maplot"],
            chrom=config["chrom_list"]

    output: html="output/report/report.html",
            rmd="output/report/report.Rmd"
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
output: 
    html_document:
        self_contained: no

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
suppressMessages(library(ade4))
suppressMessages(library(adegraphics))
suppressMessages(library(lattice))

```
        """.format(title=params.title,
                   user=params.user,
                   path=os.path.join(params.wdir, 'output', 'report'))
    
        SAMPLES = params.smp.split()

        ALL += HEAD


        #-----------------------------------------------------------------------
        # Graph (global)
        #-----------------------------------------------------------------------
        
        WF_1= """
### The template workflow

```{r , echo=FALSE, results='asis'}

find_img_and_dotable(glob="../../output/report/rulegraph.png", 
                    title="Worklow (global level)",
                    width=800,
                    ncol=1)
```         
        """
        
        ALL += WF_1

        WF_2= """
###  The complete workflow (sample-wise)

```{r , echo=FALSE, results='asis'}

find_img_and_dotable(glob="../../output/report/dag.png", 
                    title="Worklow (global level)",
                    width=800,
                    ncol=1)
```         
        """
        
        ALL += WF_2


        #-----------------------------------------------------------------------
        # Info
        #-----------------------------------------------------------------------

        INFO= """
### Selected chromosome for analysis

To avoid mapping biases, reads were mapped to all chromosomes (this may included haplotype variants unassigned sequence...). 
The subsequent analysis (after mapping) will concentrate on chromosomes:

- {chrom}
        """.format(chrom=params.chrom.replace(","," "))
        
        ALL += INFO


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


### Quality control on raw reads ({read}, {smp})

Quality control of sample: **{smp}**. Statistics were performed using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
More information on diagram interpretation can be obtained [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

```{{r, echo=FALSE, results='asis'}}
table_list <- list()
nb_table <- 1
img <- Sys.glob("../../output/fastqc_raw/{smp}/*/Images/*.png")
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
        # Mapping stats by chrom
        #-----------------------------------------------------------------------

        MAPPING_STATS_BY_CHR = """


### Mapping statistics by chromosome ({smp})

Mapping statistics for sample **{smp}**. Statistics are provided by chromosome.

```{{r, echo=FALSE, results='asis'}}
img <- Sys.glob("../../output/bam_stat_by_chrom/{smp}_bam_stats.png")
img_html <- vector2_html_img(img, pos=5, width=300)

cat(img_html)

d <- read.table("../../output/bam_stat_by_chrom/{smp}_bam_stats.txt", head=F,row=1)
colnames(d) <- c("Nb read mapped")
kable(d, align="l", digits = 2, padding=0)

```

        """

        for i in SAMPLES:
            ALL += MAPPING_STATS_BY_CHR.format(smp=i)  


        #-----------------------------------------------------------------------
        # Check transcript coverage
        #-----------------------------------------------------------------------

        COVERAGE = """


### 5' to 3' coverage plot

This diagram displays gene body coverage obtained from representative genes.
Plot was computed using [geneBody_coverage.py](http://rseqc.sourceforge.net/). See [this file]({house}) to get the gene list used.
A pdf version of the diagram can be found [here](../../output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.pdf).


```{{r, engine='bash', echo=FALSE}}
convert ../../output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.pdf ../../output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.png
``` 

```{{r, echo=FALSE, results='asis'}}
img <- "../../output/rseqc_coverage_diag/rseqc_cov.geneBodyCoverage.curves.png"
img_html <- vector2_html_img(img, pos=5, width=300)

cat(img_html)

```

        """.format(house=params.house)

        ALL += COVERAGE

        #-----------------------------------------------------------------------
        # Novel transcripts
        #-----------------------------------------------------------------------

        CUFFMERGE = """

### Transcript discovery : statistics

Cufflinks was used to search for novel transcripts. These transcripts were merged with known transcripts to create an exhaustive gene annotation track (gtf format) that can be found [here](../../output/inferred_gene_annotation/known_transcripts_and_selected_class_code_u.gtf). 
The diagram below displays the number of transcripts found by class-code. The class-code meaning is available [here](http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/).
In the current analysis the transcript with class code **{class_code}** have been selected and used for quantification at the gene level. This may be changed in the config file.  

|Priority  |  Code  |  Description|
|----------|--------|-------------|
|1|    =|    Complete match of intron chain|
|2|    c|    Contained|
|3|    j|    Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript|
|4|    e|    Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment.|
|5|    i|    A transfrag falling entirely within a reference intron|
|6|    o|    Generic exonic overlap with a reference transcript|
|7|    p|    Possible polymerase run-on fragment (within 2Kbases of a reference transcript)|
|8|    r|    Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case|
|9|    u|    Unknown, intergenic transcript|
|10|    x|    Exonic overlap with reference on the opposite strand|
|11|    s|    An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)|
|12|    .|    (.tracking file only, indicates multiple classifications)|

```{{r, echo=FALSE, results='asis'}}
img <- "../../output/cuffmerge/cuffmerge_stats.png"
img_html <- vector2_html_img(img, pos=5, width=300)

cat(img_html)

```
        """.format(class_code=params.class_code)

        ALL += CUFFMERGE

        #-----------------------------------------------------------------------
        # Various correlation spots
        #-----------------------------------------------------------------------

        CORRELATION = """



### Various correlation plots

```{r, echo=FALSE, results='asis'}

table_list <- list()
nb_table <- 1
img <- Sys.glob("../../output/corr_plot/*png")

find_img_and_dotable(glob="../../output/corr_plot/*png", 
                    title="Correlation diagram",
                    ncol=3)

img <- Sys.glob("../../output/corr_plot/Pairs_plot.png")
img_html <- vector2_html_img(img, pos=3, width=900)
cat(img_html)



```

        """


        ALL += CORRELATION
 
        #-----------------------------------------------------------------------
        # MA PLOT
        #-----------------------------------------------------------------------

        MAPLOT = """


### MA plot by sample

These plots are intented to highlight genes that may be specific of each samples. For two samples A and B, the $M$ and $A$ values are computed as follow for each gene $g$.

$$M_{g} = log(A_{g}) - log(B_{g})$$
$$A_{g} = (log(A_{g}) + log(B_{g}))/2$$

When more than two samples are available, each sample (e.g $A$) is compared to a 
representative pseudo-sample that is computed from the row (gene) medians.

**Hint**: Warn with the list of differential expressed genes when no replicates are available. When working with polyA selected RNA, a classical bias is to observed histone genes that are known to be non-polyadenylated and that underline selection issues.
**Note**: The diagram shows normalized expression data (scaling factor as proposed by DESeq2).
 
```{r maplot, echo=FALSE, results='asis'}

img_list <- list()
dir.create("../maplot", showWarnings = FALSE)

d <- read.table("../../output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt",
                sep="\t",head=T,row=1)

sf <- estimSf(d)
d <- round(sweep(d, 2, sf, "/"),0)

d <- log2(d+1)

if(ncol(d) <= 2){

    p <- maplot(d[,1], d[,2], 
                gene_names=rownames(d), 
        title=paste(colnames(d)[1], "vs" , colnames(d)[2]))
    print(p)
    dev.off()
}else{

    B <- apply(d,1,median)
    for(i in 1:ncol(d)){
        p <- maplot(d[,i], B, 
                    gene_names=rownames(d), 
                    title=paste(colnames(d)[i], "vs reference (median gene expression)"))
        nm <- paste(colnames(d)[i], "vs_reference", ".png", sep="")
        path <- file.path("../maplot/",nm)
        png(path, width = 1500, height = 1200, res=250)
        print(p)
        dev.off()

    }
}

find_img_and_dotable(glob="../maplot/*.png",
                    title="MA plot",
                    ncol=2)

```

        """


        ALL += MAPLOT 

        #-----------------------------------------------------------------------
        # PCA
        #-----------------------------------------------------------------------

        PCA = """
        
### PCA diagrams

**From Wikipedia**: "Principal component analysis (PCA) is a statistical procedure that uses an orthogonal transformation to convert a set of observations of possibly correlated variables into a set of values of linearly uncorrelated variables called principal components."

The objective is to use a this newly discovered coordinate system to represents objects through a limited number of dimension while conserving most of the intrinsic variance. Most of the time one consider that dimension choosen for representation should capture 90% of the observed variance.

The following diagrams depict:

- The samples displayed in the new coordinate system made of the two principal components. Class boundaries are depicted by ellipses.
- The original variable (genes) that show the best correlations with the two principal components.
- The cumulative sums of eigen values (expressed as a percentage of variance)
- The samples (with label) displayed in the new coordinate system made of the two principal components.
    
```{{r pca, results='asis', echo=FALSE}}


d <- read.table("../../output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt",
                sep="\t",head=T,row=1)

sf <- estimSf(d)
d <- round(sweep(d, 2, sf, "/"),0)

d <- log2(d+1)

pheno <- "{pheno}"

pheno <- strsplit(pheno, " ")[[1]]


dir.create("../pca_mds", showWarnings = FALSE)

acp <- dudi.pca(t(d),  scannf=FALSE, nf=3)



if(ncol(acp$li) >= 2){{

    png("../../output/pca_mds/PCA_ade_g.png", 
            units = 'in', res = 300, width=6, height=6)
    s.class(acp$li, fac = as.factor(pheno), 
            pellipses.lwd = 2, pellipses.border = 2:5, 
            pellipses.col = 2:5)
    dev.off()

    png("../../output/pca_mds/PCA_ade_g_label.png", 
            units = 'in', res = 300, width=6, height=6)
    s.label(acp$li, ppoints.col = "red", plabels = list(box = list(draw = FALSE), optim = TRUE), plot = T)
    dev.off()

    png("../../output/pca_mds/PCA_ade_g_eigen.png", 
            units = 'in', res = 300, width=6, height=6)
        (kip <- 100 * acp$eig/sum(acp$eig))
        barplot(cumsum(kip), xlab="PC", ylab="Cumulative sum (%)")
    dev.off()

    png("../../output/pca_mds/PCA_ade_g_corr_circle.png", 
            units = 'in', res = 300, width=6, height=6)
    i <- rownames(head(acp$co[order(acp$co$Comp1),],20))
    j <- rownames(tail(acp$co[order(acp$co$Comp1),],20))
    k <- rownames(head(acp$co[order(acp$co$Comp2),],20))
    l <- rownames(tail(acp$co[order(acp$co$Comp2),],20))
    adegpar(list(plabels.cex=0.5))
    s.corcircle(acp$co[unique(c(i,j,k,l)),])
   dev.off()
}}

```

```{{r pca_insert, results='asis', echo=FALSE}}

find_img_and_dotable(glob="../../output/pca_mds/PCA_ade_g*png", width=250, ncol=2, pos=3, title="PCA with ade4")


```

        """.format(pheno = params.classes)

        ALL += PCA

        #-----------------------------------------------------------------------
        # MDS PLOT
        #-----------------------------------------------------------------------

        MDS = """
### Multidimensionnal scaling

```{{r mds, echo=FALSE, results='asis'}}

suppressMessages(library(MASS))
suppressMessages(library(ggplot2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))

d <- read.table("../../output/quantification_known_and_novel_genes/gene_counts_known_and_novel_mini.txt",
                sep="\t",head=T,row=1)

sf <- estimSf(d)
d <- round(sweep(d, 2, sf, "/"),0)

d <- log2(d+1)
pheno <- "{pheno}"

pheno <- strsplit(pheno, ",")[[1]]



dissimilarity <- 1 - cor(d)
distance <- as.dist(dissimilarity)
    
fit <- suppressWarnings(isoMDS(distance, k=2)) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]


df <- data.frame(x, y, row.names=names(x), 
        group=as.factor(pheno),
        colour=as.character(pheno),
        name=names(x))
p <- ggplot(df, aes(x = x, y = y)) + geom_blank()
p <- p + theme_bw()
p <- p + geom_point(size=0.8, alpha=0.5,na.rm=T) 
p <- p + geom_text_repel(data=df, aes(label=name), color = 'gray25')
png("../../output/pca_mds/MDS.png", width = 1500, height = 1200, res=250)
print(p)
dev.off()

find_img_and_dotable(glob="../../output/pca_mds/MDS.png", 
                width=250, ncol=2, pos=3, title="MDS")

```
        """.format(pheno = pheno)

        ALL += MDS


        #-----------------------------------------------------------------------
        # MDS PLOT
        #-----------------------------------------------------------------------

        MAPLOTBYCLASS = """
### Differential analysis ({comp})

- **Tested samples**: {sample}
    - **class 1 contains**: {cls1}
    - **class 2 contains**: {cls2}

- **Corrected p-value threshold (EdgeR_BH)**
    - {thresh}


```{{r , echo=FALSE, results='asis'}}


img_list <- list()

d <- read.table("../../output/comparison/{comp}/DESeq2_EdgeR_pval_and_norm_count_log2.txt",
                sep="\t",head=T,row=1,
                check.names = FALSE)
d.save <- d

class1 <- unlist(strsplit("{cls1}", " "))
class2 <- unlist(strsplit("{cls2}", " "))

d <- d[,c(class1, class2)]

cat(colnames(d), file="toto")
p <- maplot_pval(rowMeans(d[,class1]),
                 rowMeans(d[,class2]), 
                 gene_names=rownames(d),
                 thresh={thresh},
                 pval=d.save$EdgeR_BH, 
                 title="{comp}",
                 cex=1.5)

path <- file.path("../../output/comparison/{comp}/maplot_{comp}_report.png")
png(path, width = 1500, height = 1200, res=250)
suppressWarnings(print(p))
dev.off()

# clustering (heatmap of samples correlation)


d.clust <- na.omit(d[d.save$EdgeR_BH <= {thresh}, ])
if(nrow(d.clust) == 0){{
        print("No significant genes were found.")
}}else{{
        pear <- cor(d.clust, method="pearson")
        palette <-colorRampPalette(c("yellow", "black","blueviolet"))
        path <- file.path("../../output/comparison/{comp}/cor_heatmap_{comp}_report.png")

        png(path, width = 1500, height = 1200, res=250)
        levelplot(pear,col.regions=palette, scales=list(cex=0.4))
        dev.off()
        
        # clustering (tree)

        pear <- as.dist((1-pear)/2)
        hp <- hclust(pear, method="average")
        path <- file.path("../../output/comparison/{comp}/hclust_{comp}_report.png")
        png(path, width = 1500, height = 1200, res=250)
        plot(hp,hang=-1, lab=colnames(d.clust), cex=0.4)
        dev.off()

}}

## Summary
a <- length(d.save$EdgeR_BH[d.save$EdgeR_BH < {thresh}])
b <- length(d.save$DESeq2_padj[d.save$DESeq2_padj < {thresh} & !is.na(d.save$DESeq2_padj)])

nbsig <- data.frame(NbSig_DESeq2_padj=b,NbSig_EdgeR_BH=a)
kable(nbsig)

```


```{{r, results='asis', echo=FALSE}}
find_img_and_dotable(glob="../../output/comparison/{comp}/*_report.png", 
                    title="Comparison plot",
                    ncol=3)
```

- Two types of differential analysis are provided (EdgeR and DESeq2). Data (log2(counts + 1)) together with pvalues, adjusted-pvalues (...) are available here:

```{{r , echo=FALSE, results='asis'}}
vector2_md_link("../../output/comparison/{comp}/DESeq2_EdgeR_pval_and_norm_count_log2.txt",
                chunknb=5, 
                insert=F)
```
        """

        for i in params.comp:
            cls = config["comparison"][i]
            cls1,cls2 = config["comparison"][i].split(" ")
            ALL += MAPLOTBYCLASS.format(comp=i,
                                        sample=config["comparison"][i].replace(","," "),
                                        cls1=cls1.replace(","," "),
                                        cls2=cls2.replace(","," "),
                                        thresh=params.threshold_maplot)          
      
        #-----------------------------------------------------------------------
        # Write
        #-----------------------------------------------------------------------

        fh = open("output/report/report.Rmd", "w")
        fh.write(ALL)
        fh.close()



        #-----------------------------------------------------------------------
        # Process Rmd to html
        #-----------------------------------------------------------------------
        
        CMD = """ cd output/report; Rscript -e "rmarkdown::render('report.Rmd')" """
        
        os.system(CMD)
        
