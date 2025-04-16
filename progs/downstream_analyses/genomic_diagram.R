library(Gviz)
library(rtracklayer)
chr <- "10"
gm <- GeneRegionTrack("/shared/projects/project_dputhier/043_irla_elife_2022/input/gtf/Mus_musculus.GRCm38.97.gtf", chromosome=chr)
gm <- gm[gene(gm) == "ENSMUSG00000000731" & feature(gm) %in% c("exon", "CDS", "UTR")]
plotTracks(gm, thinBoxFeature="UTR")

itrack <- IdeogramTrack(genome = "mm10", chromosome = "10")
axis_track <- GenomeAxisTrack()
plotTracks(list(gm, itrack, axis_track), thinBoxFeature="UTR", collapse=FALSE, transcriptAnnotation="Aire")


 
allChromosomeCoverage <- import.bw("/shared/projects/project_dputhier/043_irla_elife_2022/output/bwig/K14~mTEClo~int~SRR11000298_plus.bw",
                                   as="GRanges")
accDT <- DataTrack(allChromosomeCoverage,chomosome="10") 
plotTracks(list(gm, itrack, axis_track, accDT), thinBoxFeature="UTR", collapse=FALSE, transcriptAnnotation="Aire")
