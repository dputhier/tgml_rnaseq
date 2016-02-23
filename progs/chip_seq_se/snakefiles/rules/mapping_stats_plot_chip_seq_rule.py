rule hist_mapping_stats:
    input: "output/mapping_stats/{smp}.stats"
    output: "output/mapping_stats/{smp}.stats.png"
    threads: 1
    run: R("""
      col.pal <- c("#7fc97f","#beaed4", "#fdc086","#ffff99","#386cb0")
      f <- Sys.glob("*/mapping_stats/{wildcards.smp}*.stats")
      m <- matrix(NA, nc=length(f), nr=5)
      colnames(m) <- basename(f)
      colnames(m) <- gsub(".stats","",colnames(m))
      rownames(m) <- c("TOTAL", "TRIM", "MAPPED", "QUALITY_FILTERED", "RM_DUP")

      for(i in 1:length(f)){{
         d <- read.table(f[i], fill=NA)
         m[1,i] <- d[2,"V3"]
         m[2,i] <- d[4,"V3"]
         m[3,i] <- as.numeric(as.character(d[6,1]))
         m[4,i] <- as.numeric(as.character(d[8,1]))
         m[5,i] <- as.numeric(as.character(d[10,1]))
         
      }}
     png("{output[0]}", width = 350, hei=350)
     barplot(m, beside=T, col=col.pal, cex.names=0.7 , las=2)
     abline(h=seq(10e6,1e9, 10e6), col="gray", lty=2)
     legend("bottomleft", leg=c("Raw", "Trimmed", "Mapped", "Qual_filt", "Rm_dup"), fill=col.pal)
     dev.off()

          
    """)