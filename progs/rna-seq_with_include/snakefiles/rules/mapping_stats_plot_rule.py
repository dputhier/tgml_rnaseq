rule hist_mapping_stats:
    input: "output/mapping_stats/{smp}_R1.stats", "output/mapping_stats/{smp}_R2.stats"
    output: "output/mapping_stats/{smp}.stats.png", "output/mapping_stats/{smp}.stats.pdf"
    threads: 1
    run: R("""
      f <- Sys.glob("*/mapping_stats/{wildcards.smp}_R*.stats")
      m <- matrix(NA, nc=length(f), nr=3)
      colnames(m) <- basename(f)
      colnames(m) <- gsub(".stats","",colnames(m))
      rownames(m) <- c("TOTAL", "TRIM", "MAPPED")

      for(i in 1:length(f)){{
         d <- read.table(f[i], fill=NA)
         m[1,i] <- d[2,"V3"]
         m[2,i] <- d[4,"V3"]
         m[3,i] <- as.numeric(as.character(d[6,1]))
      }}
     png("{output[0]}", width = 350, hei=350)
     barplot(m, beside=T, col=c("#FDD369", "#F98D32", "#1BBBDA"), cex.names=0.7 , las=2)
     abline(h=seq(10e6,1e9, 10e6), col="gray", lty=2)
     legend("bottomleft", leg=c("Raw", "Trimmed", "Mapped"), fill=c("#FDD369", "#F98D32", "#1BBBDA"))
     dev.off()

     pdf("{output[1]}")
     barplot(m, beside=T, col=c("#FDD369", "#F98D32", "#1BBBDA"), border="white", cex.names=0.7 , las=2)
     abline(h=seq(10e6,1e9, 10e6), col="gray", lty=2)
     legend("bottomleft", leg=c("Raw", "Trimmed", "Mapped"), fill=c("#FDD369", "#F98D32", "#1BBBDA"))
     dev.off()
          
    """)