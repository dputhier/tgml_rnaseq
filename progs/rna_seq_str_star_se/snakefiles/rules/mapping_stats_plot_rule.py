rule barplot_stats:
    input: "output/mapping_stats/{smp}.stats"
    output: "output/mapping_stats/{smp}.stats.png", "output/mapping_stats/{smp}.stats.pdf"
    params: mem="4G"
    threads: 1
    shell: """
     module load r/4.4.1
     echo " f <- Sys.glob('output/mapping_stats/{wildcards.smp}.stats');
      m <- matrix(NA, nc=length(f), nr=3);
      colnames(m) <- basename(f);
      colnames(m) <- gsub('.stats','',colnames(m));
      rownames(m) <- c('TOTAL', 'TRIM', 'MAPPED');
      for(i in 1:length(f)){{d <- read.table(f[i], fill=NA);
         m[1,i] <- as.numeric(as.character(d[2,'V3']));
         m[2,i] <- as.numeric(as.character(d[4,'V3']));
         m[3,i] <- as.numeric(as.character(d[6,1])) }};
     png('{output[0]}', width=6, height=6, units = 'in', res = 300);
     barplot(m, beside=T, col=c('#FDD369', '#F98D32', '#1BBBDA'), cex.names=0.7 , las=2);
     abline(h=seq(10e6,1e9, 10e6), col='gray', lty=2);
     legend('bottomleft', leg=c('Raw', 'Trimmed', 'Mapped'), fill=c('#FDD369', '#F98D32', '#1BBBDA'));
     dev.off();
     pdf('{output[1]}');
     barplot(m, beside=T, col=c('#FDD369', '#F98D32', '#1BBBDA'), border='white', cex.names=0.7 , las=2);
     abline(h=seq(10e6,1e9, 10e6), col='gray', lty=2);
     legend('bottomleft', leg=c('Raw', 'Trimmed', 'Mapped'), fill=c('#FDD369', '#F98D32', '#1BBBDA'));
     dev.off();
    write.table(m, 'output/mapping_stats/{wildcards.smp}.stats.txt', sep='\\t', quote=F, col.names=NA);   
    " | R --slave
    """