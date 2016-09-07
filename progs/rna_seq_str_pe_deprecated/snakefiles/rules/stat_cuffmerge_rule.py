rule stat_cuffmerge:
    input: "output/cuffmerge/cuffmerge_stats.txt"
    output: "output/cuffmerge/cuffmerge_stats.png"
    threads: 1
    message: "Plotting stats for cuffmerge"
    run: R("""
        source("progs/scripts/diagram_fun.R")
        d <- read.table("{input}", sep="\t", header=F, row=2)
        png("{output}", width = 6, height = 6, units = 'in', res = 300)
        df_barplot(d, beside=TRUE)
        dev.off()
""")
