rule prep_stat_cuffmerge:
    input: "output/cuffmerge/merged.gtf"
    output: "output/cuffmerge/cuffmerge_stats.txt"
    threads: 1
    message: "Preparing stats for cuffmerge"
    shell:"""
perl -ne '/transcript_id "(.*?)"/; $tx_id = $1; if(/class_code "(.*?)"/){{ print $tx_id,"\t",$1,"\n"}}' {input} \
| sort | uniq | cut -f2 | sort | uniq -c| sed 's/^ *//'| perl -npe 's/ /\t/' > {output}
"""