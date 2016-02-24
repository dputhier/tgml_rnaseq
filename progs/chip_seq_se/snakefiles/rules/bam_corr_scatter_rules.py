rule bam_corr_scatter:
    input: npz="output/multiBamSummary/merge_peaks_coverage.npz"
    output: scatter="output/multiBamSummary/scatterplot_spearman.png", \
            scatter_mat="output/multiBamSummary/scatterplot_pearson.tab", \
            heatmap="output/multiBamSummary/heatmap_spearman.png", \
            heatmap_mat="output/multiBamSummary/heatmap_spearman.tab", \
            pca="output/multiBamSummary/pca.png"

    shell: """
    plotCorrelation \
    -in {input} \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation computed across merged peaks." \
    --whatToPlot scatterplot \
    -o {output.scatter}  \
    --outFileCorMatrix {output.scatter_mat}
    
    plotCorrelation \
    -in {input} \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o {output.heatmap}   \
    --outFileCorMatrix {output.heatmap_mat} 

    plotPCA -in {input} \
        -o {output.pca} \
        -T "PCA of read counts"
    """