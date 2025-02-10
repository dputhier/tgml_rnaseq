rule gtftk_get_biotypes:
  input:    gtf=config["gtf"]
  output:   biotypes="output/gene_info/biotypes.tsv", lncrna="output/gene_info/all_gene_coord.bed"
  params:   mem="40G"
  threads:  1
  message:  """--- Searching for gene biotypes."""
  shell:    """
            echo "Loading modules"
            module unload python
            module unload bedtools
            module load bedtools/2.30.0     
            module load  python/3.9
            pip install -U pygtftk

            echo "Calling gtftk (1)"
            gtftk convert_ensembl -V 3 -i {input.gtf} | \
                        gtftk select_by_key -V 3 -g | \
                        gtftk tabulate -V 3 -k gene_id,gene_name,gene_biotype -u -x -o {output.biotypes}

            echo "Calling gtftk (2)"
            gtftk convert_ensembl -V 3 -i {input.gtf} | \
                        gtftk select_by_key -V 3 -g | \
                        gtftk convert -n gene_name,gene_biotype -f bed6 -o {output.lncrna}


  """
