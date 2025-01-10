import os
def get_file_1(wildcards) :
	smp = config["comparison"][wildcards.comp].split(" ")
	smp = smp[0].split(",")
	smp = [os.path.join("output/bam" , x + ".bam") for x in smp]
	return(smp)

def get_file_2(wildcards) :
	smp = config["comparison"][wildcards.comp].split(" ")
	smp = smp[1].split(",")
	smp = [os.path.join("output/bam" , x + ".bam") for x in smp]
	return(smp)

rule cuffdiff:
  input:    a = get_file_1 , b= get_file_2
  output:   "output/cuffdiff/{comp}/{comp}_gene_exp.diff"
  params:   args=config["cuffdiff"]["args"] , gtf=config["gtf"], mem="40G" , labels = lambda wildcards: wildcards.comp.replace("vs", ",") 
  threads:  config["cuffdiff"]["threads"]
  message:  """--- counting with Cuffdiff ."""
  shell:    """
        mkdir -p output/cuffdiff/{wildcards.comp} 
        module load cufflinks/2.2.1 
        a=$(echo {input.a} | sed 's/ /,/')
        b=$(echo {input.b} | sed 's/ /,/')
        cuffdiff {params.args} -o output/cuffdiff/{wildcards.comp}/ -L {params.labels} -p {threads} {params.gtf} $a $b 
        cd output/cuffdiff/{wildcards.comp}
        for i in *; do mv $i {wildcards.comp}_$i; done
  """