import os
def get_file_1(wildcards) :
	smp = config["comparison"]["RipmOVAvsOTII"].split(" ")
	smp = smp[0].split(",")
	smp = [os.path.join("output/bam" , x + ".bam") for x in smp]
	return(smp)

def get_file_2(wildcards) :
	smp = config["comparison"]["RipmOVAvsOTII"].split(" ")
	smp = smp[1].split(",")
	smp = [os.path.join("output/bam" , x + ".bam") for x in smp]
	return(smp)




rule cuffdiff:
  input:    a = get_file_1 , b= get_file_2
  output:   "output/cuffdiff/RipmOVA_OTII.done"
  params:   args=config["cuffdiff"]["args"] , gtf=config["gtf"], mem="40G" , labels="RipmOVA,OTII"
  threads:  config["cuffdiff"]["threads"]
  message:  """--- counting with Cuffdiff ."""
  shell:    """
        mkdir -p output/cuffdiff/ 
        module load cufflinks/2.2.1 
		a=`echo {input.a} | sed 's/ /,/' `
		b=`echo {input.b} | sed 's/ /,/' `
        cuffdiff {params.args} -o output/cuffdiff/ -L {params.labels} -p {threads} {params.gtf} $a $b 
        touch {output}
            """