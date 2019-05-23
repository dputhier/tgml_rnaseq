import os

rule trimming_pe:
  input:  fwd="input/fastq/{smp}_R1.fq.gz", rev="input/fastq/{smp}_R2.fq.gz"
  output: fwd="output/trimmed/{smp}_R1_t.fq.gz", 
          rev="output/trimmed/{smp}_R2_t.fq.gz",
          single="output/trimmed/{smp}_R1_singletons.fq.gz"
  params: qual=config["sickle"]["threshold"], len=config["sickle"]["length"], qt=config["sickle"]["qualtype"]
  threads: 1
  conda: os.path.join(config["workingdir"], "conda", "rnaseq.yaml")
  message: """--- Trimming."""
  shell: """
        sickle pe -g -f {input.fwd} -r {input.rev}  -l {params.len} -q {params.qual} -t {params.qt}  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log 
  """
  
