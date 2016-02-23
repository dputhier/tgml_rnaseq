rule trimming_se:
  input:  "input/fastq/{smp}.fq.gz"
  output: "output/trimmed/{smp}_t.fq.gz"
  params: qual=config["sickle"]["threshold"], len=config["sickle"]["length"], qt=config["sickle"]["qualtype"]
  threads: 1
  message: """--- Trimming."""
  shell: """
    sickle se -g -f {input} -l {params.len} -q {params.qual} -t {params.qt}  -o {output} &> {input}.log 
  """
  
