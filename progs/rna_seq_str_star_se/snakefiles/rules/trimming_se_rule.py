rule trimming_pe:
  input:  fwd="input/fastq/{smp}.fq.gz"
  output: fwd="output/trimmed/{smp}_t.fq.gz"
  params: qual=config["sickle"]["threshold"], len=config["sickle"]["length"], qt=config["sickle"]["qualtype"]
  threads: 1
  message: """--- Trimming."""
  shell: """
        sickle se -g -f {input.fwd}  -l {params.len} -q {params.qual} -t {params.qt}  -o {output.fwd}  &> {input.fwd}.log 
  """
  
