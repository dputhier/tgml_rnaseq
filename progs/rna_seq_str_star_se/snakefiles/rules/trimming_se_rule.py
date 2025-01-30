rule trimming:
  input:  fwd="input/fastq/{smp}.fq.gz"
  output: fwd="output/trimmed/{smp}.fq.gz"
  params: qual=config["sickle"]["threshold"], len=config["sickle"]["length"], qt=config["sickle"]["qualtype"], mem="5G"
  threads: 1
  message: """--- Trimming---"""
  shell: """
        module load sickle-trim/1.33
        sickle se -g -f {input.fwd}  -l {params.len} -q {params.qual} -t {params.qt}  -o {output.fwd}  &> {input.fwd}.log 
  """
  
