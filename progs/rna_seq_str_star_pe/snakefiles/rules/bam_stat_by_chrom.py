import os
script = os.path.join(config["workingdir"], "progs", "scripts", "reads_per_chrom.py")
python_path = config["python2_path"]


rule bam_stat_by_chrom:
        input: bam="output/bam/{smp}.bam"
        output: o="output/bam_stat_by_chrom/{smp}_bam_stats.txt", d="output/bam_stat_by_chrom/{smp}_bam_stats.png"
        params: script=script, python_path=python_path
        message: """--- Performing some stats about mapping (by chrom)."""
        threads: 1
        shell:  """
         {params.python_path} {params.script} -i {input} -o {output.o} -d  {output.d} &> {output[0]}.log
                """
