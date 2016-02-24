SAMPLES = config["samples"].split()
CHIP = SAMPLES
CHIP.remove(config["input"])

rule merge_peaks:
        input: expand("output/macs/{smp}/{smp}_sorted_peaks.bed", smp=CHIP)
        output: "output/merged_peaks/merged_peaks.bed"
        message: """--- Merging peaks."""
        shell: """
        cat output/macs/*/*_sorted_peaks.bed | sortBed | mergeBed > {output}
        """