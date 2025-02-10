# Makefile

help:
	@echo "To start the analyses"
	@echo ">>> type make run to run the Workflow <<<"
	@echo ">>> type make dry for a dry run <<<"
	@echo ">>> type make clean for the cleaning of slurm.out files <<<"
	@echo ">>> type make queue to check the current jobs in queue <<<"

# Run the workflow
run:
	@bash -c "module unload snakemake; \
	module load snakemake/7.25.0 ; \
	snakemake -s snakefiles/Snakefile.py --cluster 'sbatch -c {threads} --mem {params.mem} --partition=fast --account=2427_data_master' -c 3000 -j 500 --rerun-incomplete --rerun-trigger mtime"
	
dry:
	@bash -c "module unload snakemake; \
	module load snakemake/7.25.0 ; \
	snakemake -s snakefiles/Snakefile.py --cluster 'sbatch -c {threads} --mem {params.mem} --partition=fast --account=2427_data_master' -c 3000 -j 500 --rerun-incomplete --rerun-trigger mtime -np"

	
#Clean
clean:
	rm /shared/projects/project_dputhier/043_irla_elife_2022/slurm*.out

scancel:
	scancel `squeue -u $USER | perl -npe 's/^\s+//' | grep -v JOBID | perl -npe 's/\s+.*//'`

queue:
	@squeue -u $$USER
