# demoTape
Computational demultiplexing of targeted single-cell sequencing (tapestri) data

---
To submit the pipeline on an HPC via slurm executor:

    sbatch 
        -n 1 
        --cpus-per-task=1
        --time=48:00:00
        --mem-per-cpu=1024
        --output="logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).log"
        --open-mode=truncate
        --wrap="snakemake -s workflow/Snakefile_analysis -j 500 --configfile configs/multiplex_sims.yaml --executor slurm --rerun-incomplete --drop-metadata --use-conda -k"