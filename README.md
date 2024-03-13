# demoTape
Computational demultiplexing of targeted single-cell sequencing (tapestri) data

---
To submit the pipeline on an HPC via slurm:

    sbatch 
        -n 1 
        --cpus-per-task=1
        --time=48:00:00
        --mem-per-cpu=1024
        --output="logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out"
        --error="logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).err"
        --open-mode=truncate
        --wrap="snakemake -s Snakefile -j 500 --configfile ../configs/multiplex_sims.yaml --profile <slurm_profile> --rerun-incomplete --drop-metadata -k"
