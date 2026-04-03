# DemoTape
## Computational demultiplexing of targeted single-cell sequencing (tapestri) data

DemoTape is a computational demultiplexing method for targeted single-cell DNA sequencing (scDNA-seq) data, namely MissiobIo Tapestri data, based on a distance metric between individual cells at single-nucleotide polymorphisms loci. 

The corresponding preprint can be found [here](https://www.biorxiv.org/content/10.1101/2024.12.06.627152v1).

## Requirements

- Python >=v3.X
- [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) >=8.X


All other requirements are installed automatically via Snakemake in separate conda envs


## Usage

The whole demoTape pipeline can be executed via:
```bash
snakemake 
    -s workflow/Snakefile_analysis
    -j 500
    --configfile configs/MS1_analysis.yaml
    --executor slurm
    --rerun-incomplete
    --drop-metadata
    -k 
    --use-conda
```
According to the running environment (local/HPC), the `executor` needs to be adjusted

In the config file, the following variables need to be specified:
```
analysis:
  specific:
    input-dir: <INPUT_DIR>
    output-dir: <OUTPUT_DIR>
output:
  prefix: <PREFIX>
```



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