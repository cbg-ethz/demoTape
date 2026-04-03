# DemoTape
## Computational demultiplexing of targeted single-cell sequencing (tapestri) data

DemoTape is a computational demultiplexing method for targeted single-cell DNA sequencing (scDNA-seq) data, namely MissiobIo Tapestri data, based on a distance metric between individual cells at single-nucleotide polymorphisms loci. 

The corresponding preprint can be found [here](https://www.biorxiv.org/content/10.1101/2024.12.06.627152v1).

## Requirements

### Software
- Python >=v3.X
- [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) >=8.X


All other requirements are installed automatically via Snakemake in separate conda envs.

### Resources
The following resources need to be
- The annotation file (.bed) for the used Tapestri panel
- dbsnp file (.bed or .txt) for the used reference genome (e.g., hg19)


## Usage

### Running
#### Only demultiplexing
To run only DemoTape, you can run:
```bash
python workflow/scripts/run_demoTape.py -i <VARIANTS.VCF> -n <NO_SAMPLES>
```
where `<VARIANTS.VCF>` is the .csv file produced by the MissionBio [Mosaic Pipeline](https://github.com/MissionBio/mosaic).


Alternatively, starting from the loom file, you can also first run 
```bash
python workflow/scripts/mosaic_preprocessing.py -i <INPUT.LOOM>
```
(This is what happens if the whole DemoTape pipeline is run)


#### Whole pipeline (demultiplexing, assigning sample→patient, plotting, downstream analysis)
The whole DemoTape analysis pipeline can be executed via:
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
According to the running environment (local/HPC), the `executor` needs to be adjusted.


### Config 
In the [config file](configs/MS1_analysis.yaml), the following variables need to be specified:
```
analysis:
  specific:
    input-dir: <INPUT_DIR>
    output-dir: <OUTPUT_DIR>
  general:
    panel_annotation: resources/<ANNOTATED_TAPESTRI_PANEL>.bed

output:
  prefix: <PREFIX>
```
The Tapestri panel file can be annotated (i.e., gene names assigned to loci) via [BED Annotation](https://github.com/vladsavelyev/bed_annotation).


Additionally, to run downstream analysis with [BnpC](https://github.com/cbg-ethz/BnpC) or [COMPASS](https://github.com/cbg-ethz/compass), the corresponding software needs to be downloaded and the py/exe files. specified

## Simulations

To run the simulation pipeline, execute:
```bash
snakemake 
    -s workflow/Snakefile_simulations
    -j 500
    --configfile configs/simulations.yaml
    --executor slurm
    --rerun-incomplete
    --drop-metadata
    -k 
    --use-conda
```
where `input-looms` as well as the `exe` files for souporcell and scSplit needs to be adjusted