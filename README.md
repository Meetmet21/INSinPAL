
# INSinPAL: Large Insertion Detection in Palindromic Fragile Sites

## Description:

INSinPAL is a Snakemake-based workflow designed to detect large insertion events, classified as structural variants (SVs) with a minimum size of 50 base pairs, from mapped paired-end sequencing reads. Given the high number of insertion calls generated by SV callers from whole genome sequencing (WGS) data, INSinPAL focuses specifically on palindromic fragile sites, which are known to undergo genomic rearrangment. This targeted approach enhance the detection of large insertions in WGS data relavant for rare disease diagnostics (i.e see [*DNA inverted repeats and human disease*](https://www.imrpress.com/journal/FBL/3/4/10.2741/A284)).


### Workflow steps

INSinPAL Snakemake workflow will proceed as following:
 - Mine palindromic sites from a reference genome and select fragile ones.
 - Call insertion SVs on sample BAM and merge callers output.
 - Select insertion SVs mapping to palindromic fragile sites.
 - Format VCFs to BEDs and annotate for MEIs, putative inserted sequences source and size.
 - An additional step for BAM formatting exists if mate score tags are absent, needed for INSurVeyo caller.

### Documentation

See ![DAG of jobs](./doc/dag.pdf) for workflow rules structure and/or read the [user guide](./doc/user_guide.md) detailed information about the workflow as biological concepts, input/output formats, etc...
 

## Requirements

* Python2 and 3
* Conda == 24.1.2
* Singularity == 3.2.0
* Make
* Linux-based machine

Note that, INSinPAL was not tested for requirements versions lower or higher than the ones cited above.

## Getting Started

INSinPAL is a type of metacaller using three different SV callers and requires a few external data files to be set up.

* Callers:
  * [INSurVeyor 1.1.2](https://github.com/kensung-lab/INSurVeyor)
  * [Manta 1.6.0](https://github.com/Illumina/manta)
  * [Basil 1.1.0](https://github.com/seqan/anise_basil)
  * [SCRAMble 1.0.2](https://github.com/GeneDx/scramble)
* Data:
  * [GRCh37 reference genome with repeated regions masked](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.masked.gz)
  * [GRCh37 chromosomes](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz)
  * [MEIs reference sequences](resources/data/MEI_consensus_seqs_SCRAMble_plus_MOBSTER.fa)

To set up the INSinPAL environment, one can use the `install.sh` bash file to download all required files and software if all the requirements are met.
```bash
# Clone this repo
git clone git@github.com:Meetmet21/INSinPAL.git
# Go to working directory
cd INSinPAL
# Use sudo privelge to build singularity images
sudo singularity build resources/singularity/Anis-Basil/1.2.0/anisebasil.sif resources/singularity/Anis-Basil/1.2.0/Anise_Basil_1_2_0.def
# Run remaining installation steps
bash install.sh
```
Note that, the `config/parameters.py` lists all the paths necessary for Snakemake to work. If one wants to change data sources 
such as the reference genome, they should also update this file to the correct paths. Moreover, the default cluster profile 
is set to LSF in `workflow/profile/lsf-status.py` and needs to be change if one wants to submit jobs on SLURM. Finally, the 
`workflow/profile/config.yaml` file sets Snakemake's default parameters and can be changed if needed.

## Usage

To perform an analysis with INSinPAL, one can use the `run_analysis.sh` script:
```bash
# Access the help section
bash run_analysis.sh --help
# Run analysis
bash run_analysis.sh --sample SAMPLE_ID --path SAMPLE_BAM
```
Note that, `SAMPLE_ID` should not contain any sepecial characters except underscores and `SAMPLE_PATH` should lead to the absolute path
of a BAM and respective BAI file.

For each analysis, if `run_analysis.sh` is used, logs related to rules are stored in `workflow/logs/sample_ID/` and the 
ones related to LSF jobs are stored in `logs/` with dry runs for a given sample in `logs/dry_runs/`.

## Contact

You can cantact me via my email: mehmetsehir1@gmail.com.

