
# INSinPAL: Large Insertion Detection in Palindromic Fragile Sites

## Description:

INSinPAL is a Snakemake workflow designed to detect large insertion events, classified as structural variants (SVs) with a minimum size of 50 base pairs, from mapped paired-end sequencing reads. Given the high number of insertion calls, from whole genome sequencing (WGS) data, generated by SV callers, INSinPAL focuses specifically on palindromic fragile sites, which are known to undergo genomic rearrangment. This targeted approach enhance the detection of large insertions in WGS data relavant for rare disease diagnostics (i.e see [*DNA inverted repeats and human disease*](https://www.imrpress.com/journal/FBL/3/4/10.2741/A284)).


### Workflow steps

INSinPAL Snakemake workflow will proceed as following:
 - Mine palindromic sites from a reference genome and select fragile ones.
 - Call insertion SVs on sample BAM and merge callers output.
 - Select insertion SVs mapping to palindromic fragile sites.
 - Format VCFs to BEDs and annotate for MEIs, putative inserted sequences source and size.
 - An additional step for BAM formatting exists if mate score tags are absent, needed for INSurVeyo caller.

### Documentation

See ![DAG of jobs](./doc/dag.pdf) for workflow rules structure and/or read the [user guide](./doc/user_guide.md) detailed information about the workflow as setup steps, biological concepts, etc...
 

## Requirements

* Python2 and 3
* Conda == 24.1.2
* Singularity == 3.2.0-1

Note that, INSinPAL was not tested for versions lower or higher than the ones cited above.

## Getting Started

INSinPAL is a kind of metacaller using thre different SV callers and a few external data files that need to be setup:

* Callers:
  * [INSurVeyor 1.1.2](https://github.com/kensung-lab/INSurVeyor)
  * [Manta 1.6.0](https://github.com/Illumina/manta)
  * [Basil 1.1.0](https://github.com/seqan/anise_basil)
  * [SCRAMble 1.0.2](https://github.com/GeneDx/scramble)
* Data:
  * [GRCh37 reference genome with repeated regions masked](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.masked.gz)
  * [GRCh37 chromosomes](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz)
  * [MEIs reference sequences](resources/data/MEI_consensus_seqs_SCRAMble_plus_MOBSTER.fa)

To setup INSinPAL environment, one can use the `install.sh` bash file to donwload all required files and software if all the requirements are met.
```bash
bash install.sh
```
Note that, the `config/parameters.py` lists all paths necessary for Snakemake to work. If one want to change data sources such as the reference genome, he should also change this file to the correct paths.



