#!/bin/bash

# Restrictive error handling
set -o errexit
set -o nounset
set -o pipefail

# Print error message
log_error() {
        echo -e 'Error: '"$1" >&2
}

# Print log messsage
log_stdout() {
        echo -e 'Install log: '"$1"
}

### Conda ###

# Check for conda bin
log_stdout "Checking for Conda bin."
if ! $(command -v conda &> /dev/null)
then
        log_error "Conda is not detected and need to be installed and initiated for INSinPAL setup process."
        exit 1
fi

# INSinPAL conda environment setup
log_stdout "Installing INSinPAL Conda environment."
SNK_ENV_NAME="INSinPAL"
# Snakemake environment already present
if $(conda list --name "$SNK_ENV_NAME" &> /dev/null)
then
        log_stdout ""$SNK_ENV_NAME" environment is already present."
# Build INSinPAL conda environment then
elif $(conda env create --name "$SNK_ENV_NAME" --file=resources/conda/snakemake_env.yaml)
then
        log_stdout ""$SNK_ENV_NAME" environment successfully created."
else
        log_error ""$SNK_ENV_NAME" environment coulnd't be build."
        exit 1
fi

# Activate Snakemake Conda env
log_stdout "Activating INSinPAL Conda environment."
if $(source activate "$SNK_ENV_NAME")
then
        log_stdout ""$SNK_ENV_NAME" was successfully activated."
else
        log_error ""$SNK_ENV_NAME" couldn't be activated."
        exit 1
fi

### Data ###
log_stdout "Data acquisition."

# Reference genome hg19 masked repeat region from UCSC ftp
PATH_DATA_HG19="resources/data/Genome/hg19"
HG19_FASTA_NAME="genome_PAR_masked.fasta.gz"
if $(curl -o "$PATH_DATA_HG19"/"$HG19_FASTA_NAME" "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.masked.gz")
then
        log_stdout "hg19 reference genome was successfully downloaded to "$PATH_DATA_HG19" from UCSC ftp."
else
        log_error "hg19 reference genome couldn't be downloaded from UCSC ftp."
        exit 1
fi

# Uncompress reference genome OR skip if already exists
if [[ -f "$PATH_DATA_HG19"/${HG19_FASTA_NAME%.gz} ]] || $(gunzip -d "$PATH_DATA_HG19"/"$HG19_FASTA_NAME")
then
        log_stdout "hg19 reference genome successfully uncompressed."
else
        log_error "hg19 reference genome couldn't be uncompressed by gunzip."
        exit 1
fi

# Reference genome hg19 masked repeat regions per chromosome
HG19_CHROMOSOMES="chromosomes.tar.gz"
if $(curl -o "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES" https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz)
then
        log_stdout "hg19 reference genome chromosomes have beend downloaded to "$PATH_DATA_HG19"/chromosomes from UCSC ftp."
else
        log_error "hg19 reference genome chromosomes couldn't be downdloaded."
        exit 1
fi

# Untar chromosomes
if $(tar -xvzf "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES")
then
        log_stdout "Chromosomes have been untared."
else
        log_error "Untar process has failed with tar."
fi

# Rename chromosomes files to 1.fa than chr1.fa and remove files mathcing alternative alleles or MT
for

### Progs ###


### Test ###
