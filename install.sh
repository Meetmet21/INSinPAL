#!/bin/bash

# Restrictive error handling
set -o errexit
set -o nounset
set -o pipefail

### Functions

# Print error message
log_error() {
		echo -e '######'
        echo -e 'Error: '"$1" >&2
}

# Print log messsage
log_stdout() {
		echo -e '############'
        echo -e 'Install log: '"$1"
}

### Conda ###

# Check for conda bin

log_stdout "Checking for Conda bin."

if ! $(command -v conda &> /dev/null)
then

	log_error "Conda is not detected and need to be installed."

	exit 1
else

	log_stdout "Conda bin found."

fi

# INSinPAL conda environment setup

log_stdout "Installing INSinPAL Conda environment."

SNK_ENV_NAME="INSinPAL"
# Snakemake environment already present
if $(conda list --name "$SNK_ENV_NAME" &> /dev/null)
then

        log_stdout ""$SNK_ENV_NAME" environment is already present."

# Build INSinPAL conda environment then
elif $(conda create --name "$SNK_ENV_NAME" --file=resources/conda/snakemake_env.yaml)
then

        log_stdout ""$SNK_ENV_NAME" environment was successfully created."

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

log_stdout "Data acquisition: reference genome hg19 masked regions."

# Reference genome hg19 masked repeat region from UCSC ftp
mkdir -p resources/data/Genome/hg19/chromosomes
PATH_DATA_HG19="resources/data/Genome/hg19"
HG19_FASTA_NAME="genome_PAR_masked.fasta.gz"
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.masked.gz"
if [[ -f "$PATH_DATA_HG19"/"$HG19_FASTA_NAME" ]] || $(curl -o "$PATH_DATA_HG19"/"$HG19_FASTA_NAME" "$URL")
then

        log_stdout "hg19 reference genome was successfully downloaded to "$PATH_DATA_HG19" from UCSC ftp."

else

        log_error "hg19 reference genome couldn't be downloaded from UCSC ftp."

        exit 1
fi

# Uncompress reference genome OR skip if already exists
if [[ -f "$PATH_DATA_HG19"/${HG19_FASTA_NAME%.gz} ]]
then

        log_stdout "hg19 reference genome successfully uncompressed."

        rm "$PATH_DATA_HG19"/"$HG19_FASTA_NAME"
elif $(gunzip -d "$PATH_DATA_HG19"/"$HG19_FASTA_NAME")
then

        log_stdout "hg19 reference genome was successfully uncompressed."

else

        log_error "hg19 reference genome couldn't be uncompressed by gunzip."

        exit 1
fi

log_stdout "Data acquisition: reference genome fasta per chromosome."

# Reference genome hg19 masked repeat regions per chromosome
HG19_CHROMOSOMES="chromosomes.tar.gz"
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz"
if [[ -f "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES" ]] || $(curl -o "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES" "$URL")
then

        log_stdout "hg19 reference genome chromosomes have been downloaded to "$PATH_DATA_HG19"/chromosomes from UCSC ftp."

else

        log_error "hg19 reference genome chromosomes couldn't be downdloaded."

        exit 1
fi

# Untar chromosomes
if $(tar -xvzf "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES" --directory "$PATH_DATA_HG19"/chromosomes &> /dev/null)
then

        log_stdout "Chromosomes have been untared."

        # Remove tar file
        rm "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES"
else

        log_error "Untar process has failed with tar."

        exit 1
fi

# Rename chromosomes files to 1.fa than chr1.fa and remove files mathcing alternative alleles or MT

log_stdout "Rename chromosomes files from chr1.fa.masked to 1.fa and remove unnecessary fasta files."

for file in $(ls "$PATH_DATA_HG19"/chromosomes/*); do
	# Filename
	filename=$(basename "$file")
	# Remove contif, alternative chr files and Mitochondrial file
	if [[ $(grep -cP "chr.+_" "$file") -eq 1 ]] || [[ $(grep -cP "chrM" "$file") -eq 1 ]];
	then
		rm "$file"
	else
		# Change file name from chr1.fa.masked to 1.fa
		newname=${filename#chr}
		newname=${newname%.masked}
		# Change name in directory
		mv "$file" "$PATH_DATA_HG19"/chromosomes/"$newname"
	fi
done

### Progs ###

log_stdout "Downloading needed software."

PROGS_DIR="resources/singularity"
# Check for singularity

log_stdout "Checking for singularity bin."

if [[ $(command -v singularity /dev/null) ]]
then

	log_stdout "Singularity bin found."

else

	log_error "Singularity not found, please install it."

	exit 1
fi

# INSurVeyor
INSURVEYOR_DIR="$PROGS_DIR"/"INSurVeyor/1.1.2/"
mkdir -p "$INSURVEYOR_DIR"
URL="https://github.com/kensung-lab/INSurVeyor/releases/download/1.1.2/insurveyor.sif"

log_stdout "Downloading INSurVeyor singularity image from https://github.com/kensung-lab/INSurVeyor/releases/download/1.1.2/insurveyor.sif"

if $(wget --quiet --directory-prefix "$INSURVEYOR_DIR" "$URL")
then

	log_stdout "INSurVeyor was successfully donwloaded to "$INSURVEYOR_DIR"."

else

	log_error "INSurVeyor couldn't be downloaded."

	exit 1
fi

# Manta
MANTA_DIR="$PROGS_DIR"/"manta/1.6.0/"
mkdir -p "$MANTA_DIR"
URL="https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2"

log_stdout "Downloading Manta binary distribution from https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2."

if $(wget --quiet --directory-prefix "$MANTA_DIR" "$URL")
then

	log_stdout "Manta 1.6.0 binary distribution was successfully downloaded to "$MANTA_DIR"."

else

	log_error "Manta couldn't be installed."

	exit 1
fi


log_stdout "Untar Manta distribution file."

if $(tar -xjvf "$MANTA_DIR"/$(basename "$URL") --directory "$MANTA_DIR" &> /dev/null)
then
	# Remove the compress file
	rm "$MANTA_DIR"/$(basename "$URL")

	log_stdout "Manta was successfully uncompressed in "$MANTA_DIR"."

else

	log_error "Manta couldn't be uncompressed."

	exit 1
fi

# Basil

log_stdout "Building Basil singularity image."

BASIL_DIR="$PROGS_DIR"/"Anis-Basil/1.2.0"
if $(singularity build "$BASIL_DIR"/"anisebasil.sif" "$BASIL_DIR"/"Anise_Basil_1_2_0.def")
then

	log_stdout "Basil singularity image was successfully built in "$BASIL_DIR"."

else

	log_error "Basil singularity image couldn't be build."

	exit 1
fi

# AnnotSV

URL="https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.4.2.tar.gz"
log_stdout "Downloading AnnotSV code source from "$URL"."
ANNOTSV_DIR="$PROGS_DIR"/"AnnotSV/3.4"
mkdir -p "$ANNOTSV_DIR"

if $(wget --quiet --directory-prefix "$ANNOTSV_DIR" "$URL")
then

	log_stdout "AnnotSV code source was succesfully downloaded to "$ANNOTSV_DIR"."

else
	log_error "AnnotSV couldn't be downloaded."
	exit 1
fi

log_stdout "Untaring AnnotSV."

if (tar -xvzf "$ANNOTSV_DIR"/$(basename "$URL")  --directory "$ANNOTSV_DIR" &> /dev/null)
then
	log_stdout "AnnotSV was successfully untared."
else
	log_error "AnnotSV couldn't be untared."
fi

log_stdout "Building code source for AnnotSV."
if
### Test ###
