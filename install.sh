#!/bin/bash

################################################ GLOBAL VARIABLES ################################################


# Error handling parameters
set -o errexit
set -o nounset
set -o pipefail


# Path to main directory
MAIN_DIR=$(echo $PWD)


################################################ FUNCTIONS ################################################


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


################################################ CONDA ENVIRONMENT SETTINGS ################################################


# Check for conda bin
log_stdout "Checking for Conda bin."

if ! command -v conda &> /dev/null
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
if conda list --name "$SNK_ENV_NAME" &> /dev/null
then
        log_stdout ""$SNK_ENV_NAME" environment is already present."

# Build INSinPAL conda environment then
elif conda env create --yes --file resources/conda/environment.yml
then
        log_stdout ""$SNK_ENV_NAME" environment was successfully created."
else
        log_error ""$SNK_ENV_NAME" environment coulnd't be build."
        exit 1
fi


# Activate INSinPAL Conda env
log_stdout "Activating INSinPAL Conda environment."

if source activate "$SNK_ENV_NAME"
then
        log_stdout ""$SNK_ENV_NAME" was successfully activated."
else
        log_error ""$SNK_ENV_NAME" couldn't be activated."
        exit 1
fi


################################################ DATA ACQUISITION ################################################


###################### REFERENCE GENOME hg19 FULL SEQUENCE ######################


# Reference genome hg19 masked repeat region from UCSC ftp with patches
# with the UCSC hg19 chrM "NC_001807" sequence and GenBank revised chrMT "NC_012920" sequence by Cambridge Reference Sequence.
log_stdout "Downloading reference genome hg19 (GRCh37.p13 + MT) with masked repeat regions."
# Set up downlading DIR
mkdir -p resources/data/Genome/hg19/chromosomes
PATH_DATA_HG19="resources/data/Genome/hg19"
# Final file name
HG19_FASTA_NAME="genome.masked.fasta"
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/p13.plusMT/hg19.p13.plusMT.fa.gz"

if wget --timestamping --quiet -O tmp.gz "$URL"
then
        log_stdout "hg19 reference genome was successfully downloaded to "$PATH_DATA_HG19" from UCSC ftp."
else
        log_error "hg19 reference genome couldn't be downloaded from UCSC ftp."
        exit 1
fi

# Uncompress reference genome OR skip if already exists
if gunzip --uncompress tmp.gz
then
        log_stdout "hg19 reference genome was successfully uncompressed."
else
        log_error "hg19 reference genome couldn't be uncompressed by gunzip."
        exit 1
fi

# Remove alternate haplotypes and loci, fix loci and keeps only the rCRS MT as they cause caller shut down
log_stdout "Remove alternative haplotypes from hg19."

awk 'BEGIN {RS=">"; ORS=""; FS="\n"; } (NR>1 && $1!~/hap/ && $1!~/fix/ && $1!~/alt/ && $1!="chrM") {print ">"$0}' tmp > "${PATH_DATA_HG19}"/"${HG19_FASTA_NAME}"
rm tmp

log_stdout "Rename contigs from >chr1 to >1 in reference genome."
if ! sed --in-place --regexp-extended -e "/^>/ s/chr//g" "${PATH_DATA_HG19}"/"${HG19_FASTA_NAME}"
then
	log_error "Couldn't change the contig names with sed."
	exit 1
fi

# Index fasta file
if samtools faidx "${PATH_DATA_HG19}"/"${HG19_FASTA_NAME}"
then
	log_stdout "hg19 fasta file was successfully indexed."
fi


###################### REFERENCE GENOME hg19 PER CHROMOSOME ######################


# Reference genome hg19 masked repeat regions per chromosome
log_stdout "Data acquisition: reference genome fasta per chromosome."
# Rename tar to
HG19_CHROMOSOMES="chromosomes.tar.gz"
URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.chromFa.tar.gz"

if wget --timestamping --quiet -O "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES" "$URL"
then
        log_stdout "hg19 reference genome chromosomes have been downloaded to "$PATH_DATA_HG19"/chromosomes from UCSC ftp."
else
        log_error "hg19 reference genome chromosomes couldn't be downdloaded."
        exit 1
fi


# Untar hg19 chromosomes

if tar -xvzf "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES" --directory "$PATH_DATA_HG19"/chromosomes
then
        log_stdout "Chromosomes have been untared."
        # Remove tar file
        rm "$PATH_DATA_HG19"/chromosomes/"$HG19_CHROMOSOMES"
        # Etract files from tar dir to target dir
        mv "$PATH_DATA_HG19"/chromosomes/chroms/* "$PATH_DATA_HG19"/chromosomes && rm -rf "$PATH_DATA_HG19"/chromosomes/chroms
else
        log_error "Untar process has failed with tar."
        exit 1
fi


# Rename chromosomes files to 1.fa instead of chr1.fa and remove files mathcing alternative alleles or MT
log_stdout "Rename chromosomes files from chr1.fa to 1.fa and remove unnecessary fasta files."

for file in $(ls "$PATH_DATA_HG19"/chromosomes/*); do
	# Get filename
	filename=$(basename "$file")
	# Remove contig, alternative chr files and Mitochondrial file
	if [[ $(grep --count --perl-regexp "chr.+_" "$file") -eq 1 ]] || [[ $(grep --count "chrM" "$file") -eq 1 ]]
	then
		rm "$file"
	else
		# Change file name from chr1.fa.masked to 1.fa
		newname=${filename#chr}
		# Change name in directory
		mv "$file" "$PATH_DATA_HG19"/chromosomes/"$newname"
	fi
done


################################################ SOFTWARE INSTALLATION ################################################

log_stdout "Downloading necessary software for INSinPAL environment."


###################### SINGULARITY ######################


# Check for singularity
log_stdout "Checking for singularity bin."

if command -v singularity > /dev/null
then
	log_stdout "Singularity bin found."
else
	log_error "Singularity not found, please install it."
	exit 1
fi


###################### INSURVEYOR 1.1.2 ######################


# DIrectory path for software
PROGS_DIR="resources/progs"
log_stdout "Downloading INSurVeyor singularity image from https://github.com/kensung-lab/INSurVeyor/releases/download/1.1.2/insurveyor.sif"
# DIrectory path for iNSurVeyor
INSURVEYOR_DIR="$PROGS_DIR"/"INSurVeyor/1.1.2/"
mkdir -p "$INSURVEYOR_DIR"
URL="https://github.com/kensung-lab/INSurVeyor/releases/download/1.1.2/insurveyor.sif"

if wget --timestamping --quiet --directory-prefix "$INSURVEYOR_DIR" "$URL"
then
	log_stdout "INSurVeyor was successfully donwloaded to "$INSURVEYOR_DIR"."
else
	log_error "INSurVeyor couldn't be downloaded."
	exit 1
fi


###################### MANTA 1.6.0 ######################


log_stdout "Downloading Manta binary distribution from https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2."
# Directory path for software
MANTA_DIR="$PROGS_DIR"/"Manta/1.6.0/"
mkdir -p "$MANTA_DIR"
URL="https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2"

if wget --timestamping --quiet --directory-prefix "$MANTA_DIR" "$URL"
then
	log_stdout "Manta 1.6.0 binary distribution was successfully downloaded to "$MANTA_DIR"."
else
	log_error "Manta couldn't be installed."
	exit 1
fi


log_stdout "Untar Manta distribution file."

if tar -xjvf "$MANTA_DIR"/$(basename "$URL") --directory "$MANTA_DIR"
then
	# Remove the compress file
	rm "$MANTA_DIR"/$(basename "$URL")
	log_stdout "Manta was successfully uncompressed in "$MANTA_DIR"."
else
	log_error "Manta couldn't be uncompressed."
	exit 1
fi


###################### BASIL 1.2.0 ######################


log_stdout "Checking for anisebasil.sif singularity image in "$PROGS_DIR"/Anis-Basil/1.2.0."
# Directory path to software
BASIL_DIR="$PROGS_DIR"/"Anis-Basil/1.2.0"

# Build anisebasil.sif if image not present
if [[ -f "$BASIL_DIR"/"anisebasil.sif" ]]
then
	log_stdout "Basil singularity image was successfully located in "$BASIL_DIR"."
else
	log_error "Basil singularity image is not present, please use the following command:\nsudo singularity build "$BASIL_DIR"/anisebasil.sif "$BASIL_DIR"/Anise_Basil_1_2_0.def."
	exit 1
fi


###################### SCRAMBLE 1.0.2 ######################


SCRAMBLE_DIR="${PROGS_DIR}"/"scramble/1.0.2"
log_stdout "Checking for scramble.sif singularity image in "${SCRAMBLE_DIR}"."

# Build scramble.sif if image not present
if [[ -f "${SCRAMBLE_DIR}"/"scramble.sif" ]]
then
	log_stdout "SCRAMBle singularity image was successfully located in "${SCRAMBLE_DIR}"."
else
	log_error "SCRAMBle singularity image is not present, please use the following command:\nsudo singularity build "$SCRAMBLE_DIR"/scramble.sif "$SCRAMBLE_DIR"/scramble.def."
	exit 1
fi


###################### ANNOTSV 3.4.2 ######################


log_stdout "Downloading AnnotSV code source from "$URL"."
# Directory path to software
ANNOTSV_DIR="$PROGS_DIR"/"AnnotSV/3.4"
URL="https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.4.2.tar.gz"
# Set directory tree
mkdir -p "$ANNOTSV_DIR"

if wget --timestamping --quiet --directory-prefix "$ANNOTSV_DIR" "$URL"
then
	log_stdout "AnnotSV code source was succesfully downloaded to "$ANNOTSV_DIR"."
else
	log_error "AnnotSV couldn't be downloaded."
	exit 1
fi


log_stdout "Untaring AnnotSV."

if tar -xvzf "$ANNOTSV_DIR"/$(basename "$URL")  --directory "$ANNOTSV_DIR"
then
	log_stdout "AnnotSV was successfully untared."
else
	log_error "AnnotSV couldn't be untared."
fi


log_stdout "Building code source for AnnotSV."
# Binary ditribution directory
BUILD_DIR="$ANNOTSV_DIR"/"AnnotSV"
# Rename folder
mv "$ANNOTSV_DIR"/"AnnotSV-3.4.2" "$BUILD_DIR"

# Build and capture exit code
(
	cd "$BUILD_DIR" && \
	make DESTDIR= PREFIX=. install && \
	make DESTDIR= PREFIX=. install-human-annotation
)
if [[ "$?" -eq 0 ]]
then
	log_stdout "AnnotSV was successfully built to "$BUILD_DIR"."
else
	log_error "AnnotSV couldn't be build."
	exit 1
fi


############################################### TEST ################################################


# Ask user for test or not with a test sample
while true; do
	read -p "Do you want to run a test with HG00096 sample from 1k genomes phase 3 data ? (y/n)" -n 1 ANS
	case $ANS in
		[Yy] ) log_stdout "Running test.."
			break;;
		[nN] ) log_stdout "INSinPAL is successfully installed."
			exit 0;;
		* ) echo "Invalid response";;
	esac
done

# Sample ID from GIAB
SAMPLE_ID="HG00096"
# Sample BAM path
SAMPLE_PATH="$PWD"/"$SAMPLE_ID".bam

log_stdout "Downloading HG00096.bam and bai."

# Download HG00096 BAM nd BAI from 1k genome ftp
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"

#if ! wget --timestamping --quiet -O "$SAMPLE_PATH" "$URL"
#then
#	log_error "Couldn't download HG00096.bam from 1k genome ftp."
#	exit 1
#fi

# BAI
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai"

#if ! wget -timestamping --quiet -O "$SAMPLE_PATH".bai "$URL"
#then
#	log_error "Couldn't download HG00096.bam.bai from 1k genome ftp."
#	exit 1
#fi

# Expected output from INSinPAL
OUTFILE="results/"${SAMPLE_ID}"/"${SAMPLE_ID}".xlsx"

log_stdout "Running test with sample ID: "${SAMPLE_ID}" and sample path: "${SAMPLE_PATH}"."

# Run analysis with run_analysis.sh
bash run_analysis.sh --sample "${SAMPLE_ID}" --path "${SAMPLE_PATH}"

# Check if exit code worked
if [[ $? -eq 0 ]] && [[ -f "${OUTFILE}" ]]
then
	log_stdout "Tests are successfully completed and INSinPAL is fully installed."
else
	log_error "Test didn't work well."
	exit 1
fi
