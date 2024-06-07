#!/bin/bash

############################################### GLOBAL VARIABLES ###############################################


set -o errexit
set -o nounset
set -o pipefail

# No arguments given
if [[ "$#" -eq 0 ]]
then
	echo "Usage:"
	echo -e "\t"$0" -s/--sample <sample_id> -p/--path <sample_absolute_path> [-h/--help]"
	exit 0
fi

# Help function
help() {
	echo -e ""$0"\n"
	echo -e "This script runs INSinPAL snakemake-based workflow given two paraemeters:\n"
	echo -e "\t-s/--sample:\tSample ID for the run. Note that all resulting files will contain the given sample ID."
	echo -e "\t-p/--path:\tSAMPLE BAM file absolute path for the run.\n"
	exit 0
}

SAMPLE_ID=
SAMPLE_PATH=

# Parse Command-line arguments
while [[ "$#" -gt 0 ]]; do
	case $1 in
		-s|--sample)
			SAMPLE_ID="$2"
			shift 2 # Shift to next position argument
			;;
		-p|--path)
			SAMPLE_PATH="$2"
			shift 2
			;;
		-h|--help)
			help
			exit 0
			;;
		-*|--*)
			echo "ERROR: Unknown option: "$1""
			exit 1
			;;
		*)
			echo "ERROR: Invalid parameter: "$1""
			exit 1
			;;
	esac
done

# Check sample and path are defined
if [[ -z "$SAMPLE_ID" ]] || [[ -z "$SAMPLE_PATH" ]]
then
	echo "ERROR: You must set the -s/--sample and/or -p/--path argument."
	exit 1
fi
############################################### FUNCTIONS ###############################################

# Print error message
log_error() {
	echo 'Launch-run Error: '"$1" >&2
}

# Print log message
log_stdout() {
	echo "Launch-run: ""$1""."
}


# Grep the available BAM file name and remove pattern
SAMPLE=$(basename "$DIR_BAM" | sed "s/.genome.bam//")

log_stdout "Preparing "$SAMPLE" run..."

# Snakemake config file
path_config_file="$PWD""/config/config.yaml"
# New sample name and path in snakemake config file
if ! $(sed -r -i -e "s|^[[:space:]]+[^:]+|  "$SAMPLE"|" -e "s|[^:]+$| "$DIR_BAM"|" "$path_config_file")
then
	log_error "sed command error with main config file"
	exit 1
fi

log_stdout "Modifing Snakemake config file"

# Snakemake profile config file 
lsf_config="$PWD""/workflow/profile/config.yaml"
# LSF log file name for new sample
if ! $(sed -r -i -e "s|lsf.*\.err|lsf_"$SAMPLE".err|" -e "s|lsf.*\.out|lsf_"$SAMPLE".out|" "$lsf_config")
then
	log_error "sed command error with lsf config file"
	exit 1
fi

log_stdout "Modifing workflow config file"

# Activate snakemake env
source activate snakemake

log_stdout "Activating Conda snakemake"

# Dry run pipeline for sample
snakemake --profile workflow/profile -n > logs/dry_runs/"$SAMPLE"_dry_run.txt

log_stdout "Running dry-run for "$SAMPLE""

# Pipeline running in screen session
screen -S pipeline -X stuff "\rsnakemake --profile workflow/profile --rerun-triggers mtime\n"

log_stdout "Running pipeline for "$SAMPLE" in screen session"

# Program succeed
exit 0
