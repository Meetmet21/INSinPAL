#!/bin/bash

############################################### GLOBAL VARIABLES ###############################################


set -o errexit
set -o nounset
set -o pipefail

###################### PARSE CL ARGUMENTS ######################

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
	echo "######"
	echo -e 'Error: '"$1" >&2
}

# Print log message
log_stdout() {
	echo "###########"
	echo -e "Launch-run: ""$1""."
}


############################################### MAIN PROGRAM ###############################################


log_stdout "Preparing "$SAMPLE_ID" run with "$SAMPLE_PATH" file.""

# Snakemake config file
CONFIG="config/config.yaml"

# New sample name and path in snakemake config file
if ! $(sed -r -i -e "s|^[[:space:]]+[^:]+|  "$SAMPLE_ID"|" -e "s|[^:]+$| "$SAMPLE_PATH"|" "$CONFIG")
then
	log_error "sed command error with main config file. Couldn't change "$CONFIG" file with new sample ID and path.""
	exit 1
else
	log_stdout "Changing sample ID and path in "$CONFIG"."
fi

# Snakemake profile config file 
CONFIG_LSF="workflow/profile/config.yaml"

# LSF log file name for new sample
if ! $(sed -r -i -e "s|lsf.*\.err|lsf_"$SAMPLE_ID".err|" -e "s|lsf.*\.out|lsf_"$SAMPLE_ID".out|" "$CONFIG_LSF")
then
	log_error "sed command error with lsf config file IN "$CONFIG_LSF"."
	exit 1
else
	log_stdout "Modifing Snakemake config file for LSF log file names."
fi

# Activate Conda env
ENV="INSinPAL"

if $(source activate "$ENV")
then
	log_stdout "Activaiting "$ENV" Conda environment."
else
	log_error "Couldn't acitvate "$ENV" with Conda."
	exit 1
fi

# Dry run pipeline for sample
log_stdout "Running dry-run for "$SAMPLE_ID" with Snakemake. TXT dry_run in logs/dry_runs/"$SAMPLE_ID".txt"
snakemake --profile workflow/profile -n | tee logs/dry_runs/"$SAMPLE_ID".txt

# Pipeline running in screen session
screen -S pipeline -X stuff "\rsnakemake --profile workflow/profile --rerun-triggers mtime\n"

log_stdout "Running pipeline for "$SAMPLE" in screen session"

# Program succeed
exit 0
