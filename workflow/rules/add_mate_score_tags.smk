#!/usr/bin/env python3

"""
:Author: Mehmet Sehir
:Date: 14.05.2024

Set of rules related to BAM formatting. Final purpose is to add MC tags to BAMs for INSurVeyor caller, otherwise,
it doesn't work.
"""

# MODULES
import sys
from os.path import join

# PARAMETER FILE
sys.path.append("config/")
import parameters

# Path to pipeline main directory directories
paths = parameters.WorkFlowPaths()

# CONFIG file
configfile: join(paths.config, "config.yaml")

def get_path_sample_from_config(wildcards):
    """
    Access to sample WGS data path from config file.
    :param wildcards:
    :return:
    """
    return config["samples"][wildcards.sample]

rule sort_bam_by_name:
    """
    Sort BAM by QNAMW for downstream formatting.
    """
    input:
        bam = get_path_sample_from_config
    output:
        sorted_bam = temp("{sample}.sortedname.bam")
    conda:
        join(paths.envs, "mapping.yaml")
    threads: 12
    shadow: "full"
    shell:
        "samtools sort -n -@ {threads} -O BAM -o {output.sorted_bam} {input.bam}"

rule fixmate_bam:
    """
    Add mate score tags to BAM.
    """
    input:
        names_sorted = rules.sort_bam_by_name.output.sorted_bam
    output:
        fixmate_bam = temp("{sample}.fixmate.bam")
    conda:
        join(paths.envs,"mapping.yaml")
    threads: 12
    shadow: "full"
    shell:
        "samtools fixmate -m -@ {threads} -O BAM {input.names_sorted} {output.fixmate_bam}"

rule sort_final_bam_by_coordinates:
    """
    Resort to coordinates.
    """
    input:
        fixmate_bam = rules.fixmate_bam.output.fixmate_bam
    output:
        formatted_bam = temp(join(paths.results, "{sample}/{sample}.bam"))
    conda:
        join(paths.envs, "mapping.yaml")
    threads: 12
    shadow: "full"
    shell:
        "samtools sort -@ {threads} -O BAM -o {output.formatted_bam} {input.fixmate_bam}"

rule index_final_bam:
    """
    Index the final BAM file.
    """
    input:
        final_bam=rules.sort_final_bam_by_coordinates.output.formatted_bam
    output:
        bai=temp(rules.sort_final_bam_by_coordinates.output.formatted_bam + ".bai")
    conda:
        join(paths.envs, "mapping.yaml")
    threads: 12
    shell:
        "samtools index -b -@ {threads} {input.final_bam}"
