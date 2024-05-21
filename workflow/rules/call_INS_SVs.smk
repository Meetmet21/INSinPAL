"""
:Author: Mehmet Sehir
:Date: 14.05.2024

From sample in config file, call SVs with 3 software: INSurVeyor, Manta and Basil.
Additionally, screen for MEIs in sample with SCRAMble by splitting per chromosome to enable parallelization.
"""

# MODULES
import sys
from os.path import join, realpath, dirname

# PARAMETER FILE
sys.path.append("config/")
import parameters
# Path to pipeline main directory directories
paths = parameters.WorkFlowPaths()
# PAth to callers
caller = parameters.Progs()
# Path to reference genome
data = parameters.Data()

configfile: join(paths.config, "config.yaml")

# local rules
localrules: call_SVs_manta_configure, filter_basil_vcf, split_clusters_by_chr, get_scramble_header, merge_scramble_calls

rule call_INS_insurveyor:
    """
    Call Insertions with INSurVeyor once the sample BAM is formatted with mate score tags. This rule create an inserveyor
    directory within the result/{sample_id} directory, but only temporarly. It will be deleted at the end of the pipeline
    to save space. Finally, the caller generates a compressed bcf file, but downstream analysis need uncompressed vcf.
    """
    input:
        bam=join(paths.results, "{sample}/{sample}.bam"),
        bai=join(paths.results, "{sample}/{sample}.bam.bai"),
        ref=data.hg19_fa
    output:
        vcf=join(paths.results, "{sample}/insurveyor/{sample}.insurveyor.vcf")
    threads: 12
    log:
        join(paths.workflow,"logs/{sample}/call_INS_insurveyor.log")
    params:
        working_dir = join(paths.results,"{sample}/insurveyor"),
        insurveyor_sif=caller.insurveyor,
        min_INS_size=50,
        path_ref=lambda wildcards, input: dirname(input.ref),
        oath_bam=lambda wildcards, input: dirname(input.bam)
    conda:
        join(paths.envs, "mapping.yaml")
    shell:
        "if [[ -d ${params.working_dir} ]]; then rm -rf {params.working_dir}; fi; "
        "mkdir -p {params.working_dir}; "
        "singularity run -B {params.working_dir} -B {params.path_ref} -B {params.oath_bam} "
        "{params.insurveyor_sif} "
        "--threads {threads} "
        "--min-insertion-size {params.min_INS_size} "
        "{input.bam} {params.working_dir} {input.ref} | "
        "tee {log} 2>&1; "
        "bcftools view -O v -o {output.vcf} {params.working_dir}/out.pass.vcf.gz"

# Extract path bam file
def get_path_sample_from_config(wildcards):
    return config["samples"][wildcards.sample]

# Call SVs with manta (To be filtred after)
# First set the configuration directory and files
rule call_SVs_manta_configure:
    """
    First step for Manta caller consisting in generating a config file for the workflow. Here, it's necessary to give 
    the path to a python interpreter version 2.x. Same for INSurVeyor, the manta directory is a tmp.
    """
    input:
        bam = get_path_sample_from_config,
        bai = lambda wildcards: get_path_sample_from_config(wildcards) + ".bai",
        ref = data.hg19_fa
    output:
        # No direct outptu from caller
        workflow_file=join(paths.results, "{sample}/manta", "runWorkflow.py")
    log:
        join(paths.workflow,"logs/{sample}/call_SVs_manta_configure.log")
    params:
        working_dir = lambda wildcards, output: dirname(output.workflow_file),
        manta_config_file=caller.manta_config,
        # Manta needs python 2.*
        python2_7= caller.python27
    shell:
        "if [[ -d {params.working_dir} ]]; then rm -rf {params.working_dir}; fi; "
        "mkdir -p {params.working_dir}; " 
        "{params.python2_7} {params.manta_config_file} "
        "--bam {input.bam} "
        "--referenceFasta {input.ref} "
        "--runDir {params.working_dir} | "
        "tee {log} 2>&1"

# Run manta workflow to call SVs
rule call_SVs_manta_workflow:
    """
    Second step for Manta caller using the previously generated config file. For downstream analysis, the vcf is 
    uncompressed.
    """
    input:
        workflow_script=rules.call_SVs_manta_configure.output.workflow_file
    output:
        vcf=join(paths.results, "{sample}/manta/{sample}.manta.vcf"),
    threads: 12
    log:
        join(paths.workflow,"logs/{sample}/call_SVs_manta_workflow.log")
    params:
        # Need python 2
        python2_7=caller.python27,
        working_dir=lambda wildcards, output: dirname(output.vcf)
    conda:
        join(paths.envs, "mapping.yaml")
    shell:
        "{params.python2_7} {input.workflow_script} -j {threads} | tee {log} 2>&1; "
        "bcftools view -O v -o {output.vcf} {params.working_dir}/results/variants/diploidSV.vcf.gz"

# Call SVs with Basil
rule call_SVs_basil:
    """
    First step for Basil caller. As before, the resulting directoires are tmp. No need to uncompress VCF file.
    """
    input:
        bam = get_path_sample_from_config,
        bai = lambda wildcards: get_path_sample_from_config(wildcards) + ".bai",
        ref = data.hg19_fa
    output:
        vcf=temp(join(paths.results, "{sample}/basil/{sample}_tmp.vcf"))
    params:
        basil_sif=caller.basil,
        output_dir=join(paths.results, "{sample}/basil"),
        path_ref= lambda wildcards,input: dirname(input.ref),
        oath_bam= lambda wildcards,input: dirname(input.bam)
    log:
        join(paths.workflow,"logs/{sample}/call_SVs_basil.log")
    threads: 12
    shell:
        "if [[ -d {params.output_dir} ]]; then rm -rf {params.output_dir}; fi; "
        "mkdir -p {params.output_dir}; "
        "singularity exec -B {params.output_dir} -B {params.path_ref} -B {params.oath_bam} "
        "{params.basil_sif} basil "
        "--realignment-num-threads {threads} "
        "-ir {input.ref} "
        "-im {input.bam} "
        "-ov {output.vcf} | "
        "tee {log} 2>&1"

# Filter the output vcf from basil caller.
rule filter_basil_vcf:
    """
    Second step consists to filter calls with not enough supporting signals.
    --min-oea-sum MIN_OEA_SUM 5
                        Minimal total OEA (One end anchored) coverage.
    --min-clipping-each-side MIN_CLIPPING_EACH 5
                        Minimal OEA coverage on each side.
    """
    input:
        vcf=rules.call_SVs_basil.output.vcf
    output:
        vcf=join(paths.results, "{sample}/basil/{sample}.basil.vcf")
    params:
        in_vcf_path=rules.call_SVs_basil.params.output_dir,
        filter_script=caller.filtering_basil,
        basil_sif=caller.basil,
        min_oea_sum=5,
        min_clipping_sum=5
    log:
        join(paths.workflow,"logs/{sample}/filter_basil.log")
    shell:
        "singularity exec -B {params.in_vcf_path} "
        "{params.basil_sif} python3 {params.filter_script} "
        "-i {input.vcf} "
        "-o {output.vcf} "
        "--min-oea-sum {params.min_oea_sum} "
        "--min-clipping-sum {params.min_clipping_sum} | "
        "tee {log} 2>&1"

# Find soft clipped clusters in BAM file.
rule scramble_identify_cluster:
    """
    First step of calling MEIs in sample BAM. Here, we increased the singal supporting for clusters to 10, in order to
    reduce execution time.
    """
    input:
        bam = get_path_sample_from_config,
        bai = lambda wildcards: get_path_sample_from_config(wildcards) + ".bai",
    output:
        clusters=join(paths.results, "{sample}/scramble/{sample}.clusters.txt")
    threads: 1
    params:
        scramble_sif=caller.scramble,
        working_dir=lambda wildcards, output: dirname(output.clusters),
        path_bam=lambda wildcards,input: dirname(input.bam),
        min_soft_clipped_bases=10,
        min_soft_clipped_reads=10
    shell:
        "if [[ -d {params.working_dir} ]]; then rm -rf {params.working_dir}; fi; "
        "mkdir -p {params.working_dir}; "
        "singularity exec -B {params.path_bam} -B {params.working_dir} "
        "{params.scramble_sif} cluster_identifier "
        "-m {params.min_soft_clipped_bases} "
        "-s {params.min_soft_clipped_reads} "
        "{input.bam} > {output.clusters} "

rule split_clusters_by_chr:
    """
    SCRAMble has no multi-threading. So split the cluster file by chromosome to enable parallelization, by launching the
    analysis for each chromosome at the same time and merging.
    """
    input:
        clusters=rules.scramble_identify_cluster.output.clusters
    output:
        out=temp("{sample}.clusters.{chr}.txt")
    params:
        chr=lambda wildcards, output: output.out.split(".")[2]
    shell:
        "grep -E '^{params.chr}:' {input.clusters} > {output.out}"


rule scramble_cluster_analysis:
    """
    Launch SCRAMble analysis for each chromosome on the cluster.
    """
    input:
        clusters="{sample}.clusters.{chr}.txt",
        ref=data.hg19_fa,
        # Reference MEI sequences db.
        mei_refs=data.meis_ref
    output:
        outfile=temp("{sample}_{chr}_MEIs.txt")
    threads: 1
    params:
        out_name=lambda wildcards, output: realpath(output.outfile).replace("_MEIs.txt", ""),
        scramble_sif=caller.scramble,
        path_ref=lambda wildcards, input: dirname(input.ref),
        path_meis_ref=lambda wildcards, input: dirname(input.mei_refs),
        path_cluster=lambda wildcards, input: realpath(input.clusters),
        cluster_analysis=caller.cluster_analysis
    shell:
        "singularity exec -B {params.path_ref} -B {params.path_meis_ref} -B {params.path_cluster} "
        "{params.scramble_sif} "
        "Rscript --vanilla {params.cluster_analysis} "
        "--out-name {params.out_name} "
        "--cluster-file {params.path_cluster} "
        "--install-dir /app/cluster_analysis/bin/ "
        "--mei-refs {input.mei_refs} "
        "--ref {input.ref} "
        "--eval-meis "
        "--no-vcf"


rule get_scramble_header:
    """
    Extract SCRAMble header of one file for futur merging.
    """
    input:
        "{sample}_1_MEIs.txt"
    output:
        temp("{sample}_scramble_header.txt")
    shell:
        "grep 'Insertion' {input} > {output}; "


rule merge_scramble_calls:
    """
    Merge all chromosomes calls in one file with header.
    """
    input:
        meis=expand("{{sample}}_{chr}_MEIs.txt", chr=data.chrom),
        header=rules.get_scramble_header.output
    output:
        join(paths.results, "{sample}/scramble/{sample}_MEIs.txt")
    shell:
        "cat {input.header} > {output}; "
        "for file in {input.meis}; do grep -v 'Insertion' $file >> tmp; done; "
        "sort -V tmp >> {output}; "
        "rm tmp"

