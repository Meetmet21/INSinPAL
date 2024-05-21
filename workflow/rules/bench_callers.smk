# MODULES
import sys
from os.path import join


# PARAMETER FILE
sys.path.append("config/")
import parameters
# Path to pipeline main directory directories
paths = parameters.WorkFlowPaths()

# local rules
localrules: trueset_vcf_to_bed, TP_within_BED_SV, compute_bench_metrics

rule trueset_vcf_to_bed:
    input:
        config["trueset"]["SIM_INS_SET"]
    output:
        expand(join(paths.resources, "{trueset}.bed"), trueset=config["trueset"])
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "vcf_to_bed.py")

rule TP_within_BED_SV:
    input:
        trueset=rules.trueset_vcf_to_bed.output,
        sv_bed="{sample}.{progs}.MEIS.source.bed"
    output:
        filtred_bed="{sample}.{progs}.MEIS.source.TP.bed"
    params:
        window = 5
    conda:
        join(paths.envs, "mapping.yaml")
    shell:
        "bedtools window -w {params.window} -u -a {input.sv_bed} -b {input.trueset} > {output.filtred_bed}"

rule compute_bench_metrics:
    input:
        tp_bed=rules.TP_within_BED_SV.output.filtred_bed,
        bench_bed=rules.trueset_vcf_to_bed.output
    output:
        summary=join(paths.results, "{sample}/bench/{progs}/summary.txt")
    params:
        caller="{progs}"
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "compute_bench_metrics.py")

