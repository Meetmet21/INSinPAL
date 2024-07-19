"""
:Author: Mehmet Sehir
:Date: 15.05.2024

This is the merging step to constitute a Metacaller for large Insertions in WGS data within unstable palindromic regions.
The merge is done by concating all three INS BED files (Three callers) and removing breakpoints that overlaps in a window of 50 bp.
"""

localrules: concat_SV_bed, merge_close_breakpoints

rule concat_SV_bed:
    """
    Simple concatenation of each caller INS BED file.
    fields:
        - Chr
        - Start
        - Stop
        - SV type
        - SV size
        - MEI type
        - Source
    """
    input:
        expand(join(paths.results, "{{sample}}/{{sample}}.{progs}.Pal.MEIS.source.bed"), progs=caller.names)
    output:
        temp("{sample}.merged.bed")
    shell:
        "cat {input} >> tmp; "
        "sort -V tmp > {output}; "
        "rm tmp"

rule merge_close_breakpoints:
    """
    To avoid duplicate events, after concatenating files, merge close breakpoints.
    """
    input:
        rules.concat_SV_bed.output
    output:
        temp("{sample}.merged.modup.bed")
    conda:
        join(paths.envs, "python_base.yaml")
    params:
        # Window in which breakpoints are merged
        window = 50
    script:
        join(paths.scripts, "remove_duplicates_from_merged.py")
