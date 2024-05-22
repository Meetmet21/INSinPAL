"""
:Author: Mehmet Sehir
:Date: 15.05.2024

Few rules to change a VCF into a INS BED file, where only SV INS are kept with a size of 50 bp at least. Then, a few
annotations are done by intersecting different annotation files with INS BED file which was filtred for breakpoints
overlapping with putative unstable palindromic regions.
    - MEIs annotation via the SCRAMble results
    - Putative source of insertion in ref via the extraction of discordant reads around beakpoint
"""

# MODULES
import sys
from os.path import join

# PARAMETER FILE
sys.path.append("config/")
import parameters
# Path to pipeline main directory directories
paths = parameters.WorkFlowPaths()

localrules: meis_to_bed, vcf_to_bed, annotate_MEIS_part1, annotate_MEIS_part2, INS_in_pal_regions, annotate_source_INS

rule meis_to_bed:
    """
    Extract from SCRAMble final result only interesting fields which are in this case:
        - Chr name
        - Start
        - Stop
        - MEI type
    """
    input:
        meis=join(paths.results, "{sample}/scramble/{sample}_MEIs.txt")
    output:
        bed=join(paths.results, "{sample}/scramble/{sample}.MEIS.bed")
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "meis_to_bed.py")

rule vcf_to_bed:
    """
    Extract from callers final results (VCF) interesting Insertion events.
    The final output format is a BED with the following fields:
        - Chr name
        - Start (1 based) 
        - Stop
        - SV type (Onlx INS)
        - Size (From callers)
    """
    input:
        join(paths.results, "{sample}/{progs}/{sample}.{progs}.vcf")
    output:
        temp("{sample}.{progs}.bed")
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "vcf_to_bed.py")

# Select INS breakpoints overlapping palindromic regions.
rule INS_in_pal_regions:
    """
    Keep only Insertions overlapping with palindromes considered as unstable (see mine_palindromes_in_genome.smk for
    features to consider a palindrome as unstable). This rule is the main source of selection for INS events, so new
    features can be selected to have more stringent parameters for Recombinogenic_palindromes_bysize.bed to reduce or
    increase the number of mappable sites. 
    """
    input:
        pal=join(paths.resources, "Recombinogenic_palindromes_bysize.bed"),
        ins_bed=rules.vcf_to_bed.output
    output:
        temp("{sample}.{progs}.Pal.bed")
    conda:
        join(paths.envs, "mapping.yaml")
    params:
        # From this range, considered as overlapping breakpoints.
        window=10
    shell:
        "bedtools window -w {params.window} -u -a {input.ins_bed} -b {input.pal} > {output}"

rule annotate_MEIS_part1:
    """
    Extract MEIs matching INS breakpoints in a new file.
    """
    input:
        bed=rules.INS_in_pal_regions.output,
        MEIS=rules.meis_to_bed.output
    output:
        annotated_bed=temp("{sample}.{progs}.Pal.bed.tmp")
    params:
        # From this range, considered as overlapping breakpoints.
        window=10
    conda:
        join(paths.envs, "mapping.yaml")
    shell:
        "bedtools window -w {params.window} -a {input.bed} -b {input.MEIS} | "
        "cut -f 1-3,9 > {output.annotated_bed}"

# Merge MEIs positions with INS calls.
rule annotate_MEIS_part2:
    """
    Annotate INS BED with the MEI type information, so the final result looks like:
        - Chr
        - Start
        - Stop
        - SV type
        - Size
        - MEI type
    """
    input:
        bed_w_MEIS=rules.annotate_MEIS_part1.output.annotated_bed,
        bed_SV=rules.INS_in_pal_regions.output
    output:
        final_bed=temp("{sample}.{progs}.Pal.MEIS.bed")
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "annotate_MEIS.py")

rule annotate_source_INS:
    """
    To understand where the putative Inserted sequence comes from within the genome, extract the discordant reads around
    the INS breakpoint. Here he window parameter is fixed to 300bp to capture all discordant reads realted to the given
    breakpoint. This is because the maximum length of the insert + one side read length is 300bp in average. 
    
    The extraction from BAM file follows a few criteria (see script documentation) including well mapped discordant
    reads. Feel free to modify the cutoff for MAPQ of discordant reads if you want less stringent conditions. 
    (10 -> 90% mapping quality)
    
    Once the discordant reads are extracted, the algorithm checks for the mate REF information and stores it in a data
    structure. Then, the most frequent REF are kept depnding on if their frequency reaches the halp of the most frequent
     REF (can be modified to higher cutoff through the threshold parameter). And a new field is added to INS BED:
        - Previous ones
        - Source
    
    The main idea is to capture signals of uniq Insertion sources, which can be related to something else than MEIs and 
    potential can be linked to large Insertion events.
    """
    input:
        bed=rules.annotate_MEIS_part2.output.final_bed,
        bam=join(paths.results, "{sample}/{sample}.bam"),
        bai=join(paths.results,"{sample}/{sample}.bam.bai")
    output:
        join(paths.results, "{sample}/{sample}.{progs}.Pal.MEIS.source.bed")
    params:
        # Window range on left and right of a breakpoint in which seek discordant reads.
        # Max = insert size + (left or right read length)
        window=300,
        # Frequency fraction of the most frequent source to reach to be selected as source.
        threshold=0.5,
        # Read QMAP threshold.
        cutoff=10
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "get_source_INS.py")
