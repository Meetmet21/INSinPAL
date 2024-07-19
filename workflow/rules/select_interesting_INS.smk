"""
:Author: Mehmet Sehir
:Date: 15.05.2024

Few rules to change a VCF into a SV BED file, where only SV INS are kept with a size of 50 bp at least.
"""

localrules: vcf_to_bed, INS_in_pal_regions

rule vcf_to_bed:
    """
    Extract from callers final results (VCF) interesting Insertion events.
    The final output format is a BED with the following fields:
        - Chr name
        - Start (1 based) 
        - Stop
        - SV type (Only INS)
        - Size (If estimated by callers)
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
        pal=join(paths.resources, "data/palindromes/Recombinogenic_palindromes_bysize.bed"),
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
