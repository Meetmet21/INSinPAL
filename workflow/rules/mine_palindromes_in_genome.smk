"""
:Author: Mehmet Sehir
:Date: 15.05.2024

From the hgt19 reference genome, extract all palindromic sequences considering a few conditions:
    - The spacer sequence in between both arms of a palindrome, should not be greater than 20 bp.
    - The identity of sequence (reverse complementary of one arm compared to the other), or more precisely, the
    divergence between them, should not exceed 15%.

Once found, the sequence is annotated with different stats and its start and end positions relative to the reference genome
is extracted (0 based). Then, a few step of filtering are performed to keep only unstable ones based on criteria found in
 the literature: "Long inverted repeats in eukaryotic genomes: Recombinogenic motifs determine genomic plasticity"
"""

# MODULES
import sys
from os.path import join

# PARAMETER FILE
sys.path.append("config/")
import parameters
# Path to pipeline main directory directories
paths = parameters.WorkFlowPaths()
data = parameters.Data()

localrules: select_size_pal_for_SV_filtering

rule mine_palindromes_by_chr:
    """
    This rule need to have the reference genome hg19 sequence in fasta format per chromosome.
    For each chromosome, it will mine palindromes within the genome: see the documentation in the python script for more
    details about the algorithm. But basically, the idea is to find palindrome seeds of at least 8bp (4 each side). Then 
    extend these seeds until one of the criteria within parameters section is not respected. Then, compute different
    stats to annotate each palindrome for downstream analysis.
    
    Feel free to modify the parameters, the ones here are based on the literature information about palindromes with
    potential capacity to form secondary structure and undergo SV events. 
    
    Also, note that, this step can take a lot of time, so once computed, keep the output file as a reference DB for 
    palindromic regions in hg19.
    """
    input:
        chr_fa=lambda wildcards: data.hg19_fa_chromosomes.get(wildcards.chr)
    output:
        temp("palindromes_{chr}.bed")
    params:
        # Min palindrome length to start mining
        seed_pal_len=4,
        # Miss rate between each arms of a palindrome
        miss_rate=0.15,
        # Maximum spacer length between each arm of a palindrome
        max_sp=20
    threads: 1
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "palindrome_mining_by_chr.py")

rule filter_chr_for_unstable_pals:
    """
    Filtering step based on the "d" score from the paper: 
    "Long inverted repeats in eukaryotic genomes: Recombinogenic motifs determine genomic plasticity"
    Only putative unstable palindromes remains.
    """
    input:
        chr_bed=rules.mine_palindromes_by_chr.output
    output:
        temp("tmp_dir/filtered_{chr}.bed")
    threads: 1
    conda:
        join(paths.envs, "r_base.yaml")
    script:
        join(paths.scripts, "filter_unstable_pal_by_chr.R")

rule merge_pal_chr_together:
    """
    Merge all chromosomes together.
    """
    input:
        expand("tmp_dir/filtered_{chr}.bed", chr=data.chrom)
    output:
        join(paths.resources, "data/palindromes/filtered_hg19_genome_palindromes.bed"),
    shell:
        "cat {input} >> {output}; "
        "rm -rf tmp_dir"

rule select_size_pal_for_SV_filtering:
    """
    For downstream usage, select a subset of palindromes based on their size, as longer ones are more prompt to form
    second structure in vivo.
    Here, grep is used to subset the initial dataset. To add more size group modify the parameter as follow:
    size_to_keep="group_size_level\|another_group_size_level"
    """
    input:
        rules.merge_pal_chr_together.output
    output:
        join(paths.resources, "data/palindromes/Recombinogenic_palindromes_bysize.bed")
    params:
        # Size factor levels: '0-50 bp', '51-99 bp', '100-200 bp', '>200 bp'
        size_to_keep="100-200 bp\|>200 bp"
    shell:
        "grep '{params.size_to_keep}' {input} > {output}"


