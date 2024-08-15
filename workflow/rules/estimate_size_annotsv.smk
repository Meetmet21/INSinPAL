"""
:Author: Mehmet Sehir
:Date: 15.05.2024

Another purpose of this pipeline, is to estimate the inserted sequence length through the Source annotation done before.
For this, via the breakpoint, Source field in INS BED merged and sample BAM, the aim is to extract discordant reads
mapping the edges of the inserted sequence in a given breakpoint. Then, by subtracting the reverse and forward reads
positions, we can get an idea of the putative length. This feature has been tested in simulated insertions in WGS data
and worked well. However, in real life data, sometimes very large insertion site can be proposed, the user has to check
by itself the given event. More information about the algorithm in find_insertion_length.py

finally, palindromes information is included in INS BED merged file to have more informaiton about the insertion event
contect. Then, annotation of each event is done via AnnotSV and an excel file is generated for the filtering step.

TODO: Add frequencies via the local database of Insertion events.
"""

localrules: find_insertion_size, set_user_defined_header, add_header_to_bed, AnnotSV, final_xml_file

rule find_insertion_size:
    """
    As described before, this rule try to extrapolate the insrted sequence length in the Insertion breakpoint via
    discordant reads. The algorithm process as follow: 
        - Check if a uniq source of insertion
        - Extract dscrdant reads on both side of the breakpoint
        - Separate discordant reads mate depending on if reverse or forward
        - Check if there is a consensus distribution of positions
        (This part is done by computing the standard-deviation of positions, modifs threshold in params in necessary)
        - If enough supporting signals, compute the putative length
    Note that sometimes, the length is overestimated because of lack of strong signals. But on simulated insertions in
    WGS data, it worked well.
    """
    input:
        bed="{sample}.merged.modup.bed",
        bam=join(paths.results, "{sample}/{sample}.bam"),
        bai=join(paths.results,"{sample}/{sample}.bam.bai")
    output:
        temp("{sample}.merged.nodup.size.bed")
    params:
        # Window range on left and right of a breakpoint in which seek split reads.
        window = 300,
        # Read QMAP threshold to be considered
        mapq_cutoff = 10,
        # Allowed standard deviation between positions
        # This parameter allows to select set of position containing very divergent split
        # read positions
        std_cutoff = 500
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "find_insertion_length.py")

# Add palindrome annotation
rule add_palindrome_annotation:
    """
    Add breakpoint overlapping palindrome informations for a better overview of the event.
    The final output BED file contians the following fields:
        - Chr name  
        - Start   
        - Stop    
        - SVtype  
        - SVlen   
        - MEItype 
        - Source  
        - Paltype: Perfect, Near or Spacer.
        - Pallen: length of the palindrome with the spacer if there is.
        - Splen: the spacer length between both arms.
        - Missrate%: The dissimilarity percentage in sequence betweem both arms.
        - AT%: AT content in the sequence.
        - Recomb_score: the d score for instability of a palindrome.
    """
    input:
        tmp=rules.find_insertion_size.output,
        pal_annot=join(paths.resources, "data/palindromes/Recombinogenic_palindromes_bysize.bed")
    output:
        join(paths.results, "{sample}/{sample}_full_annot.bed")
    conda:
        join(paths.envs, "mapping.yaml")
    params:
        window=10
    shell:
        "bedtools window -w {params.window} -a {input.tmp} -b {input.pal_annot} | "
        "cut -f 1-7,11-16 | "
        "bedtools merge -c 4,5,6,7,8,9,10,11,12,13 -o last > {output}"

rule set_user_defined_header:
    """
    Get BED header names for AnnotSV user defined fields.
    """
    output:
        temp("header.bed")
    params:
        header_fields="#Chrom\tStart\tStop\tSV_type\tSV_len\tMEI_type\tSource_INS"
                      "\tPal_type\tPal_len\tSp_len\tMiss%\tAT%\tRecomb_score\n"
    run:
        with open(output[0], "w") as header:
            header.write(params.header_fields)

rule add_header_to_bed:
    """
    Add header to bed for AnnotSV to keep user defined columns in bed.
    """
    input:
        bed=rules.add_palindrome_annotation.output,
        header=rules.set_user_defined_header.output
    output:
        temp("{sample}.header.bed")
    shell:
        "cat {input.header} {input.bed} > {output}"

rule AnnotSV:
    """
    AnnotSV call to annotate INS breakpoints. See config file of progs for parameters (in Progs folder).
    """
    input:
        rules.add_header_to_bed.output
    output:
        temp("{sample}.annotsv.tsv")
    params:
        annotsv = caller.AnnotSV
    shell:
        "{params.annotsv} -SVinputFile {input} -outputFile {output} -outputDir {wildcards.sample} -svtBEDcol 4 "
        "&& mv {wildcards.sample}/{output} . "
        "&& rm -rf {wildcards.sample}/"

rule clean_working_dir:
    """
    Remove callers working directory to free space.
    """
    input:
        rules.AnnotSV.output
    output:
        temp("cleaned_{sample}.flag")
    params:
        insurveyor_dir=join(paths.results, "{sample}/insurveyor"),
        basil_dir=join(paths.results, "{sample}/basil"),
        manta_dir=join(paths.results, "{sample}/manta"),
        scramble_dir=join(paths.results, "{sample}/scramble")
    shell:
        "rm -rf {params} && touch {output}"

rule final_xml_file:
    """
    From AnnotSV bed or tsv file to Excel file for better visualization.
    """
    input:
        rules.AnnotSV.output,
        rules.clean_working_dir.output
    output:
        join(paths.results, "{sample}/{sample}.xlsx"),
        touch("/data/DMGL/SMG-mesehir-data/bam_to_run/{sample}.done")
    conda:
        join(paths.envs, "python_base.yaml")
    script:
        join(paths.scripts, "from_bed_to_excel.py")
