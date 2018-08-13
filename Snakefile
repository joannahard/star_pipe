SAMPLES, = glob_wildcards("data/input/{smp}_R1_001.fastq.gz")


rule all:
    input:
        expand("data/output/{smp}/trimmed/", smp=SAMPLES),
        expand("data/output/{smp}/mapped/", smp=SAMPLES),
        expand("data/output/{smp}/fastqc_untrimmed/", smp=SAMPLES),
        expand("data/output/{smp}/fastqc_trimmed/", smp=SAMPLES),
        expand("data/output/{smp}/mapped/Aligned.sortedByCoord.out_sorted.bam", smp=SAMPLES),
        expand("data/output/{smp}/mapped/Aligned.sortedByCoord.out_sorted.bam.bai", smp=SAMPLES)



rule trimming:
    input:
        fwd_raw = "data/input/{smp}_R1_001.fastq.gz",
        rev_raw = "data/input/{smp}_R2_001.fastq.gz"
    output:
        fwd_trim = "data/output/{smp}/trimmed/{smp}_R1_001_val_1.fq.gz",
        rev_trim = "data/output/{smp}/trimmed/{smp}_R2_001_val_2.fq.gz",
        dir = "data/output/{smp}/trimmed/"
    shell:
        "trim_galore --phred33 --quality 25 --length 20 --paired {input.fwd_raw} {input.rev_raw} --output_dir {output.dir}"


rule mapping:
    input:
        fwd_trim = "data/output/{smp}/trimmed/{smp}_R1_001_val_1.fq.gz",
        rev_trim = "data/output/{smp}/trimmed/{smp}_R2_001_val_2.fq.gz",
        gtf = config["gtf"],
        ref = config ["ref_dir"]
    output:
        "data/output/{smp}/mapped/"
    shell:
        "STAR --genomeDir {input.ref} --readFilesCommand zcat --readFilesIn {input.fwd_trim} {input.rev_trim} --runThreadN 5 --outFileNamePrefix {output} --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {input.gtf}"


rule samtools_sort:
    input:
        "data/output/{smp}/mapped/Aligned.sortedByCoord.out.bam"
    output:
        "data/output/{smp}/mapped/Aligned.sortedByCoord.out_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"


rule samtools_index:
    input:
        "data/output/{smp}/mapped/Aligned.sortedByCoord.out_sorted.bam"    
    output:
        "data/output/{smp}/mapped/Aligned.sortedByCoord.out_sorted.bam.bai"
    shell:
        "samtools index {input}"



rule fastqc_untrimmed:
    input:
        fwd_raw = "data/input/{smp}_R1_001.fastq.gz",
        rev_raw = "data/input/{smp}_R2_001.fastq.gz"
    output:
        "data/output/{smp}/fastqc_untrimmed/"
    shell:
        "fastqc --quiet --outdir {output} --extract  -f fastq {input.fwd_raw} {input.rev_raw}"



rule fastqc_trimmed:
    input:
        fwd_trim = "data/output/{smp}/trimmed/{smp}_R1_001_val_1.fq.gz",
        rev_trim = "data/output/{smp}/trimmed/{smp}_R2_001_val_2.fq.gz"
    output:
        "data/output/{smp}/fastqc_trimmed/"
    shell:
        "fastqc --quiet --outdir {output} --extract  -f fastq {input.fwd_trim} {input.rev_trim}"






#rule qualimap:
#    input: 
#        "data/output/{smp}/mapped/Aligned.sortedByCoord.out.bam"
#    output:
#        "data/output/{smp}/qualimap/"
#      #report = "{filename}.qualimap/qualimapReport.html",
#      #report = "data/output/{smp}/qualimap/qualimapReport.html"
#      #gr = "{filename}.qualimap/genome_results.txt",
#      #ish = "{filename}.qualimap/raw_data_qualimapReport/insert_size_histogram.txt",
#      #ch = "{filename}.qualimap/raw_data_qualimapReport/coverage_histogram.txt",
#      #gc = "{filename}.qualimap/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"
#    #log: "{filename}.qualimap/qualimap.log"
#    params: "-sd -sdmode 0 --java-mem-size=20G -c -nw 400"
#    threads: 10
#    shell: "qualimap bamqc -nt {threads} {params} -bam {input} -outdir {output}"




# includes multimappers use htseq instead because it is on conda. or is countfeatures in qualimap. ask bioinfo support
#rule featurecounts:
#    input:
#        bam = "data/output/{smp}/mapped/Aligned.sortedByCoord.out.bam",
#        gtf = config["gtf"]
#    output:
#        "data/output/{smp}/featurecounts"
#    shell:
#        featureCounts -O -M -Q 30 -p -a {input.gtf} -o {output} {input.bam}
