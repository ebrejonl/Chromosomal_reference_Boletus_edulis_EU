
rule prep_juicer:
    input:
        ref="Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta",
        ref_fasta_index="Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta.fai"
    output:
        chrom_sizes="Results/juicer/chrom.sizes"
    container: c_popgen
    threads: 17
    shell:
        """
        bwa index {input.ref}
        cut -f1,2 {input.ref_fasta_index} > {output.chrom_sizes}
        """


rule juicer:
    input:
        fq1 = "fastq/hi_c_R1.fastq.gz",
        fq2 = "fastq/hi_c_R2.fastq.gz",
        chrom_sizes="Results/juicer/chrom.sizes",
        ref="Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta"
    params:
        dir="/home/etiennebrejon/Desktop/work/Dropbox/Reference_genome_21_02_2024/14_11/Chromosomal_reference_Boletus_edulis_EU",
        genome="BolEdBiel_h2"
    output:
        bam = "Results/juicer/contact_map.hic"
    container: c_popgen
    threads: 12
    shell:
        """
        export _JAVA_OPTIONS="-Xmx28000m -Xms28000m"
        /opt/juicer/CPU/juicer_small.sh -z {input.ref} -p {input.chrom_sizes} -d {params.dir} -g {params.genome} \
        -D /opt/juicer
        mv {params.dir}/contact_map.hic {output.bam}
        """





### add aligner 
#rule align_hic_reads:
#    input:
#        fq1 = "Data/HiC/hi_c_1.fq.gz",
#        fq2 = "Data/HiC/hi_c_2.fq.gz",
#        ref = "Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta"
#    output:
#        sam = "Results/juicer/aligned_reads.sam",
#        bam = "Results/juicer/aligned_reads.bam"
#    container: c_geno
#    threads: 17
#    shell:
#        """
#        # index reference genome
#        bwa-mem2 index {input.ref}
#        # Align Hi-C reads using BWA-MEM2 and output SAM format
#        bwa-mem2 mem -t {threads} \
#        -A 1 -B 4 -E 50 -L 0 -t 15 -M -R "@RG\\tID:HiC_Align\\tSM:HiC_Sample" \
#        {input.ref} {input.fq1} {input.fq2} > {output.sam}
#
#        # Convert SAM to BAM format using samtools
#        samtools view -bS {output.sam} > {output.bam}
#        """
#
#rule sort_and_index_bam:
#    input:
#        bam="Results/juicer/aligned_reads.bam"
#    output:
#        sorted_bam="Results/juicer/aligned_reads.sorted.bam",
#        bam_index="Results/juicer/aligned_reads.sorted.bam.bai"
#    container: c_geno
#    shell:
#        """
#        # Sort the BAM file
#        samtools sort -o {output.sorted_bam} {input.bam}
#
#        # Index the sorted BAM file
#        samtools index {output.sorted_bam}
#        """
#
#rule generate_pairs:
#    input:
#        sorted_bam="Results/juicer/aligned_reads.sorted.bam",
#        chrom_sizes="Results/juicer/chrom.sizes"
#    output:
#        pairs="Results/juicer/aligned_reads.pairs",
#        dedup_pairs="Results/juicer/aligned_reads.dedup.pairs"
#    container: c_popgen
#    shell:
#        """
#        touch {input.chrom_sizes}
#        # Parse BAM into .pairs file using the chrom.sizes file as genome
#        pairtools parse -c --chroms-path Results/juicer/chrom.sizes --output {output.pairs} {input.sorted_bam}
#
#        # Sort the .pairs file
#        pairtools header generate {output.pairs}
#        pairtools sort --output {output.pairs} {output.pairs}
#
#        # Deduplicate the .pairs file
#        pairtools dedup --output {output.dedup_pairs} {output.pairs}
#        """
#
#rule Hi_c_map:
#    input:
#        pairs="Results/juicer/aligned_reads.dedup.pairs",
#        chrom_sizes="Results/juicer/chrom.sizes"
#    output:
#        contact_map="Results/juicer/contact_map.hic"
#    params:
#        output_dir="Results/juicer/map",  # Output directory
#    container: c_popgen
#    threads: 17
#    shell:
#        """
#        # Run Juicer Tools to generate the contact map from .pairs file
#        mkdir -p {params.output_dir}
#        java -jar /opt/juicer/scripts/common/juicer_tools.1.9.9_jcuda.0.8.jar pre \
#            -p {input.pairs} {params.output_dir}/contact_map.hic \
#            -g {input.chrom_sizes}
#        """
#