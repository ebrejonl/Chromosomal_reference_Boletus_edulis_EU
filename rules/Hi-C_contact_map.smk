
rule prep_juicer:
    input:
        ref_fasta_index="Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta.fai"
    output:
        chrom_sizes="Results/juicer/chrom.sizes"
    shell:
        """
        cut -f1,2 {input.ref_fasta_index} > {output.chrom_sizes}
        """

## add aligner 
rule align_hic_reads:
    input:
        fq1 = "Data/HiC/hi_c_1.fq.gz",
        fq2 = "Data/HiC/hi_c_2.fq.gz",
        ref = "Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta"
    output:
        sam = "Results/juicer/aligned_reads.sam",
        bam = "Results/juicer/aligned_reads.bam"
    container: c_geno
    threads: 17
    shell:
        """
        # index reference genome
        bwa-mem2 index {input.ref}
        # Align Hi-C reads using BWA-MEM2 and output SAM format
        bwa-mem2 mem -t {threads} \
        -M -R "@RG\\tID:HiC_Align\\tSM:HiC_Sample" \
        {input.ref} {input.fq1} {input.fq2} > {output.sam}

        # Convert SAM to BAM format using samtools
        samtools view -bS {output.sam} > {output.bam}
        """

rule sort_and_index_bam:
    input:
        bam="Results/juicer/aligned_reads.bam"
    output:
        sorted_bam="Results/juicer/aligned_reads.sorted.bam",
        bam_index="Results/juicer/aligned_reads.sorted.bam.bai"
    container: c_geno
    shell:
        """
        # Sort the BAM file
        samtools sort -o {output.sorted_bam} {input.bam}

        # Index the sorted BAM file
        samtools index {output.sorted_bam}
        """


rule Hi_c_map:
    input:
        sorted_bam="Results/juicer/aligned_reads.sorted.bam",
        chrom_sizes="Results/juicer/chrom.sizes"
    output:
        contact_map="Results/juicer/contact_map.hic"
    params:
        output_dir="Results/juicer/map",  # Output directory
        genome_id="Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta",  # Genome identifier
    container: c_popgen
    threads: 17
    shell:
        """
        # Run Juicer Tools to generate the contact map
        mkdir -p {params.output_dir}
        java -jar /opt/juicer/scripts/common/juicer_tools.1.9.9_jcuda.0.8.jar pre \
            {input.sorted_bam} {params.output_dir}/contact_map.hic \
            -g {params.genome_id} \
            -p {input.chrom_sizes}
        """
