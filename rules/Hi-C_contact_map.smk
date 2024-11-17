
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
        # Align Hi-C reads using BWA-MEM2 and output SAM format
        bwa-mem2 mem -t {threads} \
        -M -R "@RG\\tID:HiC_Align\\tSM:HiC_Sample" \
        {input.ref} {input.fq1} {input.fq2} > {output.sam}

        # Convert SAM to BAM format using samtools
        samtools view -bS {output.sam} > {output.bam}
        """




rule Hi_c_map:
    input: 
        fq1= "Data/HiC/hi_c_1.fq.gz",
        fq2= "Data/HiC/hi_c_2.fq.gz",
        ref= "Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta",
        chrom_sizes="Results/juicer/chrom.sizes"
    output: 
        contact_map="Results/juicer/contact_map.hic"
    params:
        output_dir="Results/juicer",  # Path to output directory
        genome_id="BolEdBiel_h2",
        threads=15
    container: c_popgen
    shell:
        """
            java -jar /opt/juicer/scripts/common/juicer_tools.1.9.9_jcuda.0.8.jar pre \
            {input.fq1} {input.fq2} {output.contact_map} \
            -t {params.threads} \
            -g {params.genome_id} \
            -p {input.chrom_sizes}
        """