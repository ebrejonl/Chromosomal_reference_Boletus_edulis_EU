
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
        /opt/juicer/CPU/juicer.sh -z {input.ref} -p {input.chrom_sizes} -d {params.dir} -g {params.genome} \
        -D /opt/juicer
       # mv {params.dir}/contact_map.hic {output.bam}
        """
