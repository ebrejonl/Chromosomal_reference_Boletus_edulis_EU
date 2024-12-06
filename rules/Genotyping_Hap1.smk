## Genotyping all EU individuals on Haplotype 1 with the unplaced contigs
### genotyping on haplotype 1 with contigs ###

rule geno_hap1:
    input:
        ref= "Data/Fasta/Haplotype1.fasta"
    output:
        index_ghost="Data/Fasta/Haplotype1_geno/index_ghost.txt"
    container: c_geno
    threads: 1
    shell:
        """
        # Reference genome indexing
        bwa-mem2.avx2 index {input.ref}
        touch {output.index_ghost}
        """


rule MAPPING_EU_h1:
    input:
        index_ghost="Data/Fasta/Haplotype1_geno/index_ghost.txt",
        I1 = "Data/Fastq/{indiv}.1.trimmed.fastq.gz",
        I2 = "Data/Fastq/{indiv}.2.trimmed.fastq.gz",
        ref= "Data/Fasta/Haplotype1.fasta"
    output:
        O1="Data/MAPPING/Haplotype1_geno/{indiv}_mapping/{indiv}.bam.gz"
    container: c_geno
    threads: 40
    shell:
        """
        touch {input.index_ghost}
        # Mapping step: Save the output as a SAM file
        #bwa-mem2.avx2 mem -t 1 #
        bwa-mem2 mem -t 25 \
        -M -R "@RG\\tID:{wildcards.indiv}_1\\tSM:{wildcards.indiv}" \
        {input.ref} {input.I1} {input.I2} > {wildcards.indiv}.sam
        echo "Mapping done: SAM file generated."

        # Convert SAM to BAM
        samtools view -bS {wildcards.indiv}.sam > {output.O1}
        echo "SAM to BAM conversion done."
        """

rule Sorting_EU_h1:
    input:
        O1="Data/MAPPING/Haplotype1_geno/{indiv}_mapping/{indiv}.bam.gz"
    output:
        O2="Data/MAPPING/Haplotype1_geno/{indiv}_mapping/{indiv}.sorted.duplicates.bam.gz"
    container: c_geno
    threads: 10
    shell:
        """
        # Sorting step
        samtools sort -o {wildcards.indiv}.sorted.bam -@ 10 {input.O1}
        echo "Sorting done."

        # Mark duplicates
        picard MarkDuplicates \
        -I {wildcards.indiv}.sorted.bam \
        -M metrics_duplicates.txt \
        -O {output.O2} \
        -COMPRESSION_LEVEL 5 \
        -VALIDATION_STRINGENCY LENIENT 
        echo "Mark duplicates done."

        # Cleanup intermediate files
        rm {wildcards.indiv}.sam {wildcards.indiv}.sorted.bam
        echo "Cleanup done."

        # Indexing step
        samtools index {output.O2}
        echo "Indexing done."
        """


#### Variant Calling ### 
rule Variant_calling_Haplotype_caller_EU_h1:
    input: 
        I1="Data/MAPPING/Haplotype1_geno/{indiv}_mapping/{indiv}.sorted.duplicates.bam.gz",
        ref= "Data/Fasta/Haplotype1.fasta"
    output: 
        "Data/Variant_Calling/Haplotype1_geno/{indiv}/{indiv}_gatk.vcf.gz"
    container: c_geno
    threads: 15
    shell:
     """
    #index ref
    if [ ! -e "Data/Edulis/Haplotype2_renamed_reordered_Chr_only.dict" ]; then
    gatk CreateSequenceDictionary -R {input.ref}
    samtools faidx {input.ref}
    echo "File doesn't exist, running the command."
    fi
    samtools faidx {input.ref}
    gatk --java-options "-Xmx10g" HaplotypeCaller -ERC BP_RESOLUTION \
    -R {input.ref} \
    -I {input.I1} \
    -O {output} \
    --do-not-run-physical-phasing true
     """

rule make_scaffold_list_EU_h1:            # removing spaces in fasta seq names
    input:
        ref="Data/Fasta/Haplotype1.fasta"
    output:
        L="Data/Scaffold_hap1.list"
    shell:
        """
        # Create .list file with chromosome names
        grep "^>" {input} | sed 's/>//' > {output.L}
        """


rule DB_import_EU_h1:
     input: 
        vcf=expand( "Data/Variant_Calling/Haplotype1_geno/{indiv}/{indiv}_gatk.vcf.gz", indiv=config['individuals_EU']),
        Scaffold_list="Data/Scaffold_hap1.list"
     output: 
        ghost_EU="Data/Variant_Calling/Haplotype1_geno/ghost",
        #ghost="Variant_Calling/DB/vcfheader.vcf",
        DIR=directory("Data/Variant_Calling/Haplotype1_geno/DB")
     container: c_geno
     threads: 15
     shell:
        """
     vcf_all_samples=$(echo " {input.vcf}" | sed "s/\[//g; s/\]//g; s/,//g; s/'//g; s/ / -V /g")
     gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
        $vcf_all_samples \
        --genomicsdb-workspace-path {output.DIR} \
        --tmp-dir Variant_Calling \
        -L {input.Scaffold_list} 

     touch {output.ghost_EU}
        """

rule GenotypeGVCF_EU_h1:
    input:
        ghost_EU="Data/Variant_Calling/Haplotype1_geno/ghost",
        ref= "Data/Fasta/Haplotype1.fasta"
        #ghost="Variant_Calling/DB/vcfheader.vcf"
    output: "Data/VCF/EU_pop_all_site_unfiltered_h1.vcf.gz"
    container: c_geno
    threads: 10
    shell:
        """
    touch {input.ghost_EU}
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R {input.ref} \
        -V gendb://Haplotype1_geno/Variant_Calling/DB \
        -O {output} \
        --include-non-variant-sites true # All site vcf file
        """

rule Coverage_h1:
    input:
        I1 ="Data/VCF/EU_pop_all_site_unfiltered_h1.vcf.gz"
    output:
        H1_dp="Results/dp_h1.tsv"
    container: c_popgen
    threads: 7
    shell:
        """
        vcftools --gzvcf {input} --site-mean-depth --stdout > {output}
        """