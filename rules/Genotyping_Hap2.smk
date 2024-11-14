## Genotyping on haplotype 2 chromosomal only reference 

rule index_MAPPING_EU:
    input:
        ref= "Data/Edulis/Haplotype2_renamed_reordered_Chr_only.fasta"
    output:
        index_ghost="Data/Edulis/index_ghost.txt"
    container: c_geno
    threads: 1
    shell:
        """
        # Reference genome indexing
        bwa-mem2.avx2 index {input.ref}
        touch {output.index_ghost}
        """

rule MAPPING_EU:
    input:
        index_ghost="Data/Edulis/index_ghost.txt",
        I1 = "Data/Edulis/Keaton_EU_samples/{indiv}.1.trimmed.fastq.gz",
        I2 = "Data/Edulis/Keaton_EU_samples/{indiv}.2.trimmed.fastq.gz",
        ref= "Data/Edulis/Haplotype2_renamed_reordered_Chr_only.fasta"
    output:
        O1="Data/Edulis/Keaton_EU_samples/MAPPING/{indiv}_mapping/{indiv}.bam.gz"
    container: c_geno
    threads: 10
    shell:
        """
        touch {input.index_ghost}
        # Mapping step: Save the output as a SAM file
        #bwa-mem2.avx2 mem -t 1 #
        bwa-mem2 mem -t 1 \
        -M -R "@RG\\tID:{wildcards.indiv}_1\\tSM:{wildcards.indiv}" \
        {input.ref} {input.I1} {input.I2} > {wildcards.indiv}.sam
        echo "Mapping done: SAM file generated."

        # Convert SAM to BAM
        samtools view -bS {wildcards.indiv}.sam > {output.O1}
        echo "SAM to BAM conversion done."
        """

rule Sorting_EU:
    input:
        O1="Data/Edulis/Keaton_EU_samples/MAPPING/{indiv}_mapping/{indiv}.bam.gz",
    output:
        O2="Data/Edulis/Keaton_EU_samples/MAPPING/{indiv}_mapping/{indiv}.sorted.duplicates.bam.gz"
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
rule Variant_calling_Haplotype_caller_EU:
    input: 
        I1="Data/Edulis/Keaton_EU_samples/MAPPING/{indiv}_mapping/{indiv}.sorted.duplicates.bam.gz",
        ref="Data/Edulis/Haplotype2_renamed_reordered_Chr_only.fasta"
    output: 
        "Data/Edulis/Keaton_EU_samples/Variant_Calling/{indiv}/{indiv}_gatk.vcf.gz"
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

rule make_scaffold_list_EU:            # removing spaces in fasta seq names
    input:
        "Data/Edulis/Haplotype2_renamed_reordered_Chr_only.fasta"
    output:
        L="Data/Edulis/Keaton_EU_samples/Scaffold_hap2.list"
    shell:
        """
        # Create .list file with chromosome names
        grep "^>" {input} | sed 's/>//' > {output.L}
        """





rule DB_import_EU:
     input: 
        vcf=expand( "Data/Edulis/Keaton_EU_samples/Variant_Calling/{indiv}/{indiv}_gatk.vcf.gz", indiv=config['individuals_EU']),
        Scaffold_list="Data/Edulis/Keaton_EU_samples/Scaffold_hap2.list"
     output: 
        ghost_EU="Data/Edulis/Keaton_EU_samples/Variant_Calling/ghost",
        #ghost="Variant_Calling/DB/vcfheader.vcf",
        DIR=directory("Data/Edulis/Keaton_EU_samples/Variant_Calling/DB")
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

rule GenotypeGVCF_EU:
    input:
        ghost_EU="Data/Edulis/Keaton_EU_samples/Variant_Calling/ghost",
        ref="Data/Edulis/Haplotype2_renamed_reordered_Chr_only.fasta"
        #ghost="Variant_Calling/DB/vcfheader.vcf"
    output: "Data/Edulis/Keaton_EU_samples/Variant_Calling/EU_pop_all_site_unfiltered.vcf.gz"
    container: c_geno
    threads: 10
    shell:
        """
    touch {input.ghost_EU}
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R {input.ref} \
        -V gendb://Data/Edulis/Keaton_EU_samples/Variant_Calling/DB \
        -O {output} \
        --include-non-variant-sites true # All site vcf file
        """
