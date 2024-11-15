
rule prep_juicer:
    input:
        ref_fasta_index="Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta.fai.fai"
    output:
        chrom.sizes="Results/juicer/chrom.sizes"
    params:
        output_dir="results/juicer"
    shell:
        """
        mkdir -p {params.output_dir}
        cut -f1,2 {input.ref_fasta_index} > {output.chrom.sizes}
        """

rule Hi_c_map:
    input: 
        fq1= "Data/HiC/hi_c_1.fq.gz",
        fq2= "Data/HiC/hi_c_2.fq.gz",
        ref= "Data/Ref/Haplotype2_renamed_reordered_Chr_only.fasta",
        chrom.sizes="Results/juicer/chrom.sizes"
     output: 
        contact_map="Results/juicer/contact_map.hic"
      params:
        output_dir="Results/juicer",  # Path to output directory
        genome_id="BolEdBiel_h2",
        threads=15
     container: c_popgen
         shell:
        """
        {params.juicer_dir}/scripts/juicer.sh \
            -d {params.output_dir} \
            -z {input.ref} \
            -p {input.chrom.sizes} \
            -g {params.genome_id} \
            -t {params.threads} \
            -f {input.fq1},{input.fq2}
        """