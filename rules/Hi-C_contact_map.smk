
rule prep_juicer:
    input:
        ref_fasta="Data/Edulis/Haplotype2_renamed_reordered_Chr_only.fasta.fai"
    output:
        chrom.sizes="Results/juicer/chrom.sizes"
    params:
        output_dir="results/juicer"
    shell:
        """
        mkdir -p {params.output_dir}
        cut -f1,2 {input.ref_fasta} > {output.chrom.sizes}
        """

rule Hi_c_map:
    input: 
        fq1= "Data/hi_c_1.fq.gz",
        fq2= "Data/hi_c_2.fq.gz",
        ref= "Data/Edulis/Haplotype2_renamed_reordered_Chr_only.fasta",
        chrom.sizes="Results/juicer/chrom.sizes"
     output: 
        contact_map="Results/juicer/contact_map.hic"
      params:
        output_dir="Results/juicer",  # Path to output directory
        genome_id="hg38",  # path tp chr size file??
     container: c_popgen
     threads: 15
         shell:
        """
        {params.juicer_dir}/scripts/juicer.sh \
            -d {params.output_dir} \
            -z {input.reference} \
            -p {params.output_dir}/chrom.sizes \
            -g {params.genome_id} \
            #-s HindIII \  # Change enzyme as needed
            -t {params.threads} \
            -f {input.fastq1},{input.fastq2}
        """