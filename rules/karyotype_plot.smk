## Karyotype plot
# gather data from TE pred, this rule is for the DAG
# The data for cazymes etc... are in the annotation file in the zenodo repos.
rule Annotation_files_gather:
    input: "Data/Fasta/{genome}.softMasked.fasta"
    output: "Data/Annotation/Cazymes.tsv"
    shell:
        """
        touch {input}
        touch {output}
        """


# Gene mapping in R for Figure 3
rule Karyotype_plot:
    input:
        hap2_index="Data/Fasta/Haplotype2_renamed_reordered_with_contigs.fasta.fai",
        Caz_coor="Data/Annotation/Cazymes.tsv"
    output:
        Karyo_plot="Karyotype_plot.pdf"
    container: c_R 
    script:
        "rules/Rscripts/Karyotype_plot.R"