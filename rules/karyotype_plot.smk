rule Karyotype_plot:
    input:
        hap2_index="Data/Fasta/Haplotype2_renamed_reordered_with_contigs.fasta.fai",
        Caz_coor="Data/Annotation/Cazymes.tsv"
    output:
        Karyo_plot="Karyotype_plot.pdf"
    container: c_R 
    script:
        "rules/Rscripts/Karyotype_plot.R"