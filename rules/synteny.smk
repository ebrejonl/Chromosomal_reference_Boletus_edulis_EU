## synteny with Suillus
rule Busco_suillus:
    input:
        Suillus="Data/Fasta/Suillus_bovinus.fasta"
    output:
        ghost_busco="Results/Busco/busco_ghost_suillus.txt" # real output is in Results/Busco/Suillus_busco.tsv
    container: c_synteny
    threads: 17
    shell:
        """
        busco -i {input.Suillus} \
            -m genome \
            --force \
            -l basidiomycota \
            --augustus \
            --cpu 10
        
        """

# plotting 
rule Synteny:
    input:
    Suillus_busco="Results/Busco/Suillus_busco.tsv",
    Haplotype2_busco="Results/Busco/Haplotype2_busco.tsv"
    output:
    Circle_plot="Chord_diagram_Suillus.pdf"
    container: c_R 
    script:
        "rules/Rscripts/Suillus_synteny.R"