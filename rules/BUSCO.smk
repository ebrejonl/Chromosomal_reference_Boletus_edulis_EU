rule BUSCO:
    input:
        fasta1="Data/Fasta/Haplotype1.fasta",
        fasta2="Data/Fasta/Haplotype2_renamed_reordered_with_contigs.fasta"
    output:
        ghost_busco="Results/Busco/busco_ghost.txt"
    container: c_synteny
    threads: 17
    shell:
        """
busco -i {input.fasta1} \
    -m genome \
    --force \
    -l basidiomycota \
    --augustus \
    --out_path Results/Busco \
    --cpu 18

# haplo 2
busco -i {input.fasta2} \
    -m genome \
    --force \
    -l basidiomycota \
    --augustus \
    --out_path Results/Busco \
    --cpu 18
    touch {output.ghost_busco}
        """


rule Hap_synteny:
    input:
        ghost_busco="Results/Busco/busco_ghost.txt",
        fai="Data/Fasta/Haplotype2_renamed_reordered_Chr_only.fasta.fai"
    output:
        busco_plot="Results/Busco/Synteny_haplotypes_busco.RDS"
    container: c_R 
    script:
        "rules/Rscripts/Synteny_btw_haplo.R"


rule Gene_TE_hap1_vs_hap2:
    input:
        busco_plot="Results/Busco/Synteny_haplotypes_busco.RDS",
        gff=expand("Data/Fasta/{genome}.fasta.mod.EDTA.TEanno.gff3", genome=config["genomes"])
    output:
        te1="Results/Pgenes.RDS",
        te2="Results/Teplot.RDS",
        te3="Results/model_gene.RDS"
    container: c_R 
    script:
        "rules/Rscripts/Gene_comp.R"


rule Coveverage_R:
    input: H1_dp="Results/dp_h1.tsv"
    output: Coverage= "Results/Coverage_h1.RDS"
    container: c_R 
    script:
        "rules/Rscripts/Coverage_h1.R"


rule Plotting_hap1vshap2:
    input:
        te1="Results/Pgenes.RDS",
        te2="Results/Teplot.RDS",
        te3="Results/model_gene.RDS",
        busco_plot="Results/Busco/Synteny_haplotypes_busco.RDS",
        Coverage= "Results/Coverage_h1.RDS",
        fai="Data/Fasta/Haplotype2_renamed_reordered_Chr_only.fasta.fai"
    output:
        Fig1="Figure1.pdf"
    container: c_R 
    script:
        "rules/Rscripts/Synteny_btw_haplo.R"