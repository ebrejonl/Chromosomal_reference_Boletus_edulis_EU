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