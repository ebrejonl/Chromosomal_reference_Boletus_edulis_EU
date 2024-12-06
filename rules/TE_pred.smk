# This prediction was made by Tremble K github://
# added post-hoc as rule to the snakemake pipeline for consistency
rule TE_pred:
    input:
        genome="Data/Fasta/{genome}.fasta",
        cds="Data/Fasta/{genome}.cds.fasta",
    output:
        soft_masked_genome="Data/Fasta/{genome}.softMasked.fasta",
        funannotate_dir="Data/Fasta/{genome}_funannotate"
    params:
        species="Boletus edulis",
        isolate="bielefeld_haplotype1",
        augustus_species="boletus_edulis_bd747",
        busco_db="basidiomycota",
        funannotate_db="~/work_programs/funannotate_db_22april24/"
    threads: 32
    container: c_synteny
    shell:
        """
        # Run EDTA
        EDTA.pl --genome {input.genome} --cds {input.cds} --anno 1 --threads {threads} --sensitive 1
        # Mask the genome using bedtools
        bedtools maskfasta -fi {input.genome} -soft -fo {output.soft_masked_genome} -bed {wildcards.genome}.mod.EDTA.TEanno.gff3
        # Run funannotate predict
        funannotate predict \
            -s "{params.species}" \
            --isolate {params.isolate} \
            --augustus_species {params.augustus_species} \
            --cpus {threads} \
            -i Haplotype1_edit_onlyCH \
            -o {output.funannotate_dir} \
            --protein_evidence Data/Annotation/EVM.final.gene.gff3.pep \
            --busco_db {params.busco_db} \
            --optimize_augustus \
            -d {params.funannotate_db}
        """