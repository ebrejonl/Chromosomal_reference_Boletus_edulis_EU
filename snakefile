# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Snakemake workflow for Boletus edulis chr assembly                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

configfile: "config.yaml"
import os
c_geno = config[ 'sif_genotyping' ] # for short
c_popgen = config['sif_popgen']
c_synteny = config['sif_synteny']

rule all:
    input:
        "Results/juicer/chrom.sizes",
        "Results/juicer/aligned_reads.bam",
        "Data/VCF/EU_pop_all_site_unfiltered.vcf.gz",
        "Results/juicer/contact_map.hic",
        "Results/Busco/busco_ghost.txt",
        "Results/STRUCTURE/Whole_genomefiltered.vcf.gz",
        "Data/Annotation/TELO_telomeric_repeat_windows.tsv",
        "Results/dp_h1.tsv",
        "Karyotype_plot.pdf",
        "Data/Fasta/Haplotype2_renamed_reordered_with_contigs.fasta.fai"




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
include: "./rules/Genotyping_Hap1.smk"
include: "./rules/Genotyping_Hap2.smk"
include: "./rules/structure.smk"
include: "./rules/Telomeres_pred.smk" 
include: "./rules/Hi-C_contact_map.smk"
include: "./rules/BUSCO.smk"
include: "./rules/TE_pred.smk"
include: "./rules/karyotype_plot.smk"
