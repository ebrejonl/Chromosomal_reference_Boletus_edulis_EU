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
        #"Results/juicer/chrom.sizes",
        #"Results/juicer/aligned_reads.bam",
        #"Results/juicer/contact_map.hic"
        #"Results/Busco/busco_ghost.txt",
        "Results/STRUCTURE/Whole_genomefiltered.vcf.gz"




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#include: "./rules/Genotyping_Hap1.smk"
#include: "./rules/Genotyping_Hap2.smk"
include: "./rules/structure.smk"
#include: "./rules/Telomeres_pred.smk" 
#include: "./rules/Coverage_data.smk" 
#include: "./rules/Hi-C_contact_map.smk"
#include: "./rules/BUSCO.smk"
