### Structure analysis 
rule PCA_prep_EU_indiv:
    input:
        "Data/VCF/EU_pop_all_site_unfiltered.vcf.gz" 
    output:
        wgs_full="Chromosome_3/Whole_genomefiltered.vcf.gz"
    container: c_popgen
    threads: 20
    shell:
        """ 
        vcftools --gzvcf {input} --maf 0.05 \
           --max-missing 0.8 \
           --minQ 15 \
           --minDP 5 --maxDP 100 --recode --stdout | bgzip -c > {output.wgs_full}
        """

# Nucleotide diversity pi 
rule pi_per_pop:
    input:
        vcf="Data/VCF/EU_pop_all_site_unfiltered.vcf.gz"
    output:
        "Results/DIVERSITY/PI_{subpop}.windowed.pi"
    params:
        samples=lambda wildcards: config[wildcards.subpop] if wildcards.subpop in ["Fennoscandia", "Central", "Iceland", "Great_Britain"] else None,
        indivs=lambda wildcards: " ".join([f"--indv {ind}" for ind in config[wildcards.subpop]])
    container: c_popgen
    threads: 20
    shell:
       """
         vcftools --gzvcf {input.vcf} --minQ 15 \
            {params.indivs} \
            --minDP 5 --maxDP 100 --max-missing 0.8 \
            --window-pi 10000 --window-pi-step 1000 --stdout > {output}
       """ 


