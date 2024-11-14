### Structure analysis 

# Nucleotide diversity pi 
rule pi_per_pop:
    input:
        vcf="Data/Edulis/Keaton_EU_samples/Variant_Calling/EU_pop_all_site_unfiltered.vcf.gz"
    output:
        "RESULTS/DIVERSITY/PI_{subpop}.windowed.pi"
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

### Principial component analysis 

