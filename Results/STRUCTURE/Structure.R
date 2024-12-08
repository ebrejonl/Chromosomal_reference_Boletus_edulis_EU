library(tidyverse)
library(vcfR)
library(pcadapt)
library(jtools)
library(ggrepel)
library(maps)
library(sf)
library(rnaturalearth)
library(ggforce)
library(MetBrewer)

# Part 1) METADATA
setwd("~/Desktop/work/Dropbox/Reference_genome_21_02_2024/Proper_Chromosome_names/Chromosome3/")
metadata=read.csv("EU_metadata.csv", header=F, sep=",") |> mutate(lon=as.numeric(paste0(paste0(V5, "."), V6))) |> select(-V5, -V6) |> 
    rename(ID=V2, Pop=V7,
  lat=V4) |>  mutate(lat=as.numeric(lat))

# Part 2) genomic PCA
setwd("~/Desktop/work/Dropbox/Reference_genome_21_02_2024/14_11/Chromosomal_reference_Boletus_edulis_EU/")
wgs_full
pca_vcf <- read.vcfR("snakemake@input[['wgs_full']]")
print("VCF loaded")
vecteur_noms<-c(names(pca_vcf@gt[1,]))
vecteur_noms<-vecteur_noms[-1] # remove format
# names loaded
#create pcadapt object and add ID from VCF -> with pcadapt
respca_10_maf = pcadapt(input = read.pcadapt("Results/STRUCTURE/Whole_genomefiltered.vcf.gz", type = "vcf") ,K = 5) # maf à choisir , min.maf=0.08
scores = data.frame(respca_10_maf$scores)
rownames(scores)= vecteur_noms

print("pca object created, and vcf indiv names added")
#### Compute proportion of explained variance
PEV = paste(round((respca_10_maf$singular.values^2)*100,digits=1),"%",sep="")

# merging with metadata 
scores =scores |> mutate(ID=rownames(scores)) |> left_join(metadata, by = "ID") |> 
mutate(lon=as.numeric(str_replace(lon, ",", ".")),
lat=as.numeric(lat))
print("merged")

PEV # run first to check the values here and then add them to scores below
scores =scores|> dplyr::rename("5%"=X1,
              "3.4%"=X2,
            "3.1%"=X3,"2.7%"=X4)


saveRDS(scores, "snakemake@output[['pca']]")
