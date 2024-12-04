library(tidyverse)
library(vcfR)
library(pcadapt)
library(jtools)
library(ggrepel)


# Part 1) METADATA
setwd("~/Desktop/work/Dropbox/Reference_genome_21_02_2024/Proper_Chromosome_names/Chromosome3/")
metadata=read.csv("EU_metadata.csv", header=F, sep=",") |> mutate(lon=as.numeric(paste0(paste0(V5, "."), V6))) |> select(-V5, -V6) |> 
    rename(ID=V2, Pop=V7,
  lat=V4) |>  mutate(lat=as.numeric(lat))
  
  library(maps)
  library(sf)
  library(rnaturalearth)
  library(ggforce)
  library(MetBrewer)

# Part 2) genomic PCA
setwd("~/Desktop/work/Dropbox/Reference_genome_21_02_2024/14_11/Chromosomal_reference_Boletus_edulis_EU/")
pca_vcf <- read.vcfR("Results/STRUCTURE/Whole_genomefiltered.vcf.gz")
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

PEV

scores =scores|> dplyr::rename("5.4%"=X1,
              "4.5%"=X2,
            "4.1%"=X3,"4%"=X4)


saveRDS(scores, "Whole_genome_pca_with_pop.rds")