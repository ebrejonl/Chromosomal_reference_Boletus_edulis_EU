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
#PEV = "5%"   "3.4%" "3.1%" "2.7%" "2.6%"
scores =scores|> dplyr::rename("5%"=X1,
              "3.4%"=X2,
            "3.1%"=X3,"2.7%"=X4)


saveRDS(scores, "Results/STRUCTURE/Whole_genome_pca_with_pop.rds")



mydata=readRDS("Results/STRUCTURE/Whole_genome_pca_with_pop.rds") |> 
  mutate(Pop=case_when(
    Pop=="Scandinavia" ~ "Fennoscandia", 
    Pop=="South" ~ "Southern EU",
    Pop=="USA"~"North America",
    Pop=="GB"~"Great Britain", TRUE~as.character(Pop)))

p1=ggplot(data = mydata, aes(x=mydata[,1],y=mydata[,2], fill=Pop)) + 
      theme_apa() +
        geom_vline(xintercept = 0, lty = 3 , alpha=0.6) +
          geom_hline(yintercept = 0,lty = 3,alpha=0.6 )+
            scale_fill_manual(values =met.brewer("Johnson", n=7, direction=-1), 
          name="Population") +
      geom_point(size=5, shape=21, color="black", alpha=0.5)  +
geom_mark_ellipse(mapping=aes(x=mydata[,1],y=mydata[,2], fill=Pop))+
xlab("PC1 5%")+ ylab("PC2 3.4%") +theme_classic()+ 
  theme(legend.position = "top", legend.direction = "horizontal")+  
  guides(fill = guide_legend(nrow = 1))
gc()
p1



# script to make full figure is in PI_DIVERSITY -> Plotting_Pi.R