library(ape)
library(data.table)
library(tidyverse)
library(here)

hap1_gff <- read.gff("Boletus_edulis_bielefeld_haplotype1.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
    group_by(chr) %>% count(type) %>% mutate(Haplotype = "1") 

hap2_gff <- read.gff("Boletus_edulis_bielefeld_haplotype2.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
   group_by(chr) %>% count(type) %>% mutate(Haplotype = "2")


full_gff <- rbind(hap1_gff, hap2_gff)
unique(full_gff$chr)

read_fai <- \(file){
  read_tsv(file,
           col_names = c("name", "length", "offset", "linebases", "linewidth"))%>%
  dplyr::select(name, length) %>%
  mutate(g_end = cumsum(length),
         g_start= lag(g_end, default = 0),
         g_mid=(g_start+g_end)/2,
         idx=row_number(),
         even_odd=if_else(idx<12,
                          idx%%2,
                          2)) }

Hap_1_indec <- read_fai(file="Data/Fasta/Chr_only_hap1.fasta.fai")  %>% mutate(Haplotype = "1")
Hap_2_indec <- read_fai(file="Data/Fasta/Haplotype2_renamed_reordered_Chr_only.fasta.fai")  %>% mutate(Haplotype = "2")
index_full <- rbind(Hap_1_indec, Hap_2_indec)
comparison <- aggregate(length ~ name, index_full, function(x) ifelse(max(x) == x[1], "1", "2")) %>%
  mutate(chr=name)%>% select(-name)
full_gff <- full_gff %>% left_join(comparison, by = "chr")


longest_db= full_gff%>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% group_by(chr) %>% 
  filter(Comp=="Longest") %>% summarize(haplotype=mean(as.numeric(Haplotype), 
                                                       length=mean(as.numeric(length))))


pgenes <-full_gff %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% group_by(chr) %>% 
  filter(type=="gene")%>%
  ggplot()+geom_boxplot(mapping=aes(x=Comp,y=n), width=0.4) +
  geom_jitter(mapping=aes(x=Comp,y=n, fill=Haplotype), shape=21,
alpha=0.7,position = position_jitter(w = 0.15, h = 0),   color="black", size=3.5)+
  scale_fill_manual(values = c("#697787ff", "#c9400bff"))+
  theme_classic() +xlab("")+ ylab("Total gene number")+theme(legend.position = "none")

saveRDS(file=snakemake@output[['te1']], pgenes) 



Hap_1_indec <- Hap_1_indec %>% dplyr::rename("chr"="name")
Hap_2_indec <- Hap_2_indec %>% dplyr::rename(chr=name, length.y=length)
p_freg1 <- full_gff %>% filter(Haplotype=="1") %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% group_by(chr) %>% 
  filter(type=="gene") %>%
  left_join(Hap_1_indec, by=c("chr", "Haplotype"))
p_freg2 <- full_gff %>% filter(Haplotype=="2") %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% group_by(chr) %>% 
  filter(type=="gene") %>%
  left_join(Hap_2_indec, by=c("chr", "Haplotype"))
p_freg <- rbind(p_freg1, p_freg2) %>%
  mutate(length.y=as.integer(length.y))

library(lme4)

p_freg %>% filter(Haplotype=="2")
preg <- glmer(formula=n ~ ength.y/1000000 + (1|chr) , data = p_freg, family="poisson")
preg <- glm(formula=p_freg$n ~ p_freg$length.y , data = p_freg, family="poisson")
summary(preg)


## Plotting
req <- ggplot(p_freg)+
  geom_smooth(aes(x=length.y/1000000, y=n), method="glm",  
method.args = list(family = "poisson"),color = "black", linewidth=0.3)+
  geom_point(aes(x=length.y/1000000, y=n, fill=Haplotype), shape=21,
color="black", size =4, alpha=0.7)+
  theme_classic()+ theme(axis.title.x.bottom = element_text(vjust = -0.5))+
  scale_fill_manual(values = c("#697787ff", "#c9400bff"))+xlab("Chromosome Length (Mb)")+
  ylab("Gene number")
req
saveRDS(req, snakemake@output[['te3']])


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TE in bp length ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

hap1_gff <- read.gff("Data/Fasta/Haplotype1_onlyCHR.fasta.mod.EDTA.TEanno.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
  mutate(length=end-start)  %>%
  group_by(chr, type) %>% summarize(L=sum(length)) %>% mutate(Haplotype = "1")
# Haplotype 2
hap2_gff <- read.gff("../../Plotting_chr_synteny/Kearton send gff/Haplotype2_onlyCHR.fasta.mod.EDTA.TEanno.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
    mutate(length=end-start)  %>%
  group_by(chr, type) %>% summarize(L=sum(length))  %>% mutate(Haplotype = "2")
full_TE_gff <- rbind(hap1_gff, hap2_gff) 


full_TE_gff  <- full_TE_gff %>% left_join(comparison, by = "chr")
unique(full_TE_gff$type)

p_TE <- full_TE_gff %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% 
  group_by(chr, Haplotype) %>% 
  summarize(l=sum(L), Comp=Comp) %>% distinct() %>%
  ggplot() +geom_boxplot(mapping=aes(x=Comp,y=l/1000000), width=0.4) +
  geom_jitter(mapping=aes(x=Comp,y=l/1000000,fill=Haplotype), shape=21,
alpha=0.7,position = position_jitter(w = 0.15, h = 0),   color="black", size=3.5)+
  scale_fill_manual(values = c("#697787ff", "#c9400bff"))+
  theme_classic() +xlab("Chromosome")+ ylab("Length in TEs (Mb)")+
  theme(legend.position = "none")+theme(axis.title.x.bottom = element_text(vjust = -0.5))

saveRDS(p_TE, snakemake@output[['te2']])

