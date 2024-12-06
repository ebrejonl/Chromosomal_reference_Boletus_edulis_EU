library(ape)
library(data.table)
library(tidyverse)
setwd("~/Desktop/work/Dropbox/Reference_genome_21_02_2024/Proper_Chromosome_names/Comp_gene_TE/")

# Haplotype 1
hap1_gff <- read.gff("../../../Reference_genome_21_02_2024/Plotting_chr_synteny/Kearton send gff/Boletus_edulis_bielefeld_haplotype1.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
    group_by(chr) %>% count(type) %>% mutate(Haplotype = "1") |> 
      mutate(chr=case_when(
        chr=="group0" ~ "Chr1",
        chr=="group1" ~ "Chr2",
        chr=="group2" ~ "Chr3",
        chr=="group3" ~ "Chr4",
        chr=="group4" ~ "Chr5",
        chr=="group5" ~ "Chr6",
        chr=="group6" ~ "Chr7",
        chr=="group7" ~ "Chr8",
        chr=="group8" ~ "Chr9",
        chr=="group9" ~ "Chr10",
        chr=="group10" ~ "Chr11"))

# Haplotype 2 proper order/name of chr
hap2_gff <- read.gff("../../../Reference_genome_21_02_2024/Plotting_chr_synteny/Kearton send gff/Boletus_edulis_bielefeld_haplotype2.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
   group_by(chr) %>% count(type) %>% mutate(Haplotype = "2") |> 
  mutate(chr=case_when(
    chr=="group0" ~ "Chr7",
    chr=="group1" ~ "Chr5",
    chr=="group2" ~ "Chr1",
    chr=="group3" ~ "Chr3",
    chr=="group4" ~ "Chr6",
    chr=="group5" ~ "Chr8",
    chr=="group6" ~ "Chr10",
    chr=="group7" ~ "Chr2",
    chr=="group8" ~ "Chr4",
    chr=="group9" ~ "Chr11",
    chr=="group10" ~ "Chr9"))



# Full gff tibble
full_gff <- rbind(hap1_gff, hap2_gff)
unique(full_gff$chr)

## now I need to add a column which specify which haplotype has the shortest chromosome
# I need to read the index fai files from each haplotypes
# function to read fai files
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

Hap_1_indec <- read_fai(file="../Chr_only_hap1.fasta.fai")  %>% mutate(Haplotype = "1")
Hap_2_indec <- read_fai(file="../Haplotype2_renamed_reordered_Chr_only.fasta.fai")  %>% mutate(Haplotype = "2")
index_full <- rbind(Hap_1_indec, Hap_2_indec)
comparison <- aggregate(length ~ name, index_full, function(x) ifelse(max(x) == x[1], "1", "2")) %>%
  mutate(chr=name)%>% select(-name)
# Now let's add the longest haplotype in the gff table.
full_gff <- full_gff %>% left_join(comparison, by = "chr")


longest_db= full_gff%>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% group_by(chr) %>% 
  filter(Comp=="Longest") %>% summarize(haplotype=mean(as.numeric(Haplotype), 
                                                       length=mean(as.numeric(length))))
#write.csv2(file="Longest_chromosome_list.csv", longest_db)


###blue
pgenes <-full_gff %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% group_by(chr) %>% 
  filter(type=="gene")%>%
  ggplot()+geom_boxplot(mapping=aes(x=Comp,y=n), width=0.4) +
  geom_jitter(mapping=aes(x=Comp,y=n, fill=Haplotype), shape=21,
alpha=0.7,position = position_jitter(w = 0.15, h = 0),   color="black", size=3.5)+
  scale_fill_manual(values = c("#697787ff", "#c9400bff"))+
  theme_classic() +xlab("")+ ylab("Total gene number")+theme(legend.position = "none")#+
  #theme(aspect.ratio = 1.5)
pgenes

setwd("~/Desktop/work/Dropbox/Reference_genome_21_02_2024/Proper_Chromosome_names/MAKE_ALL_FIGURES/")
saveRDS(file="pgenes_blue.rds", pgenes) 



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
saveRDS(req, "model_genes_nonoblue.RDS")

### FIN DU PLAYGROUND

#



#~~~~~~~~~~~~~~~~~~~~~~~~~~     Same with repeats/Tes   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#
#hap1_gff <- read.gff("../../Plotting_chr_synteny/Kearton send gff/Haplotype1_onlyCHR.fasta.mod.EDTA.TEanno.gff3") %>%
#    set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
#      group_by(chr) %>% count(type) %>% mutate(Haplotype = "1") |> 
#          mutate(chr=case_when(
#            chr=="group0" ~ "Chr1",
#            chr=="group1" ~ "Chr2",
#            chr=="group2" ~ "Chr3",
#            chr=="group3" ~ "Chr4",
#            chr=="group4" ~ "Chr5",
#            chr=="group5" ~ "Chr6",
#            chr=="group6" ~ "Chr7",
#            chr=="group7" ~ "Chr8",
#            chr=="group8" ~ "Chr9",
#            chr=="group9" ~ "Chr10",
#            chr=="group10" ~ "Chr11"))
#  # Haplotype 2
#  hap2_gff <- read.gff("../../Plotting_chr_synteny/Kearton send gff/Haplotype2_onlyCHR.fasta.mod.EDTA.TEanno.gff3") %>%
#    set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
#     group_by(chr) %>% count(type) %>% mutate(Haplotype = "2") |> 
#        mutate(chr=case_when(
#          chr=="group0" ~ "Chr7",
#          chr=="group1" ~ "Chr5",
#          chr=="group2" ~ "Chr1",
#          chr=="group3" ~ "Chr3",
#          chr=="group4" ~ "Chr6",
#          chr=="group5" ~ "Chr8",
#          chr=="group6" ~ "Chr10",
#          chr=="group7" ~ "Chr2",
#          chr=="group8" ~ "Chr4",
#          chr=="group9" ~ "Chr11",
#          chr=="group10" ~ "Chr9"))
#      
#
#  full_TE_gff <- rbind(hap1_gff, hap2_gff) 
#
#
#full_TE_gff  <- full_TE_gff %>% left_join(comparison, by = "chr")
#unique(full_TE_gff$type)
#
#
#p_TE <- full_TE_gff %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% 
#  group_by(chr, Haplotype) %>% 
#  summarize(n=sum(n), Comp=Comp) %>% distinct() %>%
#  ggplot() +geom_boxplot(mapping=aes(x=Comp,y=n), width=0.4) +
#  geom_jitter(mapping=aes(x=Comp,y=n, fill=Haplotype), shape=21,
# alpha=0.5,position = position_jitter(w = 0.12, h = 0),   color="black", size=3.5)+
#  scale_fill_manual(values = c("#117b5c", "#a3429f"))+
#theme_bw() +xlab("Chromosome")+ ylab("Total TEs")
#
#
#saveRDS(p_TE, "TE_plot.RDS")
#p_TE
#
#library(patchwork)
##ppp <- pgenes/p_TE
##ggsave(file="Chromosome_comparison.pdf", ppp, width=6, height = 10)
#
### Te type per chr
#p_Type <- full_TE_gff %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% 
#group_by(Comp, type) %>% summarize(m=mean(n)) %>%
#  ggplot()+
#  geom_col(mapping = aes(x = type, y = m, fill=Comp), alpha=0.8)+
#  coord_flip()+
#  scale_fill_manual(values = c("#117b5c", "#ee9566"), name="Chromosome")+
#  xlab("")+ylab("")+
#  theme_bw()
#p_Type
#  #geom_jitter(mapping = aes(x = type, y = n, fill=Comp),alpha=0.5,
#  #position = position_jitter(w = 0.1, h = 0))
#
#
#
#
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ In bp length ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
setwd("~/Desktop/work/Dropbox/Reference_genome_21_02_2024/Proper_Chromosome_names/Comp_gene_TE/")

hap1_gff <- read.gff("../../Plotting_chr_synteny/Kearton send gff/Haplotype1_onlyCHR.fasta.mod.EDTA.TEanno.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
  mutate(length=end-start)  %>%
  group_by(chr, type) %>% summarize(L=sum(length)) %>% mutate(Haplotype = "1")
# Haplotype 2
hap2_gff <- read.gff("../../Plotting_chr_synteny/Kearton send gff/Haplotype2_onlyCHR.fasta.mod.EDTA.TEanno.gff3") %>%
  set_names("chr", "method", "type", "start", "end", "dot", "direction", "phase", "rest") %>%
    mutate(length=end-start)  %>%
  group_by(chr, type) %>% summarize(L=sum(length))  %>% mutate(Haplotype = "2")
full_TE_gff <- rbind(hap1_gff, hap2_gff) %>%
    mutate(chr=case_when(
      chr=="group0" ~ "Chr1",
      chr=="group1" ~ "Chr2",
      chr=="group2" ~ "Chr3",
      chr=="group3" ~ "Chr4",
      chr=="group4" ~ "Chr5",
      chr=="group5" ~ "Chr6",
      chr=="group6" ~ "Chr7",
      chr=="group7" ~ "Chr8",
      chr=="group8" ~ "Chr9",
      chr=="group9" ~ "Chr10",
      chr=="group10" ~ "Chr11"))


full_TE_gff  <- full_TE_gff %>% left_join(comparison, by = "chr")
unique(full_TE_gff$type)
#
#
#p_TE <- full_TE_gff %>% mutate(Comp=ifelse(length==Haplotype, "Longest", "Shortest")) %>% 
#  group_by(chr, Haplotype) %>% 
#  summarize(l=sum(L), Comp=Comp) %>% distinct() %>%
#  ggplot() +geom_boxplot(mapping=aes(x=Comp,y=l/1000000), width=0.4) +
#  geom_jitter(mapping=aes(x=Comp,y=l/1000000,fill=Haplotype), shape=21,
#alpha=0.5,position = position_jitter(w = 0.15, h = 0),   color="black", size=3.5)+
#  scale_fill_manual(values = c("#117b5c", "#a3429f"))+
#  #geom_line(mapping = aes(x=Comp,y=l/1000000,group=chr),  color="#d7dbdd", linewidth=0.3)+
#  theme_classic() +xlab("Chromosome")+ ylab("Length in TEs (Mb)")+
#  theme(legend.position = "none")+theme(axis.title.x.bottom = element_text(vjust = -0.5))#+
#   # theme(aspect.ratio = 1.5)
#p_TE
#
#
#saveRDS(p_TE, "TE_plot_length.rds")
