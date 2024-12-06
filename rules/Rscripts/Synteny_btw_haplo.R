library(tidyverse)
library(ggrastr)
library(data.table)
library(here)

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

Hap_1_indec <- read_fai(file="../../Chr_only_hap1.fasta.fai")  %>% 
  mutate(Haplotype = "1")%>%
  rename("seq_id"="name")
Hap_2_indec <- read_fai(file="snakemake@input[['fai']]")  %>% 
  mutate(Haplotype = "2")%>%
  rename("seq_id"="name") #|>
index_full <- rbind(Hap_1_indec, Hap_2_indec) 

index_full = index_full|> group_by(Haplotype) |> mutate(size=g_end-g_start,
        end_with_skip=cumsum(size+150000),
      start=lag(end_with_skip, default=0),
    end=start+size, mid=(start + end) /2) |> rename("Chr"=seq_id) |> mutate(#y=as.numeric(Haplotype),
    y=case_when(
      Haplotype=="1"~ 2,
      Haplotype=="2" ~1),
      label_sign=ifelse(Haplotype=="2", -1, 1))



flipscaf=c("Chr5", "Chr6", "Chr7", "Chr9")

# BUSCO
df1 <-  read.delim("Hap1_busco.tsv", skip=3, header=FALSE, stringsAsFactors = FALSE) |> 
  mutate(Haplotype="1") |>  select(-V6, -V7,-V8) |> 
  rename("ID"=V1, "Status"=V2, "Chr" = V3,"start" = V4,"end" = V5) |> select(-Status, -V10, -V9) |> 
  left_join(index_full |> select(Haplotype,Chr, gstart_skip=start, even_odd)) |> 
  mutate(gstart = gstart_skip + start, gend = gstart_skip + end) |>  
  select(Haplotype, Chr, gstart,  gend, ID, even_odd) |> 
  pivot_longer(gstart:gend, values_to = "x") |> 
  mutate(y=case_when(
    Haplotype=="1"~ 2)) |>    
   arrange(ID, Haplotype) 
head(df1)

df2 <- read.delim("Hap2_busco.tsv", skip =3, header = FALSE, stringsAsFactors = FALSE ) |>
  mutate(Haplotype="2")|> select(-V6, -V7,-V8) |> 
  rename("ID"=V1, "Status"=V2, "Chr" = V3,"start" = V4,"end" = V5) |> select(-Status, -V10, -V9) |>
  left_join(index_full |> select(Haplotype,Chr, gstart_skip=start, even_odd)) |> 
  mutate(gstart = gstart_skip + start, gend = gstart_skip + end) |> 
  select(Haplotype, Chr, gstart,  gend, ID, even_odd) |> 
  pivot_longer(gstart:gend, values_to = "x") |> 
  mutate(y=case_when(
    Haplotype=="2"~ 1)) |>   
  arrange(ID, Haplotype)
head(df2)


full=rbind(df1, df2) |> arrange(ID, Haplotype) 
full |>  head()

 flipscaf=c("Chr5", "Chr6", "Chr7", "Chr9")
 full=full |> left_join(index_full |> select(Chr, mid)) |> 
   group_by(Chr) %>%
   mutate(x = if_else(Chr %in% flipscaf & Haplotype=="2", 2 * mid - x, x))

## Plotting this 
library(ggforce)
library(prismatic)
library(jtools)

genome_width <- .08
skip <- .08

synt = full |> group_by(ID) |> filter(!is.na(x)) |> mutate(n=length(ID)) |> filter(n==4) |> 
  ggplot() +
  geom_diagonal_wide(aes(x = x, y = y, group = ID, color=as.factor(even_odd),
fill=as.factor(even_odd)),
                     orientation = "y", linewidth = .4)  +
  scale_fill_manual(values = c("#4c4763","#cacfd2" )) +
  scale_color_manual(values = c("#4c4763","#cacfd2" )) +
  geom_rect(data = index_full,
            aes(xmin = start, xmax = end,
                ymin = y + (label_sign * skip),
                ymax = y + (label_sign * (skip + genome_width)*0.2)),
            linewidth = .1,
            color = "gray60", fill = "gray80")+
  geom_text(data = index_full,
            aes(x = mid,
                y = y + (label_sign *  (2.2 * skip + genome_width)*0.8),
                label = Chr)) + theme_void() +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
axis.line.x = element_blank(),axis.line.y = element_blank(),
        axis.text = element_blank(), legend.position = "none")
#synt

saveRDS(synt, "snakemake@input[['busco_plot']]")