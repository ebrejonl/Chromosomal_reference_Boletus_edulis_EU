library(tidyverse)
library(ggplot2)
library(circlize)
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

Suillus_index <- read_fai(file="Data/Fasta/Suillus_bovinus.fasta.fai")  %>% 
  mutate(Haplotype = "Suillus")%>%
  rename("seq_id"="name")

Hap_2_indec <- read_fai(file="Data/Fasta/Haplotype2_renamed_reordered_Chr_only.fasta.fai")  %>% 
  mutate(Haplotype = "Edulis")%>%
  rename("seq_id"="name") 

index_full <- rbind(Suillus_index, Hap_2_indec) 


df1 <-  read.delim("snakemake@input[['Suillus_busco']]", skip=3, header=FALSE, stringsAsFactors = FALSE) |> 
  select(-V6, -V7,-V8) |> 
  rename("ID"=V1, "Status"=V2, "qid" = V3,"qs" = V4,"qe" = V5) |> select(-Status, -V10, -V9) 

df2 <- read.delim("snakemake@input[['Haplotype2_busco']]", skip =3, header = FALSE, stringsAsFactors = FALSE ) |>
 select(-V6, -V7,-V8) |> 
    rename("ID"=V1, "Status"=V2, "rid" = V3,"rs" = V4,"re" = V5) |> select(-Status, -V10, -V9) 

df= left_join(df1,df2, by=c("ID")) |> filter(!rid=="") |> filter(!qid=="")

ref_chromosomes <- unique(df$rid)
query_chromosomes <- unique(df$qid)
ref_chromosomes
segment_matrix <- matrix(0, nrow = length(ref_chromosomes), ncol = length(query_chromosomes))
segment_matrix
rownames(segment_matrix) <- ref_chromosomes
colnames(segment_matrix) <- query_chromosomes
for (i in 1:nrow(df)) {
  ref <- df$rid[i]
  query <- df$qid[i]
  start <- df$rs[i]
  end <- df$re[i]
  segment_matrix[ref, query] <- segment_matrix[ref, query] + (end - start)
}
order <- c(paste0("Chr", 1:11))
segment_matrix= segment_matrix[order,]
rownames(segment_matrix) =c(1:11)

colnames(segment_matrix) <- paste0(" ", 1:12)

#grid.col <- c(CM027197.1="darkgrey", CM027198.1="darkgrey", CM027199.1#="darkgrey", CM027200.1="darkgrey", CM027201.1="darkgrey", CM027202.1#="darkgrey", CM027203.1="darkgrey", CM027204.1="darkgrey", CM027205.1#="darkgrey", CM027206.1="darkgrey",
# C2="#fb9a99",C1="#1f78b4", C3="#fdbf6f", C4="#6a3d9a", C5="#a6cee3", 
#  C7="#C89F9C", C8="#132E32", C9="#A882DD", C10="#FF312E", C11="black" )

library(RColorBrewer)
library(MetBrewer)
circos.clear()
circos.par(start.degree = 180) #flip horizontally

grid.col <- c(met.brewer(n=11, "Johnson", direction = -1),rep("darkgrey",12)) # Johnson
pdf("snakemake@input[['Circle_plot']]")
chordDiagram(segment_matrix, grid.col = grid.col, 
             big.gap = 18,
             annotationTrack = "grid", 
             annotationTrackHeight = 0.1, 
             preAllocateTracks = 1, 
             transparency = 0.4) 
panel_fun <- function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")  
  label_y <- (ylim[1] + ylim[2]) / 2  
  circos.text(x = mean(xlim), y = label_y,
              labels = sector.name,
               facing =  "outside",
              niceFacing = TRUE,
              adj = c(-0.2, 0.5),  
              cex = 0.2)  
}
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", 
              niceFacing = TRUE, adj = c(0.5, 0), col= "black")
}, bg.border = NA)
dev.off()
