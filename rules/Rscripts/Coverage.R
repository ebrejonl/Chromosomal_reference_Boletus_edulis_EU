library(tidyverse)
library(plyranges)
library(data.table)
library(here)
gc()
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


fai <- read_fai(file="Data/Fasta/Haplotype1.fasta.fai")


window_step <- 50000
window_size <- 500000

Coverage <-  read.table(snakemake@input[['dp_h1.tsv']], header = TRUE, sep = "\t") %>%
    mutate(DEPTH = as.numeric(MEAN_DEPTH)) %>%
    select(-MEAN_DEPTH) %>%
    as_tibble() %>%
    {GRanges(seqnames = .$CHROM, ranges = IRanges(start = .$POS, end = .$POS), cov = .$DEPTH)} %>%
    as_tibble() %>%
    group_by(seqnames) %>%
    mutate(window = floor((start - 1) / (window_size - window_step))) %>%
    group_by(seqnames, window) %>%
    summarise(start = min(start),
              end = max(start + window_size - 1),
              meancov = mean(cov, na.rm = TRUE),
              sdcov=sd(cov, na.rm = TRUE),
              .groups = 'drop') %>%
    mutate(midpoint = (start + end) / 2) %>%
   select(CHROM = seqnames, start, end, meancov, sdcov, midpoint) %>%
   left_join(fai, by = "CHROM") %>%  
    filter(!is.na(CHROM)) %>%  
    mutate(POS_absolute = midpoint + g_start)


Coverage |> filter(str_detect(CHROM, "Chr")) |> summarize(m=mean(meancov,na.rm=T)) # mean coverage chromosomes
Coverage |> filter(!str_detect(CHROM, "Chr")) |> summarize(m=mean(meancov,na.rm=T))# mean coverage unscaffolded contigs


library(ggrastr)
library(ggpubr)

ppp_legend <- ggplot()+
  geom_rect(data=fai, mapping=aes(xmin=g_start, xmax=g_end, ymin=-Inf,ymax=Inf, 
                                         fill=factor(as.character(even_odd))), color="transparent") + 
  scale_fill_manual(values=c(`0`="white", `1`="#cacfd2", `2`="transparent"),guide="none")+
  scale_x_continuous(labels = function(x){sprintf("%.0f",x*1e-6)}, 
                     sec.axis = sec_axis(trans = identity, 
                                         breaks=fai$g_mid[1:11],
                                         # labels=str_c("scf_",1:11)
                                         labels=fai$CHROM[1:11]), expand = c(0, 0))+
  theme(axis.text.x.top = element_text(vjust = 0.5))+
    new_scale_fill() +
    #angle=90, vjust = 0.5))+
  rasterize(
    geom_ribbon(Coverage, mapping = aes(x=POS_absolute, 
      ymin = meancov-sdcov, 
      ymax = meancov+sdcov,fill="Standart\n deviation"), alpha=.8))+
      geom_line(data = Coverage, aes(x = POS_absolute, y = meancov, color = "Mean"), linewidth = 0.3) +
      # Custom color and fill legends
      scale_color_manual(name = "Coverage", values = c("Mean" = "#f8c471")) +
      scale_fill_manual(name = "", values = c("Standart\n deviation" = "#34495e")) +
    theme_bw()+ xlab("Genomic position (Mb)")+ ylab("Depth of coverage") +
    ylim(-5,60)+ geom_bracket(xmin = 41000000, xmax = 50000000, y.position = 50,
          vjust = -0.5,  label.size = 3.5,tip.length = 0.025, label="Unplaced scaffolds")

saveRDS(ppp_legend, snakemake@output[['Coverage']])