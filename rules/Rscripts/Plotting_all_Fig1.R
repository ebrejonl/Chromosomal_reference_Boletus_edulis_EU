library(tidyverse)
library(patchwork)
library(jtools)
library(ggforce)
library(ggrastr)
library(ggpubr)
library(ggplot2)
library(prismatic)
library(jtools)
library(here)

synt=readRDS(snakemake@input[['busco_plot']]) 
cov=readRDS(snakemake@input[['Coverage']])+ theme_bw()
te1=readRDS(snakemake@input[['te1']])+ theme(axis.line.x = element_line(linewidth=0.2), 
axis.line.y = element_line(linewidth=0.3)) + scale_x_discrete(labels=c("Longest\n chromosome", "Shortest\n chromosome"))

te2=readRDS(snakemake@input[['te2']])  + theme(axis.line.x = element_line(linewidth=0.3), 
axis.line.y = element_line(linewidth=0.3))+ xlab("")+
  scale_fill_manual(values = c("#697787ff", "#c9400bff"),
labels=c("BolEdBiel_h1","BolEdBiel_h2"))+ theme(legend.position = "right")+  scale_x_discrete(labels=c("Longest\n chromosome", "Shortest\n chromosome"))

te3=readRDS(snakemake@input[['te3']])+ theme(axis.line.x = element_line(linewidth=0.3), 
axis.line.y = element_line(linewidth=0.3), legend.position = "none")+
  xlab("Chromosome Length (Mb)")






# plot
genome_width <- .1
skip <- .08
cov =cov+ theme(axis.text.x.top = element_text(size=11, color="black"),
axis.text.x.bottom = element_text(size=11, color="black"))
t=(te3 | te1 | te2) 
pp=plot_spacer() + cov + plot_spacer() + plot_layout( nrow=1,
  widths = c(0.00001, 0.98, 0.000001))
ppp=plot_spacer() + t+ plot_spacer() + plot_layout(widths = c(0, 1.2, 0.001))
p_full= synt /plot_spacer() / pp/ plot_spacer() /ppp + plot_layout(heights = c(.55,0.05,0.60,0.05,0.65)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")& theme(plot.tag = element_text(size=20))


ggsave(snakemake@output[['Fig1']], p_full, width=12, height = 10)