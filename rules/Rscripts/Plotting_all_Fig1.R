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
axis.line.y = element_line(linewidth=0.3))
te2=readRDS(snakemake@input[['te2']]) + theme(axis.line.x = element_line(linewidth=0.3), 
axis.line.y = element_line(linewidth=0.3))+ xlab("")
te3=readRDS(snakemake@input[['te3']])+ theme(axis.line.x = element_line(linewidth=0.3), 
axis.line.y = element_line(linewidth=0.3), legend.position = "none")+
  xlab("Chromosome Length (Mbp)")


# plot
genome_width <- .1
skip <- .08
cov =cov+ theme(axis.text.x.top = element_text(size=11, color="black"),
axis.text.x.bottom = element_text(size=11, color="black"))
t=(te3 | te1 | te2) 
pp=plot_spacer() + cov + plot_spacer() + plot_layout(widths = c(0.01, 1, 0.01))
ppp=plot_spacer() + t+ plot_spacer() + plot_layout(widths = c(0.005, 1, 0.005))
p_full= synt /plot_spacer() / pp/ plot_spacer() /ppp + plot_layout(heights = c(.5,0.05,0.6,0.05,0.6)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")

ggsave(snakemake@output[['Fig1']], p_full, width=12, height = 10)