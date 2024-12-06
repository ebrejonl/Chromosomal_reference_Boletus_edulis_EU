library(tidyverse)
library(purrr)
library(here)


Central= read_tsv("PI_Central.windowed.pi", col_names = T) |> mutate(pop = "Central")
Fenno=read_tsv("PI_Fennoscandia.windowed.pi", col_names =T ) |> mutate(pop ="Fennoscandia")
Great_Britain= read_tsv("PI_Great_Britain.windowed.pi",col_names = T ) |> mutate(pop ="United Kingdom")
Iceland =read_tsv("PI_Iceland.windowed.pi", col_names  =T )  |> mutate(pop ="Iceland")

full=rbind(Central,Fenno, Great_Britain, Iceland )

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

Hap_2_indec <- read_fai(file="../Haplotype2_renamed_reordered_Chr_only.fasta.fai")  %>% 
  mutate(Haplotype = "2")%>%
  rename("CHROM"="name") 

fullfai= full |> left_join(Hap_2_indec, by = "CHROM") 
fullfai= fullfai |> mutate(BIN_MID=((BIN_END+BIN_START)/2)+g_start)

library(ggrastr)
library(ggnewscale)



PI_plot=ggplot(fullfai) +
  geom_rect(data=Hap_2_indec, mapping=aes(xmin=g_start, xmax=g_end, ymin=-Inf,ymax=Inf, 
    fill=factor(as.character(even_odd))), color="transparent") + 
  scale_fill_manual(values=c(`0`="white", `1`="#e5e8e8", `2`="transparent"),guide="none")+
  scale_x_continuous(labels = function(x){sprintf("%.1fMbp",x*1e-6)}, 
                     sec.axis = sec_axis(trans = identity, 
                                         breaks=Hap_2_indec$g_mid[1:11],
                                         labels=Hap_2_indec$CHROM[1:11]))+
  theme(axis.text.x.top = element_text(vjust = 0.5))+
    new_scale_fill() +
    geom_point(data = fullfai, aes(x=BIN_MID, y=PI, fill=pop, color=pop), shape=21, size=1.7,
     stroke=0.1, alpha=0.45)+ 
      scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
    "Central"="#f6c200","United Kingdom"="#7ec68f"))+
    scale_color_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
    "Central"="#f6c200","United Kingdom"="#7ec68f"))+
  facet_wrap(~factor(pop, levels=c("Iceland","Fennoscandia","United Kingdom","Central")),  ncol=1, strip.position = "right")+ theme_bw()+ theme(legend.position = "none")+
 xlab("Genomic position")+ ylab("") + 
  theme(axis.text.y.left = element_blank(),
axis.ticks.y.left = element_blank()) + theme(axis.text.x.top = element_text(size=11, color="black"),
axis.text.x.bottom = element_text(size=11, color="black"))

PI_plot

saveRDS(file="PI_plot.rds", PI_plot)

ggsave(PI_plot, file="Pi_plot.pdf", width = 8, height = 5)


library(patchwork)
library(jtools)
library(Cairo)
box=ggplot(fullfai) +
  geom_violin(aes(x=1, y=PI ,fill=pop), alpha=0.2, width = 0.3)+
    geom_boxplot(aes(x=1, y=PI, fill=pop), outlier.shape = NA , width=0.1) +
      scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
    "Central"="#f6c200","United Kingdom"="#7ec68f"))+
  facet_wrap(~factor(pop, levels=c("Iceland","Fennoscandia","United Kingdom","Central")),  ncol=1)+ 
  theme_apa()+ xlab(NULL)+ ylab("π")+
  theme(axis.text.x.bottom = element_blank(),axis.ticks.x.bottom = element_blank(),
axis.title.y = element_text(size = 30, vjust=.5,angle = 0, hjust=.8, family = "serif"),
strip.text.x.top = element_blank(), legend.position = "none")+theme(axis.text.y.left = element_text(size=11))
box


PImean=mean(fullfai$PI) 
PI_plotb=ggplot(fullfai) +
  geom_rect(data=Hap_2_indec, mapping=aes(xmin=g_start, xmax=g_end, ymin=-Inf,ymax=Inf, 
    fill=factor(as.character(even_odd))), color="transparent") + 
  scale_fill_manual(values=c(`0`="white", `1`="#e5e8e8", `2`="transparent"),guide="none")+
  scale_x_continuous(labels = function(x){sprintf("%.1fmb",x*1e-6)}, 
                     sec.axis = sec_axis(trans = identity, 
                                         breaks=Hap_2_indec$g_mid[1:11],
                                         # labels=str_c("scf_",1:11)
                                         labels=Hap_2_indec$CHROM[1:11]))+
  theme(axis.text.x.top = element_text(vjust = 0.5))+
    new_scale_fill() +
    geom_point(data = fullfai, aes(x=BIN_MID, y=PI, fill=pop, color=pop), shape=21, size=1.7,
     stroke=0.1, alpha=0.45)+ 
      scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
    "Central"="#f6c200","United Kingdom"="#7ec68f"))+
 scale_color_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
 "Central"="#f6c200","United Kingdom"="#7ec68f"))+
  facet_wrap(~factor(pop, levels=c("Iceland","Fennoscandia","United Kingdom","Central")), 
  ncol=1, strip.position = "right")+ theme_bw()+ 
  theme(legend.position = "none", axis.title.y = element_blank(),
    panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank(),
  axis.ticks.y.left = element_blank(),
  axis.text.y.left = element_blank())+ xlim(0,40000000)+
  geom_hline(yintercept = PImean, color="white", size =.3, linetype="dashed")+
 xlab("Genomic position")+ ylab("PI")

pi_full=box + PI_plot + plot_layout(widths=c(.1,.9))
pi_full
saveRDS(pi_full, file="pi_full.rds")
ggsave(pi_full,filename = "pi_full.pdf", device=cairo_pdf, width = 10, height = 10)

library(maps)
library(sf)
library(rnaturalearth)
library(MetBrewer)
library(ggforce)

#fae48b (light yellow)
#f2c53d (gold)
#df8a13 (orange)
#cf4e1b (red-orange)
#993232 (brownish-red)
#3b2a1e (dark brown)
#1e1e1e (almost black)

ghost=read.csv()
mydata=readRDS("Results/STRUCTURE/Whole_genome_pca_with_pop.rds") |> 
#mydata=readRDS("../MAKE_ALL_FIGURES/Whole_genome_pca_with_pop_Dec4.rds") |> 
  mutate(Pop=case_when(
    Pop=="Scandinavia" ~ "Fennoscandia", 
    Pop=="South" ~ "Southern EU",
    Pop=="USA"~"North America",
    Pop=="GB"~"United Kingdom", TRUE~as.character(Pop)))

p1=ggplot(data = mydata, aes(x=mydata[,1],y=mydata[,2], fill=Pop)) + 
      theme_apa() +
        geom_vline(xintercept = 0, lty = 3 , alpha=0.6) +
          geom_hline(yintercept = 0,lty = 3,alpha=0.6 )+
            scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
          "Central"="#f6c200","United Kingdom"="#7ec68f", "North America"="#df8a13",
          "Russia"="#d04e00", 
          "Southern EU"="#a00e00")) +
      geom_point(size=4, shape=21, color="black", alpha=0.5)  +
geom_mark_ellipse(mapping=aes(x=mydata[,1],y=mydata[,2], fill=Pop))+
xlab("PC1 5%")+ ylab("PC2 3.4%") +theme_classic()+ 
  theme(legend.position = "top", legend.direction = "horizontal",
legend.title = element_blank(), legend.text = element_text(size=13),
axis.text.x.bottom =  element_text(size=11, color="black"), 
axis.text.y.left = element_text(size=11, color="black"))+
  guides(fill = guide_legend(nrow = 1))

p1
world <- map_data("world")
 worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                          returnclass = 'sf')
p2=ggplot(data = worldmap) + geom_sf(data = worldmap) +
    geom_sf(data = worldmap |> filter(name_fr=="Islande"), fill="#132b69", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Suède"), fill="#0086a8"  , alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Norvège"), fill="#0086a8", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Finlande"), fill="#0086a8", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Estonie"), fill="#0086a8", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Danemark"), fill="#0086a8", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Allemagne"), fill="#f6c200", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Pologne"), fill="#f6c200", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Angleterre"), fill="#7ec68f", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Irlande"), fill="#7ec68f", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_en=="Wales"), fill="#7ec68f", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_en=="Northern Ireland"), fill="#7ec68f", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Écosse"), fill="#7ec68f", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Russie"), fill="#d04e00", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Bulgarie"), fill="#a00e00", alpha=0.4)+
    geom_sf(data = worldmap |> filter(name_fr=="Russie"), fill="#d04e00", alpha=0.4)+
    coord_sf(xlim = c(-25, 50), ylim = c(40, 73), expand = FALSE) +
     geom_jitter(data=mydata |> mutate(Nothing=""), aes(x=lon, y=lat),
   shape=21, size=4, alpha=0.7, color="white", fill="#1b2631",width = 0.5, height = 0.2)+ 
        theme_classic()+
    theme(legend.position = "none", axis.text.x.bottom = element_blank(),
  axis.text.y.left = element_blank(), axis.ticks.x.bottom = element_blank(), 
  axis.ticks.y.left = element_blank())+
  xlab("")+ylab("")

  pp=p2 + p1 + plot_layout(widths = c(0.55, 0.45))
  ppp= guide_area() +pp  + 
    plot_layout(guides = "collect", heights = c(.1,.5))
  pppp= ppp/ pi_full + 
    plot_layout(guides = "collect", heights = c(.5,.6))

pppp

ggsave(pppp, file="Figure4_4Dec.pdf",device = cairo_pdf, width = 11, height = 12)