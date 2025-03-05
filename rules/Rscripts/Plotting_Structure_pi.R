library(tidyverse)
library(purrr)
library(here)

# Snakemake version 
inputs <- snakemake@input[['pi_files']]
data_list <- lapply(inputs$pi_files, function(file) {
  pop_name <- str_extract(file, "PI_(.*?)\\.windowed\\.pi") %>%
    str_replace("PI_", "")
    read_tsv(file, col_names = TRUE) %>%
    mutate(pop = pop_name)
})
full<- bind_rows(data_list)

## read PI_files
Central= read_tsv("24Feb_old_filters/PI_Central.windowed.pi", col_names = T) |> mutate(pop = "Central Europe")
Fenno=read_tsv("24Feb_old_filters/PI_Fennoscandia.windowed.pi", col_names =T ) |> mutate(pop ="Fennoscandia")
Great_Britain= read_tsv("24Feb_old_filters/PI_Great_Britain.windowed.pi",col_names = T ) |> mutate(pop ="United Kingdom")
Iceland =read_tsv("24Feb_old_filters/PI_Iceland.windowed.pi", col_names  =T )  |> mutate(pop ="Iceland")

full=rbind(Central,Fenno, Great_Britain, Iceland )


# reading in the fai for haplotype 2
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

# merging 
fullfai= full |> left_join(Hap_2_indec, by = "CHROM") 
fullfai= fullfai |> mutate(BIN_MID=((BIN_END+BIN_START)/2)+g_start)

library(ggrastr)
library(ggnewscale)


PImean=mean(fullfai$PI) # 0.008343298

PI_plot=ggplot(fullfai) +
  geom_rect(data=Hap_2_indec, mapping=aes(xmin=g_start, xmax=g_end, ymin=-Inf,ymax=Inf, 
    fill=factor(as.character(even_odd))), color="transparent") + 
  scale_fill_manual(values=c(`0`="white", `1`="#e5e8e8", `2`="transparent"),guide="none")+
  scale_x_continuous(labels = function(x){sprintf("%.0f",x*1e-6)}, 
                     sec.axis = sec_axis(trans = identity, 
                                         breaks=Hap_2_indec$g_mid[1:11],
                                         # labels=str_c("scf_",1:11)
                                         labels=Hap_2_indec$CHROM[1:11]), expand = c(0, 0))+
  theme(axis.text.x.top = element_text(vjust = 0.5))+
    new_scale_fill() +
    geom_point(data = fullfai, aes(x=BIN_MID, y=PI, fill=pop, color=pop), shape=21, size=1.7,
     stroke=0.1, alpha=0.45)+ 
      scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
    "Central Europe"="#f6c200","United Kingdom"="#7ec68f"))+
    scale_color_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
    "Central Europe"="#f6c200","United Kingdom"="#7ec68f"))+
  facet_wrap(~factor(pop, levels=c("Iceland","Fennoscandia","United Kingdom","Central Europe")),  ncol=1, strip.position = "right")+ theme_bw()+ theme(legend.position = "none")+
 xlab("Genomic position (Mb)")+ ylab("") + 
  geom_hline(yintercept = PImean, color="white", size =.5, linetype="dashed")+
  theme(axis.text.y.left = element_blank(),
axis.ticks.y.left = element_blank()) + theme(axis.text.x.top = element_text(size=11, color="black"),
axis.text.x.bottom = element_text(size=11, color="black"))

PI_plot

saveRDS(file="PI_plot24Feb.rds", PI_plot)



library(patchwork)
library(jtools)
library(Cairo)
box=ggplot(fullfai) +
  geom_violin(aes(x=1, y=PI ,fill=pop), alpha=0.2, width = 0.3)+
    geom_boxplot(aes(x=1, y=PI, fill=pop), outlier.shape = NA , width=0.1) +
      scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
    "Central Europe"="#f6c200","United Kingdom"="#7ec68f"))+
  facet_wrap(~factor(pop, levels=c("Iceland","Fennoscandia","United Kingdom","Central")),  ncol=1)+ 
  theme_apa()+ xlab(NULL)+ ylab("Nucleotide Diversity (π)")+
  theme(axis.text.x.bottom = element_blank(),axis.ticks.x.bottom = element_blank(),
#axis.title.y = element_text(size = 30, vjust=.5,angle = 0, hjust=.8, family = "serif"),
strip.text.x.top = element_blank(), legend.position = "none")+theme(axis.text.y.left = element_text(size=11))
box


PImean=mean(fullfai$PI) # 0.00789434 # old pi mean
PImean=mean(fullfai$PI) # 0.008343298


pi_full=box + PI_plot + plot_layout(widths=c(.1,.9))
pi_full


library(maps)
library(sf)
library(rnaturalearth)
library(MetBrewer)
library(ggforce)


mydata=readRDS("../Structure/Whole_genome_pca_with_pop24Feb.rds") |> 
  mutate(Pop=case_when(
    Pop=="Central"~ "Central Europe",
    Pop=="Scandinavia" ~ "Fennoscandia", 
    Pop=="Southern EU" ~ "Southern Europe",
    Pop=="USA"~"North America",
    Pop=="Great Britain"~"United Kingdom", TRUE~as.character(Pop)))

p1=ggplot(data = mydata, aes(x=mydata[,1],y=mydata[,2], fill=Pop)) + 
      theme_apa() +
        geom_vline(xintercept = 0, lty = 3 , alpha=0.6) +
          geom_hline(yintercept = 0,lty = 3,alpha=0.6 )+
            scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
          "Central Europe"="#f6c200","United Kingdom"="#7ec68f", "North America"="#df8a13",
          "Russia"="#d04e00", 
          "Southern Europe"="#a00e00")) +
      geom_point(size=4, shape=21, color="black", alpha=0.5)  +
#geom_mark_ellipse(mapping=aes(x=mydata[,1],y=mydata[,2], fill=Pop))+
xlab("PC1 4.9 %")+ ylab("PC2 3.5 %") +theme_classic()+ 
  theme(legend.position = "top", legend.direction = "horizontal",
legend.title = element_blank(), legend.text = element_text(size=13),
axis.text.x.bottom =  element_text(size=11, color="black"), 
axis.text.y.left = element_text(size=11, color="black"))+
  guides(fill = guide_legend(nrow = 1))+theme(panel.border = element_rect(linewidth =0.75, fill="transparent"), 
axis.line.x.bottom = element_blank(),
axis.line.y.left = element_blank())

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
   shape=21, size=2.5, alpha=0.7, color="white", fill="black",width = 0.5, height = 0.2)+ 
      #scale_fill_manual(values = c("Iceland"="#132b69","Fennoscandia"="#0086a8",
      #"Central"="#f6c200","Great Britain"="#7ec68f")) +
        theme_classic()+
    theme(legend.position = "none", axis.text.x.bottom = element_blank(),
  axis.text.y.left = element_blank(), axis.ticks.x.bottom = element_blank(), 
  axis.ticks.y.left = element_blank(), panel.border = element_rect(linewidth =0.75, fill="transparent"), 
  axis.line.x.bottom = element_blank(),
  axis.line.y.left = element_blank())+
                             #xlab("Longitude")+ylab("Latitude")
  xlab("")+ylab("")

  
p2
pp=p2 + p1 + plot_layout(widths = c(0.55, 0.45))
  ppp= guide_area() +pp  + 
    plot_layout(guides = "collect", heights = c(.1,.5))
  pppp= ppp/ pi_full + 
    plot_layout(guides = "collect", heights = c(.5,.6))
pppp

ggsave(pppp, file="Figure4.pdf",device = cairo_pdf, width = 11, height = 12)