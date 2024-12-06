## Karyotype plot for figure 2
library(tidyverse)
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

Hap_2_indec <- read_fai(file="snakemake@input[['hap2_index']]")  


Cazymes <- fread("snakemake@input[['Caz_coor']]", col.names = c("Caz", "name", "Start", "End")) %>%
  filter(!name=="h2tg000050l",!name=="h2tg000083l") %>%
  mutate(type= "CAZymes")
                          

Cazymes |> group_by(name) |> summarise(n())
nrow(fread("Data/Caz_with_coordinates.tsv", col.names = c("Caz", "name", "Start", "End")))


Effectors <- fread("Data/Effector_with_coordinates.tsv", col.names = c("Effector", "name", "Start", "End")) %>%
  filter(!name=="h2tg000050l",!name=="h2tg000083l") %>%
  mutate(type="Effector")

Effectors|> group_by(name) |> summarise(n())

 
MAT <- data.frame(Start=c(953944,3072277 ), End=c(980720,3085748), g_end=c(32057100,4111618), g_start=c(29387091, 0),# add chr data here before merge
            g_mid=c(30722096,2055809)  ,   type=c("MATa", "MATb"), name=c("Chr9","Chr1"))

missing_columns <- setdiff(names(fai_deeper), names(MAT))
MAT[missing_columns] <- NA

MAT <- MAT[names(fai_deeper)]
 
# Telomeres 
Telomeres=read_tsv("Results/Telomeres/TELO_telomeric_repeat_windows.tsv") |> rename (name="id") |> filter(forward_repeat_number > 24 | reverse_repeat_number >24) |> filter(!forward_repeat_number==27 ) |> 
  rename(Start="window") |> mutate(End = Start + 10000) |> mutate(type="Telomere") |> 
  select(name, Start, End, type)

fai_smaller=fai_deeper |> select(name,length,g_start,g_end,g_mid) |> distinct()
Telomeres = Telomeres |> left_join(fai_smaller, by ="name") |> mutate(idx=NA, even_odd=NA, Caz=NA, Effector=NA)


first_merge <- merge(Cazymes, Effectors,  by=c("name", "Start", "End", "type"), all = T) 
fai_deeper <- left_join(Hap_2_indec, first_merge, by="name") %>% filter(!name=="NA") 
fai_deeper <- rbind(MAT, fai_deeper)
fai_deeper$type <- factor(fai_deeper$type, levels = c("CAZymes", "Effector", "MATa", "MATb"))
fai_deeper_eff <- fai_deeper %>% filter(type=="Effector")
fai_deeperM <- fai_deeper 
#fai_deeper <- fai_deeper %>% filter(!type=="MATa", !type=="MATb")
order=paste0("Chr", 1:11)
fai_deeper = fai_deeper |> mutate(name=fct_relevel(name, order))
fai_deeper = rbind(fai_deeper, Telomeres[names(fai_deeper)]) 
library(ggchicklet)
library(ggrepel)
library(MetBrewer)
library(data.table)

#met.brewer(n=11, "Johnson")
karyo <-fai_deeper |> mutate(name=fct_relevel(name, order)) |> 
ggplot(., mapping=aes(ymin = g_start-g_start+1, ymax = g_end-g_start+1, 
                               xmin = 0, xmax = 0.2))+#, col=as.factor(type)))+#,col="black", fill="#85929e") +
  ggchicklet:::geom_rrect(#col="black",
                          fill="#e9eaea", cex=0.01,##cccccc", 
                          r = unit(0.4, 'npc')) +
  geom_segment(fai_deeper |> filter(!type=="MATa", !type=="MATb", !type=="Telomere") ,
               mapping = aes(x = 0.1, xend = 0.1, y = Start-5000, yend = End+5000, color=type ), #color=type),
               linewidth=12)+
  xlab(NULL)+ylab(NULL)+
  scale_color_manual(aesthetics = "color", values =c("#03286b","#ae0000") ,
                     labels=c( "CAZymes", "Effector\nProtein", "MATa", "MATb"))+ 
  geom_segment(fai_deeper_eff, mapping = aes(x = 0.1, xend = 0.1, y = Start-1000, yend = End+1000),  color="#ae0000",
               linewidth=12)+
  geom_label_repel(fai_deeper |> filter(type=="MATb" | type=="MATa"), mapping = aes(x = 0.1, y = End, label=type),  color="#86a995",
              size=3.3 , fontface="bold")+
  geom_point(fai_deeper %>% filter(type=="Telomere" & Start>230000), mapping = aes(x  = 0.1, y  = Start),shape=24, fill="yellow",color="brown", size=4)+
  geom_point(fai_deeper %>% filter(type=="Telomere" & Start<250000), mapping = aes(x  = 0.1, y  = Start),shape=25, fill="yellow",color="brown", size=4)+
  theme_light()+
  facet_wrap(~as.numeric(substr(name,4,5)), nrow = 1, strip.position = "bottom") + # ncol=6
   # annotate("text", x=3, y=200000, label="}")+
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_blank(), panel.spacing = unit(1.5, "cm"), panel.border = element_blank(),
        axis.ticks.x = element_blank() , panel.grid.minor = element_blank(),
        strip.text = element_text(size = 20, color = "black", hjust = 0.14), axis.ticks.y.left = element_blank(),
        # strip.background = element_rect(fill="#ecf0f3"),
        strip.background = element_blank(), 
        legend.title = element_blank(),
        legend.key.size = unit(1,"line"),
        legend.key.height = unit(3,"line"),
        legend.text = element_text(size=20))  #;karyo

ggsave("snakemake@output[['Karyo_plot']]", karyo, height = 9, width =14)
