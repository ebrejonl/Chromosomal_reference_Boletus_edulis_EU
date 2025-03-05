## Karyotype plot for figure 2
library(tidyverse)
library(data.table)
library(here)

dummy_files <- snakemake@input[['dummy']]

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

# MAT loci
# 
MAT <- data.frame(Start=c(953944,3072277 ), End=c(980720,3085748), g_end=c(32057100,4111618), g_start=c(29387091, 0),# add chr data here before merge
            g_mid=c(30722096,2055809)  ,   type=c("MATa", "MATb"), name=c("Chr9","Chr1"))

# Centromeres
Centromeres=data.frame(name=paste0("Chr", rep(1:11)), 
                       Start=c(1100000,1800000,2750000,2200000,
                               2050000,1150000,1350000,1400000,1250000,2400000,1600000), 
                       End=c(1100000,1800000,2750000,2200000,
                                             2050000,1150000,1350000,1400000,1250000,2400000,1600000),
                       type="Centromere")

Telomeres=rbind(Telomeres, Centromeres)


first_merge <- merge(Cazymes, Effectors,  by=c("name", "Start", "End", "type"), all = T) 
fai_deeper <- left_join(Hap_2_indec, first_merge, by="name") %>% filter(!name=="NA") 

fai_smaller=fai_deeper |> select(name,length,g_start,g_end,g_mid) |> distinct()
Telomeres = Telomeres |> left_join(fai_smaller, by ="name") |> mutate(idx=NA, even_odd=NA, Caz=NA, Effector=NA)

missing_columns <- setdiff(names(fai_deeper), names(MAT))
MAT[missing_columns] <- NA
MAT <- MAT[names(fai_deeper)]

fai_deeper <- rbind(MAT, fai_deeper)
fai_deeper$type <- factor(fai_deeper$type, levels = c("CAZymes", "Effector", "MATa", "MATb"))
fai_deeper_eff <- fai_deeper %>% filter(type=="Effector")
fai_deeperM <- fai_deeper 
order=paste0("Chr", 1:11)
fai_deeper = fai_deeper |> mutate(name=fct_relevel(name, order))
fai_deeper = rbind(fai_deeper, Telomeres[names(fai_deeper)]) 

fai_deeper |> head()


### The features first 
data_features <-  fai_deeper # already feature format ?

centro_length <- 150000
chr_width <- .6
clr <- "#E5E5E5"
#clr_h <- "gray60"
clr_h <- "gray50"


karyo2= fai_deeper |> 
  filter(type=="Centromere") |> 
  mutate( centro_pos=Start,
          centro_start=Start - .5 *centro_length,
          centro_end=Start +.5 * centro_length) |> 
  mutate(name=as.factor(name))


karyo_poly2 <- karyo2 |>
  mutate(
    x = map(as.numeric(name), \(x, w = chr_width * .5) {
      rep(x + w * c(1, -1, -1, 1, 1), 2)
    }),
    y = map2(length, centro_pos, \(l, C, cl = centro_length * .5) {
      c(0, 0, C - cl, C - cl, 0, C + cl, C + cl, l, l, C + cl)
    }),
    y_sharp = map2(length, centro_pos, \(l, C, cl = centro_length * .5) {
      lw = (C - cl) * .5
      up = C + (l - C + cl) * .5
      c(lw, lw, C - cl, C - cl, lw, C + cl, C + cl, up, up, C + cl)
    }),
    part = list(rep(c("lower", "upper"), each = 5))
  ) |>
  unnest(cols = c(x, y, y_sharp, part))|> as.data.frame()


karyo_hour2 <- karyo2 |>
  mutate(
    x = map(as.numeric(name), \(x, w = chr_width * .5) {
      x + w * c(1, -1, 1, -1, 1)
    }),
    y = map2(length, centro_pos, \(l, C, cl = centro_length * .5) {
      c(C - cl, C - cl, C + cl, C + cl, C - cl)
    })
  ) |>
  unnest(cols = c(x, y)) |> as.data.frame()

karyo_hour2

Feature=c("CAZymes",    "Effector")


karyo2.0=karyo2 |>
  ggplot(aes(x = name)) +
  geom_polygon(
    data = karyo_poly2,
    aes(x = x, y = y_sharp, group = str_c(name, part)),
    fill = clr
  ) +
  geom_shape(
    data = karyo_poly2,
    aes(x = x, y = y, group = str_c(name, part)),
    r = unit(20, "pt"),
    fill = clr
  ) +
  geom_polygon(
    data = karyo_hour2,
    aes(x = x, y = y, group = name),
    fill = clr_h
  )+
  geom_linerange(
    data = fai_deeper |> filter(type %in% Feature),
    inherit.aes = FALSE,
    aes(
      xmin = as.numeric(name) - chr_width * .41,
      xmax = as.numeric(name) + chr_width * .41,
      y = Start,
      color = type
    ),
    linewidth = 1
  )  +    
  scale_color_manual(aesthetics = "color", values =c("#ff4c4c","#3777FF"))+

  geom_label(fai_deeper |> filter(type=="MATb" | type=="MATa"), 
      inherit.aes = FALSE,
      mapping = aes(x=as.numeric(name), y = End, label=type),  
      color="#86a995", size=3.5 , fontface="bold", position = "identity")+

  geom_point(fai_deeper %>% filter(type=="Telomere" & Start>230000), 
  mapping = aes(x  = as.numeric(name), y  = Start),
  shape=24, fill="yellow",color="brown", size=4)+

  geom_point(fai_deeper %>% filter(type=="Telomere" & Start<250000), 
  mapping = aes(x  = as.numeric(name), y  = Start),
  shape=25, fill="yellow",color="brown", size=4)+

  scale_x_continuous(
    breaks = as.numeric(unique(fai_deeper$name)),
    labels = as.character(substr(unique(fai_deeper$name), 4, 5))
  ) +
      
    theme_minimal()+ xlab(NULL)+ylab("Genomic position")+ theme(axis.title.y.left = element_text(color="black", size=17, vjust=2.5),
  axis.text.y.left = element_text(color="black", size=15),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.text.x.bottom = element_text(size=18, color="black",face = "plain"),
legend.title = element_blank(), legend.text = element_text(size=18))#;karyo

ggsave("snakemake@output[['Karyo_plot']]", karyo, height = 9, width =14)
