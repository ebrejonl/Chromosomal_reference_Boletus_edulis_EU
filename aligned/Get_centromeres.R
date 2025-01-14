library(strawr)
library(tidyverse)
library(ggrastr)

gc()
chromosomes <- 1:11  
all_hic_data <- data.frame()

for (chr1 in chromosomes) {
  for (chr2 in chromosomes) {
    # Extract Hi-C data for the chromosome pair
    hic_df <- strawr::straw(
      "KR",                       
      "aligned/inter.hic",         
      as.character(chr1),         
      as.character(chr2),          
      "BP",25000)
    
    # Add columns for chromosome pair
    hic_df$chromosome1 <- paste0("chr", chr1)
    hic_df$chromosome2 <- paste0("chr", chr2)
    
    # Combine with the main data frame
    all_hic_data <- bind_rows(all_hic_data, hic_df)
  }
}

all_hic_data$chromosome1 <- factor(all_hic_data$chromosome1, levels = paste0("chr", 1:11))
all_hic_data$chromosome2 <- factor(all_hic_data$chromosome2, levels = rev(paste0("chr", 1:11)))
# 
# Add log-transformed counts column (log10)
all_hic_data <- all_hic_data %>%
mutate(log_counts = log2(counts + 1)) |> 
    filter(as.numeric(str_sub(chromosome1,4,5)) <= as.numeric(str_sub(chromosome2,4,5)) )

all_hic_data |> filter(log_counts>10) |> 
  ggplot(aes(x=y, y=log_counts ))+
  geom_point()+ 
  facet_wrap(~chromosome1)




library(zoo)
p=all_hic_data |> filter(log_counts>10) |> ## keeping only the highest values
  group_by(chromosome1) |> 
  mutate(roll_c=rollmean(log_counts, k = 40, align="center", na.pad = TRUE))

p |> ggplot()+
  geom_point(aes(x = y, y = log_counts))+
  facet_wrap(~chromosome1)+
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5)

# per chromosome
p |> filter(chromosome1=="chr1") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts), size=2.5) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1100000)

p |> filter(chromosome1=="chr2") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1800000)

p |> filter(chromosome1=="chr3") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 2750000)

p |> filter(chromosome1=="chr4") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 2200000)

p |> filter(chromosome1=="chr5") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 2050000)

p |> filter(chromosome1=="chr6") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1150000)

p |> filter(chromosome1=="chr7") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1350000)

p |> filter(chromosome1=="chr8") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1400000)

p |> filter(chromosome1=="chr9") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1250000)

p |> filter(chromosome1=="chr10") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
  geom_vline(xintercept = 2400000)

p |> filter(chromosome1=="chr11") |> 
  ggplot()+ 
  geom_point(aes(x = y, y = log_counts)) +
  geom_vline(xintercept = 1600000)

ql=unique(p$chromosome1)  |> as.data.frame()
ql$centromere=c(1100000,1800000,2750000,2200000,2050000,1150000,1350000,1400000,1250000,2400000,1600000)
colnames(ql) <- c("chromosome1", "Centromere")

p =  p   |> left_join(ql)
supp_centromeres=p |> 
  ggplot(aes(x=y/1000000, y=log_counts ))+
  geom_point()+ xlab("Position (Mbp)")+ylab("Hi-C interaction >10")+
  geom_segment(aes(x=Centromere/1000000, y=-Inf, xend =Centromere/1000000, yend = Inf), color="red")+
  facet_wrap(~chromosome1, scales = "fixed", ncol=5)

supp_centromeres

ggsave(supp_centromeres, file = "Supplementary_figure_centromeres.pdf", width = 12, height = 8)
