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
p=all_hic_data |> filter(log_counts>10) |> 
  group_by(chromosome1) |> 
  mutate(roll_c=rollmean(log_counts, k = 45, align="center", na.pad = TRUE))

p |> ggplot()+
  geom_point(aes(x = y, y = log_counts))+
  facet_wrap(~chromosome1)+
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5)
