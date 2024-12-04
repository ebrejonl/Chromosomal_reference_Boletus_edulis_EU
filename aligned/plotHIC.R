#remotes::install_github("aidenlab/straw/R")
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


hi_c_plot=ggplot(all_hic_data, aes(x = x, y = y, fill = log_counts)) +
  geom_tile() +
  facet_grid(chromosome2 ~ chromosome1,space = "free", scales = "free")+#, scales = "free") +  # Facet by chromosome pairs
  scale_fill_gradient2(low = "#ebf5fb", mid = "#e74c3c",   high = "black", midpoint = 10) + # "#fef9e7" # choose colors
  theme_void() +
  labs(fill = "Interaction (log2)") + theme(panel.spacing = unit(0, "lines"))

ggsave("hi_c_plot2.pdf", hi_c_plot, width=8, height = 6)


#ggsave("hi_c_plot_small.pdf", hi_c_plot, width=6, height = 5)

