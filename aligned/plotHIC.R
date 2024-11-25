remotes::install_github("aidenlab/straw/R")
library(strawr)
library(tidyverse)


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
      "BP",1000)
    
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
  mutate(log_counts = log10(counts + 1))

hi_c_plot=ggplot(all_hic_data, aes(x = x, y = y, fill = log_counts)) +
  geom_tile() +
  facet_grid(chromosome2 ~ chromosome1, scales = "free") +  # Facet by chromosome pairs
  scale_fill_gradient(low = "white", high = "blue") +
  theme_void() +
  labs(fill = "Interaction (log10)")

ggsave("hi_c_plot.pdf", hi_c_plot, width=12, height = 12)