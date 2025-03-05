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

head(all_hic_data)
all_hic_data |> filter(log_counts>10) |>
  ggplot(aes(x=y, y=log_counts ))+
  geom_point()+
  facet_wrap(~chromosome1)


#~~~~~~~~~~~~~~~~~~ AT content ~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(data.table)
library(tidyverse)
library(Biostrings)
library(ggrastr)
calculate_sliding_at_content <- function(fasta_file, window_size = 1000, step_size = 100) {
  # Read the FASTA file - now keeping all sequences
  genome <- readDNAStringSet(fasta_file)
  
  # Initialize list to store results for each chromosome
  results_list <- list()
  
  # Process each chromosome
  for (i in seq_along(genome)) {
    sequence <- as.character(genome[[i]])
    chr_name <- names(genome)[i]
    # Calculate windows
    positions <- seq(1, nchar(sequence) - window_size + 1, by = step_size)
    at_content <- numeric(length(positions))
    # Calculate AT content for each window
    for (j in seq_along(positions)) {
      start <- positions[j]
      end <- start + window_size - 1
      window <- substr(sequence, start, end)
      # Count A and T occurrences
      a_count <- sum(strsplit(window, "")[[1]] == "A")
      t_count <- sum(strsplit(window, "")[[1]] == "T")
      # Calculate AT content percentage
      at_content[j] <- (a_count + t_count) * 100 / window_size
    }
    # Store results for this chromosome
    results_list[[chr_name]] <- data.frame(
      Position = positions,
      AT_Content = at_content,
      Chromosome = chr_name
    )
  }
  
  # Combine all results
  all_results <- do.call(rbind, results_list)
  return(all_results)
}

results <- calculate_sliding_at_content("Data/Fasta/Haplotype2_renamed_reordered_Chr_only.fasta")
results |> head()

gc()

results$Chromosome=factor(results$Chromosome, levels=c("Chr1",
"Chr2","Chr3","Chr4","Chr5", "Chr6","Chr7","Chr8","Chr9", "Chr10","Chr11"))


## Add in centromeres
regions <- data.frame(cbind(
  Chromosome = paste0("Chr",seq(1:11)),
  startC = as.numeric(c(1099000, 1799000, 2749000, 2199000, 2049000, 1149000, 1349000, 1399000, 1249000, 2399000, 1599000)),
  endC = as.numeric(c(1101000, 1801000, 2751000, 2201000, 2051000, 1151000, 1351000, 1401000, 1251000, 2401000, 1601000)))
)

results =results |> left_join(regions) |> 
  mutate(startC=as.numeric(startC),
          endC=as.numeric(endC),
        Chromosome=gsub("C", "c", Chromosome))

head(results)
p=ggplot(results)+ rasterize(
  geom_point(aes(x=Position, y=AT_Content), alpha=0.3, 
  color="steelblue"))+
  geom_rect(aes(xmin=startC, xmax=endC, ymin=-Inf, ymax=Inf), 
  color="black", linewidth=0.1, fill="transparent")+
  facet_wrap(~Chromosome, scales = "free_x")
p




######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~########
Centro=data.frame(results |> select(Chromosome, startC, endC)  |> filter(!duplicated(Chromosome)))

regions <- data.frame(
  Chromosome = factor(paste0("chr", 1:11), levels = paste0("chr", 1:11)),
  startC = as.numeric(c(1099000, 1799000, 2749000, 2199000, 2049000, 1149000, 1349000, 1399000, 1249000, 2399000, 1599000)),
  endC = as.numeric(c(1101000, 1801000, 2751000, 2201000, 2051000, 1151000, 1351000, 1401000, 1251000, 2401000, 1601000))
)

combined_data <- bind_rows(
  results %>% 
    mutate(
      type = "AT Content", 
      value = AT_Content, 
      Chromosome = factor(Chromosome, levels = paste0("chr", 1:11))
    ),
  all_hic_data %>% 
    filter(log_counts > 10) %>%
    mutate(
      type = "Hi-C Counts", 
      value = log_counts, 
      Position = x,
      Chromosome = factor(chromosome1, levels = paste0("chr", 1:11))
    ) |> left_join(Centro) |>   mutate(Chromosome = factor(Chromosome, levels = paste0("chr", c(1:5, 6:11))))
)

first_chroms <- paste0("chr", 1:5)
second_chroms <- factor(paste0("chr", 6:11), levels=paste0("chr", 6:11))

custom_labeller <- as_labeller(c(
  "AT Content" = "AT content (%)", 
  "Hi-C Counts" = "Hi-C interactions\n (>10 only)",
  "chr1" = "Chr1",
  "chr2" = "Chr2",
  "chr3" = "Chr3",
  "chr4" = "Chr4",
  "chr5" = "Chr5",
  "chr6" = "Chr6",
  "chr7" = "Chr7",
  "chr8" = "Chr8",
  "chr9" = "Chr9",
  "chr10" = "Chr10",
  "chr11" = "Chr11"))

p_first <- ggplot(combined_data %>% filter(Chromosome %in% first_chroms), 
                  aes(x = Position, y = value)) +
  rasterize(geom_point(aes(color = type), alpha = 0.6, size = 2)) +
  geom_rect(aes(xmin = startC, xmax = endC, ymin = -Inf, ymax = Inf), 
            color = "red", linewidth = 0.1, fill = "transparent") +
  facet_grid(rows = vars(type), cols = vars(Chromosome), scales = "free_y",
labeller = custom_labeller) +
  scale_color_manual(values = c("steelblue", "black")) +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  theme(strip.text.y = element_text(angle = 0, size=20),
strip.text.x.top = element_text(angle = 0, size=20),
        panel.spacing = unit(0, "lines"), 
        legend.position = "none")

p_second <- ggplot(combined_data %>% filter(Chromosome %in% second_chroms), 
                   aes(x = Position, y = value)) +
  rasterize(geom_point(aes(color = type), alpha = 0.6, size = 2)) +
  geom_rect(aes(xmin = startC, xmax = endC, ymin = -Inf, ymax = Inf), 
            color = "red", linewidth = 0.1, fill = "transparent") +
  facet_grid(rows = vars(type), cols = vars(Chromosome), scales = "free_y",labeller = custom_labeller) +
  scale_color_manual(values = c("steelblue", "black")) +
  theme_minimal() +
  labs(x = "Position", y = "", title = "") +
  theme(strip.text.y = element_text(angle = 0, size=20),
        strip.text.x.top = element_text(angle = 0, size=20),
        panel.spacing = unit(0, "lines"), axis.title.x.bottom = element_text(angle = 0, size=20),
        legend.position = "none")

p_full <- p_first / p_second



ggsave(p_full, file="Supplementary_panel_Feb10.pdf", width = 25, height = 15)
