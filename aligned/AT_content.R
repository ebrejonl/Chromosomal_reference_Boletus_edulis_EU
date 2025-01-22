#### get the AT content of centromeres vs rest of chromosome, potentially genome wide
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
          endC=as.numeric(endC))


p=ggplot(results)+ rasterize(
  geom_point(aes(x=Position, y=AT_Content), alpha=0.3, 
  color="steelblue"))+
  geom_rect(aes(xmin=startC, xmax=endC, ymin=-Inf, ymax=Inf), 
  color="black", linewidth=0.1, fill="transparent")+
  facet_wrap(~Chromosome, scales = "free_x")
p

ggsave("Genome_wide_AT_per_chr_vs_Centromeres.pdf", p, width=8, height = 8)



resultC=results %>%
  group_by(Chromosome) %>%
  mutate(targetMidpoint = (startC + endC)/2,
closest1 = Position[order(abs( Position - ((startC + endC)/2)))[1]],
closest2 = Position[order(abs( Position - ((startC + endC)/2)))[2]],
closest3 = Position[order(abs(Position - ((startC + endC)/2)))[3]],
closest4 = Position[order(abs(Position - ((startC + endC)/2)))[4]],
closest5 = Position[order(abs(Position - ((startC + endC)/2)))[5]],
closest6 = Position[order(abs(Position - ((startC + endC)/2)))[6]],
closest7 = Position[order(abs(Position - ((startC + endC)/2)))[7]],
closest8 = Position[order(abs(Position - ((startC + endC)/2)))[8]],
closest9 = Position[order(abs(Position - ((startC + endC)/2)))[9]],
closest10 = Position[order(abs(Position - ((startC + endC)/2)))[10]]) %>%
  ungroup()

resultC |> head()

resultCwith= resultC|> group_by(Chromosome) |>  
  filter(Position==closest1 |
         Position==closest2 |
         Position==closest3 |
         Position==closest4 |
         Position==closest5 |
         Position==closest6 |
         Position==closest7 |
         Position==closest8 |
         Position==closest9 |
         Position==closest10)

resultCwithout= resultC|> group_by(Chromosome) |>  
  filter(!Position==closest1&
         !Position==closest2 &
         !Position==closest3 &
         !Position==closest4 &
         !Position==closest5 &
         !Position==closest6 &
         !Position==closest7 &
         !Position==closest8 &
         !Position==closest9 &
         !Position==closest10)

resultCwith$Chromosome=factor(resultCwith$Chromosome, levels=c("Chr1",
"Chr2","Chr3","Chr4","Chr5", "Chr6","Chr7","Chr8","Chr9", "Chr10","Chr11"))
        
resultCwithout$Chromosome=factor(resultCwithout$Chromosome, levels=c("Chr1",
"Chr2","Chr3","Chr4","Chr5", "Chr6","Chr7","Chr8","Chr9", "Chr10","Chr11"))
      


p2=ggplot()+
  geom_boxplot(resultCwith,mapping=aes(x="Centromere", y=AT_Content), fill="purple", alpha=0.3,outliers=F)+
  geom_boxplot(resultCwithout,mapping=aes(x="WGS", y=AT_Content), fill="steelblue", alpha=0.3,outliers=F)+
  facet_wrap(~Chromosome)+
  theme_bw()+ xlab(NULL)+ylab("AT %")

p2
ggsave("AT_per_chr_vs_Centromeres.pdf", p2, width=8, height = 8)


p3=ggplot()+
  geom_boxplot(resultCwith,mapping=aes(x="Centromere", y=AT_Content), fill="purple", alpha=0.3,outliers=F)+
  #geom_jitter(resultCwith,mapping=aes(x="Centromere", y=AT_Content), size=3, fill="purple", width=0.3)+
  #geom_jitter(resultCwith,mapping=aes(x="WGS", y=AT_Content), size=3,shape=21,fill="steelblue", width=0.3)+
  geom_boxplot(resultCwithout,mapping=aes(x="WGS", y=AT_Content),shape=21, fill="steelblue", alpha=0.3,outliers=F)+
  theme_bw()+ xlab(NULL)+ylab("AT %")

p3
ggsave("WGS_vs_Centromeres_boxplot.pdf", p3, width=6, height = 6)
