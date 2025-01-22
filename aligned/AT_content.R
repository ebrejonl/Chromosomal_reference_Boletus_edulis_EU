#### get the AT content of centromeres vs rest of chromosome, potentially genome wide

library(data.table)
library(tidyverse)
library(Biostrings)

calculate_sliding_at_content <- function(fasta_file, window_size = 1000, step_size = 100) {
  # Read the FASTA file
  genome <- readDNAStringSet(fasta_file)
  sequence <- as.character(genome[[1]])
  
  # Initialize vectors to store results
  positions <- seq(1, nchar(sequence) - window_size + 1, by = step_size)
  at_content <- numeric(length(positions))
  
  # Calculate AT content for each window
  for (i in seq_along(positions)) {
    start <- positions[i]
    end <- start + window_size - 1
    window <- substr(sequence, start, end)
    
    # Count A and T occurrences
    a_count <- sum(strsplit(window, "")[[1]] == "A")
    t_count <- sum(strsplit(window, "")[[1]] == "T")
    
    # Calculate AT content percentage
    at_content[i] <- (a_count + t_count) * 100 / window_size
  }
  
  # Return results as a data frame
  return(data.frame(
    Position = positions,
    AT_Content = at_content
  ))
}

results <- calculate_sliding_at_content("Data/Fasta/Haplotype2_renamed_reordered_Chr_only.fasta", window_size=1000, step_size=100)

results |> head()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# per chromosomes:

library(Biostrings)

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