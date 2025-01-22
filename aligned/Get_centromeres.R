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




library(dplyr)
library(zoo)
p=all_hic_data |> filter(log_counts>10) |> ## keeping only the highest values
  group_by(chromosome1) |>
  mutate(roll_c=rollmean(log_counts, k = 40, align="center", na.pad = TRUE))

p |> ggplot()+
  geom_point(aes(x = y, y = log_counts))+
  facet_wrap(~chromosome1)+
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5)
# calculate AT content per 100 bp windows
df_wide[, start := floor((POS - 1) / 100) * 100 + 1, by = CHROM]
df_wide[, end := start + 99, by = CHROM]

# per chromosome
p |> filter(chromosome1=="chr1") |>
  ggplot()+
  geom_point(aes(x = y, y = log_counts), size=2.5) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1100000)
df_wide[, A_T_count := (grepl("A|T", REF, ignore.case = TRUE) +
                        grepl("A|T", ALT, ignore.case = TRUE))]

p |> filter(chromosome1=="chr2") |>
  ggplot()+
  geom_point(aes(x = y, y = log_counts)) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 1800000)
window_size <- 100

p |> filter(chromosome1=="chr3") |>
  ggplot()+
  geom_point(aes(x = y, y = log_counts)) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 2750000)
result <- df_wide |>
  lazy_dt() %>%
  group_by(CHROM) %>%
  arrange(POS) %>%
  mutate(window_position = rollapply(A_T_count, width = window_size, FUN = mean, fill = NA)) -> result_windowed

p |> filter(chromosome1=="chr4") |>
  ggplot()+
  geom_point(aes(x = y, y = log_counts)) +
#  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 2200000)
# apply the sliding window to calculate AT content
result_windowed <- result %>%
  group_by(CHROM, window_position) %>%
  summarise(mean_at_count = mean(A_T_count, na.rm = TRUE))

p |> filter(chromosome1=="chr5") |>
  ggplot()+
  geom_point(aes(x = y, y = log_counts)) +
  geom_line(aes(x = y,y = roll_c), color ='red', linewidth=1.5) +
  geom_vline(xintercept = 2050000)
AT_summary = result_windowed |> as.data.table()
AT_summary$CHROM <- gsub("^C", "c", AT_summary$CHROM)

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
ql=unique(AT_summary$CHROM) |> as.data.frame()
ql$centromere=c(1100000,1800000,2750000,2200000,2050000,1150000,1350000,1400000,1250000,2400000,1600000)
colnames(ql) <- c("chromosome1", "Centromere")

p =  p   |> left_join(ql)
supp_centromeres=p |>
  ggplot(aes(x=y/1000000, y=log_counts ))+
  geom_point()+ xlab("Position (Mbp)")+ylab("Hi-C interaction >10")+
  geom_segment(aes(x=Centromere/1000000, y=-Inf, xend =Centromere/1000000, yend = Inf), color="red")+
AT_summary=paste0(AT_summary |> left_join(ql))
supp_centromeres=AT_summary |>
  ggplot(aes(x=window_position/1000000, y = mean_at_count ))+
  geom_point(alpha=0.5, size = 1.5, shape=21) +
  geom_segment(aes(x=Centromere/1000000, y=-Inf, xend =Centromere/1000000, yend = Inf), color="red", alpha=0.15)+
  facet_wrap(~chromosome1, scales = "fixed", ncol=5)

supp_centromeres

ggsave(supp_centromeres, file = "Supplementary_figure_centromeres.pdf", width = 12, height = 8)

## through AT content?
library(dtplyr)
library(zoo)
library(data.table)
library(vcfR)
library(tidyverse)
library(genomalicious)
gc()
#df=vcf2DT(vcfFile="Results/STRUCTURE/Whole_genomefiltered.vcf.gz", dropCols = NULL, keepComments = FALSE, keepInfo = FALSE)
vcf=read.vcfR("Results/STRUCTURE/Whole_genomefiltered.vcf.gz")
gc()
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##############

ref_seqs <- getREF(vcf)
alt_seqs <- getALT(vcf)
positions <- getPOS(vcf)
chromosomes <- getCHROM(vcf)

dvf=data.frame(ref_seqs, alt_seqs,positions,chromosomes)
dvf=as.data.table(dvf)
head(dvf)

calculate_sliding_1kb <- function(data, window_size = 1000, step_size = 100) {
  positions_diff <- diff(data$positions)
  window_rows <- sum(positions_diff <= window_size)

   data %>% lazy_dt() %>%
    # Filter and sort by chromosome and position
    group_by(chromosomes) %>%
    arrange(positions) %>%
    # Count A and T occurrences in ref_seqs and alt_seqs
    mutate(A_T_count = (ref_seqs %in% c("A", "T")) + (alt_seqs %in% c("A", "T"))) %>%
    # Apply sliding window
    mutate(sliding_avg_A_T = rollapply(
        A_T_count,
        width = window_rows,
        FUN = sum,
        align = "center",
        partial = TRUE
      )
    ) %>% as.data.frame()
  #  ungroup() %>% as.data.frame()
}

result <- calculate_sliding_1kb(dvf)
gc()
#saveRDS(result, file="AT_content.RDS")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(data.table)
library(dtplyr)
library(tidyverse)
library(GenomicRanges)

result=readRDS("AT_content.RDS")


# Example usage:
regions <- data.frame(cbind(
  chromosomes = paste0("Chr",seq(1:11)),
  startC = as.numeric(c(1099000, 1799000, 2749000, 2199000, 2049000, 1149000, 1349000, 1399000, 1249000, 2399000, 1599000)),
  endC = as.numeric(c(1101000, 1801000, 2751000, 2201000, 2051000, 1151000, 1351000, 1401000, 1251000, 2401000, 1601000)))
)
regions

result[,-6] |> head()


window_step <- 1000
window_size <- 10000



AT= result[,3:5] %>% #as.data.table() |> lazy_dt()  %>%
      {GRanges(seqnames=.$chromosomes,
               ranges=IRanges(start = .$positions, end=.$positions), AT=.$A_T_count)}  |>
      as.data.table() |> lazy_dt()  %>%  
  #as_tibble() |> 
      group_by(seqnames) |> 
  mutate(window=floor((start-1)/(window_size - window_step))) |> 
  mutate(end=start + window_size - 1) |> 
  group_by(seqnames, window) |> 
  summarize(start=min(start),
            end=max(end),
            ATcontent = sum(AT),
            .groups = 'drop') %>%
              mutate(midpoint = (start + end) / 2) %>%
              select(chromosomes = seqnames, start, end,  ATcontent, midpoint) |> as.data.table()


  
library(ggrastr)
library(ggpubr)
AT = AT |> left_join(regions, by="chromosomes") |> 
  mutate(startC=as.numeric(startC),
          endC=as.numeric(endC))

AT |> head()
gc()
AT$chromosomes=factor(AT$chromosomes, levels=
  c("Chr1",
"Chr2",
"Chr3",
    "Chr4","Chr5", "Chr6","Chr7","Chr8","Chr9", "Chr10","Chr11"
  ))
pppp <- ggplot(AT)+
  rasterize(
        geom_point(mapping=aes(x=midpoint, y=ATcontent), 
        size=2.5, alpha=0.5, color="steelblue"))+
  geom_rect(mapping=aes(xmin=startC-50000 ,xmax=endC+50000, 
    ymin=-Inf,ymax=Inf), 
      color="red", fill="transparent", linewidth=0.1) + 
facet_wrap(~chromosomes, nrow=3)+
  theme_bw()
  
pppp

ggsave("AT_per_Chr.pdf", pppp, width = 8, height = 5)

pppp






##############################################""


















vcf@fix
geno=vcfR::extract.gt(df)


## careful with sites with more than 2 alleles, need to load the genotype file drom vcftools too and compare
head(df, n=10)
gc()
df= df |> lazy_dt() |> filter(!GT=="./.") |> as.data.table()
df$GTs="Na"

# GTs in geno format
df[, GTs := fcase( GT %in% c("0|0", "0/0"), "0",
                   GT %in% c("1|0", "1/0", "0|1", "0/1"), "1",
                   GT %in% c("1|1", "1/1"), "2", default = "9")]

df_wide= df  |> lazy_dt() |>
          pivot_wider(
            names_from = SAMPLE,
            values_from = GTs,
            id_cols = c(LOCUS, REF, ALT, CHROM, POS)) |> as.data.table()

head(df_wide, n = 10)


## calculate AT content per 100 bp windows

# initiate windows taking chroms into account
df_wide[, start := floor((POS - 1) / 100) * 100 + 1, by = CHROM]
df_wide[, end := start + 99, by = CHROM]

# Calculate AT content
df_wide[, AT_count := (grepl("A|T", REF, ignore.case = TRUE) +
                        grepl("A|T", ALT, ignore.case = TRUE))]

df_wide %>% head()

AT_summary =df_wide |> lazy_dt() |> group_by(CHROM, start) |> summarise(AT=sum(AT_count)) |> as.data.table()

AT_summary$CHROM <- gsub("^C", "c", AT_summary$CHROM)

AT <- AT_summary |> lazy_dt() |> rename(chromosome1=CHROM) |> left_join(ql) |>
  mutate(
    centro = case_when(
      start > Centromere - 2000 & start < Centromere + 2000 ~ "centro",
      TRUE ~ "other")) |>
  as_tibble() |>
  #filter(AT>50) |>
ggplot(aes(x = start/1000000, y = AT, fill=centro, color=centro )) + geom_point(alpha=0.5, size = 1.5, shape=21) +
  geom_segment(aes(x=Centromere/1000000, y=-Inf, xend =Centromere/1000000, yend = Inf), color="red", alpha=0.15)+
  scale_fill_manual(values = c("#9999CC","#FF9933")) +
  scale_color_manual(values = c("#9999CC","#FF9933"))+
    facet_wrap(~chromosome1, ncol=5)

ggsave("AT_summary.pdf", plot = AT, width=12, height = 6)
ggsave("Supplementary_figure_centromeres.pdf", plot = supp_centromeres, width = 12, height = 8)
ql=unique(AT_summary$CHROM) |> as.data.frame()
df_wide[, start := floor((POS - 1) / 100) * 100 + 1, by = CHROM]
AT <- AT_summary |> lazy_dt() |> rename(chromosome1=CHROM) |> left_join(ql) |>
  mutate(
    centro = case_when(
      start > Centromere - 2000 & start < Centromere + 2000 ~ "centro",
      TRUE ~ "other")) |>
  as_tibble() |>
  #filter(AT>50) |>
ggplot(aes(x = start/1000000, y = AT, fill=centro, color=centro )) + geom_point(alpha=0.5, size = 1.5, shape=21) +
df_wide[, end := start + 99, by = CHROM]
ql$centromere=c(1100000,1800000,2750000,2200000,2050000,1150000,1350000,1400000,1250000,2400000,1600000)
colnames(ql) <- c("chromosome1", "Centromere")

AT_summary=paste0(AT_summary |> left_join(ql))
supp_centromeres=AT_summary |>
  ggplot(aes(x=window_position/1000000, y = mean_at_count ))+
  geom_point(alpha=0.5, size = 1.5, shape=21) +

  scale_fill_manual(values = c("#9999CC","#FF9933")) +
  scale_color_manual(values = c("#9999CC","#FF9933"))+
    facet_wrap(~chromosome1, ncol=5)

ggsave("AT_summary.pdf", plot = AT, width=12, height = 6)

  facet_wrap(~chromosome1, scales = "fixed", ncol=5)

ggsave("Supplementary_figure_centromeres.pdf", plot = supp_centromeres, width = 12, height = 8)
      TRUE ~ "other")) |> 
  as_tibble() |> 
  #filter(AT>50) |> 
ggplot(aes(x = start/1000000, y = AT, fill=centro, color=centro )) + geom_point(alpha=0.5, size = 1.5, shape=21) +
  geom_segment(aes(x=Centromere/1000000, y=-Inf, xend =Centromere/1000000, yend = Inf), color="red", alpha=0.15)+
  scale_fill_manual(values = c("#9999CC","#FF9933")) +
  scale_color_manual(values = c("#9999CC","#FF9933"))+
    facet_wrap(~chromosome1, ncol=5)

ggsave("AT_summary.pdf", plot = AT, width=12, height = 6)
