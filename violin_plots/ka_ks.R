# load library
library(dplyr)
library(readr)
library(ape)
library(tibble)
library(tidyr)
library(ggplot2)

# read in what should have been the input file from the KAKS python script

# grep '^#\|missense_variant\|synonymous_variant' annotated_vcf.vcf > mis_syn.txt

setwd("C:/Users/chfal/OneDrive/Desktop/Distichus_stuff/sliding_window")

# read data
miss_syn <- read.delim("miss_syn.txt", sep = c("\t"))

# extract variant type from the string that has all the info in it
miss_syn$variant_type <- sapply(strsplit(miss_syn$INFO, split="\\|"), '[', 2)

# extract the gene id from the string that has all the info in it
miss_syn$gene_id <- sapply(strsplit(miss_syn$INFO, split="\\|"), '[', 5)

# data cleaning and widening - group by gene ID, then count the number of each variant. pivot this table wider so that each column has its own name and type of variant.
miss_syn_ka_ks <- miss_syn %>%
  group_by(gene_id) %>%
  count(variant_type) %>%
  pivot_wider(names_from=variant_type, values_from=n)

# data calculating. sum up each synonymous or nonsynonymous variant.
miss_syn_ka_ks2 <- miss_syn_ka_ks %>%
  mutate(sum_synonymous = sum(synonymous_variant, `splice_region_variant&synonymous_variant`, na.rm=T)) %>%
  mutate( sum_nonsynonymous = sum(missense_variant, `missense_variant&splice_region_variant`, stop_gained, start_lost , `stop_lost&splice_region_variant`, na.rm=T))


# calculate KA_KS
miss_syn_ka_ks2$ka_ks <- miss_syn_ka_ks2$sum_nonsynonymous / miss_syn_ka_ks2$sum_synonymous


# get distinct gene IDs
miss_syn_gene_id <- miss_syn %>%
  select(gene_id,X.CHROM) %>%
  distinct()

# left join this with gene IDs
final_miss_syn <- left_join(miss_syn_ka_ks2,miss_syn_gene_id)

# set colors
colors <- c("#006838","#EE2A7B","#A4509F","#00AEEF")

# set types of each chromosome
final_miss_syn$chr_type <- "Autosome"
final_miss_syn[final_miss_syn$X.CHROM %in% c("scaffold_8"),]$chr_type <- "X1"
final_miss_syn[final_miss_syn$X.CHROM %in% c("scaffold_12"),]$chr_type <- "X2"
final_miss_syn[final_miss_syn$X.CHROM %in% c("scaffold_13","scaffold_16"),]$chr_type <- "Y"

# filter just the autosomes and sex chromosomes you want
final_miss_syn_filtered <-final_miss_syn[final_miss_syn$X.CHROM %in% c("scaffold_8","scaffold_9","scaffold_10","scaffold_11","scaffold_12","scaffold_13","scaffold_16"),]

# plot the data
ka_ks_plot <- ggplot(final_miss_syn_filtered,aes(y=chr_type,x=ka_ks, fill=chr_type)) +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
    geom_violin() +
  labs(x="Nonsynonymous to Synonymous Substitution Ratio Distribution", y=NULL) +
  scale_y_discrete(limits=c("Autosome","X1","X2","Y")) +
  theme(legend.position = "none") +
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

ka_ks_plot

## Statistical tests Ka/Ks ####

# get vectors of X1, X2, and Y
ancient_x <- final_miss_syn_filtered %>% filter(chr_type=="X1")
neo_x <- final_miss_syn_filtered %>% filter(chr_type=="X2")
y <- final_miss_syn_filtered %>% filter(chr_type=="Y")
autosomes <- final_miss_syn_filtered %>% filter(chr_type=="Autosome")

# compare X1 to autosomes, x1 is higher than autosomes
wilcox.test(ancient_x$ka_ks,autosomes$ka_ks,"greater")

# compare X2 to autosomes, x2 is higher than autosomes
wilcox.test(neo_x$ka_ks,autosomes$ka_ks,"greater")

# compare y to autosomes, y is higher than autosomes
wilcox.test(y$ka_ks,autosomes$ka_ks,"greater")

