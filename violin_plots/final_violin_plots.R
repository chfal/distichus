# load libraries
library(dplyr)
library(readr)
library(ape)
library(tibble)
library(tidyr)
library(ggplot2)
library(forcats)
library(cowplot)

## Reading and Managing Gene / TE density ####

# Read in data for gene density
AnoDis_gene_1mb <-read.delim("gene_density_1mb.txt", sep = "\t", header = FALSE) 

colnames(AnoDis_gene_1mb)<-c("Scaffold", "start", "end", "bp_count")                        


# Read in data for transposable element density
final_te_1mb <- read.delim("final_1mb_repeat_density.txt", sep = "\t", header = FALSE)

colnames(final_te_1mb)<-c("Scaffold", "start", "end", "bp_count")            

# Get cumulative count dataframe, which says how big each scaffold is
cumulative_count  <- AnoDis_gene_1mb %>% 
  group_by(Scaffold) %>% 
  summarise(maxbp= max(start)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(maxbp)), default = 0)) %>% 
  select(Scaffold, bp_add)

# add cumulative count dataframe to gene dataframe
AnoDis_gene_1mb <- AnoDis_gene_1mb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

# add cumulative count dataframe to te dataframe
final_te_1mb <- final_te_1mb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

# sets the axis for violin plots
axis_set_1mb <- AnoDis_gene_1mb %>% 
  group_by(Scaffold) %>% 
  summarize(center = mean(bp_cum))


## TE Density Plot ####

# add a new factor, which is the chromosome type

final_te_1mb$chr_type <- "Autosome"
final_te_1mb[final_te_1mb$Scaffold %in% c("8"),]$chr_type <- "X1"
final_te_1mb[final_te_1mb$Scaffold %in% c("12"),]$chr_type <- "X2"
final_te_1mb[final_te_1mb$Scaffold %in% c("13","16"),]$chr_type <- "Y"

# add custom vector of colors so it matches the rest of the plots
colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")

# select only the autosomes and sex chromosomes we want
final_te_1mb_filtered <- final_te_1mb %>%
  filter(Scaffold %in% c(8,9,10,11,12,13,16))


# plot a violin plot
final_te_violin <- ggplot(final_te_1mb_filtered,
                      aes(x= bp_count, y=as_factor(chr_type),
                          fill = as_factor(chr_type))) +  
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  scale_y_discrete(limits=c("Autosome","X1","X2","Y")) +
  theme_minimal() +
  labs(x="Transposable Element Density",y=NULL) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

# display plot in viewer
plot(final_te_violin)

ggsave("final_te_violin.pdf",final_te_violin, width=5, height=3)

## Gene Density Plot ####

# add factor, chromosome type, to each
AnoDis_gene_1mb$chr_type <- "Autosome"
AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("8"),]$chr_type <- "X1"
AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("12"),]$chr_type <- "X2"
AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("13","16"),]$chr_type <- "Y"

# select only the autosomes and sex chromosomes we want

AnoDis_gene_1mb_filtered <- AnoDis_gene_1mb %>%
  filter(Scaffold %in% c(8,9,10,11,12,13,16))

# plot gene density violin plot
gene_violin <- ggplot(AnoDis_gene_1mb_filtered,
                      aes(x= bp_count, y=as_factor(chr_type),
                          fill = as_factor(chr_type))) +  
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  scale_y_discrete(limits=c("Autosome","X1","X2","Y")) +
  theme_minimal() +
  labs(x="Gene Density",y=NULL) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

# display plot in viewer
plot(gene_violin)

# save plot to working directory
ggsave("gene_violin.pdf",gene_violin, width=5, height=3)

## SNPEFF Read in Data ####

# get high effects in 1mb windows
mb_high <-read.delim("1mb_high.txt", sep = "\t", header = FALSE) 

# get moderate effects in 1mb windows
mb_moderate <-read.delim("1mb_moderate.txt", sep = "\t", header = FALSE) 

# add column names to high and moderate effects
colnames(mb_high)<-c("Scaffold", "start", "end", "bp_count")            

colnames(mb_moderate)<-c("Scaffold", "start", "end", "bp_count")            

# create cumulative counts for this script
cumulative_count  <- mb_high %>% 
  group_by(Scaffold) %>% 
  summarise(maxbp= max(start)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(maxbp)), default = 0)) %>% 
  select(Scaffold, bp_add)

# get cumulative count column added with a join
mb_high <- mb_high %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

mb_moderate <- mb_moderate %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

# set axis for this graph
axis_set_1mb <- mb_high %>% 
  group_by(Scaffold) %>% 
  summarize(center = mean(bp_cum))

## Violin Plots SNPEFF ####

# add chromosome types
mb_high$chr_type <- "Autosome"
mb_high[mb_high$Scaffold %in% c("8"),]$chr_type <- "X1"
mb_high[mb_high$Scaffold %in% c("12"),]$chr_type <- "X2"
mb_high[mb_high$Scaffold %in% c("13","16"),]$chr_type <- "Y"

mb_moderate$chr_type <- "Autosome"
mb_moderate[mb_moderate$Scaffold %in% c("8"),]$chr_type <- "X1"
mb_moderate[mb_moderate$Scaffold %in% c("12"),]$chr_type <- "X2"
mb_moderate[mb_moderate$Scaffold %in% c("13","16"),]$chr_type <- "Y"

# set colors
colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")

# select high effect snps

high_snpeff_filtered <- mb_high %>%
  filter(Scaffold %in% c(8,9,10,11,12,13,16))


moderate_snpeff_filtered <- mb_moderate %>%
  filter(Scaffold %in% c(8,9,10,11,12,13,16))


# get plot of HIGH and MODERATE effects

high_for_join <- high_snpeff_filtered %>%
  select(start, end, bp_count, Scaffold) %>% 
  rename(bp_count_high = bp_count)

high_moderate_final <- left_join(high_for_join,moderate_snpeff_filtered, by=c("start","end","Scaffold"))

high_moderate_final$bp_count_final <- high_moderate_final$bp_count + high_moderate_final$bp_count_high


# filter ALL GENE WINDOWS WHERE IT'S 0
zeroes <-AnoDis_gene_1mb %>%
  filter(bp_count==0) %>%
  select(start, end, bp_count, Scaffold)

# then an anti join so we are only plotting the windows in which there ARE genes - because it will add unnecessary zeroes if we plot snpeff in ALL windows, including those with genes
high_moderate_no_zeroes <- anti_join(high_moderate_final, zeroes, by=c("start","end","Scaffold"))


final_snpeff_violin <- ggplot(high_moderate_no_zeroes,aes(x=bp_count,y=as_factor(chr_type),fill=as_factor(chr_type))) +
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  labs(x="High and Moderate Effect Variants",y=NULL) +
  scale_y_discrete(limits=c("Autosome","X1","X2","Y")) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

plot(final_snpeff_violin)


# save plot to working directory
ggsave("snpeff_violin.pdf",final_snpeff_violin, width=5, height=3)


## Statistical Tests ####

# GENES

gene_autosome <- filter(AnoDis_gene_1mb_filtered, chr_type=="Autosome") %>%
  select(bp_count)

gene_x1 <- filter(AnoDis_gene_1mb_filtered, chr_type=="X1") %>%
  select(bp_count)

gene_x2 <- filter(AnoDis_gene_1mb_filtered, chr_type=="X2") %>%
  select(bp_count)

gene_y <- filter(AnoDis_gene_1mb_filtered, chr_type=="Y") %>%
  select(bp_count)


# wilcox test:wilcox.test(the one you think is less, the one that you think is greater, "less")
# p-value will be in support of it being less

# wilcox.test(the one you think is greater, the one you think is less, 'greater')
# p-value will be in support of it being greater

# there are less genes on the Y chromosome than there are on autosomes
wilcox.test(gene_y$bp_count,gene_autosome$bp_count, "less")

# there are more genes on the ancient X chromosome than there are on autosomes but this is not significantly different
wilcox.test(gene_x1$bp_count,gene_autosome$bp_count, "greater")

# there are less genes on the neo-x chromosome than there are on autosomes
wilcox.test(gene_x2$bp_count,gene_autosome$bp_count, "less")


# TEs

te_autosome <- filter(final_te_1mb_filtered, chr_type=="Autosome") %>%
  select(bp_count)

te_x1 <- filter(final_te_1mb_filtered, chr_type=="X1") %>%
  select(bp_count)

te_x2 <- filter(final_te_1mb_filtered, chr_type=="X2") %>%
  select(bp_count)

te_y <- filter(final_te_1mb_filtered, chr_type=="Y") %>%
  select(bp_count)

# there are less TEs on the y than there are on Autosomes
wilcox.test(te_y$bp_count,te_autosome$bp_count,"less")


# there are less tes on the ancient x than there are on autosomes
wilcox.test(te_x1$bp_count, te_autosome$bp_count, "less")

# there are more TEs on the neo X than there are on autosomes
wilcox.test(te_x2$bp_count, te_autosome$bp_count, "greater")


# SNP EFF Effects 

snp_eff_autosome <- filter(high_moderate_no_zeroes, chr_type=="Autosome") %>%
  select(bp_count)

snp_eff_x1 <- filter(high_moderate_no_zeroes, chr_type=="X1") %>%
  select(bp_count)

snp_eff_x2 <- filter(high_moderate_no_zeroes, chr_type=="X2") %>%
  select(bp_count)

snp_eff_y <- filter(high_moderate_no_zeroes, chr_type=="Y") %>%
  select(bp_count)

# altogether there are more high and medium effects on sex chromosomes than autosomes

wilcox.test(snp_eff_x1$bp_count,snp_eff_autosome$bp_count, "less")

wilcox.test(snp_eff_x2$bp_count,snp_eff_autosome$bp_count, "less")

wilcox.test(snp_eff_y$bp_count,snp_eff_autosome$bp_count, "less")
