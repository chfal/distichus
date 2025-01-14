library(tidyverse)
library(cowplot)
# gotta set our colors EARLY!
colors <- c("#006838","#EE2A7B", "#A4509F","#00AEEF")

## Chromosome Coverage ####


# grep out scaffold with sed -i 
depth_by_chrm <- read.delim("sequencing_depth_by_chromosome.txt",sep=" ",header=F)

# set column names
colnames(depth_by_chrm)<-c("Scaffold", "scaffold_length", "bases", "coverage")                        
# set chromosome types
depth_by_chrm$chr_type <- "Autosome"
depth_by_chrm[depth_by_chrm$Scaffold %in% c(8),]$chr_type <- "X1"
depth_by_chrm[depth_by_chrm$Scaffold %in% c(12),]$chr_type <- "X2"
depth_by_chrm[depth_by_chrm$Scaffold %in% c(13,16),]$chr_type <- "Y"

# set the scaffold ID to be factor
depth_by_chrm$Scaffold <- as.factor(depth_by_chrm$Scaffold)

# plot coverage by each chromosome
coverage_by_chromosome <- ggplot(depth_by_chrm[depth_by_chrm$Scaffold %in% c(8,9,10,11,12,13,16),], aes(x=Scaffold, y=coverage, fill=as.factor(chr_type))) +
  geom_col() +
  scale_fill_manual(values=colors) +
  geom_text(aes(label = round(coverage, digits=2)), vjust = -0.2) +
  guides(fill=guide_legend(title="Chromosome Type")) +
  theme_bw() +
  theme(legend.text.align = 0) +
  ylim(0,60)

# display plot
coverage_by_chromosome

# save plot
ggsave("coverage_by_chromosome.pdf", coverage_by_chromosome, width=10, height=8)

## Window Coverage 10kb and 50kb ####

# read in all data
bedcov_10kb <- read.delim("bedcov_10kb.txt",sep="\t",header=F)

bedcov_50kb <- read.delim("bedcov_50kb.txt",sep="\t",header=F)

colnames(bedcov_10kb)<-c("Scaffold", "start", "end", "bedcov")                        

colnames(bedcov_50kb)<-c("Scaffold", "start", "end", "bedcov")                        


## Cumulative count, add to each

cumulative_count  <- bedcov_10kb %>% 
  group_by(Scaffold) %>% 
  summarise(maxbp= max(start)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(maxbp)), default = 0)) %>% 
  select(Scaffold, bp_add)

bedcov_10kb <- bedcov_10kb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add) %>%
  mutate(coverage = bedcov/10000)

bedcov_50kb <- bedcov_50kb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add) %>%
  mutate(coverage=bedcov/50000)


# set each chromosome type
bedcov_10kb$chr_type <- "Autosome"
bedcov_10kb[bedcov_10kb$Scaffold %in% c(8),]$chr_type <- "Ancient X"
bedcov_10kb[bedcov_10kb$Scaffold %in% c(12),]$chr_type <- "New X"
bedcov_10kb[bedcov_10kb$Scaffold %in% c(13,16),]$chr_type <- "Y"


bedcov_50kb$chr_type <- "Autosome"
bedcov_50kb[bedcov_50kb$Scaffold %in% c(8),]$chr_type <- "Ancient X"
bedcov_50kb[bedcov_50kb$Scaffold %in% c(12),]$chr_type <- "New X"
bedcov_50kb[bedcov_50kb$Scaffold %in% c(13,16),]$chr_type <- "Y"


# plot coverages in 10kb and 50kb windows

ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(7,8,9,10,11,12,13,14,15,16),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values=colors) +
  labs(title="10 kb window coverage") +
  guides(col= guide_legend(title= "Chromosome Type"))


ggplot(bedcov_50kb[bedcov_50kb$Scaffold %in% c(7,8,9,10,11,12,13,14,15,16),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values=colors) +
  labs(title="50 kb window coverage") +
  guides(col= guide_legend(title= "Chromosome Type"))


# plot individual chromosome coverage in 10kb windows

x1 <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(8),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#EE2A7B") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 8 - Ancient X", x="Cumulative Base Pairs") +
  theme_bw() +
  guides(colour="none")



y_13 <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(13),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#00AEEF") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 13 - Y", x="Cumulative Base Pairs") +
  theme_bw() +
  guides(colour="none")



y_16 <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(16),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#00AEEF") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 16 - Y", x="Cumulative Base Pairs") +
  theme_bw() +
  guides(colour="none")


x2 <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(12),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#A4509F") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 12 - Neo - X", x="Cumulative Base Pairs") +
  theme_bw() +
  guides(colour="none")

# plot together

coverage_windows <- plot_grid(x1, x2, y_13, y_16, ncol=2)

coverage_windows


## Statistical Tests ####

autosome_coverage <- filter(bedcov_50kb, chr_type=="Autosome") %>%
  select(coverage,Scaffold) %>%
  filter(Scaffold==9 | Scaffold==10 | Scaffold==11)

test_autosome_coverage <- filter(bedcov_50kb, chr_type=="Autosome") %>%
  select(coverage,Scaffold) %>%
  filter(Scaffold==14 | Scaffold==15 | Scaffold==7)

scaffold_7 <- filter(bedcov_50kb, Scaffold==7)

wilcox.test(scaffold_7$coverage,autosome_coverage$coverage)

y_coverage <- filter(bedcov_50kb, chr_type=="Y")

scaffold_13 <- filter(bedcov_50kb,Scaffold=="13")

scaffold_16 <- filter(bedcov_50kb,Scaffold=="16")

x2_coverage <- filter(bedcov_50kb, chr_type=="New X")

x1_coverage <- filter(bedcov_50kb, chr_type=="Ancient X")


wilcox.test(y_coverage$coverage, autosome_coverage$coverage,"less")

wilcox.test(x2_coverage$coverage, autosome_coverage$coverage,"less")

wilcox.test(x1_coverage$coverage,autosome_coverage$coverage,"less")


wilcox.test(scaffold_13$coverage,autosome_coverage$coverage,"less")

wilcox.test(scaffold_16$coverage,autosome_coverage$coverage,"less")

