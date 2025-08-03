library(tidyverse)

# get high effects in 1mb windows
mb_high <-read.delim("1mb_high.txt", sep = "\t", header = FALSE) 

# get moderate effects in 1mb windows
mb_moderate <-read.delim("1mb_moderate.txt", sep = "\t", header = FALSE) 

# get low effects in 1mb windows
mb_low <-read.delim("1mb_low.txt", sep = "\t", header = FALSE) 


# add column names
colnames(mb_high)<-c("Scaffold", "start", "end", "count_high")            
colnames(mb_moderate)<-c("Scaffold", "start", "end", "count_moderate")

colnames(mb_low)<-c("Scaffold", "start", "end", "count_low")      

# join together data
join_1 <- left_join(mb_high,mb_moderate, by=c("start","end","Scaffold"))

join_2  <- left_join(join_1, mb_low, by=c("start","end","Scaffold"))

final_df <- join_2 %>%
  mutate(denominator = count_high + count_low + count_moderate) %>%
  mutate(numerator = count_high + count_moderate) %>%
  mutate(proportion = numerator / denominator)
  

## filter regions with no genes in them

AnoDis_gene_1mb <-read.delim("gene_density_1mb.txt", sep = "\t", header = FALSE) 

colnames(AnoDis_gene_1mb)<-c("Scaffold", "start", "end", "bp_count")                        

zeroes <-AnoDis_gene_1mb %>%
  filter(bp_count==0) %>%
  select(start, end, bp_count, Scaffold)

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


# then an anti join so we are only plotting the windows in which there ARE genes - because it will add unnecessary zeroes if we plot snpeff in ALL windows, including those with genes
no_zeroes_df <- anti_join(final_df, zeroes, by=c("start","end","Scaffold"))

## prepare to plot

# set axis for this graph
axis_set_1mb <- AnoDis_gene_1mb %>% 
  group_by(Scaffold) %>% 
  summarize(center = mean(bp_cum))

## Violin Plots SNPEFF ####

# add chromosome types
no_zeroes_df$chr_type <- "Autosome"
no_zeroes_df[no_zeroes_df$Scaffold %in% c("8"),]$chr_type <- "X1"
no_zeroes_df[no_zeroes_df$Scaffold %in% c("12"),]$chr_type <- "X2"
no_zeroes_df[no_zeroes_df$Scaffold %in% c("13","16"),]$chr_type <- "Y"


colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")

no_zeroes_df_filtered <- no_zeroes_df %>%
  filter(Scaffold %in% c(8,9,10,11,12,13,16))

proportion_violin <- ggplot(no_zeroes_df_filtered,aes(x=proportion,y=as_factor(chr_type),fill=as_factor(chr_type))) +
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  labs(x="Proportion of High and Moderate Effects to All Effects",y=NULL) +
  scale_y_discrete(limits=c("Autosome","X1","X2","Y")) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

proportion_violin


ggsave("proportion_violin.png",proportion_violin, width=5, height=3)


## Filtering for statistics

snp_eff_autosome_prop <- filter(no_zeroes_df_filtered, chr_type=="Autosome") %>%
  select(proportion)

snp_eff_x1_prop <- filter(no_zeroes_df_filtered, chr_type=="X1") %>%
  select(proportion)

snp_eff_x2_prop <- filter(no_zeroes_df_filtered, chr_type=="X2") %>%
  select(proportion)

snp_eff_y_prop <- filter(no_zeroes_df_filtered, chr_type=="Y") %>%
  select(proportion)

median(snp_eff_autosome_prop$proportion, na.rm=T)

median(snp_eff_x1_prop$proportion, na.rm=T)

median(snp_eff_x2_prop$proportion, na.rm=T)

median(snp_eff_y_prop$proportion, na.rm=T)
## Tests

# X1 has a higher proportion of effects than autosomes does
wilcox.test(snp_eff_x1_prop$proportion,snp_eff_autosome_prop$proportion, "greater")

#X2 has a higher proportion of effects than autosomes does
wilcox.test(snp_eff_x2_prop$proportion,snp_eff_autosome_prop$proportion, "greater")

# Y also has higher proportion of effects than autosomes does
wilcox.test(snp_eff_y_prop$proportion,snp_eff_autosome_prop$proportion, "greater")


