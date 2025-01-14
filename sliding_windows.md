

```
# Get lengths of fasta entries
samtools faidx ${SAMPLE}.fa
cut -f1-2 ${SAMPLE}.fa.fai > ${SAMPLE}.chrom.sizes

# make bed file of windows
bedtools makewindows -g ${SAMPLE}.chrom.sizes -w 10000 > 10kb_windows.bed
bedtools makewindows -g ${SAMPLE}.chrom.sizes -w 1000000 > 1Mb_windows.bed

# get counts for each window of genes, repeats, and snp-eff outputs
bedtools coverage -a 1Mb_windows.bed -b genes.gff -counts > gene_density.txt

bedtools coverage -a 10kb_windows.bed -b repeats.gff  -counts > repeat_density.txt
# for repeat density, did 2: one of jody's Repeatmasker and one of Dovetail's repeatmasker

# did this also for snp-eff outputs

# had to make a file with and without the header and join them together bc bedtools doesn't like vcf headers - then had to add back in
```


<details><summary>density_lanscape.R</summary>

```

##########################################################
# Density landscape of genome features
# Author: Jody M. Taft
##########################################################
# Description: 
# Here I present the landscape of genome features for the worst lizard in the world, Anolis Distichus
##########################################################

# Load packages 

library(tidyverse)

# Import dataset

setwd("C:/Users/chfal/OneDrive/Desktop/sliding_window")

# Remember to remove "scaffold_" from file before importing I used sed.

# read in data
# AnoDis_gene_10kb <-read.delim("gene_density_10kb.txt", sep = "\t", header = FALSE) 
AnoDis_gene_1mb <-read.delim("gene_density_1mb.txt", sep = "\t", header = FALSE) 

# Dovetail_te_10kb <- read.delim("dovetail_10kb_repeat_density.txt", sep = "\t", header = FALSE)

Dovetail_te_1mb <- read.delim("dovetail_1mb_repeat_density.txt", sep = "\t", header = FALSE)

# jody_te_10kb <- read.delim("jody_10kb_repeat_density.txt", sep = "\t", header = FALSE)

jody_te_1mb <- read.delim("jody_1mb_repeat_density.txt", sep = "\t", header = FALSE)


# add heading names to dataframe    
# colnames(AnoDis_gene_10kb)<-c("Scaffold", "start", "end", "bp_count")                        
colnames(AnoDis_gene_1mb)<-c("Scaffold", "start", "end", "bp_count")                        
# colnames(Dovetail_te_10kb)<-c("Scaffold", "start", "end", "bp_count")           

colnames(Dovetail_te_1mb)<-c("Scaffold", "start", "end", "bp_count")            

# colnames(jody_te_10kb)<-c("Scaffold", "start", "end", "bp_count")           

colnames(jody_te_1mb)<-c("Scaffold", "start", "end", "bp_count")            


# Add a cumulative column for the x-axis
cumulative_count  <- AnoDis_gene_1mb %>% 
  group_by(Scaffold) %>% 
  summarise(maxbp= max(start)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(maxbp)), default = 0)) %>% 
  select(Scaffold, bp_add)


# AnoDis_gene_10kb <- AnoDis_gene_10kb %>% 
#   inner_join(cumulative_count, by = "Scaffold") %>% 
#   mutate(bp_cum = end + bp_add)

AnoDis_gene_1mb <- AnoDis_gene_1mb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)


# Dovetail_te_10kb <- Dovetail_te_10kb %>% 
#   inner_join(cumulative_count, by = "Scaffold") %>% 
#   mutate(bp_cum = end + bp_add)

Dovetail_te_1mb <- Dovetail_te_1mb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

# jody_te_10kb <- jody_te_10kb %>% 
#   inner_join(cumulative_count, by = "Scaffold") %>% 
#   mutate(bp_cum = end + bp_add)

jody_te_1mb <- jody_te_1mb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)


# axis_set_10kb <- AnoDis_gene_10kb %>% 
#   group_by(Scaffold) %>% 
#   summarize(center = mean(bp_cum))

axis_set_1mb <- AnoDis_gene_1mb %>% 
  group_by(Scaffold) %>% 
  summarize(center = mean(bp_cum))


# genes_10kb <- ggplot(AnoDis_gene_10kb[AnoDis_gene_10kb$Scaffold %in% c("8","12","13","16"),], 
#                 aes(x = bp_cum, y = bp_count, 
#                 color = as_factor(Scaffold))) +  
#   geom_point(alpha = 0.75) +
#   scale_x_continuous(label = axis_set_10kb$Scaffold, breaks = axis_set_10kb$center) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,50)) + 
#   # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = NULL, 
#        y = "Gene Density") + 
#   theme_minimal() +
#   theme( 
#     legend.position = "none",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     ) +
#   labs(title="10kb gene density")
#   # geom_smooth(method = "gam", formula = y~s(x))
#   
# plot(genes_10kb) 


genes_1mb <- ggplot(AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("8","12","13","16"),], 
                aes(x = bp_cum, y = bp_count, 
                    color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,50)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb window gene density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(genes_1mb)


# dovetail_10kb <- ggplot(Dovetail_te_10kb[Dovetail_te_10kb$Scaffold %in% c("8","12","13","16"),], 
#                      aes(x = bp_cum, y = bp_count, 
#                          color = as_factor(Scaffold))) +  
#   geom_point(alpha = 0.75) +
#   scale_x_continuous(label = axis_set_10kb$Scaffold, breaks = axis_set_10kb$center) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,50)) + 
#   # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = NULL, 
#        y = "Repeat Density") + 
#   theme_minimal() +
#   theme( 
#     legend.position = "none",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#   ) +
#   labs(title="10kb dovetail TE density")
# # geom_smooth(method = "gam", formula = y~s(x))
# 
# plot(dovetail_10kb) 


dovetail_1mb <- ggplot(Dovetail_te_1mb[Dovetail_te_1mb$Scaffold %in% c("8","12","13","16"),], 
                        aes(x = bp_cum, y = bp_count, 
                            color = as_factor(Scaffold))) +  
  geom_point() +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,15000)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb dovetail TE density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(dovetail_1mb)


# jody_10kb <- ggplot(jody_te_10kb[jody_te_10kb$Scaffold %in% c("8","12","13","16"),], 
#                         aes(x = bp_cum, y = bp_count, 
#                             color = as_factor(Scaffold))) +  
#   geom_point(alpha = 0.75) +
#   scale_x_continuous(label = axis_set_10kb$Scaffold, breaks = axis_set_10kb$center) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,50)) + 
#   # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = NULL, 
#        y = "Repeat Density") + 
#   theme_minimal() +
#   theme( 
#     legend.position = "none",
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#   ) +
#   labs(title="10kb repeatmasker jody TE density")
# # geom_smooth(method = "gam", formula = y~s(x))
# 
# plot(jody_10kb) 


jody_1mb <- ggplot(jody_te_1mb[jody_te_1mb$Scaffold %in% c("8","12","13","16"),],
                       aes(x = bp_cum, y = bp_count, 
                           color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,20000)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb repeatmasker jody TE density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(jody_1mb)



### VIOLIN PLOTS


jody_violin <- ggplot(jody_te_1mb[jody_te_1mb$Scaffold %in% c("8","12","13","16"),],
                   aes(x= bp_count, y=as_factor(Scaffold),
                       fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb repeatmasker jody TE violin plot")

dovetail_violin <- ggplot(Dovetail_te_1mb[Dovetail_te_1mb$Scaffold %in% c("8","12","13","16"),],
       aes(x= bp_count, y=as_factor(Scaffold),
           fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb repeatmasker dovetail TE violin plot")

plot(dovetail_violin)

gene_violin <- ggplot(AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("8","12","13","16"),],
       aes(x= bp_count, y=as_factor(Scaffold),
           fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb gene density violin plot")

plot(gene_violin)
```
</details>


<details><summary>Snp_Eff_Landscape.R</summary>

```

##########################################################
# Density landscape of genome features
# Author: Jody M. Taft
##########################################################
# Description: 
# Here I present the landscape of genome features for the worst lizard in the world, Anolis Distichus
##########################################################

# Load packages 

library(tidyverse)

# Import dataset

setwd("C:/Users/chfal/OneDrive/Desktop/sliding_window")

# Remember to remove "scaffold_" from file before importing I used sed.

# read in data
kb_high <- read.delim("10kb_high.txt", sep = "\t", header = FALSE)

mb_high <-read.delim("1mb_high.txt", sep = "\t", header = FALSE) 

kb_moderate <- read.delim("10kb_moderate.txt", sep = "\t", header = FALSE)

mb_moderate <-read.delim("1mb_moderate.txt", sep = "\t", header = FALSE) 

kb_low <- read.delim("10kb_low.txt", sep = "\t", header = FALSE)

mb_low <-read.delim("1mb_low.txt", sep = "\t", header = FALSE) 

             
colnames(kb_high)<-c("Scaffold", "start", "end", "bp_count")                        
colnames(mb_high)<-c("Scaffold", "start", "end", "bp_count")            

colnames(kb_moderate)<-c("Scaffold", "start", "end", "bp_count")            

colnames(mb_moderate)<-c("Scaffold", "start", "end", "bp_count")            

colnames(kb_low)<-c("Scaffold", "start", "end", "bp_count")            

colnames(mb_low)<-c("Scaffold", "start", "end", "bp_count")            


cumulative_count  <- kb_high %>% 
  group_by(Scaffold) %>% 
  summarise(maxbp= max(start)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(maxbp)), default = 0)) %>% 
  select(Scaffold, bp_add)


kb_high <- kb_high %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

mb_high <- mb_high %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

kb_moderate <- kb_moderate %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

mb_moderate <- mb_moderate %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

kb_low <- kb_low %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

mb_low <- mb_low %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

axis_set_10kb <- kb_high %>% 
   group_by(Scaffold) %>% 
   summarize(center = mean(bp_cum))

axis_set_1mb <- mb_high %>% 
  group_by(Scaffold) %>% 
  summarize(center = mean(bp_cum))



plot_kb_high <- ggplot(kb_high[kb_high$Scaffold %in% c("8","12","13","16"),], 
                    aes(x = bp_cum, y = bp_count, 
                        color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,10)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="10kb window high snp-eff density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(plot_kb_high)

plot_mb_high <- ggplot(mb_high[mb_high$Scaffold %in% c("8","12","13","16"),], 
                       aes(x = bp_cum, y = bp_count, 
                           color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,7)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb window high snp-eff density")

plot(plot_mb_high)

plot_kb_moderate <- ggplot(kb_moderate[kb_moderate$Scaffold %in% c("8","12","13","16"),], 
                       aes(x = bp_cum, y = bp_count, 
                           color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,20)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="10 kb window moderate snp-eff density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(plot_kb_moderate)


plot_mb_moderate <- ggplot(mb_moderate[mb_moderate$Scaffold %in% c("8","12","13","16"),], 
                           aes(x = bp_cum, y = bp_count, 
                               color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1 mb window moderate snp-eff density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(plot_mb_moderate)


plot_kb_low <- ggplot(kb_low[kb_low$Scaffold %in% c("8","12","13","16"),], 
                           aes(x = bp_cum, y = bp_count, 
                               color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,30)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="10 kb window low snp-eff density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(plot_kb_low)


plot_mb_low <- ggplot(mb_low[mb_low$Scaffold %in% c("8","12","13","16"),], 
                           aes(x = bp_cum, y = bp_count, 
                               color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,120)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1 mb window low snp-eff density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(plot_mb_low)


### VIOLIN PLOTS


violin_mb_high <- ggplot(mb_high[mb_high$Scaffold %in% c("8","12","13","16"),],
                      aes(x= bp_count, y=as_factor(Scaffold),
                          fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb high snpeff violin plot")

plot(violin_mb_high)


violin_mb_moderate <- ggplot(mb_moderate[mb_moderate$Scaffold %in% c("8","12","13","16"),],
                         aes(x= bp_count, y=as_factor(Scaffold),
                             fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb moderate snpeff violin plot")

plot(violin_mb_moderate)


violin_mb_low <- ggplot(mb_low[mb_low$Scaffold %in% c("8","12","13","16"),],
                             aes(x= bp_count, y=as_factor(Scaffold),
                                 fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb low snpeff violin plot")

plot(violin_mb_low)

```
</details>

  


<details><summary>Final_R_Script.R</summary>

  ```
library(tidyverse)

setwd("C:/Users/chfal/OneDrive/Desktop/sliding_window")

# Remember to remove "scaffold_" from file before importing I used sed.

## Gene & TE Density ####

AnoDis_gene_1mb <-read.delim("gene_density_1mb.txt", sep = "\t", header = FALSE) 

colnames(AnoDis_gene_1mb)<-c("Scaffold", "start", "end", "bp_count")                        


# te density
jody_te_1mb <- read.delim("jody_1mb_repeat_density.txt", sep = "\t", header = FALSE)

colnames(jody_te_1mb)<-c("Scaffold", "start", "end", "bp_count")            

cumulative_count  <- AnoDis_gene_1mb %>% 
  group_by(Scaffold) %>% 
  summarise(maxbp= max(start)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(maxbp)), default = 0)) %>% 
  select(Scaffold, bp_add)

AnoDis_gene_1mb <- AnoDis_gene_1mb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

jody_te_1mb <- jody_te_1mb %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

axis_set_1mb <- AnoDis_gene_1mb %>% 
  group_by(Scaffold) %>% 
  summarize(center = mean(bp_cum))

genes_1mb <- ggplot(AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("8","12","13","16"),], 
                    aes(x = bp_cum, y = bp_count, 
                        color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,50)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb window gene density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(genes_1mb)

jody_1mb <- ggplot(jody_te_1mb[jody_te_1mb$Scaffold %in% c("8","12","13","16"),],
                   aes(x = bp_cum, y = bp_count, 
                       color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,20000)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb repeatmasker jody TE density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(jody_1mb)


## Violin Plots ####

# add a new factor
jody_te_1mb$chr_type <- "Autosome"
jody_te_1mb[jody_te_1mb$Scaffold %in% c("8"),]$chr_type <- "Ancient X"
jody_te_1mb[jody_te_1mb$Scaffold %in% c("12"),]$chr_type <- "New X"
jody_te_1mb[jody_te_1mb$Scaffold %in% c("13","16"),]$chr_type <- "Y"

colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")

jody_te_1mb_selected <- jody_te_1mb[jody_te_1mb$Scaffold %in% c("8","9","12","13","16"),]

jody_violin <- ggplot(jody_te_1mb_selected,
                      aes(x= bp_count, y=as_factor(chr_type),
                          fill = as_factor(chr_type))) +  
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  scale_y_discrete(limits=c("Autosome","Ancient X","New X","Y"),labels=c("Autosome","Ancient X","New X","Y")) +
  theme_minimal() +
  labs(x=NULL,y=NULL) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

plot(jody_violin)

ggsave("C:/Users/chfal/OneDrive/Desktop/sliding_window/jody_te.pdf", jody_violin)


AnoDis_gene_1mb$chr_type <- "Autosome"
AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("8"),]$chr_type <- "Ancient X"
AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("12"),]$chr_type <- "New X"
AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("13","16"),]$chr_type <- "Y"

AnoDis_gene_1mb_selected <- AnoDis_gene_1mb[AnoDis_gene_1mb$Scaffold %in% c("8","9","10","11","12","13","16"),]


gene_violin <- ggplot(AnoDis_gene_1mb_selected,
                      aes(x= bp_count, y=as_factor(chr_type),
                          fill = as_factor(chr_type))) +  
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  scale_y_discrete(limits=c("Autosome","Ancient X","New X","Y"),labels=c("Autosome","Ancient X","New X","Y")) +
  theme_minimal() +
  labs(x=NULL,y=NULL) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

plot(gene_violin)

ggsave("C:/Users/chfal/OneDrive/Desktop/sliding_window/gene_violin.pdf", gene_violin)


## SNPEFF ####

# read in data

mb_high <-read.delim("1mb_high.txt", sep = "\t", header = FALSE) 

mb_moderate <-read.delim("1mb_moderate.txt", sep = "\t", header = FALSE) 

mb_low <-read.delim("1mb_low.txt", sep = "\t", header = FALSE) 


colnames(mb_high)<-c("Scaffold", "start", "end", "bp_count")            


colnames(mb_moderate)<-c("Scaffold", "start", "end", "bp_count")            


colnames(mb_low)<-c("Scaffold", "start", "end", "bp_count")            


cumulative_count  <- mb_high %>% 
  group_by(Scaffold) %>% 
  summarise(maxbp= max(start)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(maxbp)), default = 0)) %>% 
  select(Scaffold, bp_add)


mb_high <- mb_high %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)


mb_moderate <- mb_moderate %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

mb_low <- mb_low %>% 
  inner_join(cumulative_count, by = "Scaffold") %>% 
  mutate(bp_cum = end + bp_add)

axis_set_1mb <- mb_high %>% 
  group_by(Scaffold) %>% 
  summarize(center = mean(bp_cum))

plot_mb_high <- ggplot(mb_high[mb_high$Scaffold %in% c("8","12","13","16"),], 
                       aes(x = bp_cum, y = bp_count, 
                           color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,7)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb window high snp-eff density")

plot(plot_mb_high)

plot_mb_moderate <- ggplot(mb_moderate[mb_moderate$Scaffold %in% c("8","12","13","16"),], 
                           aes(x = bp_cum, y = bp_count, 
                               color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1 mb window moderate snp-eff density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(plot_mb_moderate)

plot_mb_low <- ggplot(mb_low[mb_low$Scaffold %in% c("8","12","13","16"),], 
                      aes(x = bp_cum, y = bp_count, 
                          color = as_factor(Scaffold))) +  
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set_1mb$Scaffold, breaks = axis_set_1mb$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,120)) + 
  # scale_color_manual(values = rep(c("#00d6bbff", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Gene Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1 mb window low snp-eff density")
# geom_smooth(method = "gam", formula = y~s(x))

plot(plot_mb_low)


## Violin Plots SNPEFF ####

mb_high$chr_type <- "Autosome"
mb_high[mb_high$Scaffold %in% c("8"),]$chr_type <- "Ancient X"
mb_high[mb_high$Scaffold %in% c("12"),]$chr_type <- "New X"
mb_high[mb_high$Scaffold %in% c("13","16"),]$chr_type <- "Y"

mb_moderate$chr_type <- "Autosome"
mb_moderate[mb_moderate$Scaffold %in% c("8"),]$chr_type <- "Ancient X"
mb_moderate[mb_moderate$Scaffold %in% c("12"),]$chr_type <- "New X"
mb_moderate[mb_moderate$Scaffold %in% c("13","16"),]$chr_type <- "Y"


mb_low$chr_type <- "Autosome"
mb_low[mb_low$Scaffold %in% c("8"),]$chr_type <- "Ancient X"
mb_low[mb_low$Scaffold %in% c("12"),]$chr_type <- "New X"
mb_low[mb_low$Scaffold %in% c("13","16"),]$chr_type <- "Y"


violin_mb_moderate <- ggplot(mb_moderate[mb_moderate$Scaffold %in% c("8","12","13","16"),],
                             aes(x= bp_count, y=as_factor(Scaffold),
                                 fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Mutations with High Effects") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb moderate snpeff violin plot") 

plot(violin_mb_moderate)


violin_mb_low <- ggplot(mb_low[mb_low$Scaffold %in% c("8","12","13","16"),],
                        aes(x= bp_count, y=as_factor(Scaffold),
                            fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb low snpeff violin plot")

plot(violin_mb_low)


colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")


high_snpeff <- mb_high[mb_high$Scaffold %in% c("8", "9","10","11,","12","13","16"),]

high_snpeff  <- high_snpeff %>%
  filter(bp_count>0)


violin_mb_high <- ggplot(high_snpeff,aes(x=bp_count,y=as_factor(chr_type),fill=as_factor(chr_type))) +
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  labs(x=NULL,y=NULL) +
  scale_y_discrete(limits=c("Autosome","Ancient X","New X","Y"),labels=c("Autosome","Ancient X","New X","Y")) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

plot(violin_mb_high)
ggsave("C:/Users/chfal/OneDrive/Desktop/sliding_window/violin_mb_high.pdf", violin_mb_high)


violin_mb_moderate <- ggplot(mb_moderate[mb_moderate$Scaffold %in% c("8","12","13","16"),],
                             aes(x= bp_count, y=as_factor(Scaffold),
                                 fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Mutations with High Effects") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb moderate snpeff violin plot") 

plot(violin_mb_moderate)


violin_mb_low <- ggplot(mb_low[mb_low$Scaffold %in% c("8","12","13","16"),],
                        aes(x= bp_count, y=as_factor(Scaffold),
                            fill = as_factor(Scaffold))) +  
  geom_violin(alpha = 0.75) +
  labs(x = NULL, 
       y = "Repeat Density") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(title="1mb low snpeff violin plot")

plot(violin_mb_low)



## MAKING MORE PLOTS ###


# get plot of HIGH and MODERATE effects

high_for_join <- mb_high %>%
  select(start, end, bp_count, Scaffold) %>% 
  rename(bp_count_high = bp_count)

high_moderate <- left_join(high_for_join,mb_moderate, by=c("start","end","Scaffold"))

high_moderate$bp_count_final <- high_moderate$bp_count + high_moderate$bp_count_high


# filter ALL GENE WINDOWS WHERE IT'S 0
zeroes <-AnoDis_gene_1mb %>%
  filter(bp_count==0) %>%
  select(start, end, bp_count, Scaffold)

# then an anti join
high_moderate_no_zeroes <- anti_join(high_moderate, zeroes, by=c("start","end","Scaffold"))


high_moderate_no_zeroes <- high_moderate_no_zeroes[new_test$Scaffold %in% c("8", "9","10","11,","12","13","16"),]


final_violin <- ggplot(high_moderate_no_zeroes,aes(x=bp_count,y=as_factor(chr_type),fill=as_factor(chr_type))) +
  geom_violin() +
  stat_summary(fun="mean",geom="crossbar",color="black") +
  scale_fill_manual(values=colors) +
  theme_minimal() +
  labs(x=NULL,y=NULL) +
  scale_y_discrete(limits=c("Autosome","Ancient X","New X","Y"),labels=c("Autosome","Ancient X","New X","Y")) +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

plot(final_violin)
ggsave("C:/Users/chfal/OneDrive/Desktop/sliding_window/no_zeroes_violin_moderate_high.pdf", final_violin)

## Statistical Tests ####

# GENES

gene_autosome <- filter(AnoDis_gene_1mb, chr_type=="Autosome") %>%
  select(bp_count)

gene_ancientx <- filter(AnoDis_gene_1mb, chr_type=="Ancient X") %>%
  select(bp_count)

gene_neox <- filter(AnoDis_gene_1mb, chr_type=="New X") %>%
  select(bp_count)

gene_y <- filter(AnoDis_gene_1mb, chr_type=="Y") %>%
  select(bp_count)


# wilcox test:wilcox.test(the one you think is less, the one that you think is greater, "less")
# p-value will be in support of it being less

# wilcox.test(the one you think is greater, the one you think is less, 'greater')
# p-value will be in support of it being greater

# there are less genes on the Y chromosome than there are on autosomes
wilcox.test(gene_y$bp_count,gene_autosome$bp_count, "less")

# there are more genes on the ancient X chromosome than there are on autosomes
wilcox.test(gene_ancientx$bp_count,gene_autosome$bp_count, "greater")

# there are less genes on the neo-x chromosome than there are on autosomes, but this is not significant
wilcox.test(gene_neox$bp_count,gene_autosome$bp_count, "less")


# TEs

te_autosome <- filter(jody_te_1mb, chr_type=="Autosome") %>%
  select(bp_count)

te_ancientx <- filter(jody_te_1mb, chr_type=="Ancient X") %>%
  select(bp_count)

te_neox <- filter(jody_te_1mb, chr_type=="New X") %>%
  select(bp_count)

te_y <- filter(jody_te_1mb, chr_type=="Y") %>%
  select(bp_count)


# there are less TEs on the y than there are on Autosomes
wilcox.test(te_y$bp_count,te_autosome$bp_count,"less")


# there are more TEs on the neo X than there are on autosomes
wilcox.test(te_neox$bp_count,te_autosome$bp_count, "greater")

# there are more tes on the ancient X than there are on autosomes
wilcox.test(te_ancientx$bp_count, te_autosome$bp_count, "greater")


# SNP EFF Effects 

snp_eff_autosome <- filter(high_moderate_no_zeroes, chr_type=="Autosome") %>%
  select(bp_count)

snp_eff_ancientx <- filter(high_moderate_no_zeroes, chr_type=="Ancient X") %>%
  select(bp_count)

snp_eff_newx <- filter(high_moderate_no_zeroes, chr_type=="New X") %>%
  select(bp_count)

snp_eff_y <- filter(high_moderate_no_zeroes, chr_type=="Y") %>%
  select(bp_count)



wilcox.test(snp_eff_ancientx$bp_count,snp_eff_autosome$bp_count, "less")

wilcox.test(snp_eff_newx$bp_count,snp_eff_autosome$bp_count, "less")

wilcox.test(snp_eff_y$bp_count,snp_eff_autosome$bp_count, "less")

```
</details>
