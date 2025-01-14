Attempting to understand coverage across both sex chromosomes and autosomes


This link links to the samtools bedcov function which adds up the read depth at each base, then we divide by window size to get per-window coverage.

https://www.htslib.org/doc/samtools-bedcov.html

<details><summary>coverage.sh</summary>

  ```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=coverage                      # job name for listing in queue
#SBATCH --mem=50G                              # memory to allocate in Mb
#SBATCH -n 1                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs

echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module load java
module load bedtools2

echo "Bash commands for the analysis you are going to run"

echo "#...SAMTOOLS Summary stats"

# from manual, does not work - produced depth but at every position
# samtools depth AnolisDistichus.addGP.bam -b 10kb_windows.bed -o GATK_bam_read_depth_10kb.txt
# samtools depth AnolisDistichus.addGP.bam -b 50kb_windows.bed -o GATK_bam_read_depth_50kb.txt

# works by chromosome, but does not actually have fine-scale windows
# samtools depth AnolisDistichus.addGP.bam | awk '{c[$1]++;s[$1]+=$3} END {for (i in c) print i, c[i], s[i], s[i]/c[i]}' > sequencing_depth_by_chromosome.txt

# does not work
# samtools depth AnolisDistichus.addGP.bam | bedtools map -a 10kb_windows.bed -b - -c 3 -o mean > sequencing_depth_in_windows_10kb.txt
# samtools depth AnolisDistichus.addGP.bam | bedtools map -a 50kb_windows.bed -b - -c 3 -o mean > sequencing_depth_in_windows_50kb.txt

# from manual, keeps failing due to memory issues
# bedtools coverage -a 10kb_windows.bed -b AnolisDistichus.addGP.bam > read_depth_bedtools_10kb.txt
# bedtools coverage -a 50kb_windows.bed -b AnolisDistichus.addGP.bam > read_depth_bedtools_50kb.txt

# new attempt with samtools, this sums up genome coverage ACROSS window size, so now need to divide by window size - 10000 for 10kb, 50000 for 50kb
# samtools bedcov 10kb_windows.bed AnolisDistichus.addGP.bam > bedcov_10kb.txt
# samtools bedcov 50kb_windows.bed AnolisDistichus.addGP.bam > bedcov_50kb.txt

# use from vcf
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' -R 10kb_windows.bed /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_allsites2.vcf.gz > bcftools_10kb.txt
bcftools query -f '%CHROM\t%POS\t%INFO/DP\n' -R 50kb_windows.bed /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_allsites2.vcf.gz > bcftools_50kb.txt
```
</details>

<details><summary>coverage.R</summary>
  
```
library(tidyverse)
library(cowplot)
setwd("C:/Users/chfal/OneDrive/Desktop/sliding_window/")
# gotta set our colors EARLY!
colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")


## Chromosome Coverage ####


# grep out scaffold with sed -i 
depth_by_chrm <- read.delim("sequencing_depth_by_chromosome.txt",sep=" ",header=F)

colnames(depth_by_chrm)<-c("Scaffold", "scaffold_length", "bases", "coverage")                        


depth_by_chrm$chr_type <- "Autosome"
depth_by_chrm[depth_by_chrm$Scaffold %in% c(8),]$chr_type <- "Ancient X"
depth_by_chrm[depth_by_chrm$Scaffold %in% c(12),]$chr_type <- "New X"
depth_by_chrm[depth_by_chrm$Scaffold %in% c(13,16),]$chr_type <- "Y"

ggplot(depth_by_chrm[depth_by_chrm$Scaffold %in% c(6,7,8,9,10,11,12,13,14,15,16),], aes(x=Scaffold, y=coverage, fill=as.factor(chr_type))) +
  geom_col() +
  scale_fill_manual(values=colors) +
  scale_x_continuous(n.breaks=17) +
  geom_text(aes(label = coverage), vjust = -0.2) +
  guides(fill=guide_legend(title="Chromosome Type"))


## Window Coverage 10kb and 50kb ####

bedcov_10kb <- read.delim("bedcov_10kb.txt",sep="\t",header=F)

bedcov_50kb <- read.delim("bedcov_50kb.txt",sep="\t",header=F)

colnames(bedcov_10kb)<-c("Scaffold", "start", "end", "bedcov")                        

colnames(bedcov_50kb)<-c("Scaffold", "start", "end", "bedcov")                        


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


bedcov_10kb$chr_type <- "Autosome"
bedcov_10kb[bedcov_10kb$Scaffold %in% c(8),]$chr_type <- "Ancient X"
bedcov_10kb[bedcov_10kb$Scaffold %in% c(12),]$chr_type <- "New X"
bedcov_10kb[bedcov_10kb$Scaffold %in% c(13,16),]$chr_type <- "Y"


bedcov_50kb$chr_type <- "Autosome"
bedcov_50kb[bedcov_50kb$Scaffold %in% c(8),]$chr_type <- "Ancient X"
bedcov_50kb[bedcov_50kb$Scaffold %in% c(12),]$chr_type <- "New X"
bedcov_50kb[bedcov_50kb$Scaffold %in% c(13,16),]$chr_type <- "Y"



ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(7,8,9,10,11,12,13,14,15,16),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values=colors) +
  labs(title="10 kb window coverage") +
  guides(col= guide_legend(title= "Chromosome Type"))



ggplot(bedcov_50kb[bedcov_50kb$Scaffold %in% c(7,8,9,10,11,12,13,14,15,16),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values=colors) +
  labs(title="50 kb window coverage") +
  # facet_grid("chr_type") +
  guides(col= guide_legend(title= "Chromosome Type"))


ancient_x <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(8),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#EE2A7B") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 8 - Ancient X") +
  guides(col= guide_legend(title= "Chromosome Type")) +
  theme_bw()


y_13 <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(13),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#00AEEF") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 13 - Y") +
  guides(col= guide_legend(title= "Chromosome Type")) +
  theme_bw()


y_16 <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(16),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#00AEEF") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 16 - Y") +
  guides(col= guide_legend(title= "Chromosome Type")) +
  theme_bw()

neo_x <- ggplot(bedcov_10kb[bedcov_10kb$Scaffold %in% c(12),], aes(x=bp_cum, coverage, color=chr_type)) +
  geom_point() +
  scale_color_manual(values="#A4509F") +
  geom_smooth(method="loess", color="black") +
  labs(title="Scaffold 12 - Neo - X") +
  guides(col= guide_legend(title= "Chromosome Type")) +
  theme_bw()


plot_grid(ancient_x, neo_x, y_13, y_16, ncol=2)

ggplot(cumulative_count, aes(x=Scaffold, y=bp_add)) +
  geom_point()

```
</details>
