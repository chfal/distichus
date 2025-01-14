# Satsuma
Used Satsuma to align genome to itself to look for synteny within sex chromosomes to each other.

## Installed Satsuma using Alyssa github
```
git clone https://git.code.sf.net/p/satsuma/code satsuma-code

cd satsuma-code
make
```

## Select scaffold using faidx
```
samtools faidx AnoDis1.0.fasta scaffold_12 > anodis_12.fasta
samtools faidx AnoDis1.0.fasta scaffold_16 > anodis_16.fasta
```

## use seqtk to get list of all of them BUT the chromosome that is the query

```
seqtk subseq anodis_1_18.fa without_8.txt > without_8.fasta

```

### satsuma.sh

Notes on running: Used high-memory node (mem). But you can only submit for 3 days. So what I did was I chose one scaffold at a time to align. Started with smallest scaffold of interest to test (16).
Pick one scaffold, select it as above, then compare it to scaffolds 1-18 (rest of genome of interest). Wow. Look what reading the manual can do for you.

```
#!/bin/bash
#SBATCH --partition=mem
#SBATCH --account=general
#SBATCH --exclude=gpuc001,gpuc002,halc068
#SBATCH --job-name=satsuma
#SBATCH --mem=700G
#SBATCH --cpus-per-task=50
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --no-requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

work_folder=/projects/f_geneva_1/chfal/distichus/satsuma/satsuma3
out_folder=/projects/f_geneva_1/chfal/distichus/satsuma/satsuma3/out
satsuma_folder=/projects/f_geneva_1/chfal/distichus/satsuma-code

time ${satsuma_folder}/SatsumaSynteny \
-t ${work_folder}/anodis_16.fasta \
-q ${work_folder}/AnoDis1_18.fasta \
-o ${out_folder} -n 50
```

reran satsuma with different parameters - created file with everything BUT that chromosome in question, used that chromosome as the target, and the others as the query, dropped probability to lower with min-prob

```
work_folder=/projects/f_geneva_1/chfal/distichus/satsuma/satsuma3
out_folder=/projects/f_geneva_1/chfal/distichus/satsuma/satsuma3/14_without_14
satsuma_folder=/projects/f_geneva_1/chfal/distichus/satsuma-code

time ${satsuma_folder}/SatsumaSynteny \
-t ${work_folder}/anodis_14.fasta \ 
-q ${work_folder}/without_14.fasta \
-o ${out_folder} -n 50 \
-min_prob 0.999
```


for outputs of satsuma, to count which chromosomes were syntenic hits

 ```
 cut -f 4 satsuma_summary.chained.out | sort | uniq -c
```


# Satsuma R analyses

<details>

<summary>cleos_cursed_satsuma.R</summary>


```
# Cleo's Cursed Satsuma
# this script plots synteny of sex chromosomes (8 - ancient X, 12 - new X, 13 - part of Y, 16 - part of Y)
# adapted from Raul's code

# Beginning Stuff ------------------------------------------------------------

# load libraries
library(ggplot2)
library(dplyr)

setwd("C:/Users/chfal/Downloads/")


# read data

chrom.sizes <- read.table("AnoDis.chrom.sizes", header=FALSE)

chrom.sizes <- chrom.sizes %>%
  mutate(bp_cum = cumsum(as.numeric(V2)))

sex_data <- read.table("sex_out.out", header = FALSE)
head(sex_data)


colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")

# target is what you started with in satsuma, query is what the "target" matched with
names <- c("target_sequence_name","first_target_base","last_target_base","query_sequence_name","first_query_base","last_query_base","identity","orientation")
sex_data <- setNames(sex_data, names)

# no synteny between 8 and 12 which is good, so do not plot

# X to Y Chromosome ------------------------------------------------------------

eight_and_13 <- sex_data[sex_data$target_sequence_name == "scaffold_8" & sex_data$query_sequence_name == "scaffold_13",] 


plot_8_13 <- as.data.frame(c(eight_and_13$first_target_base,(eight_and_13$first_query_base)))

plot_8_13$comp <- c(rep("Scaffold 8", length(eight_and_13$first_target_base)),rep("Scaffold 13", length(eight_and_13$first_query_base)))

plot_8_13$x <- c(1:length(eight_and_13$first_target_base),1:length(eight_and_13$first_query_base))

names(plot_8_13) <- c("position", "comp", "x")

save_8_13 <- ggplot(plot_8_13, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.05, col = "#00AEEF") + 
  ylim(0,max(51939674 )) +
  theme_void() +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177), col="#00AEEF", linewidth=3) +
  geom_segment(aes(x="Scaffold 8", xend="Scaffold 8",
                   y = 0, yend = 51939674), col="#EE2A7B", linewidth=3)


save_8_13

ggsave("8_13.svg", save_8_13)  

# plot 13 and 8

thirteen_and_8 <- sex_data[sex_data$target_sequence_name == "scaffold_13" & sex_data$query_sequence_name == "scaffold_8",] 


plot_13_8 <- as.data.frame(c(thirteen_and_8$first_target_base,(thirteen_and_8$first_query_base)))

plot_13_8$comp <- c(rep("Scaffold 13", length(thirteen_and_8$first_target_base)),rep("Scaffold 8", length(thirteen_and_8$first_query_base)))

plot_13_8$x <- c(1:length(thirteen_and_8$first_target_base),1:length(thirteen_and_8$first_query_base))

names(plot_13_8) <- c("position", "comp", "x")

save_13_8 <- ggplot(plot_13_8, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.05, col = "#00AEEF") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177), col="#00AEEF", linewidth=3) +
  geom_segment(aes(x="Scaffold 8", xend="Scaffold 8",
                   y = 0, yend = 51939674), col="#EE2A7B", linewidth=3) +
  theme_void()

ggsave("13_8.svg", save_13_8)  

# plot 8 and 16

eight_and_16 <- sex_data[sex_data$target_sequence_name == "scaffold_8" & sex_data$query_sequence_name == "scaffold_16",] 

plot_8_16 <- as.data.frame(c(eight_and_16$first_target_base,(eight_and_16$first_query_base)))

plot_8_16$comp <- c(rep("Scaffold 8", length(eight_and_16$first_target_base)), rep("Scaffold 16", length(eight_and_16$first_query_base)))

plot_8_16$x <- c(1:length(eight_and_16$first_target_base),1:length(eight_and_16$first_query_base))

names(plot_8_16) <- c("position", "comp", "x")


save_8_16 <- ggplot(plot_8_16, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.05, col = "blue") + 
  ylim(0,max(51939674)) +
  theme_void() +
  geom_segment(aes(x="Scaffold 8", xend="Scaffold 8",
                   y = 0, yend = 51939674), col="#EE2A7B", linewidth=3) +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528), col="blue", linewidth=3)
  

ggsave("save_8_16.svg",save_8_16)

# plot 16 and 8

sixteen_and_8 <- sex_data[sex_data$target_sequence_name == "scaffold_16" & sex_data$query_sequence_name == "scaffold_8",] 

plot_16_8 <- as.data.frame(c(sixteen_and_8$first_target_base,(sixteen_and_8$first_query_base)))

plot_16_8$comp <- c(rep("Scaffold 16", length(sixteen_and_8$first_target_base)), rep("Scaffold 8", length(sixteen_and_8$first_query_base)))

plot_16_8$x <- c(1:length(sixteen_and_8$first_target_base),1:length(sixteen_and_8$first_query_base))

names(plot_16_8) <- c("position", "comp", "x")

save_16_8 <- ggplot(plot_16_8, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.05, col = "blue") + 
  ylim(0,max(51939674)) +
  theme_void() +
  geom_segment(aes(x="Scaffold 8", xend="Scaffold 8",
                   y = 0, yend = 51939674), col="#EE2A7B", linewidth=3) +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528), col="blue", linewidth=3)

save_16_8

ggsave("save_16_8.svg",save_16_8)


# plot 12 and 13

twelve_and_13 <- sex_data[sex_data$target_sequence_name == "scaffold_12" & sex_data$query_sequence_name == "scaffold_13",] 

plot_12_13 <- as.data.frame(c(twelve_and_13$first_target_base,(twelve_and_13$first_query_base)))

plot_12_13$comp <- c(rep("Scaffold 12", length(twelve_and_13$first_target_base)),rep("Scaffold 13", length(twelve_and_13$first_query_base)))

plot_12_13$x <- c(1:length(twelve_and_13$first_target_base),1:length(twelve_and_13$first_query_base))

names(plot_12_13) <- c("position", "comp", "x")

save_12_13 <- ggplot(plot_12_13, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.1, col = "#00AEEF") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 12", xend="Scaffold 12",
                   y = 0, yend = 33619221), col="#A4509F", linewidth=3) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177), col="#00AEEF", linewidth=3) +
  theme_void()

ggsave("save_12_13.svg",save_12_13)

# plot 13 and 12

thirteen_and_12 <- sex_data[sex_data$target_sequence_name == "scaffold_13" & sex_data$query_sequence_name == "scaffold_12",] 

plot_13_12 <- as.data.frame(c(thirteen_and_12$first_target_base,(thirteen_and_12$first_query_base)))

plot_13_12$comp <- c(rep("Scaffold 13", length(thirteen_and_12$first_target_base)),rep("Scaffold 12", length(thirteen_and_12$first_query_base)))

plot_13_12$x <- c(1:length(thirteen_and_12$first_target_base),1:length(thirteen_and_12$first_query_base))

names(plot_13_12) <- c("position", "comp", "x")

save_13_12 <- ggplot(plot_13_12, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.1, col = "#00AEEF") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 12", xend="Scaffold 12",
                   y = 0, yend = 33619221), col="#A4509F", linewidth=3) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177), col="#00AEEF", linewidth=3) +
  theme_void() 


ggsave("save_13_12.svg",save_13_12)

# plot 12 and 16

twelve_and_16 <- sex_data[sex_data$target_sequence_name == "scaffold_12" & sex_data$query_sequence_name == "scaffold_16",] 


plot_12_16 <- as.data.frame(c(twelve_and_16$first_target_base,(twelve_and_16$first_query_base)))

plot_12_16$comp <- c(rep("Scaffold 12", length(twelve_and_16$first_target_base)),rep("Scaffold 16", length(twelve_and_16$first_query_base)))

plot_12_16$x <- c(1:length(twelve_and_16$first_target_base),1:length(twelve_and_16$first_query_base))

names(plot_12_16) <- c("position", "comp", "x")

save_12_16 <- ggplot(plot_12_16, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.05, col = "blue") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528), col="blue", linewidth=3) +
  geom_segment(aes(x="Scaffold 12", xend="Scaffold 12",
                   y = 0, yend = 33619221), col="#A4509F", linewidth=3) +
  theme_void()

ggsave("save_12_16.svg", save_12_16)
  

# plot 16 and 12

sixteen_and_12 <- sex_data[sex_data$target_sequence_name == "scaffold_16" & sex_data$query_sequence_name == "scaffold_12",] 

plot_16_12 <- as.data.frame(c(sixteen_and_12$first_target_base,(sixteen_and_12$first_query_base)))

plot_16_12$comp <- c(rep("Scaffold 16", length(sixteen_and_12$first_target_base)),rep("Scaffold 12", length(sixteen_and_12$first_query_base)))

plot_16_12$x <- c(1:length(sixteen_and_12$first_target_base),1:length(sixteen_and_12$first_query_base))

names(plot_16_12) <- c("position", "comp", "x")

save_16_12 <- ggplot(plot_16_12, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.05, col = "blue") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528), col="blue", linewidth=3) +
  geom_segment(aes(x="Scaffold 12", xend="Scaffold 12",
                   y = 0, yend = 33619221), col="#A4509F", linewidth=3) +
  theme_void()

ggsave("save_16_12.svg", save_16_12)


# Y to Y Chromosome ------------------------------------------------------------

# plot 13 and 16
thirteen_and_16 <- sex_data[sex_data$target_sequence_name == "scaffold_13" & sex_data$query_sequence_name == "scaffold_16",] 

plot_13_16 <- as.data.frame(c(thirteen_and_16$first_target_base,(thirteen_and_16$first_query_base)))

plot_13_16$comp <- c(rep("Scaffold 13", length(thirteen_and_16$first_target_base)),rep("Scaffold 16", length(thirteen_and_16$first_query_base)))

plot_13_16$x <- c(1:length(thirteen_and_16$first_target_base),1:length(thirteen_and_16$first_query_base))

names(plot_13_16) <- c("position", "comp", "x")

save_13_16 <- ggplot(plot_13_16, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.1, col = "blue") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177)) +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528)) +
  theme_void()

save_13_16 

ggsave("13_16.svg", save_13_16)


# plot 16 and 13

sixteen_and_13 <- sex_data[sex_data$target_sequence_name == "scaffold_16" & sex_data$query_sequence_name == "scaffold_13",] 

plot_16_13 <- as.data.frame(c(sixteen_and_13$first_target_base,(sixteen_and_13$first_query_base)))

plot_16_13$comp <- c(rep("Scaffold 16", length(sixteen_and_13$first_target_base)),rep("Scaffold 13", length(sixteen_and_13$first_query_base)))

plot_16_13$x <- c(1:length(sixteen_and_13$first_target_base),1:length(sixteen_and_13$first_query_base))

names(plot_16_13) <- c("position", "comp", "x")

save_16_13 <- ggplot(plot_16_13, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.1, col = "darkblue") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                 y = 0, yend = 26556177)) +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528)) +
  theme_minimal()

save_16_13

ggsave("16_13.svg", save_16_13)


# SANDBOX -------------------------------------------------------------------------

### THIS DOESN'T ACTUALLY WORK BECAUSE THERE IS NO WAY TO TELL IT WHICH LINES ARE DRAWN FROM WHICH POINTS

test_13_16 <- sex_data[sex_data$target_sequence_name == "scaffold_13" & sex_data$query_sequence_name == "scaffold_16" | sex_data$target_sequence_name == "scaffold_16" & sex_data$query_sequence_name == "scaffold_13",] 

test_2 <- as.data.frame(c(test_13_16$first_target_base,(test_13_16$first_query_base)))

test_2$comp <- c(rep("Scaffold 13", length(test_13_16$first_target_base)),rep("Scaffold 16", length(test_13_16$first_query_base)))

test_2$x <- c(1:length(test_13_16$first_target_base),1:length(test_13_16$first_query_base))

names(test_2) <- c("position", "comp", "x")

ggplot(test_2, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.1, col = "blue") + 
  ylim(0,max(26556177)) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177), col="orange") +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528))

```

</details>


### Reran satsuma with masked genome file to see if there is any difference
This reduces the noise.

Soft masked genome has lower case a, t, c, g to identify where repeat elements are.
Hard masked genome replaces lower case a,t,c,g (regions of repeat elements) with uppercase N.

Reran Satsuma Synteny on what I thought was the hard masked genome but wasn't.


The code below takes the Satsuma synteny output had 8 syntenic to 12 (wasn't before) and gets the regions of scaffold 8 that were syntenic to scaffold 12.
```
# to use samtools faidx, we need to get it into the format
chr:start-stop
for example:
scaffold_8:1-10

# get lists of scaffold 8 and also where it starts
got scaffold_8, got start using cut f-1 and cut f-2 respectively

#now add the start through pasting
paste --delimiters=':' scaffold_8.txt 8_start.txt > 8_start_2.txt

# now add the end 
paste --delimiters='-' 8_start_2.txt 8_end.txt > 8_final.txt

# samtools command
samtools faidx AnoDis1.0.simple_mask.soft.complex_mask.hard.fasta -r 8_final.txt > 8_portions.fasta

```


## How to Make the Figure

Best to do 8/16 and 8/12 in one, then 12/13 and 12/16 separately.

Using the black mouse, put all 8/16 and 16/8 in one place in Adobe Illustrator, move them on top of each other, and lock it. Make sure to reflect across X axis if needed. Then, lock that object.

Highlight all and release ALL clipping masks.

Then, switch to the white mouse and one by one move the paths for 8/13 and 13/8 on top of the previous image.

Then, lock and repeat process for 12/13 and 12/16.

