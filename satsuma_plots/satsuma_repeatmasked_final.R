# Satsuma analysis script
# this script plots synteny of sex chromosomes (8 - ancient X, 12 - new X, 13 - part of Y, 16 - part of Y)
# adapted from Raul Araya Donoso's code

# Beginning ------------------------------------------------------------

# load libraries
library(ggplot2)
library(dplyr)


# read data in

# starting with chromosome sizes 
chrom.sizes <- read.table("AnoDis.chrom.sizes", header=FALSE)

chrom.sizes <- chrom.sizes %>%
  mutate(bp_cum = cumsum(as.numeric(V2)))

# this is the sex chromosome data
sex_data <- read.table("sex_chromosomes.out", header = FALSE)
head(sex_data)


# set colors for all chromosomes
colors <- c("#EE2A7B","#006838","#A4509F","#00AEEF")

# target is what you started with in satsuma, query is what the "target" matched with
names <- c("target_sequence_name","first_target_base","last_target_base","query_sequence_name","first_query_base","last_query_base","identity","orientation")
sex_data <- setNames(sex_data, names)

# no synteny between 8 and 12 which is good, so do not need to plot

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

ggsave("save_8_13_RM.svg", save_8_13)  

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

save_13_8

ggsave("save_13_8_RM.svg", save_13_8)  

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


ggsave("save_8_16_RM.svg",save_8_16)

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

ggsave("save_16_8_RM.svg",save_16_8)


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

ggsave("save_12_13_RM.svg",save_12_13)

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


ggsave("save_13_12_RM.svg",save_13_12)

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

ggsave("save_12_16_RM.svg", save_12_16)


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

save_16_12

ggsave("save_16_12_RM.svg", save_16_12)


# Y to Y Chromosome ------------------------------------------------------------

# plot 13 and 16
thirteen_and_16 <- sex_data[sex_data$target_sequence_name == "scaffold_13" & sex_data$query_sequence_name == "scaffold_16",] 

plot_13_16 <- as.data.frame(c(thirteen_and_16$first_target_base,(thirteen_and_16$first_query_base)))

plot_13_16$comp <- c(rep("Scaffold 13", length(thirteen_and_16$first_target_base)),rep("Scaffold 16", length(thirteen_and_16$first_query_base)))

plot_13_16$x <- c(1:length(thirteen_and_16$first_target_base),1:length(thirteen_and_16$first_query_base))

names(plot_13_16) <- c("position", "comp", "x")

save_13_16 <- ggplot(plot_13_16, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.5, col = "blue") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177), linewidth=3, col="#00AEEF") +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528),linewidth=3, col="blue") +
  theme_void()

save_13_16 

ggsave("13_16_RM.svg", save_13_16)


# plot 16 and 13

sixteen_and_13 <- sex_data[sex_data$target_sequence_name == "scaffold_16" & sex_data$query_sequence_name == "scaffold_13",] 

plot_16_13 <- as.data.frame(c(sixteen_and_13$first_target_base,(sixteen_and_13$first_query_base)))

plot_16_13$comp <- c(rep("Scaffold 16", length(sixteen_and_13$first_target_base)),rep("Scaffold 13", length(sixteen_and_13$first_query_base)))

plot_16_13$x <- c(1:length(sixteen_and_13$first_target_base),1:length(sixteen_and_13$first_query_base))

names(plot_16_13) <- c("position", "comp", "x")

save_16_13 <- ggplot(plot_16_13, aes(x = comp, y = position, group = x)) + geom_line(alpha = 0.1, col = "blue") + 
  ylim(0,max(51939674)) +
  geom_segment(aes(x="Scaffold 13", xend="Scaffold 13",
                   y = 0, yend = 26556177), linewidth=3, col="#00AEEF") +
  geom_segment(aes(x="Scaffold 16", xend="Scaffold 16",
                   y = 0, yend = 15407528), linewidth=3, col="blue") +
  theme_void()

save_16_13

ggsave("save_16_13_RM.svg", save_16_13)
