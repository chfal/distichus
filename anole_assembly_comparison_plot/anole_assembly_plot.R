library(readr)
library(ggplot2)

# read data
anole_assemblies <- read_csv("anole_assemblies_n50s_buscos.csv")

# select only anolis distichus assemblies
ano_dis <- anole_assemblies[1:9,]

# select non-anolis distichus assemblies
all_others <- anole_assemblies[10:18,]

# plot anolis distichus to see if it improved iteratively, which it did
ggplot(ano_dis, aes(x=main_genome_bases, y=non_missing, color=species)) +
  geom_point() +
  # geom_text(hjust=0, vjust=0) +
  xlim(0,350000000) +
  labs(x="N50", y="Non-Missing BUSCOs")

# select anolis distichus final assembly
anodis1.0 <- subset(anole_assemblies, species=="AnoDis1.0")

# plot anolis distichus final assembly in comparison to the other assemblies
anole_assembly_plot <- ggplot(all_others, aes(x=main_genome_bases, y=non_missing)) +
  geom_point(size=4) +
  geom_point(data=anodis1.0, aes(x=main_genome_bases, y=non_missing), color='magenta', size=4) +
  ggrepel::geom_text_repel(data= anodis1.0, aes(label = anodis1.0$species, fontface = "italic"), color='magenta', nudge_y=50) +
  xlim(0,350000000) +
  labs(x="N50", y="Non-Missing BUSCOs") +
  ggtitle("Anole Assembly Statistics") +
  ggrepel::geom_text_repel(aes(label = all_others$species, fontface = "italic"), color='black', nudge_y=-5) +
  ylim(2500,3500) +
  theme_bw()

# view plot
anole_assembly_plot

# save plot
ggsave("anole_assembly_stats.pdf", anole_assembly_plot,width = 5.27, height = 3.5, units = "in")
