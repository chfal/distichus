library(strawr)
library(ggplot2)
library(dplyr)

s13_14 <- strawr::straw("NONE", "inter_30.hic", "scaffold_13", "scaffold_14", "BP", 1000000, "observed")


r13_14 <- c(13,14,sum(s13_14$counts))

s13_15 <- strawr::straw("NONE", "inter_30.hic", "scaffold_13", "scaffold_15", "BP", 1000000, "observed")

r13_15 <- c(13,15,sum(s13_15$counts))


s13_16 <- strawr::straw("NONE", "inter_30.hic", "scaffold_13", "scaffold_16", "BP", 1000000, "observed")

r13_16 <- c(13,16,sum(s13_16$counts))


s13_17 <- strawr::straw("NONE", "inter_30.hic", "scaffold_13", "scaffold_17", "BP", 1000000, "observed")

r13_17 <- c(16,17,sum(s13_17$counts))


s13_18 <- strawr::straw("NONE", "inter_30.hic", "scaffold_13", "scaffold_18", "BP", 1000000, "observed")

r13_18 <- c(13,18,sum(s13_18$counts))


compare_13 <- rbind(r13_18, r13_14,r13_15,r13_16,r13_17)

compare_13<- as.data.frame.array(compare_13)


contacts_plot <- ggplot(compare_13, aes(x=as.factor(compare_13$V2), y=compare_13$V3)) +
  geom_col(fill=c("#006838","#006838","#006838","#00AEEF","#006838")) +
  labs(x="Scaffold Name", y="Chromatin Contacts with Scaffold 13") +
  theme_bw()

contacts_plot 
ggsave("hic_contacts_plot.pdf", contacts_plot, width=7, height=5)
