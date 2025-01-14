# load libraries
library(dplyr)
library(readr)
library(ape)
library(tibble)

#read gff
gff <- read.gff("PO3115_Anolis_distichus.annotation.gff")
gff$ANNID <- sapply(strsplit(gff$attributes,"="), '[', 2)
gff$ANNID <- sapply(strsplit(gff$ANNID, '-'),'[',1)
gff$ANNID <- sapply(strsplit(gff$ANNID, ';'),'[',1)

# read in a list of already identified genes

list_identified <- tibble(read_lines("maker_identified_list.txt"))

list_identified$ANNID<-list_identified$`read_lines("maker_identified_list.txt")`

# get whatever isn't in the identified genes list
not_identified <- anti_join(gff,list_identified)

# all annotation IDs regardless of Blast search that have more than 3 exons
MultiExons <- not_identified %>%
  filter(not_identified$type %in% c("exon","gene")) %>%
  group_by(ANNID) %>%
  count("exon")

# this count was off by 1, so remove 1 from all counts and return actual count
MultiExons$Actual_Count <- MultiExons$n - 1


# get a list of all of the nonidentified ones that have 3 exons
ANNIDS_3_exons<- MultiExons %>%
  filter(Actual_Count>=3) %>%
  select(ANNID)

# get a list of all of the identified ones
ANNID_identified <- list_identified %>%
  select(ANNID)

# get a final list and write it out
final_list <- full_join(ANNIDS_3_exons,ANNID_identified) %>%
  distinct()

write_csv(final_list, "final_list.csv")
