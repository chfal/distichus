# SnpEff

---

## Make Your Own Database 

Installed snpEFF in home directory.

```
cd

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

```
Added database to config.file

```
nano snpEff.config

AnoDis.genome : AnoDis
```

Made folder in data file:

```
mkdir data
cd data
mkdir AnoDis
cd AnoDis
```

Copied over data files: YOU MUST RENAME THEM EXACTLY, OR THE DATABASE WILL FAIL TO BUILD.

```
sequences.fa # this is the assembly fasta file, previously called AnoDis1.0.fasta (but you have to rename it or it will not run)

# this is the CDS that is from gffread in nucleotide sequence, we get this from this line of gffread:
gffread -x Distichus_CDS_Final.fasta -g /projects/f_geneva_1/chfal/distichus/AnoDis1.0.fasta
/projects/f_geneva_1/chfal/distichus/annotation_analyses/Final_Renamed.gff

# we rename Distichus_CDS_Final.fasta to cds.fa or else it will not run 
cds.fa 


# get final protein sequences from filtered GFF, we get this from this line of gffread:

gffread -y Distichus_Protein_Final.fasta -g /projects/f_geneva_1/chfal/distichus/AnoDis1.0.fasta /projects/f_geneva_1/chfal/distichus/annotation_analyses/Final_Renamed.gff

# we rename Distichus_Protein_Final.fasta to protein.fa or else it will not run
protein.fa

genes.gff # this is the gff from dovetail that has been renamed/filtered. we have put SO MUCH WORK into filtering as seen in the filtering scripts. This used to be called Final_Gff but as with the other files you must rename it or it will not run.

```

All names must be the same (ANODIS, not ANNID) - must be named consistently.


Go back to snpeff directory where the java.jar file is; run this line of code in interactive job.

```
interactive # my alias starting an interacative job on cluster
module load java

java -jar snpEff.jar build -gff3 -v AnoDis
```

Database will build. Then you will have lots of stuff like this.

```
-rw-rw-r-- 1 chf29 chf29       19M Oct 12 09:26 snpEffectPredictor.bin
-rw-rw-r-- 1 chf29 chf29       26M Oct 12 09:26 sequence.scaffold_1.bin
-rw-rw-r-- 1 chf29 chf29      2.4M Oct 12 09:27 sequence.scaffold_10.bin
-rw-rw-r-- 1 chf29 chf29      2.1M Oct 12 09:27 sequence.scaffold_11.bin
-rw-rw-r-- 1 chf29 chf29      1.6M Oct 12 09:27 sequence.scaffold_12.bin
-rw-rw-r-- 1 chf29 chf29      876K Oct 12 09:27 sequence.scaffold_14.bin
-rw-rw-r-- 1 chf29 chf29      972K Oct 12 09:27 sequence.scaffold_15.bin
-rw-rw-r-- 1 chf29 chf29      794K Oct 12 09:27 sequence.scaffold_17.bin
-rw-rw-r-- 1 chf29 chf29      342K Oct 12 09:27 sequence.scaffold_18.bin
-rw-rw-r-- 1 chf29 chf29       22M Oct 12 09:27 sequence.scaffold_2.bin
-rw-rw-r-- 1 chf29 chf29       18M Oct 12 09:28 sequence.scaffold_3.bin
-rw-rw-r-- 1 chf29 chf29       17M Oct 12 09:28 sequence.scaffold_4.bin
-rw-rw-r-- 1 chf29 chf29       14M Oct 12 09:28 sequence.scaffold_5.bin
-rw-rw-r-- 1 chf29 chf29      6.3M Oct 12 09:29 sequence.scaffold_6.bin
-rw-rw-r-- 1 chf29 chf29      5.2M Oct 12 09:29 sequence.scaffold_7.bin
-rw-rw-r-- 1 chf29 chf29      2.9M Oct 12 09:29 sequence.scaffold_8.bin
-rw-rw-r-- 1 chf29 chf29      2.4M Oct 12 09:29 sequence.scaffold_9.bin
-rw-rw-r-- 1 chf29 chf29      1.2M Oct 12 09:29 sequence.bin
```

---

## Run SnpEff

Move SNP file into the home directory, run this in an interacive job.


```
java -Xmx8g -jar snpEff.jar AnoDis AnoDis.vcf > AnoDis_out.vcf
```

### KA KS CALCULATION

SNPEFF vcf file goes brr!!!


<details><summary>KA_KS.R</summary>
  
```
# load library
library(dplyr)
library(readr)
library(ape)
library(tibble)
library(tidyr)

# read in what should have been the input file from the KAKS python script

# grep '^#\|missense_variant\|synonymous_variant' AnoDis_out.vcf > mis_syn.txt

setwd("C:/Users/chfal/OneDrive/Desktop/sliding_window")

miss_syn <- read.delim("miss_syn.txt", sep = c("\t"))

# extract variant type from the string that has all the info in it
miss_syn$variant_type <- sapply(strsplit(miss_syn$INFO, split="\\|"), '[', 2)

# extract the gene id from the string that has all the info in it
miss_syn$gene_id <- sapply(strsplit(miss_syn$INFO, split="\\|"), '[', 5)

# data cleaning and widening - group by gene ID, then count the number of each variant. pivot this table wider so that each column has its own name and type of variant.
miss_syn_ka_ks <- miss_syn %>%
  group_by(gene_id) %>%
  count(variant_type) %>%
  pivot_wider(names_from=variant_type, values_from=n)

# data calculating. sum up each synonymous or nonsynonymous variant.
miss_syn_ka_ks2 <- miss_syn_ka_ks %>%
  mutate(sum_synonymous = sum(synonymous_variant, `splice_region_variant&synonymous_variant`, na.rm=T)) %>%
  mutate( sum_nonsynonymous = sum(missense_variant, `missense_variant&splice_region_variant`, `splice_region_variant&intron_variant`, `splice_donor_variant&intron_variant`, stop_gained, `splice_acceptor_variant&intron_variant`, start_lost , `stop_lost&splice_region_variant`, na.rm=T))


# calculate KA_KS
miss_syn_ka_ks2$ka_ks <- miss_syn_ka_ks2$sum_synonymous / miss_syn_ka_ks2$sum_nonsynonymous

```

</details>
