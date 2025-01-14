# Annotation Analysis

The annotation came basically 90% done by Dovetail. However, to get the final list of genes, we still had to filter some of the gene models that were created. Then, we still had to exclude some annotated elements and then rename the remaining ones in a sequential way. 

```
# install gff3sort
conda create -n gff3sort
conda activate gff3sort
conda install -c bioconda gff3sort
```

<summary>make_annotation.sh</summary>

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002,halc068               # exclude CCIB GPUs
#SBATCH --account=general
#SBATCH --job-name=R                          # job name for listing in queue
#SBATCH --mem=30G                              # memory to allocate in Mb
#SBATCH -n 1                                   # number of cores tao use
#SBATCH -N 2                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=24:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs


# make list of genes that have already been identified by Maker
grep "\sgene\s" PO3115_Anolis_distichus.annotation.gff | grep "=Similar\s" | cut -f 9 | cut -f 2 -d "=" | cut -f 1 -d ";" | cut -f 1 -d "-" | sort | uniq > maker_identified_list.txt


# R CODE NOT RUNNING BECAUSE PACKAGES ARE NOT INSTALLING PROPERLY, BUT SCRIPT LINKED AND CAN BE RUN LOCALLY
# module purge
# module load R
# srun Rscript /projects/f_geneva_1/chfal/distichus/annotation_analyses/r_manipulations.R

# get what is in the final list (made by R code above) and select which annotations we are going to keep
grep -f final_list.csv PO3115_Anolis_distichus.annotation.gff > kept_annotations.gff

# replaces _AED with AED
sed -i 's/_AED/AED/g' kept_annotations.gff

# replaces _eAED with eAED
sed -i 's/_eAED/eAED/g' kept_annotations.gff

# replaces OX%#3D with NCBI:txid
sed -i 's/OX%3D/NCBI:txid/g' kept_annotations.gff

# we are sorting the output gff using gff3sort package

# clear modules and activate the right one
module purge
eval "$(conda shell.bash hook)" 
conda activate gff3sort

# when sorting with gff3sort we will use the chromosome order and then use the natural flag to make sure they are sorted correctly
gff3sort.pl --precise --chr_order natural kept_annotations.gff > gff3sort_sorted.gff

# now we are going to get the list of geneIDs that were kept
grep "maker\sgene\s" gff3sort_sorted.gff | cut -f 9 | cut -f 2 -d "=" | cut -f 1 -d "-" | cut -f 1 -d ";" | cut -f 1 -d "," | uniq > gene_ids_list.txt

# take the previous list, add line numbers and the prefix ANODIS
nl -n rz gene_ids_list.txt | cut -c 2- | sed '1,$ s/^/AnoDis/' > old_new_names.txt

# copy over output file in case of errors
cp gff3sort_sorted.gff gff3sort_renamed.gff

# for loop to replace and renumber
while IFS= read -r line; do
  find_text=$(echo "$line" | awk '{print $2}')
  replace_text=$(echo "$line" | awk '{print $1}')
  sed -i "s/$find_text/$replace_text/g" gff3sort_renamed.gff
done < old_new_names.txt


# copy over sorted file to again to introduce NAME field but to make sure no files get eaten by amarel
cp gff3sort_renamed.gff gff3sort_geneIDs.gff


# make list of names of gene codes
grep "maker\sgene\s" gff3sort_geneIDs.gff | cut -f 9 | grep "Similar to " | cut -f 1 -d":" | sed 's/ID=//g' | sed 's/;Note=Similar to/\t/g' > gff3sort_renamed_ids.txt

# original
# ID=ANN29722;Note=Similar to APOA1

# new, ideally
# ID=ANN29722;Name=APOA1,Note=Similar to APOA1

# for loop to add "NAME" field to each line that has an identified gene code
while IFS= read -r line; do
  ann_name=$(echo "$line" | awk '{print $1}')
  gene_code=$(echo "$line" | awk '{print $2}')
  sed -i "s/$ann_name;Note=Similar to $gene_code/$ann_name;Name=$gene_code,Note=Similar to $gene_code/g" gff3sort_geneIDs.gff
done < gff3sort_renamed_ids.txt

```


<summary>R_manipulations.R</summary>

```
setwd("/projects/f_geneva_1/chfal/distichus/annotation_analyses")

# load library
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
```


<summary>test_gffread.sh</summary>

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --account=general
#SBATCH --job-name=gffread                          # job name for listing in queue
#SBATCH --mem=15G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons


gffread -w gff3sort_renamed_transcripts.fasta -g /projects/f_geneva_1/chfal/distichus/AnoDis1.0.fasta /projects/f_geneva_1/chfal/distichus/annotation_analyses/gff3sort_geneIDs.gff

gffread -x gff3sort_renamed_cds.fasta -g /projects/f_geneva_1/chfal/distichus/AnoDis1.0.fasta /projects/f_geneva_1/chfal/distichus/annotation_analyses/gff3sort_renamed.gff

gffread -y gff3sort_renamed_protein.fasta -g /projects/f_geneva_1/chfal/distichus/AnoDis1.0.fasta /projects/f_geneva_1/chfal/distichus/annotation_analyses/gff3sort_geneIDs.gff
```

# Look for x-linked gene locations

had file of x-linked ACAR genes.

had to add Name= to beginning of line and ',' to end of line because or else i was getting matches of 'ran' in 'transposable elements' and it was too messy and not actually sorting.

```
sed -e 's/^/Name=/' AcarX_genes_list.txt > AcarX_edited.txt
sed -i 's/$/,/' AcarX_edited.txt
grep -f AcarX_edited.txt Final_Sorted_Gff.gff > ACAR_Genes.gff
```


## Used GAG to get information for supplement

Downloaded GAG program from here to my desktop, added it to Amarel using ondemand in the annotation_analyses folder, and unzipped it

```
https://genomeannotation.github.io/GAG/
```

Ran this script. It makes a new output folder in the genomeannotation-GAG-997e384 folder with the output name chosen.

-f is the fasta
-g is the gff
--fix_start_stop fixes the start and stop codons so they actually show up
-o is the output folder

```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=gag
#SBATCH --mem=10G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load python/2.7.17-gc563

cd /projects/f_geneva_1/chfal/distichus/annotation_analyses/genomeannotation-GAG-997e384/

python2.7 gag.py \
-f /projects/f_geneva_1/chfal/distichus/annotation_analyses/AnoDis1.0.fasta \
-g /projects/f_geneva_1/chfal/distichus/annotation_analyses/Final_Sorted_Gff.gff \
--fix_start_stop \
-o gag_distichus2
```

## Blast search for Y-linked sequences

Got Y-linked sequences from FASTA file from the paper published by Marin et al., 2017

Used the ANODIS1.0.FASTA as my database, and searched for where the Y-linked sequences might be

```
makeblastdb -in AnoDis1.0.fasta -out AnoDisDatabase -dbtype nucl -parse_seqids


module purge
module load blast/2.10.1-zz109
blastn -query /projects/f_geneva_1/chfal/distichus/blast/A_carolinensis_Y_transcripts.fa -db /projects/f_geneva_1/chfal/distichus/blast/AnoDisDatabase -evalue 1e-10 -num_threads 32 -out /projects/f_geneva_1/chfal/distichus/blast/Y_transcripts_Blast_out_max_target_1.txt -max_target_seqs 1 -outfmt "6 sscinames qseqid sseqid pident length mismatch evalue bitscore"
```

- sscinames: scientific names if available
- qseqid: query sequence
- sseqid: target sequence
- pident: % identical positions
- length: length of alignment of sequences
- mismatch: # of mismatches
- e-value: expect value, smaller e-value = better match
- bitscore: higher bitscore = more similarity





<details><summary>A Bunch of Stuff That Didn't Work </summary>
   
1. Get a list of unique Maker Annotation IDs that were identified (had a matching BLAST search).
   This is actually two steps. First, we need to grep for the words "gene" and "similar" in the original annotation, and send that to a .txt file (first line of code in code block). This .txt file will then be filtered because we only want the ANNIDs that match, which is accomplished by the one-two punch of cut, sort, and uniq commands. (second line)

   This step was accomplished by this line of code:
   ```
    # this code gets the already identified genes by MAKER pipeline
    grep "\sgene\s" PO3115_Anolis_distichus.annotation.gff | grep "=Similar\s" > identified_genes_maker.txt

    # this code gets the list of ANNIDs from the above output
   cut -f 9 identified_genes_maker.txt | cut -f 2 -d "=" | cut -f 1 -d ";" | cut -f 1 -d "-" | sort | uniq > maker_identified_list2.txt

   ```
2. Get a list of the unidentified genes, but that which have 3 or more exons. Merge it with the before list so that there is a long list of distinct Annotation IDs. This was accomplished by some of the worst R code I have ever written.

  ```
  library(tidyverse)
  library(ape)
  
  setwd("C:/Users/chfal/Downloads/")
  
  # get whole annotation
  gff <- read.gff("PO3115_Anolis_distichus.annotation.gff")
  gff$ANNID <- sapply(strsplit(gff$attributes,"="), '[', 2)
  gff$ANNID <- sapply(strsplit(gff$ANNID, '-'),'[',1)
  gff$ANNID <- sapply(strsplit(gff$ANNID, ';'),'[',1)
  
  # read in a list of already identified genes
  
  ## this list was made using this code in BASH
  
  # grep "\sgene\s" PO3115_Anolis_distichus.annotation.gff | grep "=Similar\s" > identified_genes_maker.txt
  
  list_identified <- tibble(read_lines("maker_identified_list2.txt"))
  
  list_identified$ANNID<-list_identified$`read_lines("maker_identified_list2.txt")`
  
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
  
  final_list <- full_join(ANNIDS_3_exons,ANNID_identified) %>%
     distinct()
  
  write_csv(final_list, "final_list.csv")  
  
  ```
  This final_list.csv has a length of 26373, which is the number of gene models we want (yay!)

3. Use the list of matching Annotation IDs to filter the original gff file to get only the selected ones we want using GREP.
   ```
   grep -f final_list.csv PO3115_Anolis_distichus.annotation.gff > kept_annotations.gff
   ```
   However, these are not sequential because we left out some, so we need to rename them, which brings us to step 4.

4. Make the MAKER MAP file and rename the selected annotations.

  SCRIPT 1:
  This script we have the maker_map_ids function, with the ANODIS prefix. The justify argument adds that amount of leading zeros, in this case 5 so we could have enough for our 26,373 genes. This makes the map file which maps the ANNID to the new renamed and renumbered ID. 

  ```
  #!/bin/bash
  #SBATCH --partition=cmain
  #SBATCH --exclude=gpuc001,gpuc002
  #SBATCH --job-name=maker_map_id
  #SBATCH --mem=20G
  #SBATCH -n 10
  #SBATCH -N 1
  #SBATCH --time=3-00:00:00
  #SBATCH --requeue
  #SBATCH --mail-user=chf29@scarletmail.rutgers.edu
  #SBATCH --mail-type=FAIL
  
  
  cd /projects/f_geneva_1/chfal/distichus/annotation_analyses
  
  echo "load modules"
  module purge
  module load singularity/3.1.0
  
  MAKER_IMAGE=/projects/f_geneva_1/programs/maker:2.31.11-repbase.sif
  
  echo ""
  echo "map maker IDs to numerical IDs with a specified prefix"
  singularity exec $MAKER_IMAGE  maker_map_ids --prefix ANODIS_ --justify 5 \
  kept_annotations.gff > test.maker.map
  
  echo ""
  echo "done"
  ```

  SCRIPT 2:
  This script takes the MAKER MAP file and renames the headers in the gff and the fasta files of the transcripts (both protein and nucleotide fasta).

  ```
  #!/bin/bash
  #SBATCH --partition=cmain
  #SBATCH --exclude=gpuc001,gpuc002
  #SBATCH --job-name=maker_change_ids
  #SBATCH --mem=20G
  #SBATCH -n 10
  #SBATCH -N 1
  #SBATCH --time=3-00:00:00
  #SBATCH --requeue
  #SBATCH --mail-user=chf29@scarletmail.rutgers.edu
  #SBATCH --mail-type=FAIL
  
  
  cd /projects/f_geneva_1/chfal/distichus/annotation_analyses
  
  echo "load modules"
  module purge
  module load singularity/3.1.0
  
  MAKER_IMAGE=/projects/f_geneva_1/programs/maker:2.31.11-repbase.sif
  
  echo ""
  echo "change IDs of a given GFF/FASTA file based off of the map file created"
  
  
  echo "Master GFF"
  singularity exec $MAKER_IMAGE map_gff_ids \
  test.maker.map kept_annotations2.gff
  
  
  echo "Protein FASTA"
  singularity exec $MAKER_IMAGE map_fasta_ids \
  test.maker.map distichus_transcripts_renamed.fasta
  
  
  echo "Transcript fasta"
  singularity exec $MAKER_IMAGE map_fasta_ids \
  test.maker.map distichus_transcripts_protein_renamed.fasta
  
  echo ""
  echo "done"
  
  ```

### Sed and grep final edits

Copy final file in case you screw up. (to test.gff in this case). These are the operations I did on it to make it better. These are simple remove character operations.

```
# replaces ANODIS_ with ANODIS
sed -i 's/ANODIS_/ANODIS/g' test_3.gff

# replaces the word Alias up until the semicolon with nothing (removes)
sed 's/Alias[^;]*;//' test_2.gff > test_3.gff

# replaces _AED with AED
sed -i 's/_AED/AED/g' test_3.gff

# replaces _eAED with eAED
sed -i 's/_eAED/eAED/g' test_3.gff

# replaces OX%#3D with NCBI:txid
sed -i 's/OX%3D/NCBI:txid/g' test_3.gff

```

## Inserting name into name field - this comes from maker pipeline, so again, a bunch of stuff that didn't work 

The previous iteration of the files looked like this.

```
scaffold_11     maker   gene    1345255 1346681 .       -       .       ID=AnoDis24241;Name=;Note=Similar to APOA1: Apolipoprotein A-I (Anas platyrhynchos OX%3D8839);
```

However, we wanted the name field to be filled out with the actual name of the protein like so.

```
scaffold_11     maker   gene    1345255 1346681 .       -       .       ID=AnoDis24241;Name=APOA1;Note=Similar to APOA1: Apolipoprotein A-I (Anas platyrhynchos OX%3D8839);

```

This turned out to be a little more work intensive than we thought.

```
# This grep command gets the ANODIS ID and then the protein code if it's all one letter.
grep "maker\sgene\s" test_3.gff | cut -f 9 | grep "Similar to " | cut -f 1 -d":" | sed 's/ID=//g' | sed 's/;Name=;Note=Similar to /\t/g' | grep -v " " > geneIDs_Names_list.txt 
```

This code runs a bunch of sed commands in R. It is very cursed. It is very slow. It is most certainly not the optimal way to do it.
```
list_to_loop <- read.delim("geneIDs_Names_list.txt",header = F)


for (i in 1:length(list_to_loop[[1]])){
  cmd_string <- NULL
  cmd_string <- paste("sed -i 's/", list_to_loop[i,1], ";Name=;/", list_to_loop[i,1], ";Name=", list_to_loop[i,2], ";/g' test_3.gff" , sep="")
  system(cmd_string)
}
```


# Don't actually run this, but this would rename based on non-taxid

```
# this code gets all of the things that didn't match from the output gff of the above step  and still have a remaining name issue
grep "maker\sgene\s" test_4.gff | grep "Name=;" | grep "Similar to" > no_match.gff

# this code creates another text file to fix the unmatched ones
cut -f 9 no_match.gff | grep "AnoDis" | cut -f2- -d "=" | cut -f 1  -d "(" > fix_unmatched.txt

# selects only the part we want to create a similar file to the above file (geneIDS_names_list)
sed -i 's/;Name=;Note=Similar to /\t/g' fix_unmatched.txt

# had an issue with empty whitespade, this removes
sed -i 's/\ *$//g' fix_unmatched.txt

# grep "Similar to[^:(]*" no_match.gff

```

Reran the R loop again

```
new_list_to_loop <- read.delim("fix_unmatched.txt", header=F)


for (i in 1:length(new_list_to_loop[[1]])){
  cmd_string <- NULL
  cmd_string <- paste("sed -i 's/", new_list_to_loop[i,1], ";Name=;/", new_list_to_loop[i,1], ";Name=", new_list_to_loop[i,2], ";/g' test_5.gff" , sep="")
  system(cmd_string)
}

```

There were six errors that I encountered and they were all because of issues with the "/" character. I will do these by hand because at least there aren't 432 of them.

```
scaffold_2      maker   gene    10701274        10712531        .       -       .       ID=AnoDis05152;Name=;Note=Similar to DLA class I histocompatibility antigen%2C A9/A9 alpha chain (Canis lupus familiaris OX%3D9615);
scaffold_11     maker   gene    21966035        21973515        .       -       .       ID=AnoDis24584;Name=;Note=Similar to BPTI/Kunitz domain-containing protein (Fragment) (Haliotis asinina OX%3D109174);
scaffold_10     maker   gene    26759742        26809609        .       -       .       ID=AnoDis24168;Name=;Note=Similar to CNK3/IPCEF1: CNK3/IPCEF1 fusion protein (Homo sapiens OX%3D9606);
scaffold_3      maker   gene    99345808        99371610        .       +       .       ID=AnoDis11059;Name=;Note=Similar to Putative lysosomal acid lipase/cholesteryl ester hydrolase (Crotalus adamanteus OX%3D8729);
scaffold_6      maker   gene    57869819        57875268        .       +       .       ID=AnoDis20460;Name=;Note=Similar to Gastrin/cholecystokinin-like peptide (Gallus gallus OX%3D9031);
scaffold_4      maker   gene    243742487       243778167       .       +       .       ID=AnoDis16561;Name=;Note=Similar to Sodium/hydrogen exchanger 8 (Gallus gallus OX%3D9031);
scaffold_6      maker   gene    956494  988323  .       -       .       ID=AnoDis19728;Name=;Note=Similar to Casein kinase II subunit alpha' (Gallus gallus OX%3D9031)
```


can verify it worked with this
```
grep "maker\sgene\s" test_5.gff | grep "Name=;" | grep "Similar to"

and this

grep "maker\smRNA\s" test_5.gff | grep "Name=;" | grep "Similar to"

```

# LETS MCFREAKING LOSE IT

I sorted it with this:

```
module load singularity          # Load singularity module

#Then:

# Get the chosen AGAT container version
singularity pull docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0

# run the container
singularity run agat_1.0.0--pl5321hdfd78af_0.sif

agat_convert_sp_gxf2gxf.pl --gff test_4.gff -o sorted.gff
```


```
# the GFF is sorted. now we are going to take those names IN THE EXACT ORDER WE HAD THEM.

cut sorted_gff.gff -f 9 | cut -f 2 -d "=" | cut -f 1 -d "-" | cut -f 1 -d ";" | cut -f 1 -d "," | uniq > new_order.txt 

# then cut off the top of the new order because it has a .gff file name which would cause the numbering to be off - this is just a nano command

# this command adds line numbers with leading zeros to our new_order field
nl -n rz new_order.txt > new_order2.txt

# remove first trailing zero because that was still in there, and it was annoying me
 cut -c 2- new_order2.txt > new_order_3.txt

# append the prefix NEWNAME - will change this back to AnoDis once it actually works (if it does ever work).
sed '2,$ s/^/NEWNAME/' new_order_3.txt > new_order_4.txt
```

This makes a .txt dataframe with
NEW ID | OLD ID

Weirdly enough, the first few are okay, but if you run tail on this file, you can see the differences.

```
NEWNAME26364    AnoDis26360
NEWNAME26365    AnoDis26361
NEWNAME26366    AnoDis26351
NEWNAME26367    AnoDis26321
NEWNAME26368    AnoDis26343
NEWNAME26369    AnoDis26367
NEWNAME26370    AnoDis26373
NEWNAME26371    AnoDis26365
NEWNAME26372    AnoDis26362
NEWNAME26373    AnoDis26357
```

Code (from chatgpt) to go through and rename sequentially.

```
while IFS= read -r line; do
  find_text=$(echo "$line" | awk '{print $2}')
  replace_text=$(echo "$line" | awk '{print $1}')
  sed -i "s/$find_text/$replace_text/g" test_maker_rename.gff
done < new_order_4.txt
```

There are 809303 instances of the word AnoDis in the original file. There should be 809303 instances of NEWNAME in the new file. Which there are :)
```
grep -o -i "AnoDis" sorted_gff.gff | wc -l

grep -o -i "NEWNAME" more_time.gff | wc -l
```

# What Fresh Hell Is This?

Is it possible to mess with a GFF too much? Possibly. What happened was after the GFF was sorted, the gffread -w output of the file was INSANELY large like 252MB and had added 37K genes. 

This is wrong obviously.

So I HAD TO GO BACK AND RESORT IT. I think I maybe messed with it too much and that caused it? So I tried to find the issue. I wrote a bunch of test scripts.

### Tried gffsort, which didn't work. 
```
conda create -n gff3sort

conda install -c bioconda gff3sort

conda activate gff3sort

# make a FIFTH TEST FILE
cp test_4.gff test_5.gff

gff3sort.pl --precise --chr_order natural test_5.gff > gff3sort_sort.gff

# get order made in AGAT
cut sorted_gff.gff -f 9 | cut -f 2 -d "=" | cut -f 1 -d "-" | cut -f 1 -d ";" | cut -f 1 -d "," | uniq > agat_order.txt 

# get order made in gff3sort
cut  gff3sort_sort.gff.gff -f 9 | cut -f 2 -d "=" | cut -f 1 -d "-" | cut -f 1 -d ";" | cut -f 1 -d "," | uniq > gff3sort_order.txt 

# for some reason there were duplicates in this
 awk '!x[$0]++' gff3sort_order.txt > gff3sort_ordered_dupes_removed.txt

diff gff3sort_order.txt gff3sort_ordered_dupes_removed.txt
```

So I resorted it from the step of kept_annotations.gff, and then renamed. We'll see what happens next.

```
# loaded up agat, sorted kept_annotations.gff, tested it with gffread, the output was okay.

# repeated above steps to rename ANNIDs to NEWNAME IDs sequentially.

# used sed and grep to change out _AED, _eAED, and rename txid

# used sed and grep script to add "name" field to those that had it

grep "maker\sgene\s" kept_annotations_sorted_renamed.gff | cut -f 9 | grep "Similar to " | cut -f 1 -d":" | sed 's/ID=//g' | sed 's/;Note=Similar to/\t/g' > kept_annotations_sorted_renamed_geneids_names.txt
```
</details>


BLAST scaffold 8 and scaffold 12 transcripts against the rest of the genome (8, 12, 13 and 16)

```
makeblastdb -in without_sex_scaffolds.fasta -out without_sex_scaffolds -dbtype nucl -parse_seqids
```


## Blast search of genes in 8 and 12 against the whole Y scaffolds

```
grep '\<scaffold_8\>' Final_Sorted_Gff.gff > scaffold_8.gff
gffread -w scaffold_8_transcripts.fasta -g anodis_8.fasta scaffold_8.gff

grep '\<scaffold_12\>' Final_Sorted_Gff.gff > scaffold_12.gff
gffread -w scaffold_12_transcripts.fasta -g anodis_12.fasta scaffold_12.gff

```

Then also get scaffold 13 and 16 from the above subseq results


```
makeblastdb -in 13_16_blastdb.fasta -out blastdb1316 -dbtype nucl -parse_seqids

```

```
module purge
module load blast/2.10.1-zz109
blastn -query /projects/f_geneva_1/chfal/distichus/blast2/scaffold_8_transcripts.fasta -db /projects/f_geneva_1/chfal/distichus/blast2/blastdb1316 -evalue 1e-10 -num_threads 32 -out /projects/f_geneva_1/chfal/distichus/blast2/blast_8_against_13_16.txt -max_target_seqs 5 -outfmt "6 sscinames qseqid sseqid pident length mismatch evalue bitscore"

blastn -query /projects/f_geneva_1/chfal/distichus/blast2/scaffold_12_transcripts.fasta -db /projects/f_geneva_1/chfal/distichus/blast2/blastdb1316 -evalue 1e-10 -num_threads 32 -out /projects/f_geneva_1/chfal/distichus/blast2/blast_12_against_13_16.txt -max_target_seqs 5 -outfmt "6 sscinames qseqid sseqid pident length mismatch evalue bitscore"

```
