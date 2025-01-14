## GATK Pipeline

Genome Analysis ToolKit. This has also been run by many others in lab - Alyssa, Devon, and Jody - these scripts are modified for Anolis distichus.

GATK is based off of the Illumina Hi-C / Omni-C reads. To start, you need a .bam file. Minimap makes a bam file from PacBio reads, this is **not** the .bam file you are looking for. It is necessary to generate a **new** .bam file.

### FastQC, Trimmomatic, FastQC, BWA

The first step is to use Fastqc, Trimmomatic, and Fastqc again to trim reads and make sure that the read quality is good. Lastly, we use BWA to map the reads to the reference genome to make a BAM file.

<details><summary>fastqc_trim_bwa.sh</summary>

```
!/bin/bash
#SBATCH --partition=cmain
#SBATCH --account=general
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=trim_bwa
#SBATCH --mem=100G
#SBATCH -n 15
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL,END

echo "load any Amarel modules that script requires"
module purge
module load java
module load FastQC
module load samtools
module load bwa

NAME=$1
DATA_DIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/data"
READS_OUTDIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/trimmed_reads"
FASTQC_OUTDIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/fastqc"
BAM_OUTDIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/bam"
BASE_DIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/"

echo "Bash commands for the analysis you are going to run"

echo "##################### fastqc initial quality analysis"
fastqc -t 20 \
${DATA_DIR}/${NAME}_R1_001.fastq.gz \
${DATA_DIR}/${NAME}_R2_001.fastq.gz \
-o ${FASTQC_OUTDIR}

cd ${READS_OUTDIR}

echo ""
echo "##################### trimmomatic"
java -jar /projects/f_geneva_1/programs/trimmomatic/trimmomatic-0.39.jar PE \
-threads 20 -phred33 -trimlog ${readset}_trim.log \
${DATA_DIR}/${NAME}_R1_001.fastq.gz ${DATA_DIR}/${NAME}_R2_001.fastq.gz \
${NAME}_filtered.R1.fq.gz ${NAME}_filtered.unpaired.R1.fq.gz \
${NAME}_filtered.R2.fq.gz ${NAME}_filtered.unpaired.R2.fq.gz \
ILLUMINACLIP:/projects/f_geneva_1/programs/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:4 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:13:20 MINLEN:23

cd ${BASE_DIR}

echo ""
echo "##################### fastqc trimmomatic quality analysis"
fastqc -t 20 \
${READS_OUTDIR}/${NAME}_filtered.R1.fq.gz \
${READS_OUTDIR}/${NAME}_filtered.R2.fq.gz \
-o ${FASTQC_OUTDIR}


echo ""
echo "##################### index and align with BWA"
bwa index /projects/f_geneva_1/chfal/AnoDis1.0.fasta

bwa mem -t 10 /projects/f_geneva_1/chfal/AnoDis1.0.fasta \
${READS_OUTDIR}/${NAME}_filtered.R1.fq.gz \
${READS_OUTDIR}/${NAME}_filtered.R2.fq.gz \
| samtools sort -@10 -o ${BAM_OUTDIR}/${NAME}_bwa_aligned.bam -

echo ""
echo "##################### index genome with samtools - only needs to be done once"
samtools faidx /projects/f_geneva_1/chfal/AnoDis1.0.fasta


echo ""
echo "##################### index all bam files"
samtools index -b ${BAM_OUTDIR}/${NAME}_bwa_aligned.bam


echo ""
echo "done"
```
</details>


### Add and Replace

The next step of the GATK pipeline is to add and replace read groups. This requires the PU, which is sample information that comes from reading the head of the Fastq file. Others in the lab have more complex code because their samples came from multiple lanes; my sample came off of one lane and there is only one sample. The "CORE" name is just the name of the sample. Again, for more complex reads and lanes, it would be important to make sure all the sample names are correct, but since my sample there is only one, it is fine to just hard code it as "AnolisDistichus", which is what I did.

<details><summary>fastqc_trim_bwa.sh</summary>

```
#!/bin/bash
#SBATCH --partition=p_ccib_1
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=addorreplace
#SBATCH --mem=50G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759


echo "load variables"
SAMPLE=AnolisDistichus
BAM=/projects/f_geneva_1/chfal/distichus/gatk_pipeline/bam/DTG-OmniC-2_bwa_aligned.bam
PU=GW2105063837th:3:1101
READ=DTG-OmniC-2_filtered.R1.fq.gz
CORE=AnolisDistichus
OUTDIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline"
BAM_DIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/bam"

echo "run add or replace groups"
gatk AddOrReplaceReadGroups \
-I ${BAM} \
-O ${OUTDIR}/${SAMPLE}.addGP.bam \
-LB library1 -PL illumina -PU ${PU} -SM ${CORE}


echo "index reads"
samtools index ${OUTDIR}/${SAMPLE}.addGP.bam

echo "done"
```

</details>

### Mark Duplicates

The next step of GATK is to mark duplicates. This takes about 60 minutes.


<details><summary>mark_duplicates.sh</summary>

```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=MarkDup
#SBATCH --mem=150G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@rutgers.edu
#SBATCH --mail-type=FAIL,END


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759

echo ""
echo "load variables"

echo ""
echo "run mark duplicates"
gatk MarkDuplicates \
-I /projects/f_geneva_1/chfal/distichus/gatk_pipeline/AnolisDistichus.addGP.bam \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/AnolisDistichus.addGP.marked.bam \
-M /projects/f_geneva_1/chfal/distichus/gatk_pipeline/AnolisDistichus.addGP.metrics.txt \
--REMOVE_DUPLICATES false --ASSUME_SORTED true --CREATE_INDEX true

echo ""
echo "done"

```
</details>


### Make Dictionary

Before we run Haplotype Caller, we need a dictionary made from the reference Fasta file. Haplotype Caller requires this as an input.

<details><summary>make_dictionary.sh</summary>

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=gatk_dict                      # job name for listing in queue
#SBATCH --mem=10G                                # memory to allocate in Mb
#SBATCH -n 2                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,FAIL,END,REQUEUE      # email for the following reasons


echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759

gatk CreateSequenceDictionary -R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/AnoDis1.0.fasta

echo "This is a run"
echo "Now it is done"
```
</details>

### Haplotype Caller

Haplotype caller...calls haplotypes. This is the longest stage in the GATK pipeline. To get this done in 5 days before maintenance, I split up the runs by scaffold list. I did scaffold 1, 2, 3, 4, 5, and 6 all separately. Then I did scaffolds 7-241 (I knew I had 241 scaffolds because of stats.sh).

I did not run HC in the correct mode at first, which led to a problem downstream of not retaining my invariant sites. The invariant sites file was quite literally empty and had no sites in it. There is a different mode that HC should be run in (BP_RESOLUTION), as well as the flag -EMIT-ALL-CONFIDENT-SITES, which makes sure that every site is put into the file. Elsewise, the homozygous regions get put into the header and then they get lost during the combineGVCF step. This is because GATK is typically used for population genomics, to find the SNP and indel sites so it doesn't really matter about the invariant sites.


<details><summary>hap_call_1.sh</summary>

```
#!/bin/bash
#SBATCH --partition=p_ccib_1
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=gatk_HC_split_1
#SBATCH --mem=50G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --time=5-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL,END


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759
module load samtools


echo ""
echo "load variables"
SAMPLE=$1
INDIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller"
GEN_DIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller"
OUTDIR="/projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller"

echo ""
echo "run haplotype caller"
gatk --java-options "-Xms50G -Xmx50g -XX:ParallelGCThreads=2" HaplotypeCaller --native-pair-hmm-threads 2 \
-I ${INDIR}/${SAMPLE}.marked.bam \
-O ${OUTDIR}/${SAMPLE}.1.g.vcf.gz \
-R ${GEN_DIR}/AnoDis1.0.fasta \
-ERC BP_RESOLUTION \
--output-mode EMIT_ALL_CONFIDENT_SITES \
--max-reads-per-alignment-start 0 \
-RF NotDuplicateReadFilter \
-L ${GEN_DIR}/scaffold_1.list


#--exclude=halc068
#--exclude=gpuc001,gpuc002

echo ""
echo "done"
```
</details>

To break the scaffolds up, I had to make files called ```scaffold_1.list```, ```scaffold_2.list```, ```scaffold_3.list```, and so on. Each one was a .list file which had the names of the scaffold. That means for scaffold_1.list through scaffold_6.list, it was just the 
name of that single scaffold. For scaffold_7.list, I listed all 7-241 scaffolds.

Example of scaffold_7.list file is below:
```
scaffold_7
scaffold_8
scaffold_9
...
scaffold_241
```

This also meant I needed to have 7 different hap_call.sh scripts, which I made but are all identical to the first hap_call_1.sh file. Be careful to change the output file (should be x.g.vcf.gz) and sacfold_x.list to whichever one you want to submit.


When Haplotype Caller finished (~2.5 days), I combined all of the GVCFs that the program made.


<details><summary>combine_gvcf.sh</summary>

``` 
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=combineGVCF                      # job name for listing in queue
#SBATCH --mem=30G                                # memory to allocate in Mb
#SBATCH -n 2                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=2-00:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=END,FAIL    # email for the following reasons


echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759

cd /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller

echo ""
echo "#...run combineGVCF"

gatk CombineGVCFs \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis1.0.fasta \
--variant AnolisDistichus.addGP.1.g.vcf.gz \
--variant AnolisDistichus.addGP.2.g.vcf.gz \
--variant AnolisDistichus.addGP.3.g.vcf.gz \
--variant AnolisDistichus.addGP.4.g.vcf.gz \
--variant AnolisDistichus.addGP.5.g.vcf.gz \
--variant AnolisDistichus.addGP.6.g.vcf.gz \
--variant AnolisDistichus.addGP.7.g.vcf.gz \
-O combined_AnoDis.g.vcf.gz


echo ""
echo "index gVCF file"
gatk \
IndexFeatureFile \
-I combined_AnoDis.g.vcf.gz
echo "This is a run"
echo "Now it is done"
```

</details>


We are going to genotype the GVCF file that we just made. Again with HaplotypeCaller, I need to put a flag in GenotypeGVCF to emit all of the sites or they will get erased/condensed whilst genotyping. That flag is called -all-sites TRUE and it returns all of the sites. This will make a large genotype.gvcf file, but in the next step where we filter variants, it will help us by actually creating the invariant sites file lol.


<details><summary>genotypeGVCF.sh</summary>

```
#!/bin/bash

#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=gatk_GC                      # job name for listing in queue
#SBATCH --mem=40G                                # memory to allocate in Mb
#SBATCH -n 2                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=17:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=FAIL,END    # email for the following reasons


echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759

gatk --java-options "-Xmx40g" GenotypeGVCFs \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/combined_AnoDis.g.vcf.gz \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis.genotype.g.vcf.gz
-all-sites TRUE

echo "This is a run"
echo "Now it is done"
```
</details>

We are now going to select all of the variants and filter them by type (SNP, INDEL, INVARIANT).

<details><summary>select_variants.sh</summary>

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=gatk_selectvariants                      # job name for listing in queue
#SBATCH --mem=10G                                # memory to allocate in Mb
#SBATCH -n 2                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=17:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=FAIL,END    # email for the following reasons

echo "load any Amarel modules that script requires"
module purge     # clears out any pre-existing modules
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759


#Select SNPs from GVCF

gatk SelectVariants \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis.genotype.g.vcf.gz \
--select-type-to-include SNP \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/AnoDis.snps.vcf.gz

#Select indels from GVCF

gatk SelectVariants \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis.genotype.g.vcf.gz \
--select-type-to-include INDEL \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/AnoDis.indels.vcf.gz

gatk SelectVariants \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/select_variants/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/select_variants/Test5-6-Combine.genotype.g.vcf.gz \
--select-type-to-include NO_VARIATION \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/select_variants/Test5-6.invariants.vcf.gz


echo "This is a run"
echo "Now it is done"
```
</details>


The code below takes the variants we selected and puts it in a nice table form. This makes it easy to graph these things in R.

<details><summary>variant_table.sh</summary>


```

  GNU nano 2.3.1                                                       File: variant_table.sh

#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --constraint=oarc
#SBATCH --job-name=variant_table
#SBATCH --mem=10G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --time=17:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759


echo ""
echo "load variables"

echo ""
echo "Variant Table SNP"
gatk --java-options "-Xmx10g" \
VariantsToTable \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis.snps.vcf.gz \
-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis_snps.table

echo ""
echo "Variant Table Indel"
gatk --java-options "-Xmx10g" \
VariantsToTable \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis.indels.vcf.gz \
-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/haplotype_caller/AnoDis_indels.table

echo ""
echo "done"
```

</details>

R code makes a PDF of all of these variants and their graphs of distribution based on different properties. 


<details><summary>plot_gatk.R</summary>

```

#Code used here adapted from: https://evodify.com/gatk-in-non-model-organism/
#Code written by Jody Taft and adapted by Alyssa Vanerelli for sagrei data

library('gridExtra')
library('ggplot2')

# Set working directory
setwd("c:/Users/chfal/Downloads/")


#... Generating plots of GATK VariantsToTable Output...#
VCFsnps <- read.csv('AnoDis_snps.table', header = T, na.strings=c("","NA"), sep = "\t") 

## checking filtered table (need to do variant table again on filtered snp file)

VCFindel <- read.csv('AnoDis_indels.table', header = T, na.strings=c("","NA"), sep = "\t")

#retrieve dimensions of object
dim(VCFsnps)
dim(VCFindel)

VCF <- rbind(VCFsnps, VCFindel)
#VCF <- rbind(VCFsnps_asag)

VCF$Variant <- factor(c(rep("SNPs", dim(VCFsnps)[1]), rep("Indels", dim(VCFindel)[1])))
#VCF$Variant <- factor(c(rep("SNPs", dim(VCFsnps_asag)[1])))

snps <- '#F4CCCA'
indels <- '#A9E2E4' 

DP <- ggplot(VCF, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3) +  
  geom_vline(xintercept=c(5,60), color="#FF7F50") + xlim(0, 200)

QD <- ggplot(VCF, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, size=0.7, color="#FF7F50") 

FS <- ggplot(VCF, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) + xlim(0, 100) + geom_vline(xintercept=c(60), size=0.7, color="#FF7F50")

MQ <- ggplot(VCF, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=40, size=0.7, color="#FF7F50")
 + xlim(0, 75)

MQRankSum <- ggplot(VCF, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-12.5, size=0.7, color="#FF7F50") + xlim(-25, 25)

# Some functions are commented out as they incorporate the indel code as well - adapt as needed

SOR <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=3, size=1, color="#FF7F50")

#SOR <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
#  geom_vline(xintercept=c(3, 10), size=1)
#, colour = c(snps))

ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=-8, size=1, color="#FF7F50") + xlim(-30, 30)

#ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) + xlim(-30, 30) +
#  geom_vline(xintercept=-8, size=1)

pdf("AnoDis.pdf")
plot(DP)
plot(QD)
plot(FS)
plot(MQ)
plot(MQRankSum)
plot(SOR)
plot(ReadPosRankSum)
dev.off()
```
</details>


Based on the results of these outputted plots, we will filter variant and invariant sites. We used the standard settings for indels and SNPS, except for the dpeth filter, which we changed to be a cutoff at any value greater than 83.3, which was 2 standard deviations away from the mean.

<details><summary>filter_variants_gatk.sh</summary>

```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=filter_variants
#SBATCH --mem=20G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@rutgers.edu
#SBATCH --mail-type=FAIL,END


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load GATK/4.2.2.0-yc759

echo ""
echo "load variables"


echo ""
echo "Variant Filtration SNPs"
gatk --java-options "-Xmx40g" \
VariantFiltration \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis.snps.vcf.gz \
--filter-expression "QUAL < 0.00 || MQ < 40.00 || SOR > 3.00 || QD < 2.000 || FS > 60.000 || MQRankSum < -12.50 || ReadPosRankSum < -8.00 || ReadPosRankSum > 8.00" \
--filter-name "my_snp_filter" \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered.vcf.gz

echo ""
echo "Extract passing SNPs"
zcat /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered.vcf.gz | grep -E '^#|PASS' > /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_passed.vcf


echo ""
echo "Variant Filtration indels"
gatk --java-options "-Xmx40g" \
VariantFiltration \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis.indels.vcf.gz \
--filter-expression "QUAL < 0.00 || QD < 2.000 || FS > 60.000 || ReadPosRankSum < -8.00 || ReadPosRankSum > 8.00" \
--filter-name "my_indel_filter" \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_indels_filtered.vcf.gz

echo ""
echo "Extract passing indels"
zcat /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_indels_filtered.vcf.gz | grep -E '^#|PASS' > /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_indels_filtered_passed.vcf


echo ""
echo "Mark quality filtered SNPs"
gatk --java-options "-Xmx40g" \
VariantFiltration \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_passed.vcf \
-G-filter "DP < 3 || DP > 83.83211" \
-G-filter-name "depth_filter" \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth.vcf.gz

echo ""
echo "Extract passing SNPs"
zcat /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth.vcf.gz | grep -E '^#|PASS' > /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf


echo ""
echo "Mark quality filtered indels"
gatk --java-options "-Xmx40g" \
VariantFiltration \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_indels_filtered_passed.vcf \
-G-filter "DP < 3 || DP > 83.83211" \
-G-filter-name "depth_filter" \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_indels_filtered_depth.vcf.gz

echo ""
echo "Extract passing indels"
zcat /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_indels_filtered_depth.vcf.gz | grep -E '^#|PASS' > /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_indels_filtered_depth_passed.vcf


echo ""
echo "Mark quality filtered invariants"
gatk --java-options "-Xmx40g" \
VariantFiltration \
-R /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis1.0.fasta \
-V /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis.invariants.vcf.gz \
-G-filter "DP < 3 || DP > 83.83211" \
-G-filter-name "depth_filter" \
-O /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_invariants_filtered_depth.vcf.gz

echo ""
echo "Extract passing invariants"
zcat /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_invariants_filtered_depth.vcf.gz | grep -E '^#|PASS' > /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_invariants_filtered_depth_passed.vcf

echo ""
echo "done"
```
</details>


Next, we ran vcftools_table.sh to generate more statistics in a table format.

<details><summary>vcftools_table.sh</summary>

```

#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=vcftools_table
#SBATCH --mem=10G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load VCFtools/vcftools-v0.1.16-13-yc759

echo ""
echo "load variables"

$SNPS=/projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf
$OUTFILE_SNP=/projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools

echo ""
echo "run vcftools snps"
# vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --freq2 --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools --max-alleles 2           # Calculate allele frequency for each variant. --freq2 just outputs the frequencies without information about the alleles. Max-alleles 2 excludes sites that have more than two alleles

# vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --depth --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools                          # Calculate mean depth of coverage per individual

# vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --site-mean-depth --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools                 # Calculate mean depth per site

# vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --site-quality --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools                    # Calculate site quality score for each site

# vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --missing-indv --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools                    # Calculate proportion of missing data per sample

# vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --missing-site --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools                    # Calculate missing data per site

# vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --het --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools                             # Calculate heterozygosity and inbreeding coefficient per individual


vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf --site-depth --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_vcftools                             # Calculate heterozygosity and inbreeding coefficient per individual


echo ""
echo "done"


```
</details>



<details><summary>plot_vcftools.R</summary>

```
# Checking what parameters to filter by in vcftools
# Code used here was adopted from: https://speciationgenomics.github.io/filtering_vcfs/


# load tidyverse package
library(tidyverse)
library('ggplot2')

#wd <- "/projects/f_geneva_1/alyssa/sagrei/genotype_gvcf/vcftools_tables/"
wd <- "C://Users/chfal/Downloads/"
setwd(wd)


# Alut --------------------------------------------------------------------

# read in data
var_qual <- read_delim("AnoDis_snps_vcftools.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

var_depth <- read_delim("AnoDis_snps_vcftools.ldepth", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

var_miss <- read_delim("AnoDis_snps_vcftools.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

var_freq <- read_delim("AnoDis_snps_vcftools.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

ind_depth <- read_delim("AnoDis_snps_vcftools.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

ind_miss  <- read_delim("AnoDis_snps_vcftools.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

ind_het <- read_delim("AnoDis_snps_vcftools.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)



##### Plots
# Variant Quality 
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#a + theme_light() + xlim(0, 3000) + geom_vline(xintercept=30, size=0.7) + ggtitle("Variant Quality - SNPs")

# Variant mean depth 
b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#b + theme_light()  + ggtitle("Variant Mean Depth - SNPs") + xlim(0,400)

# The mean depth might be misleading from above plot as few variants may be with extremely high coverage
# look closely at mean depth
# summary(var_depth$mean_depth)

# Variant missingness
c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#c + theme_light() + xlim(-0.2, 1) + ggtitle("Variant Missingness - SNPs")

# Check summary data 
summary(var_miss$fmiss)

# Minor allele frequency 
# To find minor allele frequency at each site, we need to use a bit of dplyr based code 
# find minor allele frequency
d <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#d + theme_light() + xlim(-0.05, 0.6) + ggtitle("Minor Allele Freq - SNPs") 

# Check distribution in more detail 
summary(var_freq$maf)

# Mean depth per individual 
e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#e + theme_light() + ggtitle("Mean Depth per Ind - SNPs")

# Proportion of missing data per individual 
f <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#f + theme_light() + ggtitle("Proportion of Missing Data per Ind - SNPs")

# Heterozygosity and inbreeding coefficient per individual 
g <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
#g + theme_light() + ggtitle("Heterozygosity and inbreeding coefficient per Ind - SNPs")

# Once you have assessed these plots - configure your vcftools to filter for variants of interest

# Saving plots as PDF
pdf("AnoDis_vcftools.pdf")
a + theme_light() + xlim(0, 3000) + geom_vline(xintercept=30, size=0.7) + ggtitle("Variant Quality - SNPs")
b + theme_light()  + ggtitle("Variant Mean Depth - SNPs") + xlim(0,400)
c + theme_light() + xlim(-0.2, 1) + ggtitle("Variant Missingness - SNPs")
d + theme_light() + xlim(-0.05, 0.6) + ggtitle("Minor Allele Freq - SNPs") 
e + theme_light() + ggtitle("Mean Depth per Ind - SNPs")
f + theme_light() + ggtitle("Proportion of Missing Data per Ind - SNPs")
g + theme_light() + ggtitle("Heterozygosity and inbreeding coefficient per Ind - SNPs")
dev.off()

```
</details>

Information about each statistic and what it means can be found on Alyssa's sagrei GitHub page.


Next, we run VCFtools filter based on our results of the above step.


<details><summary>vcftools_filter.sh</summary>


```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=vcftools_filter_snp
#SBATCH --mem=10G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load VCFtools/vcftools-v0.1.16-13-yc759

echo ""
echo "load variables"

echo ""
echo "Set filters for vcftools"
MAF=0                                    # set Minor Allele Frequency
QUAL=30                                  # minimum quality score for a site to pass filtering threshold
MIN_DEPTH=3                              # minimum mean depth and minimum depth allowed for a genotype
MAX_DEPTH=83.83211                           # maximum mean depth and maximum depth allowed for a genotype

echo ""
echo "Run vcftools"
# ======
# --remove-indels                       # I left this here for just incase downstream

vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snps_filtered_depth_passed.vcf \
--maf ${MAF} --minQ ${QUAL} \
--min-meanDP ${MIN_DEPTH} --max-meanDP ${MAX_DEPTH} \
--minDP ${MIN_DEPTH} --maxDP ${MAX_DEPTH} --recode --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_snp_vcftools_filtered.vcf

echo ""
echo "done"

vcftools --vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_invariants_filtered_depth_passed.vcf \
--maf ${MAF} --minQ ${QUAL} \
--min-meanDP ${MIN_DEPTH} --max-meanDP ${MAX_DEPTH} \
--minDP ${MIN_DEPTH} --maxDP ${MAX_DEPTH} --recode --out /projects/f_geneva_1/chfal/distichus/gatk_pipeline/variant_filtering/AnoDis_invariant_vcftools_filtered3.vcf

echo ""
echo "done"
```
</details>

Okay, finally we are at the end of the GATK pipeline. We should now be able to combine all of the data to make an all-sites.vcf file. This is the last step!

<details><summary>combine_allsites.sh</summary>

```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=combine_allsites
#SBATCH --mem=10G
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL


echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load VCFtools/vcftools-v0.1.16-13-yc759


echo ""
echo "copy files"
cp /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_snp_vcftools_filtered2.vcf.recode.vcf
cp /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_indels_filtered_depth_passed.vcf
cp /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_invariant_vcftools_filtered3.vcf.recode.vcf

echo ""
echo "bgzip files"
bgzip /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_snp_vcftools_filtered2.vcf.recode.vcf
bgzip /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_indels_filtered_depth_passed.vcf
bgzip /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_invariant_vcftools_filtered3.vcf.recode.vcf


echo ""
echo "index files using tabix"
tabix /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_snp_vcftools_filtered2.vcf.recode.vcf
tabix /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_indels_filtered_depth_passed.vcf
tabix /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_invariant_vcftools_filtered3.vcf.recode.vcf


echo ""
echo "combine using bcftools"
bcftools concat \
--allow-overlaps \
/projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_snp_vcftools_filtered2.vcf.recode.vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_indels_filtered_depth_passed.vcf /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_invariant_vcftools_filtered3.vcf.recode.vcf
-O z -o /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_allsites.vcf.gz


echo ""
echo "index files"
bcftools index -t /projects/f_geneva_1/chfal/distichus/gatk_pipeline/all_sites/AnoDis_allsites.vcf.gz

echo ""
echo "done"

```
