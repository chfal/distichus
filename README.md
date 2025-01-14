# Assembling Anolis distichus genome
## With PACBIO long reads and Omni-C Hi-C Data
## Using the Rutgers OARC Amarel Computing Cluster
### Cleo Falvey
---

### Assembling with HiFi 
First, we need to make a conda environment in the home directory for the assembler HifiAsm and then we download the package into the conda environment. We can do this by:

```
conda create --name hifiasm
conda activate hifiasm
conda install -c bioconda hifiasm
```

The next step is to run the assembly for HiFiAsm. This is done with a Bash script.

  <details><summary>run_hifiasm.sh</summary>
  
  ```
  #!/bin/sh
  #SBATCH --partition=p_ccib_1   # which partition to run the job
  #SBATCH --job-name=hifiasm # job name for listing in queue
  #SBATCH --mem=192000 # memory to allocate in Mb
  #SBATCH -n 32 # number of cores to use
  #SBATCH -N 1 # number of nodes the cores should be on, 1 means all cores on same node
  #SBATCH --time=7-00:00:00 # maximum run time days-hours:minutes:seconds
  #SBATCH --mail-user=chf29@scarletmail.rutgers.edu # email address to send status updates
  #SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE # email for the following reasons       
  
  
  # now load any Amarel modules that your script will require
  module purge # clear existing modules
  eval "$(conda shell.bash hook)" # enable slurm to see conda environments
  conda activate hifiasm # load conda environment (which is hifiasm)                                                                                                    
  
  # Bash commands for the analysis you are going to run (we are calling two files - both fastqs from the two smart cells that were run)
  hifiasm -o distichus_hifi_only_32core.asm -t 32 /projects/f_geneva_1/data/denovo_genomes/distichus/m64086e_220711_172529.hifi_reads.fastq.gz /projects/f_geneva_1/data/denovo_genomes/distichus/m64086e_220726_054912.hifi_reads.fastq.gz
  ```
 The command syntax is: hifiasm -o [output file prefix name] , -t [num threads], and the Fastq reads.
  
 </details>

### Make Fasta Files from HifiAsm output
  
  HifiAsm makes a lot of output files, the way to interpret them is [here](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html). However, we need to make the files that come from .gfa format into .fa format. We can do this by running an interactive job. The call for the interactive job is this:
  ```
  srun -p cmain --job-name "hifi" --cpus-per-task 1 --mem-per-cpu 100g --time 01:00:00 --pty bash
  ```
  
  This "logs into" a node for you that you can do simple jobs on. You will know you are in the interactive job because the left side of your screen will say that you are not on amarel1 but on a different node (for example, halc015). Once in the interactive job, we changed all the files that ended in .gfa to end in .fa using this command found [here](https://hifiasm.readthedocs.io/en/latest/faq.html#how-do-i-get-contigs-in-fasta).
  
  ```
 awk '/^S/{print ">"$2;print $3}' file_name.gfa > file_name.fa

  ```
  
  ### Run Statistics
  
  You can either run statistics in an interactive job (using same command above to create the interactive job) with the addition of these commands, or you can submit a batch script (see Alyssa's script in grahami to learn how to do this).
  
  ```
  module purge
  module load java
  stats.sh Xmx8g in=/projects/f_geneva_1/chfal/distichus/hifiasm/distichus_hifi_only.asm.bp.p_ctg.gfa.fa
  ```
  In the last line of code where stats is actually called, the first part is calling the script stats.sh, the second part (Xmx8g) is asking for additional memory in Java (or else you will get the Java Heap Error - Running out of Memory on the virtual machine), and the last is the input file which is the direct path to the data.
  
  
  ### Renaming contigs
  
  I tried to run Busco, but this was unsuccessful because apparently some of the contigs were named the same thing and that was throwing an error. I tried running Devon's version of the script where he had all of the files named correctly but because mine were weird and I could not get it to replicate, I ended up hard coding everything. :( Here is the version of the script I used, and I ran it three times, once for the whole haplotype, and once for haplotype 1 and haplotype 2.
  

<details><summary>run_rename.sh</summary>
  
  ```
  #!/bin/bash
  #SBATCH --partition=p_ccib_1                    # which partition to run the job, options are in the Amarel guide
  #SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
  #SBATCH --job-name=sort                         # job name for listing in queue
  #SBATCH --output=/projects/f_geneva_1/chfal/distichus/rename/slurm-%j-%x.out
  #SBATCH --mem=100G                               # memory to allocate in Mb
  #SBATCH -n 1                                    # number of cores to use
  #SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
  #SBATCH --time=40:00:00                         # maximum run time days-hours:minutes:seconds
  #SBATCH --requeue                               # restart and paused or superseeded jobs
  #SBATCH --mail-user=chf29@scarletmail.rutgers.edu          # email address to send status updates
  #SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons

  module purge
  module load java

  echo "Bash commands for the analysis you are going to run"
  echo ""
  echo "##### sort by size and rename sequence"
  sortbyname.sh -Xmx100g in=/projects/f_geneva_1/chfal/distichus/hifiasm/distichus_hifi_only.asm.bp.hap2.p_ctg.gfa.fa out=distichus_sorted_hap2.fa length descending

  # add dummy to beginning
  printf ">dummy\nNNN\n" | cat - distichus_sorted_hap2.fa > temp && mv temp distichus_sorted_hap2.fa


  #rename
  rename.sh in=distichus_sorted_hap2.fa out=distichus_sorted_renamed_hap2.fa prefix=scaffold -Xmx10g fastawrap=500000000

  # remove single sequence entry from multifasta
  cat distichus_sorted_renamed_hap2.fa | awk '{if (substr($0,1) == ">scaffold_0") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}' >   distichus_fixed_hap2.fa
  rm distichus_sorted_hap2.fa
  rm distichus_sorted_renamed_hap2.fa
  ```
  
  </details>

    
  ### Run Busco
  
  To run Busco, I followed Alyssa's method of creating a new conda environment in my home directory for Busco, downloading the program, and activating it. However, when I tried to run the program some of my contigs were named the same thing so I had to rename the scaffolds. I also could not figure out how to run the script with variables so I hard coded all three of the busco runs and submitted them. Luckily I only had to run this three times. In the future, I will work on automating this.

<details><summary>run_busco.sh</summary>

  ```
  #!/bin/bash


  #SBATCH --partition=p_ccib_1                    # which partition to run the job, options are in the Amarel guide
  #SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
  #SBATCH --job-name=busco                        # job name for listing in queue
  #SBATCH --output=/projects/f_geneva_1/chfal/distichus/busco/slurm-%j-%x.out
  #SBATCH --mem=150G                               # memory to allocate in Mb
  #SBATCH -n 16                                   # number of cores to use
  #SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
  #SBATCH --time=7-00:00:00                         # maximum run time days-hours:minutes:seconds
  #SBATCH --requeue                               # restart and paused or superseeded jobs
  #SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
  #SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons

  echo "load any Amarel modules that script requires"
  module purge
  eval "$(conda shell.bash hook)"
  conda activate busco

  echo "bash commands for the analysis to run"
  busco -i /projects/f_geneva_1/chfal/distichus/rename/distichus_fixed_hap2.fasta -c 16 -l vertebrata_odb10 -o distichus_hap2 -m genome
  ```
</details>

### Assembly with Hi-C

I used HifiAsm to make a genome assembly using our Omni-C data. This is a separate assembly. I copied the same script as the Hifiasm script but changed the last line to this, below. -h1 and -h2 were each of the Omni-C data files respectively and then the last two calls are paths to the file for the two smart cell Pacbio Data. I then repeated the process of making the files into Fasta files from the HifiAsm output, running Stats and Busco on this assembly.

```
hifiasm -o distichus_hifi_hic.asm -t 32 --h1 DTG-OmniC-2_R1_001.fastq.gz --h2 DTG-OmniC-2_R2_001.fastq.gz /projects/f_geneva_1/data/denovo_genomes/distichus/m64086e_220711_172529.hifi_reads.fastq.gz /projects/f_geneva_1/data/denovo_genomes/distichus/m64086e_220726_054912.hifi_reads.fastq.gz
```
### Using Falcon

I downloaded the bioconda package pbassembly into its own conda environment to attempt to make Falcon run. Then I started trying to run the pipeline. It was a bit frustrating because part of the deprecated tutorial was still live so I got confused.

I am following this Github Repository:

https://github.com/PacificBiosciences/pb-assembly#general-overview (overall)
https://github.com/PacificBiosciences/pbbioconda/wiki/Assembling-HiFi-data:-FALCON-Unzip3 (unzipping)

https://github.com/phasegenomics/FALCON-Phase/blob/master/README.md (falcon phase)
https://github.com/PacificBiosciences/pb-assembly#configuration


```
cat {movie1}.fastq {movie2}.fastq ... > {falcon_unzip_input}.fastq
samtools fqidx {falcon_unzip_input}.fastq
```


### Using Salsa2

To use Salsa2 to scaffold you need a few things. First you need BWA Alignment where you align your Ommi-C or Hi-C reads to the current assembly. Then you need to turn the output from BWA into a .BED file. Then you need an indexed length .fasta.fai file from Samtools.

### BWA
I first had to use BWA to align my Hi-C reads to my assembly (I was using the HifiAsm assembly with Hi-C because that's the one I had in hand.) The BWA Script is here. The problem was that I was missing the -t 10 command which uses all the threads, so my BWA ran for a straight week on one core which meant really nothing. The -5SPM parameter disables paired-end reads (comes from [this](https://bioinformatics.stackexchange.com/questions/5076/what-is-the-correct-way-to-map-hi-c-data-with-bwa-mem) Stack Exchange page).


<details><summary>run_bwa.sh</summary>

```
#!/bin/bash
#SBATCH --partition=p_ccib_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --account=general
#SBATCH --job-name=bwa                          # job name for listing in queue
#SBATCH --mem=100G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons

echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module load samtools                            # load any modules needed
module load bwa

echo "Bash commands for the analysis you are going to run"

echo "##################### index and align with BWA"
# bwa index /projects/f_geneva_1/chfal/distichus/rename/distichus_fixed_hifi_hic.fasta

bwa mem -t 10 -5SPM /projects/f_geneva_1/chfal/distichus/rename/distichus_fixed_hifi_hic.fasta \
/projects/f_geneva_1/data/denovo_genomes/distichus/KU_data/raw_reads/DTG-OmniC-2_R2_001.fastq.gz \
/projects/f_geneva_1/data/denovo_genomes/distichus/KU_data/raw_reads/DTG-OmniC-2_R1_001.fastq.gz | \
samtools sort -@10 -o /projects/f_geneva_1/chfal/distichus/salsa2/bwa/mapped_hifi_hic.bam -

echo "change user group of files created"
chgrp -R g_geneva_1 /projects/f_geneva_1/chfal/distichus/salsa2/bwa             # changes group of all files in listed directory

```
</details>

### Bam To Bed Output:

<details><summary>run_bam2bed.sh</summary>

```
#!/bin/bash
#SBATCH --partition=p_geneva_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --job-name=bam2bed                       # job name for listing in queue
#SBATCH --mem=512G                               # memory to allocate in Mb
#SBATCH -n 1                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=7-00:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons

module purge
module load bedtools2/2.25.0

bamToBed -i /projects/f_geneva_1/chfal/distichus/salsa2/bwa/mapped_hifi_hic.bam > output_alignment.bed
sort -k 4 output_alignment.bed > tmp && mv tmp output_alignment.bed

```

</details>

### Samtools
Next, I used Samtools to index my contig file which did work. This code could be run in an interactive job but I ended up writing a slurm script for it. Here is the command for that.
```
samtools faidx /projects/f_geneva_1/chfal/distichus/rename/distichus_fixed_hifi_hic.fasta
```

### Actually Running Salsa

The things you need to run Salsa are:
1) assembly
2) .fasta.fai file of the indexed assembly (from Samtools)
3) aligned OMNI-C reads to assembly as a Bed File (BWA -> Bam To Bed)

In the script, I called the full path of the run_pipeline.py file (gotten from the  ```which run_pipeline.py``` command) because even though the Conda environment was activated, it could not just automatically find the pipeline. Then I put in the assembly, the fasta.fai file, the alignment as a .bed file, the enzyme is DNASE because of the fact it is an OMNI-C reads, and we are making scaffolds.


<details><summary>run_salsa.sh</summary>

```
#!/bin/sh
#SBATCH --partition=p_ccib_1   # which partition to run the job
#SBATCH --job-name=salsa2 # job name for listing in queue
#SBATCH --mem=192000 # memory to allocate in Mb
#SBATCH -n 1 # number of cores to use
#SBATCH -N 1 # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=7-00:00:00 # maximum run time days-hours:minutes:seconds
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE # email for the following reasons

module purge # clear existing modules
eval "$(conda shell.bash hook)" # enable slurm to see conda environments
conda activate salsa2 # load conda environment (which is hifiasm)


# actual salsa commands

python ~/.conda/envs/salsa2/bin/run_pipeline.py -a distichus_fixed_hifi_hic.fasta -l distichus_fixed_hifi_hic.fasta.fai -b output_alignment.bed -e DNASE -o scaffolds -m yes
```

</details>

## JUICER 

I ran and set up Juicer according to both the Juicer github and also Alyssa's github. We also had issues with the spack command, which is detailed in Jody's GitHub. However, I had to use Juicer1.6, so I had to run it specifically from the Aiden Lab's GitHub old module. This is the code that I had that is specific to my genome.

This is the code that was specific to setting up Juicer using Juicer 1.6.

```
mkdir juicer1.6
wget https://github.com/aidenlab/juicer/archive/refs/tags/1.6.zip
unzip *.zip

ln -s juicer/SLURM/scripts/ scripts                       #this will make a directory `scripts/` in my original `juicerdir/` directory containing only the SLURM scripts

cd scripts
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar       #downloading the most recent version of juicer tools jar
ln -s juicer_tools_1.22.01.jar juicer_tools.jar           #creating a file juicer_tools.jar that links to this version

cd ..
mkdir distichus1.6                                              #making my working directory for Anolis grahami
mkdir fastq                                               #directory where fastq HiC reads will go

```

```
screen
module load cuda/8.0
module load java/1.8.0_252


./scripts/juicer.sh -g distichus_hifi_hic -d /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_1.6 -p /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_1.6/reference/distichus_fixed_hifi_hic.chrom.sizes -y none -z /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_1.6/reference/distichus_fixed_hifi_hic.fasta -D /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_1.6 -t 20 -q p_ccib_1 -l p_ccib_1 
```

### Running 3DDNA
I downloaded 3DDNA with wget and unzipped it. I made a Conda environment called 3DDNA, which I installed all the packages needed to run the software (except Matplotlib because it was giving a fun weird error). I currently wrote a second script to call the wrapper script which looks like this but guess what it's not running. I do not think I need to load gawk or coreutils (but I did put them in the Conda environment). The only thing not in the Conda environment is Java. And also stupid Matplotlib. Oh, and Bash. I don't think I needed to load Bash. We are literally already writing everything in Bash. 

I ran into some errors running 3DDNA because I did not make it into an executable file. So this is the code to do that.
```
chmod 755 run-3ddna-pipeline.sh 
```
```
  GNU nano 2.3.1              File: run_3ddna.sh                                    
#!/bin/sh
#SBATCH --partition=p_ccib_1   # which partition to run the job
#SBATCH --job-name=3ddna # job name for listing in queue
#SBATCH --mem=192000 # memory to allocate in Mb
#SBATCH -n 32 # number of cores to use
#SBATCH -N 1 # number of nodes the cores should be on, 1 means all cores on same no$#SBATCH --time=7-00:00:00 # maximum run time days-hours:minutes:seconds
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu # email address to send status up$#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE # email for the following reasons

# load amarel modules
module purge
module load java/1.8.0_252
eval "$(conda shell.bash hook)"
conda activate 3ddna

# bash commands for analysis
./run-asm-pipeline.sh distichus_fixed_hifi_hic.fasta merged_nodups.txt
```
Somehow this became worse than the other ones when running Stats/Busco. Rip. We don't know why. We tried a new scaffolder.

### YAHS
Yahs is another scaffolder. Made new directory with Yahs. Downloaded Yahs with this command
 ```
 git clone https://github.com/c-zhou/yahs.git
 ```
 Then typed ```make``` in the folder to unpack everything.
 
 Yahs needs a .Bed file which I had to make earlier to use in Salsa2, so copied that over. Also need the assembled distichus_hifi_hic.fasta file.
 
 Ran Yahs with bash script like this:
 
 ```
   GNU nano 2.3.1                                            File: run_yahs.sh                                                                                               
#!/bin/sh
#SBATCH --partition=p_ccib_1   # which partition to run the job
#SBATCH --job-name=3ddna # job name for listing in queue
#SBATCH --mem=192000 # memory to allocate in Mb
#SBATCH -n 32 # number of cores to use
#SBATCH -N 1 # number of nodes the cores should be on, 1 means all cores on same no
#SBATCH --time=7-00:00:00 # maximum run time days-hours:minutes:seconds
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu # email address to send status up
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE # email for the following reasons

/projects/f_geneva_1/chfal/distichus/yahs/yahs/yahs distichus_fixed_hifi_hic.fasta output_alignment.bed

 ```

Had to include full file path to program because I did not add the program to my bin. Can either do that or include full file path.

I ran Busco and Stats on the YAHS output, and it was the best by far. We renamed and resorted the Yahs output using the rename script and Samtools/BWA. Then we will run the output through Juicer again. There are a lot of steps on BWA / Samtools that I just put in one ginormous script so I could have all the output files for Juicer overnight.

```
  GNU nano 2.3.1                                            File: run_bwa.sh                                                                                                
#!/bin/bash
#SBATCH --partition=p_geneva_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --account=general
#SBATCH --job-name=bwa                          # job name for listing in queue
#SBATCH --mem=192G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

echo "load any Amarel modules that script requires"
module purge # clears out any pre-existing modules
module load samtools # load any modules needed
module load bwa

echo "Bash commands for the analysis you are going to run"

echo "##################### index and align with BWA"
bwa index /projects/f_geneva_1/chfal/distichus/rename/distichus_fixed_yahs.fasta

bwa mem -t 10 -5SPM /projects/f_geneva_1/chfal/distichus/rename/distichus_fixed_yahs.fasta \
/projects/f_geneva_1/data/denovo_genomes/distichus/KU_data/raw_reads/DTG-OmniC-2_R2_001.fastq.gz \
/projects/f_geneva_1/data/denovo_genomes/distichus/KU_data/raw_reads/DTG-OmniC-2_R1_001.fastq.gz | \
samtools sort -@10 -o /projects/f_geneva_1/chfal/distichus/bwa/mapped_hifi_hic_yahs.bam -
samtools faidx /projects/f_geneva_1/chfal/distichus/rename/distichus_fixed_yahs.fasta

echo "change user group of files created"
chgrp -R g_geneva_1 /projects/f_geneva_1/chfal/distichus/bwa             # changes group of all files in listed directory


```

### Juicer on YAHS

## Set up directories (copied from Alyssa's GitHub)
You should have Juicer directory containing `scripts/`, `references/` (and optionally `restriction_sites/`), and a different working directory `AnoGra/` containing `fastq/`

[Juicer Tools jar](https://github.com/aidenlab/juicer/wiki/Download) should be installed in your `scripts/` directory

### Directory explanations
- `scripts/`
  - this folder contains SLURM scripts downloaded with a link provided below
- `references/`
  - this folder contains your reference genome and the BWA index files
- `AnoGra/`
  - this is your working directory
  - `fastq/`
    - this contains your sequence HiC reads and can remained zipped
  - `references`
    - copy over references file here as well
  - `scripts`
    - soft link scripts folder here as well
    - `ln -s ../juicer/SLURM/scripts/ scripts`

Also have to make chrom.sizes file.
```
ln -s ../juicer-1.6/SLURM/scripts scripts

cut -f1-2 distichus_fixed_yahs.fasta.fai > distichus_fixed_yahs.chrom.sizes

```
```
screen
module load cuda/8.0
module load java/1.8.0_252


./scripts/juicer.sh -g distichus_yahs1.6 -d /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_yahs1.6 -p /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_yahs1.6/reference/distichus_fixed_yahs.chrom.sizes -y none -z /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_yahs1.6/reference/distichus_fixed_yahs.fasta -D /projects/f_geneva_1/chfal/distichus/juicer1.6/distichus_yahs1.6 -t 20 -q p_geneva_1 -l p_geneva_1 
```

Juicer results can be put into JuiceBox online. We will download both inter.hic and inter_30.hic files. The inter_30.hic file indicates a MAPQ score of greater than 30. A MAPQ score shows the probability that a read is aligned correctly, with a higher score being a higher probability of the read being aligned. The MAPQ is a logarithmic scale so that means that a MAPQ score of 30 indicates a .999% or greater chance of a correct match.


### Jellyfish and Genomescope

Need to run Jellyfish to get .histo file for genomescope.

Variable is the k-mer length.

```

#!/bin/sh
#SBATCH --partition=p_ccib_1   # which partition to run the job, options are in the Amarel guide
#SBATCH --job-name=jellyfish # job name for listing in queue
#SBATCH --mem=512000 # memory to allocate in Mb
#SBATCH -n 10 # number of cores to use
#SBATCH -N 1 # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=2-00:00:00 # maximum run time days-hours:minutes:seconds, in this case 2 days
#SBATCH --mail-user=chf29@rutgers.edu # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE # email for the following reasons


# now load any Amarel modules that your script will require
module purge # clear existing modules
eval "$(conda shell.bash hook)" # enable slurm to see conda environments
conda activate jellyfish # load conda environment (which is jellyfish)

# Bash commands for the analysis you are going to run

jellyfish count -C -m ${1} -s 1000000000 -t 10 /projects/f_geneva_1/chfal/distichus/jellyfish/*.fastq -o reads${1}.jf

jellyfish histo -t 10 reads${1}.jf > reads${1}.histo
```
