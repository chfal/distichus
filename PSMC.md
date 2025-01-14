### PSMC

I ran PSMC according to Jody's GitHub. This is the script used for this. At the bottom is the code for the bootstrapping, which I ran at 100 times.

<details><summary>run_psmc.sh</summary>

```
#!/bin/bash
#SBATCH --partition=p_geneva_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=psmc                      # job name for listing in queue
#SBATCH --mem=20G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=5-12:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=FAIL,END,REQUEUE      # email for the following reasons

echo "load modules for psmc"
module purge                                    # clears out any pre-existing modules
module load java
module load samtools
module load bcftools
module load gcc/8.2/gnuplot

PSMC=/home/chf29/psmc

echo "Check depth statistics first for bam file"

INDIR=/projects/f_geneva_1/chfal/distichus/psmc
OUTDIR=/projects/f_geneva_1/chfal/distichus/psmc

BAM=AnoDis1.0_sorted
REFDIR=/projects/f_geneva_1/chfal/distichus/psmc
REF=AnoDis1.0.fasta

echo "run steps for psmc"

# bcftools mpileup -C 50 -f ${REFDIR}/${REF} ${INDIR}/${BAM}.bam | bcftools call -c | vcfutils.pl vcf2fq -d 13 -D 90 - | gzip > ${OUTDIR}/${BAM}_diploid.fq.gz

# ${PSMC}/utils/fq2psmcfa -q20 ${OUTDIR}/${BAM}_diploid.fq.gz > ${OUTDIR}/${BAM}_diploid.psmcfa

# ${PSMC}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${OUTDIR}/${BAM}_diploid.psmc ${OUTDIR}/${BAM}_diploid.psmcfa

# ${PSMC}/utils/psmc2history.pl ${OUTDIR}/${BAM}_diploid.psmc | ${PSMC}/utils/history2ms.pl > ${OUTDIR}/ven1.1_ms-cmd.sh

${PSMC}/utils/psmc_plot.pl -u 7.17e-09 -g 1 -R -p AnoDis1.0_plot ${OUTDIR}/${BAM}_diploid.psmc

echo "run bootstrap psmc"

# #${PSMC}/utils/fq2psmcfa -q20 ${OUTDIR}/${BAM}_diploid.fq.gz > ${OUTDIR}/${BAM}_diploid.psmcfa            # This is the same as above - if file exists, no need to rerun

#${PSMC}/utils/splitfa ${OUTDIR}/${BAM}_diploid.psmcfa > ${OUTDIR}/${BAM}_split.psmcfa             # This splits your file into bins to make it easier to process

#${PSMC}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${OUTDIR}/${BAM}_diploid.psmc ${OUTDIR}/${BAM}_diploid.psmcfa     # This runs as above - need to do it on your splits.psmcfa file

#seq 100 | xargs -P 30 -i ${PSMC}/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ${OUTDIR}/round-{}.psmc ${OUTDIR}/${BAM}_split.psmcfa | sh     # This is the bootstrap step for 100 runs

#cat ${OUTDIR}/${BAM}_diploid.psmc ${OUTDIR}/round-*.psmc > ${OUTDIR}/${BAM}_combined.psmc           # Combine all the runs into a single file - need this for Rplot

# ${PSMC}/utils/psmc_plot.pl -pY50000 -u 7.17e-09 -g 1 -R combined ${OUTDIR}/${BAM}_combined.psmc    # No need to run this cause we'll plot in R, but you can run manually here if need be
```

</details>


Additionally, two scripts are needed to generate plots for PMSC. The first is the plot_psmc.R script, which comes from Emily Humble, copied below.

<details><summary>plot_psmc.R</summary>

```
# Function to plot PSMC using .psmc output
# Code from https://figshare.com/articles/Plot_PSMC_results/3996156/1
# Also available at: https://datadryad.org/stash/dataset/doi:10.5061/dryad.0618v
# which is from this paper: https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12606

psmc.result<-function(file,i.iteration=25,mu=7.17e-9,s=100,g=1)
{
  X<-scan(file=file,what="",sep="\n",quiet=TRUE)
  
  # extract data for each iteration (30)
  
  START<-grep("^RD",X) # line numbers
  END<-grep("^//",X) # line numbers
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE) # \theta_0 \rho_0
  RS<-grep("^RS",X,value=TRUE) # k t_k \lambda_k \pi_k \sum_{l\not=k}A_{kl} A_{kk}
  
  write(TR,"temp.psmc.result")
  theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-theta0/4/mu/s # scale 

  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-as.numeric(2*N0*a[,3])
  Ne<-as.numeric(N0*a[,4])
  
  a$Generation<-as.numeric(2*N0*a[,3])
  a$Ne<-as.numeric(N0*a[,4])
  
  file.remove("temp.psmc.result")
  
  n.points<-length(Ne)
  
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne,mu,g)
  #plot(Ne~YearsAgo)
}
```
</details>



Additionally, this R script takes use of the plot_psmc.R function and generates the bootstrap plot.


<details><summary>plot_psmc_bootstrap.R</summary>
  
```
#### Fig 4 - Main PSMC plot ####
# Script adapted for this data set from Emily Humble's ("3.1_psmc.R" from https://github.com/elhumble/SHO_analysis_2020)
# Calls the "psmc.result" function from the "plot_psmc.R" script from https://figshare.com/articles/Plot_PSMC_results/3996156/1 

#~~ Load requried packages
library(ggplot2)
source("C:/Users/chfal/OneDrive/Desktop/psmc/plot_psmc.R") # Load plot_psmc function. File must be in current working directory, else give path to file.
library(data.table)
library(plyr)
library(tidyr)
library(scales)
options(scipen=999) # Disable scientific notation of large numbers- see "Scipen" under ?options

#~~ Specify variables for plotting
i.iteration=25 # Number of iterations to use in file
s=100 # bin size

#~~ Read in main PSMC files. 
dis_psmc_files <- paste("C:/Users/chfal/OneDrive/Desktop/psmc/psmc_100/", list.files(path = "C:/Users/chfal/OneDrive/Desktop/psmc/psmc_100/"), sep = "")


#~~ Set mu and g (for plot_psmc function):
#mu <- 7.17e-9      # Average mutation rate from Bergeron et al 2023 - Evolution of the germline mutation rate across vertebrates (Closest relative - Pogona)
#mu <- 9.45e-9      # Min mutation rate for Pogona (from above citation)
mu <- 7.17e-9      # Max mutation rate for Pogona (from above citation)
g <- 1

#~~ Run "psmc.result" from the "plot_psmc" function
dis_psmc_buf <- lapply(dis_psmc_files, psmc.result, i.iteration, mu = mu, s, g = g)

#~~ Set names in psmc_buf to sample names
dataset_names <-list.files(path = "C:/Users/chfal/OneDrive/Desktop/psmc/psmc_100/", pattern="*.psmc")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(dis_psmc_buf) <- dataset_names

#~~ Transform list into dataframe
dis_psmc_buf <- ldply(dis_psmc_buf, .id = "Sample")  

#~~ Get bootstraps
dis_boot_files <- paste("C:/Users/chfal/OneDrive/Desktop/psmc/psmc_100/", list.files(path = "C:/Users/chfal/OneDrive/Desktop/psmc/psmc_100/", pattern="*round*"), sep = "")


#~~ Run "psmc.result" from the "plot_psmc" function
dis_boot <- lapply(dis_boot_files, psmc.result, i.iteration, mu, s, g)


#~~ Set names in psmc_buf to sample names
dataset_names <-list.files(path = "C:/Users/chfal/OneDrive/Desktop/psmc/psmc_100/", pattern="*round*")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(dis_boot) <- dataset_names

#~~ Transform list into dataframe
dis_boot_df <- ldply(dis_boot, .id = "ID") %>%
  separate(ID, c("Sample", "Boot"), sep = "round-", remove = F) %>%
  mutate(Boot = gsub("-", "", Boot))

dis_boot_df$Sample= sub("","AnoDis", dis_boot_df$Sample)

#~~ Import individual run for each species 

wd<-"C:/Users/chfal/OneDrive/Desktop/psmc/"
list.files(wd)
setwd(wd)

AnoDis<-read.table("AnoDis1.0_plot.0.txt", sep="\t", header=FALSE)       # Read in the .txt file generated from the -R flag in psmc 

library(reshape)  

# I'm just adding in coloumn headings here so that I can see what I'm calling in ggplot2

AnoDis<-rename(AnoDis,c(V1="YBP", V2="Ne"))     # YBP - years before present Ne - Effective population size (x10^4) 
head(AnoDis)

Ne_dis=c(AnoDis$Ne*10e3)
AnoDis_df<- cbind(AnoDis, Ne_dis)

#~~ Define colours
cbPalette <- c("#DA9101","black")

#~~ Make plot and store as object
psmc_main<-ggplot(mapping=aes(dis_boot_df$YearsAgo, dis_boot_df$Ne)) +
  geom_line(aes(dis_boot_df$YearsAgo, dis_boot_df$Ne, group = dis_boot_df$Boot), 
            alpha = 0.2, col = cbPalette[1]) +
  geom_line(aes(AnoDis_df$YBP, AnoDis_df$Ne_dis, col = cbPalette[1]), 
              col = cbPalette[2]) +
  theme_classic() +
  scale_color_manual(values = cbPalette,
                      labels =c(1),
                      name = expression(paste(bold("Species")))) +
  theme(legend.position = c(0.9,0.9)) +
  scale_x_log10(label = comma,
                breaks = c(1000,10000,20000,50000,100000,200000,500000,1000000,2000000)) +
  annotation_logticks(sides = "b", alpha = 0.5) +
  scale_y_continuous(label = comma) +
  xlab(expression(paste("Years before present (g = 1, ", mu," = 7.17e-9",")"))) + #Average mutation rate
  #xlab(expression(paste("Years before present (g = 1, ", mu," = 9.45e-9",")"))) + # Min mutation rate 
  #xlab(expression(paste("Years before present (g = 1, ", mu," = 3.76e-9",")"))) + # Max mutation rate 
  #ylab (expression(paste("IICR (scaled in units of 4",italic(N[e]),mu,")"))) +
  ylab ("Effective population size") +
  geom_vline(xintercept=10000, linetype = "dashed",  colour = "grey") +
  annotate("text", x = 8700, y = 356000, label = "Holocene", angle = 90, size = 3) +
  geom_vline(xintercept=22000, linetype = "dashed", colour = "grey") +
  annotate("text", x = 19000, y = 352500, label = "LGM ~22 ka", angle = 90, size = 3) + 
  geom_vline(xintercept=120000, linetype = "dashed", colour = "grey") +
  annotate("text", x = 105000, y = 347500, label = "LIG ~120-140ka", angle = 90, size = 3) +
  annotate("rect", xmin = 115000, xmax = 130000, ymin = -Inf, ymax = Inf,
           alpha = .1) +
  labs(title = "100 Bootstrap PSMC Results: A. distichus")

psmc_main
#~~ Save as pdf:
#pdf(file="figs/psmc_figure.pdf", 
#    height=8.27, width=11.69)
#psmc_main
#dev.off()

#~~ Save as svg:
ggsave("C:/Users/chfal/OneDrive/Desktop/psmc/psmc_figure.svg",
        plot = psmc_main,
        width = 11.69,
        height = 8.27,
        units = "in")
```
</details>

