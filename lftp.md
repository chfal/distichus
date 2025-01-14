
# How to put files on NCBI SRA using ftp

```
conda create -n ftp

conda activate ftp

conda install conda-forge::lftp
```

# Setup ftp file on NCBI website

First, go to the Amarel directory you want to copy over. You should put all the files you want to upload in one folder.

I am uploading four files:
2 SMRT Cells of PacBio data
Forward and reverse reads for Illumina OmniC sequencing

These files are called:

```
DTG-OmniC-2_R1_001.fastq.gz
DTG-OmniC-2_R2_001.fastq.gz
m64086e_220711_172529.hifi_reads.fastq.gz
m64086e_220726_054912.hifi_reads.fastq.gz
```

I made a directory with all four files in it called ```upload_ncbi``` and navigated to that directory.

```

/projects/f_geneva_1/data/denovo_genomes/distichus/upload_NCBI

```

Now it is time to activate lftp and ftp into the NCBI website.

```
conda activate lftp

lftp ftp://ftp-private.ncbi.nlm.nih.gov # type this all at once to get into ftp

login subftp iodrac9octOckEpp # login, username, and secret password from NCBI which may change

cd uploads/chf29_scarletmail.rutgers.edu_O6S2m9rE # move to this remote directory which was supplied from NCBI

mkdir AnoDis # make upload directory

cd AnoDis # go into the folder you just made

mput * # puts all the files that are in the directory you were in (in my case, upload_ncbi) and puts those files into the NCBI website
```
