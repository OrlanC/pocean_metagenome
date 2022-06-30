# Metagenome analysis of ocean samples

Bioinformatic contribution for analysis of Pacific Ocean Metagenome.  
This reposiory is for metagenomic analysis of pacific ocean data.  

## Data

Metagenomic sequencing reads from oceans were taken from:

Biller, S., Berube, P., Dooley, K., et al. Marine microbial  
metagenomes sampled across space and time. Sci Data 5, 180176(2018).  
https://doi.org/10.1038/sdata.2018.176  

Marine sample collected by GEOTRACES crusie ship in the Pacific Ocean (32.49567 S 164.99233 W)  

BioProject:  
PRJNA385854  

Metadata:  
Samples/Deep/ Zone/ Material/ Collection_Date  
SRR5788415 /100m/  Epipelagic/    Water/     2011-06-13 T22:40:00  
SRR5788416  /75m/   Epipelagic/    Water/     2011-06-13 T22:40:00  
SRR5788417 /50m/   Epipelagic/    Water/     2011-06-13 T22:40:00  
SRR5788420  /5100m/ Mesopelagic/   Water/     2011-06-14 T04:11:00  
SRR5788421  /1023m/ Abyssopelagic/ Water/     2011-06-14 T04:11:00  
SRR5788422  /204m/  Abyssopelagic/ Water/     2011-06-13 T22:40:00  

## Downloading the data

The accession numbers to download were chosen according to .... FIXME .... and 
are listed in `SRA_Acc_List.txt`. 

To download all the accessions inside the server MAZORKA first you need to *obtain the paths
 of the listed accessions* with the script `make_paths.sh` and run it in the server
 with the command :

:corn:
~~~
qsub make_paths.sh
~~~
{: .language-bash}
This step should be very quick.

When you have your `path_file.txt` file you are ready to *download* using the script
`download_ocean_sra.sh` and run it with:

:corn:
~~~
qsub download_ocean_sra.sh
~~~
{: .language-bash}
This step will take several hours. It is recommended that you leave it overnight.

Unfortunately, MAZORKA is not able to extract the downloaded files to obtain the `.fastq.gz` files.
So we need to * move the downloaded files to the local computer* and from the local computer to
 the server BETTERLAB. For this, run:

:computer:
~~~
scp czirion@mazorka.langebio.cinvestav.mx:/LUSTRE/usuario/czirion/SRR* .
~~~
{: .language-bash}
This step may take a few hours. Make sure you have a stable internet conection.

When this is done now *move the files to BETTERLAB*:

:computer:
~~~
scp SRR* betterlab@132.248.196.38:/home/betterlab/hackatonMetagenomica/ocean_data/raw_data/
~~~
{: .language-bash}
This step may take a few hours. Make sure you have a stable internet connection.

## Rawreads processinig
This guide is for processing raw-reads of shuntgun metagenomic libraries of Pacific Ocean.
Before starting is important installing FastQC and Trimmomatic for quality control analysis and read filtering in your computer (in mazorka these programs already installed):

1. [FastQC](https://github.com/s-andrews/FastQC/blob/master/INSTALL.txt)
2. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)


### Flow chart for processing raw-reads in each library
:corn:
It's important to use a high-performance computing cluster.

1. Quality control analysis of raw-reads to know the quality of each library using fastqc:
~~~
qsub fastqc_po.sh
~~~
{: .language-bash}
2. Filter of reads and clipping Nextera Transposase adapters (NexTranspSeq-PE.fa) using trimmomatic:
~~~
qsub trimming_pocean.sh
~~~
{: .language-bash}

3. Evaluation of filtered reads and removal adapters by trimmomatic in each library to know number of filtered reads: 
~~~
qsub fastqc_po.sh
~~~
{: .language-bash} 


## Overrepresentation analysis
[Here](https://orlanc.github.io/pocean_metagenome/Overrepresentation_analysis/overrepresentation_analysis.html) you
can see the R script made for the overrepresentation analysis. 
