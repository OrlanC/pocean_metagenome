# Metagenome analysis of ocean samples

Bioinformatic contribution for analysis of Pacific Ocean Metagenome.  
This reposiory is for metagenomic analysis of pacific ocean data.  

To use this tutorial:
- :computer: Command run in the local computer
- :corn: Command run in Mazorka server
- :microscope: Command run in Betterlab server
- :hourglass: Command that takes a lot of time

- Make sure you adapt the scripts and commands to **the filesystem of your computer and/or your username**. 
For example, some scripts run in Mazorka have`username` but when you use them you should put your username.

## Data

Metagenomic sequencing reads from oceans were taken from:

Biller, S., Berube, P., Dooley, K., et al. Marine microbial  
metagenomes sampled across space and time. Sci Data 5, 180176(2018).  
https://doi.org/10.1038/sdata.2018.176  

Marine sample collected by GEOTRACES crusie ship in the Pacific Ocean (32.49567 S 164.99233 W)  

BioProject:  
PRJNA385854  

Metadata  

|Samples    |Deep   | Zone       |Material | Collection date     | Latitud | Longitude |  
|:---------:|:-----:|:---------:|:-------:|:-------------------:|:--------:|:---------:|  
|SRR5788417 |50m    |Epipelagic  | Water   | 2011-06-13 T22:40:00| -32.49567|-164.99233|
|SRR5788416 |75m    |Epipelagic  | Water   | 2011-06-13 T22:40:00| -32.49567|-164.99233|
|SRR5788415 |100m   |Epipelagic  | Water   | 2011-06-13 T22:40:00| -32.49567|-164.99233|
|SRR5788422 |204m   |Mesopelagic | Water   | 2011-06-13 T22:40:00| -32.49567|-164.99233|
|SRR5788421 |1023m  |Abyssopelagic| Water  | 2011-06-13 T22:40:00| -32.49567|-164.99233|
|SRR5788420 |5100m  |Abyssopelagic| Water  | 2011-06-13 T22:40:00| -32.49567|-164.99233|

## Downloading the data

The accession numbers to download were chosen according to .... FIXME .... and 
are listed in `scripts_download/SRA_Acc_List.txt`. 

For this sections the used software is:  
1. [sratoolkit](https://github.com/ncbi/sra-tools) version 2.9.6  
2. [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) version 2.11.3  

To download all the accessions inside the server Mazorka first you need to **obtain the paths
 of the listed accessions** with the script `scripts_download/make_paths.sh` and run it in the server
 with the command :

:corn:
~~~
qsub make_paths.sh
~~~
{: .language-bash}  
This step should be very quick.

When you have your `path_file.txt` file you are ready to **download the SRA files** using the script
`scripts_download/download_ocean_sra.sh` and run it with:

:corn: :hourglass:
~~~
qsub download_ocean_sra.sh
~~~
{: .language-bash}  

This step will take several hours. It is recommended that you leave it overnight.

Unfortunately, Mazorka is not able to convert the downloaded files into the `.fastq.gz` files.
So we need to **move the downloaded files to the local computer** and from the local computer to
 the server Betterlab. For this, run:

:computer: :hourglass:
~~~
scp username@mazorka.langebio.cinvestav.mx:/LUSTRE/usuario/username/SRR* .
~~~
{: .language-bash}
This step may take a few hours. Make sure you have a stable internet conection.

When this is done now **move the files to Betterlab**:

:computer: :hourglass:
~~~
scp SRR* betterlab@132.248.196.38:/home/betterlab/hackatonMetagenomica/ocean_data/raw_data/
~~~
{: .language-bash}  

This step may take a few hours. Make sure you have a stable internet connection.

Enter `raw_data/`. Now you can **convert the SRA files into `.fastq.gz` files** with:

:microscope:
~~~
ls | while read line; do fasterq-dump $line -S -p -e12; gzip $line_*.fastq; done
~~~
{: .language-bash}

Once you have the `.fastq.gz` yo can **transfer them again to Mazorka** to process them there.

## Raw-reads processing
This guide is for processing the raw-reads of shutgun metagenomic libraries of Pacific Ocean.
Before starting is important installing FastQC and Trimmomatic for quality control analysis and read filtering in your computer (in Mazorka these programs are already installed):

1. [FastQC](https://github.com/s-andrews/FastQC/blob/master/INSTALL.txt)
2. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)


### Flow chart for processing raw-reads in each library

It's important to use a high-performance computing cluster.

1. **Quality control analysis** of raw-reads to know the quality of each library using FastQC:
:corn:
~~~
qsub fastqc_po.sh
~~~
{: .language-bash}  

2. **Filter of reads** and clipping Nextera Transposase adapters (NexTranspSeq-PE.fa) using Trimmomatic:
:corn:
~~~
qsub trimming_pocean.sh
~~~
{: .language-bash}
  
3. **Evaluation of filtered reads** in each library to know the quality and number of filtered reads: 
:corn:
~~~
qsub fastqc_po.sh
~~~
{: .language-bash}
  
## Overrepresentation analysis
[Here](https://orlanc.github.io/pocean_metagenome/Overrepresentation_analysis/overrepresentation_analysis.html) you
can see the R script made for the overrepresentation analysis. 

## Co-ocurrence analysis

[Here](https://orlanc.github.io/pocean_metagenome/Co-ocurrence_analysis/coocurrence.html) you can see the R script made for the co-ocurrence analysis.

## NetCoMi analysis
[Here](https://orlanc.github.io/pocean_metagenome//tutorial_NetCoMi/script_redes_NetCoMi.html) you can see the R script made for the NetCoMi analysis.
