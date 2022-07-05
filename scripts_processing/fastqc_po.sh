#PBS -N Fastqc_pocean-trimmed_1
#PBS -l nodes=1:ppn=4,vmem=8gb,walltime=100:00:00
#PBS -V
#PBS -q default

#change path to the folder with the raw_reads
cd /LUSTRE/usuario/username/


#load required module

module load  FastQC/0.11.2

#Run fastqc

fastqc -t24  *.fastq*
