#PBS -N Fastqc_pocean-trimmed_1
#PBS -l nodes=1:ppn=4,vmem=8gb,walltime=100:00:00
#PBS -V
#PBS -q default

#change path
cd /LUSTRE/usuario/mcamargo/metag_paocean_071221/metgenome_lib_opS0331-36/Trimmedseq_op


#load required module

module load  FastQC/0.11.2

#Run fastqc

fastqc -t24  *.fastq*
