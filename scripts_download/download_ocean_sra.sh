#PBS -N sraToolkit
#PBS -q default
#PBS -l nodes=1:ppn=16,mem=40g,vmem=40g,walltime=100:00:00
#PBS -e /LUSTRE/usuario/username/sradownload.error
#PBS -o /LUSTRE/usuario/username/sradownload.output
#PBS -V


module load sratoolkit/2.9.6

cd /LUSTRE/usuario/username/

cat path_file.txt | while read line; do wget $line; done
