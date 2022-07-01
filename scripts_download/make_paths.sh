#PBS -N sraToolkit
#PBS -q default
#PBS -l nodes=1:ppn=16,mem=40g,vmem=40g,walltime=100:00:00
#PBS -e /LUSTRE/usuario/username/srapath.error
#PBS -o /LUSTRE/usuario/username/srapath.output
#PBS -V


module load sratoolkit/2.9.6

cd /LUSTRE/usuario/username/

cat SRA_Acc_List.txt | while read line; do srapath $line >> path_file.txt; done
