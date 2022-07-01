#PBS -N Trimmomatic_op_1
#PBS -l nodes=1;ppn=8;mem=8gb,vmem=8gb,walltime=100:00:00
#PBS -V
#PBS -q default

cd $PBS_O_WORKDIR

for i in *_1.fastq.gz
do
        base=$(basename ${i} _1.fastq.gz)
        java -jar $TRIMMOMATIC PE -phred33 ${i} ${base}_2.fastq.gz \
        Trimmedseq_op/${base}_1.trimmed.fastq.gz Trimmedseq_op/${base}_1un.trimmed.fastq.gz \
        Trimmedseq_op/${base}_2.trimmed.fastq.gz Trimmedseq_op/${base}_2un.trimmed.fastq.gz \
        ILLUMINACLIP:NexTranspSeq-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:35
done

