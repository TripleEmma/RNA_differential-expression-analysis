#! /bin/bash

#SBATCH --job-name=salmon_second_batch
#SBATCH --output=/work/jiangc/grasshopper/salmon/%x-%A-%a.log
#SBATCH --mail-user=chongyi.jiang@uni-jena.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=10

# salmon
module load GCC/10.2.0 OpenMPI/4.0.5 Salmon/1.4.0
cd /data/idiv_schielzeth/Jiangchongyi/grasshopper/Trinity_full/trinity_out_dir/ || exit
salmon index -t Trinity.fasta -i Trinity.fasta.salmon.idx
cd /data/idiv_schielzeth/Jiangchongyi/grasshopper/second_batch/scripts/ || exit
fq1=$(awk "NR==$SLURM_ARRAY_TASK_ID" 02_calculate-expression_array.list)
sample="$(basename "$fq1" _1.fq)"
dir_path=$(dirname "$fq1")
fq2="$dir_path/${sample}_2.fq"
output_dir="/data/idiv_schielzeth/Jiangchongyi/grasshopper/second_batch/salmon_trinity_full/"
sample_quant="$output_dir/${sample}_quant"
ref="/data/idiv_schielzeth/Jiangchongyi/grasshopper/Trinity_full/trinity_out_dir/Trinity.fasta.salmon.idx"
salmon quant -i $ref -l A -1 $fq1 -2 $fq2 -p 10 --gcBias -o $sample_quant
# --gcBias: Enable salmon to learn and correct for fragment-level GC biases in the input data. Specifically, this model will attempt to correct for biases in how likely a sequence is to be observed based on its internal GC content. This is recommended by DEseq2.

# RSEM v1.3.1
fq1=$(awk "NR==$SLURM_ARRAY_TASK_ID" bowtie2_RSEM.array.list)
sample="$(basename "$fq1" _1.fq)"
dir_path=$(dirname "$fq1")
fq2="$dir_path/${sample}_2.fq"
reference="/data/idiv_schielzeth/Jiangchongyi/grasshopper/Trinity_full/trinity_out_dir/Trinity"
cd /data/idiv_schielzeth/Jiangchongyi/grasshopper/Trinity_full/RSEM/expression/ || exit
/home/jiangc/RSEM/usr/local/bin/rsem-calculate-expression --bowtie2 -p 12 --calc-ci --calc-pme --paired-end $fq1 $fq2 $reference $sample
# --bowtie2: Use Bowtie 2 instead of Bowtie to align reads. Since currently RSEM does not handle indel, local and discordant alignments, the Bowtie2 parameters are set in a way to avoid those alignments. In particular, we use options '--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1' by default. The last parameter of '--score-min', '-0.1', is the negative of maximum mismatch rate. This rate can be set by option '--bowtie2-mismatch-rate'. If reads are paired-end, we additionally use options '--no-mixed' and '--no-discordant'.

