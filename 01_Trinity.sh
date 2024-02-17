#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=50G
#SBATCH --time=30-00:00:00
#SBATCH -J Trinity_Assembly
#SBATCH --output=/work/%u/Trinity/Trinity-%x-%j.o
#SBATCH --error=/work/%u/Trinity/Trinity-%x-%j.e
#SBATCH --mail-user=chongyi.jiang@uni-jena.de
#SBATCH --mail-type=ALL

module load GCC/10.2.0 OpenMPI/4.0.5 Trinity/2.9.1-Python-3.8.6    # note the version of Trinity
Trinity --seqType fq --samples_file samples_trinity.txt --CPU 12 --no_normalize_reads --max_memory 600G --min_contig_length 200 --full_cleanup
# --full_cleanup: only retain the Trinity fasta file
# potential problem is that there too many sequences

