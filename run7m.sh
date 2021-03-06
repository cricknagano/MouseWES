#!/bin/bash
#SBATCH --time=3-00:00:0
#SBATCH -e /camp/lab/swantonc/working/naganoa/EgfrMouse/output/DN21018/logrun7m.txt
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16

module purge
module load R/3.6.2-foss-2019b
export R_LIBS_USER="/camp/lab/swantonc/working/naganoa/Application/Rlibrary/3.6"

srun -c 8 -n 1 -J rcode7m time Rscript /camp/lab/swantonc/working/naganoa/EgfrMouse/code/7m_SINGLE_ANNOTATION_TRACERX_SCRIPT.R
