module purge
module load R/3.6.2-foss-2019b
#export R_LIBS_USER="/camp/lab/swantonc/working/naganoa/Application/Rlibrary/3.6.2-foss-2019b"
export R_LIBS_USER="/camp/lab/swantonc/working/naganoa/Application/Rlibrary/3.6"


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/camp/apps/eb/software/GCCcore/5.4.0/lib64
module purge
module load foss/2016b
module load GLib/2.47.5-foss-2016b
module load GCCcore/5.4.0
