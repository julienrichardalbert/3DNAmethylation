#!/bin/bash

#SBATCH --ntasks=30                      # number of MPI processes
#SBATCH --mem-per-cpu=1G               # memory; default unit is megabytes
#SBATCH --mail-user=jrichardalbert@gmail.com
#SBATCH --mail-type=ALL


########################
# TAKES LONG. RUN ONCE #
########################
THREADS=10

REF_POINTS="ctcf_peaks_3kb.bed"
COOL1="cool/ctcf_d4_TKO_merge_1000_matrix.cool"
COOL2="cool/ctcf_d4_WT_merge_1000_matrix.cool"

chicViewpointBackgroundModel \
    -m $COOL1 $COOL2 \
    --fixateRange 500000 \
    --averageContactBin 4 \
    -t $THREADS \
    -rp $REF_POINTS \
    -o background_model_d0_vs_d4_ctcf_peaks3kb_1kbb_500kb.txt


