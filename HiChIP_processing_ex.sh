# make merged bigwigs for 1D visualization
samtools merge -@ 30 \
E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
/data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1_R1_trimV1_mm10.bam \
/data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1_R2_trimV1_mm10.bam \
/data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep2_R1_trimV1_mm10.bam \
/data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep2_R2_trimV1_mm10.bam

samtools merge -@ 30 \
E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
/data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1_R1_trimV1_mm10.bam \
/data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1_R2_trimV1_mm10.bam \
/data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep2_R1_trimV1_mm10.bam \
/data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep2_R2_trimV1_mm10.bam

samtools merge -@ 30 \
E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
/data/bams/E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep1_R1_trimV1_mm10.bam \
/data/bams/E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep1_R2_trimV1_mm10.bam \
/data/bams/E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep2_R1_trimV1_mm10.bam \
/data/bams/E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep2_R2_trimV1_mm10.bam

samtools merge -@ 30 \
E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
/data/bams/E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep1_R1_trimV1_mm10.bam \
/data/bams/E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep1_R2_trimV1_mm10.bam \
/data/bams/E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep2_R1_trimV1_mm10.bam \
/data/bams/E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep2_R2_trimV1_mm10.bam

for x in *rep1-2_R1-2*bam; do samtools index $x; done

bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
  -p 20 --normalizeUsing CPM --outFileFormat bigwig \
  --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
  --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
  -b E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
  --outFileName E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw

bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
  -p 20 --normalizeUsing CPM --outFileFormat bigwig \
  --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
  --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
  -b E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
  --outFileName E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw

bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
  -p 20 --normalizeUsing CPM --outFileFormat bigwig \
  --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
  --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
  -b E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
  --outFileName E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw

bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
  -p 20 --normalizeUsing CPM --outFileFormat bigwig \
  --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
  --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
  -b E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10.bam \
  --outFileName E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw

samtools merge -@ 30 \
  E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1_R1-2_trimV1_mm10.bam \
  /data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1_R1_trimV1_mm10.bam \
  /data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1_R2_trimV1_mm10.bam

samtools merge -@ 30 \
  E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep2_R1-2_trimV1_mm10.bam \
  /data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep2_R1_trimV1_mm10.bam \
  /data/bams/E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep2_R2_trimV1_mm10.bam

samtools merge -@ 30 \
  E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1_R1-2_trimV1_mm10.bam \
  /data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1_R1_trimV1_mm10.bam \
  /data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1_R2_trimV1_mm10.bam

samtools merge -@ 30 \
  E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep2_R1-2_trimV1_mm10.bam \
  /data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep2_R1_trimV1_mm10.bam \
  /data/bams/E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep2_R2_trimV1_mm10.bam

for x in E14_D0*rep1_R1-2*bam; do samtools index $x; done
for x in E14_D4*rep2_R1-2*bam; do samtools index $x; done

bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
  -p 20 --normalizeUsing CPM --outFileFormat bigwig \
  --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
  --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
  -b E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1_R1-2_trimV1_mm10.bam \
  --outFileName E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep1_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw

  bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
    -p 20 --normalizeUsing CPM --outFileFormat bigwig \
    --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
    --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
    -b E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1_R1-2_trimV1_mm10.bam \
    --outFileName E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep1_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw

bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
  -p 20 --normalizeUsing CPM --outFileFormat bigwig \
  --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
  --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
  -b E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep2_R1-2_trimV1_mm10.bam \
  --outFileName E14_D0_WT_H3K27ac_HiChIP_AMS052022_rep2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw


bamCoverage --binSize 1 --smoothLength 0 --minMappingQuality 1 \
  -p 20 --normalizeUsing CPM --outFileFormat bigwig \
  --blackListFileName /data/reference_genomes/mm10/mm10_blackList_ENCFF547MET.bed \
  --ignoreForNormalization chrX chrM chrY --ignoreDuplicates \
  -b E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep2_R1-2_trimV1_mm10.bam \
  --outFileName E14_D0_TKO_H3K27ac_HiChIP_AMS052022_rep2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw



# call peaks THIS ACTUALLY DIDNT WORK WELL AT ALL
MACS2_ACTIVATE="source activate macs2"
MACS2="macs2"
$MACS2 callpeak -t E14_D0-4_WT-TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_all.bam  -g mm --broad
# mv NA_peaks.narrowPeak E14_D0-4_WT-TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_all_trimv1_mm10_macs2_narrow_peaks.bed
mv NA_peaks.broadPeak E14_D0-4_WT-TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_all_trimv1_mm10_macs2_broad_peaks.bed
MACS2_DEACTIVATE="source activate base"
# COMBINED ESC AND EPILC PEAKS WERE SHIT.


# CALL ESC AND EPILC PEAKS SEPARATELY AND MERGE THEM WORKED MUCH BETTER
cat E14_D0_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimv1_mm10_macs2_broad_peaks.bed \
  E14_D4_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimv1_mm10_macs2_broad_peaks.bed | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - \
  > E14_D0-4_WT-TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_all_trimv1_mm10_macs2_broad_peaks_bedtools_merge.bed
# THIS WORKED FINE AS INPUT TO HICHIPPER!!





################################################################################
#                           H3K27AC LOOP CALLING                               #
################################################################################

# Digest reference fasta genome using Arima restriction cut sites
~/bin/FUCK/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r ^GATC G^ANTC -o mm10_arima.bed mm10.fa

#run hic-pro
pwd=/scratch/hicpro_all_samples_full_depth/
nohup HiC-Pro -c config-hicpro-first-run.txt -i FASTQ/ -o hicpro_results &



# correlate replicates
# convert hic-pro matrices to h5 for correlating
# change 2500 to whatever bin size you used/want to convert
for x in *; do echo $x; hicConvertFormat --matrices $x/raw/2500/$x"_2500.matrix" \
  --outFileName $x"_2500_matrix.h5" --bedFileHicpro $x/raw/2500/$x"_2500_abs.bed" \
  --inputFormat hicpro --outputFormat h5; done

for x in *; do echo $x; hicConvertFormat --matrices $x/raw/5000/$x"_5000.matrix" \
  --outFileName $x"_5000_matrix.h5" --bedFileHicpro $x/raw/5000/$x"_5000_abs.bed" \
  --inputFormat hicpro --outputFormat h5; done




# move replicate matrix files (validPairs) to shared folders. All files in a folder are considered one sample.
# did this in FASTQ and /hicpro_results/data folders for this to work
# from Julie SEGUENI's comments and here: https://github.com/nservant/HiC-Pro/issues/121
nohup HiC-Pro -c config-hicpro-first-run.txt -i FASTQ/ -o hicpro_results -s merge_persample -s build_contact_maps -s ice_norm &


pwd
/scratch/hicpro_all_samples_full_depth/HICEXPLORER/h5_files
hicCorrelate --matrices *rep*5000_matrix.h5 --outFileNameHeatmap H3K27ac_heatmap_pearson.pdf --outFileNameScatter H3K27ac_scatterplot_pearson.pdf

# Got a really nice heatmap that shows >0.9 pearson correlation between D4 datasets and
# ~0.6-0.9 pearson correlation between D0 datasets

# run hichipper
cat D0-4_prepeak_bedtoolsmerge.yaml
peaks:
  - /scratch/hicpro_all_samples_full_depth/E14_D0-4_WT-TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_all_trimv1_mm10_macs2_broad_peaks_bedtools_merge.bed
resfrags:
  - /scratch/hicpro_all_samples_full_depth/mm10_arima.bed.gz
hicpro_output:
  - hicpro_results

# Finally, I padded the peaks by 1kb. 2.5kb could have probably been even better.
    
# Run mango to find significant loops
pwd
robin@epipax:/scratch/hicpro_all_samples_full_depth/HICHIPPER/D0-4_prepeaks_all_pad1000/DIFFLOOPS$
ln -s ../*filt.intra* .
Rscript ./mango.R

# Run diffloops to find significant loops that differ between CONDITIONS
Rscript ./diffloops.R

# Download output files to MBP14 and process in french-disco
scp robin@epipax:/scratch/hicpro_all_samples_full_depth/HICHIPPER/D0-4_prepeaks_all_pad1000/DIFFLOOPS/*csv .
for x in *mango.csv; do ./mangoToLonginteract.py $x; done
for x in *Assoc.csv; do ./diffloopToLonginteract.py $x; done
for x in *bed; do bed2bbmm10 $x; done
for x in pad1000*.bed.bb; do mv $x ${x//.bed.bb/.bb}; done
scp *bb robin@epipax:/data/track_hubs/CTCF/mm10/

function bed2bbmm10 {
  FILE=$1
  echo $FILE
  sort -k1,1 -k2,2n $FILE > tmp
  bedToBigBed tmp ~/bin/mm10.sizes tmp2 -type=bed5+13
  mv tmp2 $FILE".bb"
}
robin@epipax:/data/track_hubs$ cat group_script_bb.sh
for x in *.bb; do echo "track ${x//.bb/}"; echo "bigDataUrl $x"; echo "shortLabel ${x//.bb/}"; echo "longLabel ${x//.bb/}"; echo "parent "; echo "type bigInteract"; echo "visibility full"; echo ""; done

