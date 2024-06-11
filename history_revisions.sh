USEFUL PATHS

robin@epipax:/scratch/NAR_revisions/CTCF_HICHIP/CTCF_rep1-2_hichipper_allPeaks_1000
robin@epipax:/scratch/NAR_revisions/HICPRO_H3K27ac/HICHIPPER/D0-4_prepeaks_all_pad1000
for CUT&RUN peaks:
# jra@GREENBERG14:~/Dropbox/3dname/31_jan_2023


#
# peak  calling
# for differential loops, adjacent peaks are merged within a defined range (I used 1kb)
# so I left the peaks as-is
# However, for comparing peak calls between replicates or for measuing CTCF enrichment, I should really merge adjacent peaks.
#
#


#/Users/jra/Dropbox/3dname/REVISIONS/krebs
x="E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks_1e325_MA0139.1_peaks_data_named_limma_index_NaND4TKO.bed"
x="E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks_1e325_MA0139.1_peaks_data_named_limma_index_UPD4TKO.bed"
echo $x
wc -l $x
bedtools intersect -a $x -b CTCF_motif_CpG5_mcinfo.bed -u | wc -l
bedtools intersect -a $x -b CTCF_motif_CpG7_mcinfo.bed -u | wc -l
bedtools intersect -a $x -b CTCF_motif_CpG15_mcinfo.bed -u | wc -l
bedtools intersect -a $x -b CTCF_motif_CpG5_d4_70DNAme.bed -u | wc -l
bedtools intersect -a $x -b CTCF_motif_CpG7_d4_70DNAme.bed -u | wc -l
bedtools intersect -a $x -b CTCF_motif_CpG15_d4_70DNAme.bed -u | wc -l
bedtools intersect -a $x -b CTCF_motif_CpG5_mcinfo.bed CTCF_motif_CpG7_mcinfo.bed CTCF_motif_CpG15_mcinfo.bed -u | wc -l

130107 E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks_1e325_MA0139.1_peaks_data_named_limma_index_NaND4TKO.bed
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG5_mcinfo.bed -u | wc -l
 11654
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG7_mcinfo.bed -u | wc -l
 11021
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG15_mcinfo.bed -u | wc -l
 28708
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG5_d4_70DNAme.bed -u | wc -l
  2921
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG7_d4_70DNAme.bed -u | wc -l
  5971
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG15_d4_70DNAme.bed -u | wc -l
  5837
  bedtools intersect -a $x -b CTCF_motif_CpG5_mcinfo.bed CTCF_motif_CpG7_mcinfo.bed CTCF_motif_CpG15_mcinfo.bed -u | wc -l
     40455


  2683 E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks_1e325_MA0139.1_peaks_data_named_limma_index_UPD4TKO.bed
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG5_mcinfo.bed -u | wc -l
   955
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG7_mcinfo.bed -u | wc -l
   221
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG15_mcinfo.bed -u | wc -l
   625
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG5_d4_70DNAme.bed -u | wc -l
   807
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG7_d4_70DNAme.bed -u | wc -l
   126
(base) jra@GREENBERG14:~/Dropbox/3dname/REVISIONS/krebs$ bedtools intersect -a $x -b CTCF_motif_CpG15_d4_70DNAme.bed -u | wc -l
   520
   bedtools intersect -a $x -b CTCF_motif_CpG5_mcinfo.bed CTCF_motif_CpG7_mcinfo.bed CTCF_motif_CpG15_mcinfo.bed -u | wc -l
       1483



# loops
Can loops be called outside anchors? ****NO****
# anchors
wc -l E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks.bed
  156179 E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks.bed
# mango output of single loop file
cut -f1-3 -d ',' D4_WT_FDR0.05_mango.csv | sed 's/,/\t/g' | uniq > D4_WT_FDR0.05_mango_first_regions.bed
wc -l D4_WT_FDR0.05_mango_first_regions.bed
  101321 D4_WT_FDR0.05_mango_first_regions.bed
bedtools intersect -a D4_WT_FDR0.05_mango_first_regions.bed -b E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks.bed -u | wc -l
  101321
# output of diffloops
cut -f1-3 -d ',' D4_FDR0.05_quickAssoc.csv | sed 's/,/\t/g' | uniq > D4_FDR0.05_quickAssoc_first_region.bed
wc -l D4_FDR0.05_quickAssoc_first_region.bed
  103103 D4_FDR0.05_quickAssoc_first_region.bed


  Either I combine H3K27ac and CTCF peaks or I change software...


#
# screwing around with hichip-peaks, but I don't  think it'll be useful (same output as hichipper, plus differential peak enrichment calc)
# https://github.com/ChenfuShi/HiChIP_peaks
#
# screwing around with cLoops2, but it doesn't do "anchor vs all" comparisons, so gave up.
# robin@epipax:/scratch/NAR_revisions/CTCF_HICHIP/cloops2_testing
#
# screwing around with fithichip now
# https://hichip.readthedocs.io/en/latest/loops.html
# looks promising!
# robin@epipax:/scratch/NAR_revisions/fitchip
# robin@epipax:/home/robin/FitHiChIP
#
#
#
#
# count the number of bins/loops that overlap a peak
# max will be 50% of the diffloops count
INPUT="DiffLoops_ALL.bed"
PEAKS="E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks_1e325_MA0139.1_peaks.bed"

INPUT_COUNT=$(wc -l $INPUT | wc -l | cut -f1 -d ' ')
PEAKS_COUNT=$(wc -l $PEAKS | wc -l | cut -f1 -d ' ')

cut -f1-3 $INPUT > bin1.bed
cut -f4-6 $INPUT > bin2.bed
headRest 1 bin1.bed > bin1.1.bed
headRest 1 bin2.bed > bin2.1.bed
cat bin1.1.bed bin2.1.bed | sort -k1,1 -k2,2n > allbins.bed
BINS=$(wc -l allbins.bed | cut -f1 -d ' ')
PEAK_BINS=$(bedtools intersect -a allbins.bed -b $PEAKS -u | wc -l | cut -f1 -d ' ')
echo "LOOP BINS $INPUT overlapping peaks $PEAKS"
echo "BINS: $BINS"
echo "PEAK_BINS: $PEAK_BINS"
echo "LOOP BINS: $INPUT_COUNT"
echo "PEAKS: $PEAKS_COUNT"

#e.g.
#BINS: 80
#PEAK_BINS: 75
#would equate to 5 loops have one bin that doesn't intersect a peak; the rest of the loops overlap a peak at both ends.


INPUT="E14_D0_TKO_CTCF_HiChIP_rep1_b5000_q005.interactions_FitHiC_no11.bed"

awk '{OFS=FS="\t"}{ if ($9==1) print $1, $2, $3, "peak1"}' $INPUT >> tmp
awk '{OFS=FS="\t"}{ if ($15==1) print $4, $5, $6, "peak2"}' $INPUT >> tmp
sort -k1,1 -k2,2n tmp bedtools merge -i - -c 4 -o distinct > tmp2

awk '{OFS=FS="\t"}{ print $1, $2, $3, $9, $4, $5, $6, $15, FILENAME}' E14_D0_TKO_CTCF_HiChIP_rep1_b5000_q005.interactions_FitHiC_no11.bed






# MAking heatmaps on anchors
# Start from a csv list of loops generated by hichipper
awk '{FS=",";OFS="\t"}{ if ($18<0.05 && $14<-1) print $0}' CTCF_EpiLC_WT_vs_TKO_FDR0.05_quickAssoc_NAR_revisions_1750.csv > UPD4TKO_CTCF_loops_1750.csv

awk '{OFS=FS="\t"}{print $1, $2, $3, "loop"NR-1}' UPD4TKO_CTCF_loops_1750.bed > 1750_bin1.bed
awk '{OFS=FS="\t"}{print $4, $5, $6, "loop"NR-1}' UPD4TKO_CTCF_loops_1750.bed > 1750_bin2.bed
#  calculate enrichment over anchors using VisR

awk '{OFS=FS="\t"}{print $1, $2, $3, "1750_loop"NR-1}' UPD4TKO_CTCF_loops_1750.bed > 1750_bin1.bed
awk '{OFS=FS="\t"}{print $4, $5, $6, "1750_loop"NR-1}' UPD4TKO_CTCF_loops_1750.bed > 1750_bin2.bed

join -1 4 -2 4 -t $'\t' 1750_bin1_data.bed 1750_bin2_data.bed > 1750_bin_data.bed
awk '{OFS=FS="\t"}{ if ($6/$5 > $11/$10) print $2, $3, $4, $1}' 1750_bin_data.bed >> 1750_bin_L.bed
awk '{OFS=FS="\t"}{ if ($6/$5 < $11/$10) print $7, $8, $9, $1}' 1750_bin_data.bed >> 1750_bin_L.bed
awk '{OFS=FS="\t"}{ if ($6/$5 < $11/$10) print $2, $3, $4, $1}' 1750_bin_data.bed >> 1750_bin_R.bed
awk '{OFS=FS="\t"}{ if ($6/$5 > $11/$10) print $7, $8, $9, $1}' 1750_bin_data.bed >> 1750_bin_R.bed
sort -k1,1 -k2,2n 1750_bin_L.bed > 1750_bin_L_sort.bed
sort -k1,1 -k2,2n 1750_bin_R.bed > 1750_bin_R_sort.bed

bigwigCompare --operation log2 -b1 E14_D4_TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw -b2 E14_D4_WT_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw -bs 100 -o CTCF_log2 -p max/2
bigwigCompare --operation log2 -b1 E14_D4_TKO_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw -b2 E14_D4_WT_H3K27ac_HiChIP_AMS052022_rep1-2_R1-2_trimV1_mm10_rmDup_q10_b1_s0_CPM.bw -bs 100 -o H3K27ac_log2 -p max/2
bigwigCompare --operation log2 -b1 159_mESC_DnmtTKO_PROseq_Kreibich2023_rep1-3_trimV5_mm10_kpDup_q10_b1_s0_CPM.bw -b2 159_mESC_WT_PROseq_Kreibich2023_rep1-3_trimV5_mm10_kpDup_q10_b1_s0_CPM.bw -bs 100 -o PRO_log2 -p max/2





ac_d0_WT_merge.allValidPairs
ac_d0_TKO_merge.allValidPairs
ac_d4_WT_merge.allValidPairs
ac_d4_TKO_merge.allValidPairs
ctcf_d0_WT_merge.allValidPairs
ctcf_d0_TKO_merge.allValidPairs
ctcf_d4_WT_merge.allValidPairs
ctcf_d4_TKO_merge.allValidPairs

ctcf_d0_merge.allValidPairs@
ac_d0_merge.allValidPairs@






# Process 4C samples:
c4ctus_Demult.sh -f D4_EpiLC_4C_scrambled.fastq.gz -p barcode_NAR_revisions.fa -g EpiLC_scramble
c4ctus_Demult.sh -f D4_EpiLC_4C_epited.fastq.gz -p barcode_NAR_revisions.fa -g EpiLC_edited
for x in Demultiplexing/*fastq.gz; do echo $x; c4ctus_Mapping.sh --fastq $x  --genome mm10; done
c4ctus_Multi4Cseq.sh --lib4C 4CLib/Library_mm10_DpnII_NlaIII_30bps_Rmsk_segmentsInfos.bed --primer barcode_NAR_revisions.fa --group EpiLC_scramble --genome mm10 --Ncpu 4 --window 11
c4ctus_Multi4Cseq.sh --lib4C 4CLib/Library_mm10_DpnII_NlaIII_30bps_Rmsk_segmentsInfos.bed --primer barcode_NAR_revisions.fa --group EpiLC_edited --genome mm10 --Ncpu 4 --window 11
# Deeptools
bigwigCompare --operation subtract -b1 EpiLC_edited_Zdbf2_mm10_rep43_chr1_smoothed_11FragPerWin_Norm.bw -b2 EpiLC_scramble_Zdbf2_mm10_rep39_chr1_smoothed_11FragPerWin_Norm.bw -bs 1 -o Zdbf2_delta.bw  -p max/2
bigwigCompare --operation subtract -b1 EpiLC_edited_Nrp2_mm10_rep40_chr1_smoothed_11FragPerWin_Norm.bw -b2 EpiLC_scramble_Nrp2_mm10_rep36_chr1_smoothed_11FragPerWin_Norm.bw -bs 1 -o Nrp2_delta.bw  -p max/2
bigwigCompare --operation subtract -b1 EpiLC_edited_Csf1_mm10_rep42_chr3_smoothed_11FragPerWin_Norm.bw -b2 EpiLC_scramble_Csf1_mm10_rep38_chr3_smoothed_11FragPerWin_Norm.bw -bs 1 -o Csf1_delta.bw  -p max/2
bigwigCompare --operation subtract -b1 EpiLC_edited_Mob3b_mm10_rep41_chr4_smoothed_11FragPerWin_Norm.bw -b2 EpiLC_scramble_Mob3b_mm10_rep37_chr4_smoothed_11FragPerWin_Norm.bw -bs 1 -o Mob3b_delta.bw  -p max/2
