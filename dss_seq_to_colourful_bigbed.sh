#!/bin/bash
# JRA 2021

# Check the number of command line arguments
if [ $# -ne 2 ]; then
	script_name=$(basename $0)
	echo "Usage: $script_name DSS_seq_output.csv reference.sizes"
	exit 1
fi

INPUT=$1
SIZES=$2

headRest 1 $INPUT > tmp # remove header line, will print a new one
sed 's/"//g' tmp > tmpe # remove quotes added by R
sed $'s/,/\t/g' tmpe > tmpee # convert comma separators to tab separators
sed 's/\.[0-9]*//g' tmpee > tmpeee # remove decimal points  (fuck em)

#head tmpee
# chr     start   end     length  nCG     meanMethy1      meanMethy2      diff.Methy      areaStat
#63646   chr14   92112898        92119726  6829  123     0.688953066452432 0.0643410972397561    0.624611969212676 727.510005756673
#92902   chr18   4326686 4327986 1301    110     0.572669886173698 0.0486574925831731    0.524012393590525 571.607445868239

# BED9 format
#chr7  127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
#chr7  127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0

# get areaStat and convert to RGB.
# negative areaStat corresponds to PATERNALLY methylated DMRs.
# this only applies to how I ran DSS-seq (i.e. mat allele was "data1")

#255,0,0

echo "#chrom chromStart  chromEnd  name  score strand  thickStart  thickEnd  itemRgb" > head

awk '{OFS="\t";FS="\t"}{
  if ($10>255)
    print $2, $3, $4, ".", "0", ".", $3, $4, "255,0,0";
  else if ($10>0)
    print $2, $3, $4, ".", "0", ".", $3, $4, $10",0,0";
  else if ($10<-255)
    print $2, $3, $4, ".", "0", ".", $3, $4, "0,0,255";
  else if ($10<0)
    print $2, $3, $4, ".", "0", ".", $3, $4, "0,0,"255+$10;
}' tmpeee | sort -k1,1 -k2,2n | cat head - > tmpeeee

#head tmpeeee
#tail tmpeeee

bedToBigBed tmpeeee $SIZES ${INPUT//.csv/.bb}
