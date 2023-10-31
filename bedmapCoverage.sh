#!/bin/bash
# JRA 2018

# Check the number of command line arguments
if [ $# -ne 3 ]; then
	script_name=$(basename $0)
	echo "Usage: $script_name coverage.bed5 coordinates.bed operation"
	echo "Operations: mean, sum, count, stdev, echo-map-range"
	exit 1
fi


COVERAGE=$1
COORDINATES=$2
OPERATION=$3

# preserve original user-input interval file
DATE=$(date '+%y-%m-%d')
cp "$COORDINATES" "$COORDINATES"_"$DATE"_backup

# keep coverage file name (column 4) for data integration
COLUMN_FILENAME=$(head -n 10 $COVERAGE | cut -f4 | tail -n 1)
COLUMN_NAME="$COLUMN_FILENAME"_"$OPERATION"
echo "$COLUMN_NAME" > tmperoo

# calculate coverage
echo "calculating coverage over $COORDINATES using $COVERAGE and $OPERATION"
bedmap --header --delim "\t" --$OPERATION "$COORDINATES" $COVERAGE >> tmperoo
sed -i sedtmp 's/NAN/NaN/g' tmperoo
#for --mean & sum  (use COUNT), there are no values "0", instead you must sed NAN -> 0 #

if [ $OPERATION == "echo-map-range" ]; then
    awk '{OFS="\t";FS="\t"}{ print $3-$2 }' tmperoo > tmpe
    headRest 1 tmpe > tmpee && rm tmpe
    echo "$COLUMN_FILENAME"_range | cat - tmpee > tmperoo && rm tmpee
fi
if [ $OPERATION == "mean" ]; then
    sed -i sedtmp 's/-0.000000/0.0/g' tmperoo
fi
if [ $OPERATION == "stdev" ]; then
    sed -i sedtmp 's/-nan/NaN/g' tmperoo
fi


# generate text histogram with counts (0-10), does not work for DNAme
sed 's/NaN/0/g' tmperoo > tmpero
echo "$OPERATION frequency"
textHistogram -skip=1 -maxBinCount=11 tmpero
rm tmpero

# add column with coverage values
paste "$COORDINATES" tmperoo > tmperoooo
mv tmperoooo "$COORDINATES"
rm tmperoo

echo ""


# Edit 16 Feb 2018: add text histogram
# Edit 28 Jun 2018: add range of overlap operation
