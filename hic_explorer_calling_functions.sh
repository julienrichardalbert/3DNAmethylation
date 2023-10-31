#!/bin/bash

function plotdigital {
    # set up directory
    mkdir /scratch/high_res_hichip/explorer/plots/$OUTPUT_FOLDER
    cd /scratch/high_res_hichip/explorer/plots/$OUTPUT_FOLDER
    ln -s $COOL1 ./
    ln -s $COOL2 ./
    ln -s $REF_POINTS ./
    ln -s $BACKGROUND ./
    grep -e "$SINGLE_REFPOINT" $REF_POINTS > reftmpe.bed

    # takes a long time. Fuck me.
    chicViewpoint \
        -m $COOL1 $COOL2 \
        --fixateRange 5000000 \
        -t $THREADS \
        -rp reftmpe.bed \
        -bmf $BACKGROUND \
        --averageContactBin 4 \
        --range 500000 500000

    # produces significant.hdf5 target.hdf5 and, sometimes, errorLog.txt
    chicSignificantInteractions \
        --interactionFile chic_files.hdf5 \
        -bmf $BACKGROUND \
        --range 500000 500000 \
        --pValue 0.2 \
        --xFoldBackground 2 \
        --outFileNameSignificant significant.hdf5 \
        --outFileNameTarget target.hdf5 \
        --combinationMode dual

    # for some reason this does not save as tmp.txt... ugh.
    chicExportData --file target.hdf5 -t $THREADS --outputMode all -o tmp.tar.gz
    tar -xOf tmp.tar.gz > out

    # creates aggregate.hdf5
    chicAggregateStatistic \
        --interactionFile chic_files.hdf5 \
        --targetFile out \
        --outFileName aggregate.hdf5 \
        -t $THREADS

    # creates differential.hdf5. Only significant interactions are tested, so be lenient above.
    chicDifferentialTest \
        --aggregatedFile aggregate.hdf5 \
        --alpha 0.25 \
        --statisticTest chi2 \
        --outFileName differential.hdf5 \
        -t $THREADS

    # start keeping files now
    DATE=$(date '+%y-%m-%d-%H-%M') 

    # I need to extract the rejected null hypothesis from the differential test and calculate the logFC in order to colour the significant areas!
    chicExportData --file differential.hdf5 -t $THREADS --outputMode all -o tmp2.tar.gz
    tar -xf tmp2.tar.gz
    SIG_DIFFS=$(ls -Art *rejected_differential.txt | tail -n 1) # it's dirty but it'll do for now
    cat $SIG_DIFFS | \
        awk '{print $0, log($7/$9)/log(2)}' | \
        sed 's/-nan/log2FC(TKO\/WT)/g' > $SINGLE_REFPOINT"_digital4C_sigDiffs_"$DATE".txt"
    # rm tmp.tar.gz

    # plot all significant loops and differential loops
    chicPlotViewpoint \
        --interactionFile chic_files.hdf5 \
        --combinationMode dual \
        --range $UPSTREAM $DOWNSTREAM \
        --pValue \
        --differentialTestResult differential.hdf5 \
        --backgroundModelFile $BACKGROUND \
        --outputFormat pdf \
        --pValueSignificanceLevels 0.1 \
        --plotSignificantInteractions \
        --significantInteractions significant.hdf5 \
        --outFileName 4cplots.tar.gz
    tar -xOf 4cplots.tar.gz > $SINGLE_REFPOINT"_digital4C_"$DATE".pdf" # This only works if there is one file in the targets.tar.gz, I think
    rm 4cplots.tar.gz
    # plot to show which interactions are significant (ignore differential loops)
    chicPlotViewpoint \
        --interactionFile chic_files.hdf5 \
        --combinationMode dual \
        --range $UPSTREAM $DOWNSTREAM \
        --pValue \
        --backgroundModelFile $BACKGROUND \
        --outputFormat pdf \
        --pValueSignificanceLevels 0.1 \
        --plotSignificantInteractions \
        --significantInteractions significant.hdf5 \
        --outFileName 4cplots.tar.gz
    tar -xOf 4cplots.tar.gz > $SINGLE_REFPOINT"_digital4C_onlySig_$DATE.pdf" # This only works if there is one file in the targets.tar.gz, I think
    rm 4cplots.tar.gz
    # plot to show only which loops show a difference between samples (maybe keep this one)
    chicPlotViewpoint \
        --interactionFile chic_files.hdf5 \
        --differentialTestResult differential.hdf5 \
        --combinationMode dual \
        --range $UPSTREAM $DOWNSTREAM \
        --pValue \
        --backgroundModelFile $BACKGROUND \
        --outputFormat pdf \
        --pValueSignificanceLevels 0.1 \
        --outFileName 4cplots.tar.gz
    tar -xOf 4cplots.tar.gz > $SINGLE_REFPOINT"_digital4C_onlyDiff_$DATE.pdf" # This only works if there is one file in the targets.tar.gz, I think
    rm 4cplots.tar.gz

    # I need to include the coordinates in the file names...
    #cat $REF_POINTS | awk '{OFS=FS="\t"} { print $1":"$3-3000-10000"-"$3+75000+3000}'
    UCSC_COORDS=$(cat reftmpe.bed | awk -v up="$UPSTREAM" -v down="$DOWNSTREAM" -v pad="$REF_POINT_PADDING" '{OFS=FS="\t"} { print $1":"$3-pad-up"-"$3+down+pad}')
    echo $UCSC_COORDS > $SINGLE_REFPOINT"_UCSC_coordinates_"$DATE".txt"

    COOL1_NAME=$(basename $COOL1)
    COOL2_NAME=$(basename $COOL2)
    hicPlotMatrix \
        --log1p \
        -m $COOL1 \
        -o ${COOL1_NAME//.cool/log1p}"_"$SINGLE_REFPOINT"_"$UCSC_COORDS"_"$DATE".pdf" \
        --region $UCSC_COORDS

    hicPlotMatrix \
        --log1p \
        -m $COOL2 \
        -o ${COOL2_NAME//.cool/log1p}"_"$SINGLE_REFPOINT"_"$UCSC_COORDS"_"$DATE".pdf" \
        --region $UCSC_COORDS

    rm tmp.tar.gz tmp2.tar.gz
    rm *"$SINGLE_REFPOINT"*differential.txt
    rm out reftmpe.bed
    # rm significant.hdf5 target.hdf5 aggregate.hdf5 differential.hdf5
    cd -
}


function plotcool {
    COOL1_NAME=$(basename $COOL1)
    COOL2_NAME=$(basename $COOL2)
    hicPlotMatrix \
        --log1p \
        -m $COOL1 \
        -o ${COOL1_NAME//.cool/log1p}"_"$SINGLE_REFPOINT"_"$UCSC_COORDS"_"$DATE".pdf" \
        --region $UCSC_COORDS \
        --bigwig /scratch/high_res_hichip/explorer/mm10_TSS1kb.bw \
        --vMinBigwig 0 \
        --vMaxBigwig 1 \
        --increaseFigureHeight 1.3
    
    hicPlotMatrix \
        --log1p \
        -m $COOL2 \
        -o ${COOL2_NAME//.cool/log1p}"_"$SINGLE_REFPOINT"_"$UCSC_COORDS"_"$DATE".pdf" \
        --region $UCSC_COORDS \
        --bigwig /scratch/high_res_hichip/explorer/mm10_TSS1kb.bw \
        --vMinBigwig 0 \
        --vMaxBigwig 1 \
        --increaseFigureHeight 1.3
    
}


function plotcoollog2 {
    COOL1_NAME=$(basename $COOL1)
    hicPlotMatrix \
        -m $COOL1 \
        -o ${COOL1_NAME//.cool/}"_"$SINGLE_REFPOINT"_"$UCSC_COORDS"_"$DATE".pdf" \
        --region $UCSC_COORDS \
        --colorMap PuOr_r \
        --vMin -4 \
        --vMax 4 \
        --bigwig /scratch/high_res_hichip/explorer/mm10_TSS1kb.bw \
        --vMinBigwig 0 \
        --vMaxBigwig 1 \
        --increaseFigureHeight 1.3 \
        --dpi 1000

}
