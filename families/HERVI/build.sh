#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="LTR10A LTR10B LTR10B1 LTR10C LTR10D LTR10E LTR10F LTR10G"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERVI
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERVI-int";' > ucsc/HERVI-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR10B:  522
# LTR10B1: 547

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 1000 --minpct 0.1  ucsc/HERVI-int.hg19.gtf "ucsc/LTR10*.hg19.gtf" HERVI

### Create sequence files
build_annotation --get_sequence --flank 1000 --minpct 0.1  ucsc/HERVI-int.hg19.gtf "ucsc/LTR10*.hg19.gtf" HERVI

# Output:
# prototype           33
# oneside             25
# internal            3
# oneside*            1
# Rejected            11