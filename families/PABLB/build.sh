#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="PABL_B"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal PABL_B
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "PABL_B-int";' > ucsc/PABL_B-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# PABL_B: 667

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 700 --minpct 0.1  ucsc/PABL_B-int.hg19.gtf "ucsc/PABL_B.hg19.gtf" PABLB

### Create sequence files
build_annotation --get_sequence --flank 700 --minpct 0.1  ucsc/PABL_B-int.hg19.gtf "ucsc/PABL_B.hg19.gtf" PABLB

# Output:
# internal            51
# prototype           7
# oneside             3
# Rejected            44