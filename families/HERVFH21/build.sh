#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="LTR21A LTR21B"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERVFH21
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERVFH21-int";' > ucsc/HERVFH21-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR21A: 505
# LTR21B: 438

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 500 --minpct 0.1  ucsc/HERVFH21-int.hg19.gtf "ucsc/LTR21*.hg19.gtf" HERVFH21

### Create sequence files
build_annotation --get_sequence --flank 500 --minpct 0.1  ucsc/HERVFH21-int.hg19.gtf "ucsc/LTR21*.hg19.gtf" HERVFH21

# Output:
# oneside             16
# prototype           16
# internal            1
# Rejected            9
