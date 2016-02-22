#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="LTR6A LTR6B"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERVS71-int
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERVS71-int";' > ucsc/HERVS71-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR6A:  564
# LTR6B:  558

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 600 --minpct 0.1  ucsc/HERVS71-int.hg19.gtf "ucsc/LTR6*.hg19.gtf" HERVS71

### Create sequence files
build_annotation --get_sequence --flank 600 --minpct 0.1  ucsc/HERVS71-int.hg19.gtf "ucsc/LTR6*.hg19.gtf" HERVS71

# Output:
# prototype           40
# oneside             9
# Rejected            13