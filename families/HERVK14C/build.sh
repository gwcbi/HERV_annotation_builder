#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="LTR14C"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERVK14-int
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERVK14C-int";' > ucsc/HERVK14C-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR14C: 587

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 600 --minpct 0.1  ucsc/HERVK14C-int.hg19.gtf "ucsc/LTR14*.hg19.gtf" HERVK14C

### Create sequence files
build_annotation --get_sequence --flank 600 --minpct 0.1  ucsc/HERVK14C-int.hg19.gtf "ucsc/LTR14*.hg19.gtf" HERVK14C

# Output:
# prototype           15
# oneside             10
# internal            5
# Rejected            8
