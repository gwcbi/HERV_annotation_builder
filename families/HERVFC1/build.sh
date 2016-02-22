#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERV-Fc1-int";' > ucsc/HERV-Fc1-int.hg19.txt

# LTR
models="HERV-Fc1_LTR1 HERV-Fc1_LTR2 HERV-Fc1_LTR3"
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

# The LTR model lengths are:
# HERV-Fc1_LTR1: 346
# HERV-Fc1_LTR2: 425
# HERV-Fc1_LTR3: 476

### First run
# build_annotation --igv_preview --flank 400 --minpct 0.1  ucsc/HERV-Fc1-int.hg19.gtf "ucsc/HERV-Fc1_LTR*.hg19.gtf" HERVFC1

### Create sequence files
build_annotation --get_sequence --flank 400 --minpct 0.1  ucsc/HERV-Fc1-int.hg19.gtf "ucsc/HERV-Fc1_LTR*.hg19.gtf" HERVFC1

# Output:
# prototype           4
# oneside             2
# Rejected            1
