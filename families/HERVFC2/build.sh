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
  -e 'SELECT * FROM rmsk WHERE repName = "HERV-Fc2-int";' > ucsc/HERV-Fc2-int.hg19.txt

# LTR
models="HERV-Fc2_LTR"
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
# HERV-Fc2_LTR: 374

### First run
# build_annotation --igv_preview --flank 350 --minpct 0.1  ucsc/HERV-Fc2-int.hg19.gtf "ucsc/HERV-Fc2_LTR*.hg19.gtf" HERVFC2

### Create sequence files
build_annotation --get_sequence --flank 350 --minpct 0.1  ucsc/HERV-Fc2-int.hg19.gtf "ucsc/HERV-Fc2_LTR*.hg19.gtf" HERVFC2

# Output:
# prototype           1
# internal            1
# Rejected            0


