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
  -e 'SELECT * FROM rmsk WHERE repName = "HERVK9-int";' > ucsc/HERVK9-int.hg19.txt

# LTR
models="MER9a1 MER9a2 MER9a3 MER9B"
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
# MER9a1: 513
# MER9a2: 513
# MER9a3: 512
# MER9B:  499

### First run
# build_annotation --igv_preview --flank 500 --minpct 0.1  ucsc/HERVK9-int.hg19.gtf "ucsc/MER9*.hg19.gtf" HML3

### Create sequence files
build_annotation --get_sequence --flank 500 --minpct 0.1  ucsc/HERVK9-int.hg19.gtf "ucsc/MER9*.hg19.gtf" HML3

# Output:
# prototype           150
# oneside             82
# internal            7
# prototype*          1
# Rejected            21
