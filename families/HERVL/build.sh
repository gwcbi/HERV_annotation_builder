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
  -e 'SELECT * FROM rmsk WHERE repName = "HERVL-int";' > ucsc/HERVL-int.hg19.txt

# LTR
models="MLT2A1 MLT2A2"
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
# MLT2A1: 444
# MLT2A2: 549

### First run
# build_annotation --igv_preview --flank 500 --minpct 0.1  ucsc/HERVL-int.hg19.gtf "ucsc/MLT2A*.hg19.gtf" HERVL

### Create sequence files
build_annotation --get_sequence --flank 500 --minpct 0.1  ucsc/HERVL-int.hg19.gtf "ucsc/MLT2A*.hg19.gtf" HERVL <<EOF
IGNORE
REJECT HERVL_0067
REJECT HERVL_0122
REJECT HERVL_0127
REJECT HERVL_0835
REJECT HERVL_0812
REJECT HERVL_0955,HERVL_0227
REJECT HERVL_0291
REJECT HERVL_1035
REJECT HERVL_1165
REJECT HERVL_0597,HERVL_1285
EOF

# Output:
# prototype           483
# oneside             405
# internal            155
# oneside*            5
# prototype*          4
# internal*           1
# Rejected            331