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
  -e 'SELECT * FROM rmsk WHERE repName = "HERVH-int";' > ucsc/HERVH-int.hg19.txt

# LTR
models="LTR7 LTR7A LTR7B LTR7C LTR7Y"
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
# LTR7:  450
# LTR7A: 450
# LTR7B: 464
# LTR7C: 471
# LTR7Y: 472

### First run
# build_annotation --igv_preview --flank 500 --minpct 0.1  ucsc/HERVH-int.hg19.gtf "ucsc/LTR7*.hg19.gtf" HERVH

### Create sequence files
build_annotation --get_sequence --flank 500 --minpct 0.1  ucsc/HERVH-int.hg19.gtf "ucsc/LTR7*.hg19.gtf" HERVH <<EOF
IGNORE
TANDEM HERVH_0089+HERVH_0090
TANDEM HERVH_0829+HERVH_0830
IGNORE
IGNORE
REJECT HERVH_0175
TANDEM HERVH_0223+HERVH_0224
TANDEM HERVH_0305+HERVH_0306
IGNORE
TANDEM HERVH_0440+HERVH_0441
TANDEM HERVH_0532+HERVH_0533
TANDEM HERVH_1194+HERVH_1195
TANDEM HERVH_1181+HERVH_1182
REJECT HERVH_1289
EOF

# Output:
# prototype           967
# oneside             199
# tandem              8
# internal            8
# oneside*            2
# prototype*          2
# Rejected            105
