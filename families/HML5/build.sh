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
  -e 'SELECT * FROM rmsk WHERE repName = "HERVK22-int";' > ucsc/HERVK22-int.hg19.txt

# LTR
models="LTR22 LTR22A LTR22B LTR22C"
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
# LTR22:  571 
# LTR22A: 454
# LTR22B: 542
# LTR22C: 509

### First run
# build_annotation --igv_preview --flank 600 --minpct 0.1  ucsc/HERVK22-int.hg19.gtf "ucsc/LTR22*.hg19.gtf" HML5

### Create sequence files
build_annotation --get_sequence --flank 600 --minpct 0.1  ucsc/HERVK22-int.hg19.gtf "ucsc/LTR22*.hg19.gtf" HML5 <<EOF
REJECT HML5_0106
EOF

# Output:
# prototype           86
# oneside             35
# internal            10
# prototype*          1
# Rejected            28
