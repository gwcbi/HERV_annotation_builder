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
  -e 'SELECT * FROM rmsk WHERE repName = "HERV17-int";' > ucsc/HERV17-int.hg19.txt

# LTR
models="LTR17"
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
# LTR17: 780

### First run
# build_annotation --igv_preview --flank 750 --minpct 0.1  ucsc/HERV17-int.hg19.gtf "ucsc/LTR17*.hg19.gtf" HERVW

### Create sequence files
build_annotation --get_sequence --flank 750 --minpct 0.1  ucsc/HERV17-int.hg19.gtf "ucsc/LTR17*.hg19.gtf" HERVW <<EOF
IGNORE
EOF

# Output:
# prototype           102
# oneside             84
# internal            28
# oneside*            1
# Rejected            83
