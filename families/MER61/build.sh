#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="MER61A MER61B MER61C"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERV18
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "MER61-int";' > ucsc/MER61-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# MER61A: 341
# MER61B: 427
# MER61C: 428

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 500 --minpct 0.1  ucsc/MER61-int.hg19.gtf "ucsc/MER61?.hg19.gtf" MER61

### Create sequence files
build_annotation --get_sequence --flank 500 --minpct 0.1  ucsc/MER61-int.hg19.gtf "ucsc/MER61?.hg19.gtf" MER61 <<EOF
REJECT MER61_0028
EOF

# Output:
# prototype           151
# oneside             101
# internal            89
# oneside*            1
# prototype*          1
# Rejected            91
