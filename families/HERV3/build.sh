#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="LTR4 LTR61"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERV3-int
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERV3-int";' > ucsc/HERV3-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR4:  596
# LTR61: 559

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 600 --minpct 0.1  ucsc/HERV3-int.hg19.gtf "ucsc/LTR*.hg19.gtf" HERV3

### Create sequence files
build_annotation --get_sequence --flank 600 --minpct 0.1  ucsc/HERV3-int.hg19.gtf "ucsc/LTR*.hg19.gtf" HERV3 <<EOF
REJECT HERV3_0063
EOF

# Output:
# internal            34
# prototype           30
# oneside             15
# Rejected            17
