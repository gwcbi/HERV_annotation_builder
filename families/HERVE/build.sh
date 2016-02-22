#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="LTR2 LTR2B LTR2C"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERVE
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERVE-int";' > ucsc/HERVE-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR2:  463
# LTR2B: 483
# LTR2C: 501

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 500 --minpct 0.1  ucsc/HERVE-int.hg19.gtf "ucsc/LTR2*.hg19.gtf" HERVE

### Create sequence files
build_annotation --get_sequence --flank 500 --minpct 0.1  ucsc/HERVE-int.hg19.gtf "ucsc/LTR2*.hg19.gtf" HERVE

# Output:
# internal            54
# oneside             29
# prototype           26
# oneside*            1
# Rejected            57