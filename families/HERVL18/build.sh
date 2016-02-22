#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="LTR18A LTR18B"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal HERV18
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "HERVL18-int";' > ucsc/HERVL18-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR18A: 357
# LTR18B: 614

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 600 --minpct 0.1  ucsc/HERVL18-int.hg19.gtf "ucsc/LTR18*.hg19.gtf" HERVL18

### Create sequence files
build_annotation --get_sequence --flank 600 --minpct 0.1  ucsc/HERVL18-int.hg19.gtf "ucsc/LTR18*.hg19.gtf" HERVL18

# Output:
# prototype           86
# oneside             35
# internal            6
# prototype*          1
# oneside*            1
# Rejected            32