#! /bin/bash

### Set environment variables ############################################################
which table2gtf
if [[ "$?" != 0 ]]; then
    echo "Setting PATH"
    export PATH=$(dirname $(dirname $PWD))/bin:$PATH
    export PYTHONPATH=$(dirname $(dirname $PWD))/python:$PYTHONPATH
fi

### Model specific variables #############################################################

models="MER50 MER50B MER50C"

### Download RepeatMasker tracks from UCSC ################################################
mkdir -p ucsc

# Internal MER50-int
mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
  -e 'SELECT * FROM rmsk WHERE repName = "MER50-int";' > ucsc/MER50-int.hg19.txt

# LTR
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*.txt

### LTR model lengths ####################################################################
# The LTR model lengths are:
# MER50:  734
# MER50B: 697
# MER50C: 782

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

# Need to change file name for MER50.hg19.gtf so that globbing works
mv ucsc/MER50.hg19.gtf ucsc/MER50a.hg19.gtf

### First run
# build_annotation --flank 700 --minpct 0.1  ucsc/MER50-int.hg19.gtf "ucsc/MER50?.hg19.gtf" HERVFRD

### Create sequence files
build_annotation --get_sequence --flank 700 --minpct 0.1  ucsc/MER50-int.hg19.gtf "ucsc/MER50?.hg19.gtf" HERVFRD


# Output:
# oneside             52
# prototype           48
# internal            15
# oneside*            2
# prototype*          2
# Rejected            105
