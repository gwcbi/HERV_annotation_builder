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
  -e 'SELECT * FROM rmsk WHERE repName = "HERV9-int";' > ucsc/HERV9-int.hg19.txt

# LTR
models="LTR12 LTR12B LTR12C LTR12D LTR12E LTR12F LTR12_"
for model in $models; do
    echo "Query: $model"
    mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A \
        -e "SELECT * FROM rmsk WHERE repName = '$model';" > ucsc/$model.hg19.txt
done

wc -l ucsc/*

### LTR model lengths ####################################################################
# The LTR model lengths are:
# LTR12:   826
# LTR12B:  667
# LTR12C: 1577
# LTR12D: 1254
# LTR12E: 1322
# LTR12F:  519
# LTR12_:  688

### Format UCSC tables as GTF ############################################################
for f in ucsc/*.hg19.txt; do
    table2gtf < $f > ${f%.*}.gtf
done

### First run
# build_annotation --igv_preview --flank 1000 --minpct 0.1  ucsc/HERV9-int.hg19.gtf "ucsc/LTR12*.hg19.gtf" HERV9

### Create sequence files
build_annotation --get_sequence --flank 1000 --minpct 0.1  ucsc/HERV9-int.hg19.gtf "ucsc/LTR12*.hg19.gtf" HERV9

# Output:
# prototype           258
# oneside             90
# internal            11
# prototype*          1
# Rejected            95
