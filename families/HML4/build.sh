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
  -e 'SELECT * FROM rmsk WHERE repName = "HERVK13-int";' > ucsc/HERVK13-int.hg19.txt

# LTR
models="LTR13 LTR13A"
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
# LTR13:  1007
# LTR13A: 966

### First run
# build_annotation --igv_preview --flank 600 --minpct 0.1  ucsc/HERVK13-int.hg19.gtf "ucsc/LTR13*.hg19.gtf" HML4

### Create sequence files
build_annotation --get_sequence --flank 600 --minpct 0.1  ucsc/HERVK13-int.hg19.gtf "ucsc/LTR13*.hg19.gtf" HML4 <<EOF
TANDEM HML4_0020+HML4_0021+HML4_0022+HML4_0023
EOF

# Output:
# oneside             7
# prototype           7
# tandem              1
# Rejected            7
