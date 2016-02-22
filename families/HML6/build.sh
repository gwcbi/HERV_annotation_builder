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
  -e 'SELECT * FROM rmsk WHERE repName = "HERVK3-int";' > ucsc/HERVK3-int.hg19.txt

# LTR
models="LTR3 LTR3A LTR3B LTR3B_"
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
# LTR3:   432
# LTR3A:  484
# LTR3B:  500
# LTR3B_: 500

### First run
# build_annotation --igv_preview --flank 500 --minpct 0.1  ucsc/HERVK3-int.hg19.gtf "ucsc/LTR3*.hg19.gtf" HML6

### Create sequence files
build_annotation --get_sequence --flank 500 --minpct 0.1  ucsc/HERVK3-int.hg19.gtf "ucsc/LTR3*.hg19.gtf" HML6 <<EOF
TANDEM HML6_0099+HML6_0100
EOF

# Output:
# prototype           50
# oneside             39
# internal            2
# tandem              1
# internal*           1
# Rejected            26
