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
  -e 'SELECT * FROM rmsk WHERE repName = "HERVK-int";' > ucsc/HERVK-int.hg19.txt

# LTR
models="LTR5 LTR5A LTR5B LTR5_Hs"
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

### Setup GTF for compare ################################################################
mkdir -p tmp
cat ../other_sources/subramanianT*.hg19.gtf | sed 's/^/chr/' | sortgtf > tmp/subtables.gtf

# The LTR model lengths are:
# LTR5:    969
# LTR5A:   1033
# LTR5B:   1002
# LTR5_Hs: 968

### First run
# build_annotation --igv_preview --flank 1000 --minpct 0.1 ucsc/HERVK-int.hg19.gtf "ucsc/LTR5*.hg19.gtf" HML2

### Create sequence files
build_annotation --get_sequence --flank 1000 --minpct 0.1 ucsc/HERVK-int.hg19.gtf "ucsc/LTR5*.hg19.gtf" HML2 <<EOF
TANDEM HML2_0088+HML2_0089
EOF

# Output:
# prototype           51
# oneside             25
# internal            9
# prototype*          2
# tandem              1
# Rejected            40

### NOTE:
# The annotations from Subramanian et al. that are not included are:
# 10p12.1  - chr10:27182399-27183380   This is a solo LTR
# 12q13.2  - chr12:55727215-55728183   This is a solo LTR
# 12q24.11 - chr12:111007843-111009325 Too short, only contains 514bp of internal sequence
# 19p12b   - chr19:21841536-21841542   This is polymorphic, does not occur in hg19
