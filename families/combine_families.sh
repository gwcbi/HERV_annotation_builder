#! /bin/bash

# Clean family directories
while read fam; do
    find ./$fam -type f | grep -v build.sh | xargs rm -f
done < families.txt

# Build family annotations
while read fam; do
    echo $fam
    # Clean out any existing files
    cd $fam && . build.sh
    cd ..
done < families.txt

# Combine the annotations
mkdir -p ../combined
cat families.txt | while read fam; do cat $fam/${fam}_span.gtf; done | sortgtf > ../combined/ALLHERV_span.gtf
cat families.txt | while read fam; do cat $fam/${fam}.gtf; done | sortgtf > ../combined/ALLHERV.gtf

python <(cat <<HERE
from collections import defaultdict
import utils
families = [l.strip('\n') for l in open('families.txt','rU')]
bychrom = defaultdict(list)
for fam in families:
    tablefile = '%s/%s_table.txt' % (fam,fam)
    titer = utils.tab_line_gen(open(tablefile,'rU'))
    header = titer.next()
    for l in titer:
        l.append(fam)
        bychrom[l[3]].append(l)

sortedtable = []
for chrom in utils.CHROMNAMES:
    if chrom in bychrom:
        sortedtable.extend(sorted(bychrom[chrom],key=lambda x:int(x[4])))

with open('../combined/ALLHERV_table.txt','w') as outh:
    print >>outh, '\t'.join(header + ['family'])
    print >>outh, '\n'.join('\t'.join(l) for l in sortedtable)

HERE
)
