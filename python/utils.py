#! /bin/bash
import sys
import re
from collections import defaultdict
from subprocess import Popen,PIPE

CHROMNAMES = ['chr%d' % d for d in range(1,23)] + ['chrX','chrY']

# Chromosome sizes for hg19
CHROMSIZE = {
             'chr1':249250621,
             'chr2':243199373,
             'chr3':198022430,
             'chr4':191154276,
             'chr5':180915260,
             'chr6':171115067,
             'chr7':159138663,
             'chr8':146364022,
             'chr9':141213431,
             'chr10':135534747,
             'chr11':135006516,
             'chr12':133851895,
             'chr13':115169878,
             'chr14':107349540,
             'chr15':102531392,
             'chr16':90354753,
             'chr17':81195210,
             'chr18':78077248,
             'chr19':59128983,
             'chr20':63025520,
             'chr21':48129895,
             'chr22':51304566,
             'chrX':155270560,
             'chrY':59373566,
}

def sort_gtf(liter):
    # Organize original lines by chromosome
    alllines = defaultdict(list)  
    for l in liter:
        alllines[l[0]].append(l)
    
    for cchrom in CHROMNAMES:
        if cchrom not in alllines: continue 
        for l in sorted(alllines[cchrom],key=lambda x:int(x[3])):
            yield l

def cluster_gtf(liter):
    # Organize original lines by final column
    # bedtools cluster outputs the cluster number in the final column
    clustered = defaultdict(list)
    for l in liter:
        clustered[l[-1]].append(l)
    
    return clustered

def tab_line_gen(infile):
    return (l.strip('\n').split('\t') for l in infile if not l.startswith('#'))

def simplify_list(l):
    cur = [l[0]]
    for v in l:
        if v != cur[-1]:
            cur.append(v)
    return cur

def by_attribute(cdict, attrname):
    ret = {}
    for cname,c in cdict.iteritems():
        c.sort(key=lambda x:int(x[3])) # Sort by start position
        attrs = [dict(re.findall('(\S+)\s+"([\s\S]+?)";',l[8])) for l in c] 
        ret[cname] = simplify_list([d[attrname] for d in attrs])
    return ret

def covered_len(clust):
    ''' Calculate the covered length of region
            If multiple regions, use bedtools to merge overlaps then add lengths
    '''
    if len(clust)==1:
        return int(clust[0][4]) - int(clust[0][3])
    else:
        ret = 0
        p = Popen(['bedtools','merge','-i','-'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        o,e = p.communicate(input='\n'.join('\t'.join(l[:8]) for l in clust))
        for l in o.strip('\n').split('\n'):
            c,s,e = l.strip('\n').split('\t')
            ret += (int(e)-int(s))
        return ret
