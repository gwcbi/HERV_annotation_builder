#! /usr/bin/env python


def convert(infile,outfile,source):
    from collections import defaultdict
    import math
    chrnames = ['chr%d' % d for d in range(1,23)] + ['chrX','chrY']

    alllines = defaultdict(list)  
    lines = (l.strip('\n').split('\t') for l in infile)
    header = lines.next()

    # Find column indices
    cidx = header.index('genoName')
    sidx = header.index('genoStart')
    eidx = header.index('genoEnd')
    stidx = header.index('strand')
 
    # Organize original lines by chromosome
    for l in lines:
        alllines[ l[cidx] ].append(l)

    # Sorting and formatting lines
    counter = 1
    for cchrom in chrnames:
        if cchrom not in alllines: continue 
        for l in sorted(alllines[cchrom],key=lambda x:int(x[sidx])):
            spos = l[sidx] # if l[stidx]=='+' else l[eidx]
            epos = l[eidx] # if l[stidx]=='+' else l[sidx]
            assert int(spos)<int(epos), "ERROR start >= end %s" % l      
      
            score = l[header.index('swScore')]
            
            # Attributes
            attrs = {}
            # Total length of alignment on reference
            # attrs['aln_len'] = '%d' % (int(epos) - int(spos))
            
            # Percent of repeat sequence covered by alignment
            # if l[header.index('strand')] == '+':
            #     replen = int(l[header.index('repEnd')]) - int(l[header.index('repLeft')])
            #     attrs['pct_cov'] = '%.1f' % ((float(int(l[header.index('repEnd')]) - int(l[header.index('repStart')])) / replen) *100)
            #     repstart = l[header.index('repStart')]
            #     
            # else:
            #     replen = int(l[header.index('repEnd')]) - int(l[header.index('repStart')])
            #     attrs['pct_cov'] = '%.1f' % ((float(int(l[header.index('repEnd')]) - int(l[header.index('repLeft')])) / replen) *100)
            #     repstart = l[header.index('repLeft')]
            #
            # repend = l[header.index('repEnd')]
            # attrs['repAln'] = '%s-%s' % (repstart,repend)

            attrs['repStart'] = l[header.index('repStart')]
            attrs['repEnd']   = l[header.index('repEnd')]
            attrs['repLeft']  = l[header.index('repLeft')]

            attrs['id'] = '%s_%d' % (l[header.index('repName')], counter)
            attrs['repName'] = l[header.index('repName')]
            attrs['repType'] = 'internal' if l[header.index('repName')].endswith('-int') else 'ltr'
            attrs['repClass'] = l[header.index('repClass')]
            attrs['repFamily'] = l[header.index('repFamily')]
            # attrs['similarity'] = '%s:%s:%s' % (l[header.index('milliDiv')],l[header.index('milliDel')],l[header.index('milliIns')])
            # attrs['bin'] = l[header.index('bin')]
            attr = ' '.join('%s "%s";' % (k,v) for k,v in attrs.iteritems())
            print >>outfile, '\t'.join([cchrom,source,'exon',spos,epos,score,l[stidx],'.',attr])
            counter += 1

def main(parser):
    args = parser.parse_args()
    convert(args.infile, args.outfile, args.source)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='convert dfam hits to GTF')
    parser.add_argument('--source', default='rmsk')  
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser)
