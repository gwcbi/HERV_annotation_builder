#! /usr/bin/env python

from HERVAnnotationTool import *
import utils
from collections import Counter
from Bio import SeqIO
import sys
import os
from IGV import *
from itertools import chain
from string import letters as someletters
from glob import glob

def main(args):
    ### Read the GTF file ################################################################
    print >>sys.stderr, 'Loading GTF: %s' % args.internal_file
    gtf = [GTFLine(l) for l in utils.tab_line_gen(open(args.internal_file,'rU'))]
    
    ### Get model lengths
    # mlen = calculate_model_lengths(gtf)
    # print mlen
    mlen = calculate_model_lengths2(gtf)
    print >>sys.stderr, 'Model lengths: %s' %  mlen

    ### Correct the model coordinates ####################################################
    correct_model_coordinates(gtf,mlen)
    for g in gtf:
        if g.strand == '+':
            trueend = mlen[g.attr['repName']] + g.attr['repLeft']
        else:
            trueend = mlen[g.attr['repName']] + g.attr['repStart']
        assert trueend == g.attr['repEnd']

    ### Organize hits by chromosome ######################################################
    bychrom = defaultdict(list)
    for g in gtf:
        bychrom[g.chrom].append(g)

    ### List of HERV loci ################################################################
    print >>sys.stderr, 'Assembling HERV loci'
    all_locs = []

    ### Create HERV loci for plus strand #################################################
    for chrom in utils.CHROMNAMES:
        if chrom in bychrom:
            plus = [h for h in bychrom[chrom] if h.strand == '+']
            if not plus: continue
            plus.sort(key=lambda x: x.start)
            cur = HERVLocus(id='%s_%04d' % (args.prefix, len(all_locs)+1))
            cur.internal.append(plus[0])
            for p1 in plus[1:]:
                p0 = cur.internal[-1]
                # Genomic distance between hits
                gdist = p1.start - p0.end
                # Determine whether p1 is in sequence with locus
                if gdist <= 10: ## Overlapping (or nearly) in genome
                    insequence = True
                else:
                    ## Hits are in sequence and genomic distance is not extreme 
                    insequence = p0.attr['repLeft'] < p1.attr['repLeft']
                    insequence &= gdist < args.longdist
                if insequence:
                    cur.internal.append(p1)
                else:
                    all_locs.append(cur)
                    cur = HERVLocus(id='%s_%04d' % (args.prefix, len(all_locs)+1))
                    cur.internal.append(p1)
            all_locs.append(cur)

    ### Create HERV loci for minus strand ################################################
    for chrom in utils.CHROMNAMES:
        if chrom in bychrom:
            minus = [h for h in bychrom[chrom] if h.strand == '-']
            if not minus: continue
            minus.sort(key=lambda x: x.end, reverse=True) # Sort in reverse order
            cur = HERVLocus(id='%s_%04d' % (args.prefix, len(all_locs)+1))
            cur.internal.append(minus[0])
            for p1 in minus[1:]:
                p0 = cur.internal[-1]
                # Genomic distance between hits
                gdist = p0.start - p1.end
                # Determine whether p1 is in sequence with locus
                if gdist <= 10: ## Overlapping (or nearly) in genome
                    insequence = True
                else:
                    ## Hits are in sequence and genomic distance is not extreme 
                    insequence = p0.attr['repStart'] < p1.attr['repStart']
                    insequence &= gdist < args.longdist
                if insequence:
                    cur.internal.append(p1)
                else:
                    all_locs.append(cur)
                    cur = HERVLocus(id='%s_%04d' % (args.prefix, len(all_locs)+1))
                    cur.internal.append(p1)
            all_locs.append(cur)

    ### Add LTRs to HERV loci ############################################################
    print >>sys.stderr, 'Finding flanking LTRs'    
    for loc in all_locs:
        loc.find_ltr(args.ltr_files, args.flank)
        loc.adjust_overlaps()
    
    print >>sys.stderr, "Initial counts:"
    print >>sys.stderr, '\n'.join('%s%d' % (cat.ljust(20,' '),count) for cat,count in Counter(c.category() for c in all_locs).most_common())

    ### Filtering ########################################################################
    reject = set()
    if args.minpct > 0 or args.mincov > 0:
        print >>sys.stderr, "Removing loci with less than %d percent or %dbp model coverage" % (int(args.minpct*100), args.mincov)
        for loc in all_locs:
            if loc.model_cov() < (mlen[loc.internal_name()] * args.minpct) or loc.model_cov() < args.mincov:
                print >>sys.stderr, '%s\t%d\t%s' % (loc.id, loc.model_cov(), loc.category())
                reject.add(loc)
        
        for rloc in reject:
            all_locs.remove(rloc)
        
        print >>sys.stderr, "After filtering:"
        print >>sys.stderr, '\n'.join('%s%d' % (cat.ljust(20,' '),count) for cat,count in Counter(c.category() for c in all_locs).most_common())
        print >>sys.stderr, '%s%d' % ('Rejected'.ljust(20,' '), len(reject))


    ### Deal with overlapping loci #######################################################
    # Create GTF with all_locs
    with open('tmp.gtf','w') as outh:
        for g in utils.sort_gtf(loc.span_gtf() for loc in all_locs):
            print >>outh, '\t'.join(g)

    # Cluster overlapping and bookended using bedtools
    p1 = Popen('bedtools cluster -i tmp.gtf', shell=True, stdout=PIPE, stderr=PIPE)
    out,err = p1.communicate()
    os.remove('tmp.gtf')

    # Parse bedtools output
    overlap_groups = defaultdict(list)
    for ll in out.strip('\n').split('\n'):
        f = ll.split('\t')
        overlap_groups[f[-1]].append(GTFLine(f[:9]))
    
    # Remove clusters with one
    for k in overlap_groups.keys():
        if len(overlap_groups[k]) == 1:
            del overlap_groups[k]
    
    print >>sys.stderr, "%d overlap groups" % len(overlap_groups)
    
    if args.igv_preview and len(overlap_groups)>0:
        print >>sys.stderr, "Loading IGV"
        # Create file for IGV viewing
        with open('tmp.gtf','w') as outh:
            liter = utils.sort_gtf(chain.from_iterable(loc.each_gtf() for loc in all_locs))
            print >>outh, '\n'.join('\t'.join(_) for _ in liter)
        igv = IGV()
        igv.new()
        igv.genome('hg19')
        igv.load(os.path.join(os.getcwd(),'../other_sources/rmsk_LTR.hg19.gtf'))
        igv.load(os.path.join(os.getcwd(),'tmp.gtf'))

    tandem = []
    for k in sorted(overlap_groups.keys(), key=lambda x:int(x)):
        ogroup = overlap_groups[k]
        if args.igv_preview:
            locus_str = '%s:%s-%s' % (ogroup[0].chrom, min(gl.start for gl in ogroup)-5000, max(gl.end for gl in ogroup)+5000)
            igv.goto(locus_str)
            igv.expand()
        
        # Get locus for each member of overlap group
        og_locus = {}
        for o in ogroup:
            tmp = [c for c in all_locs if c.id == o.attr['name']]
            assert len(tmp)==1
            og_locus[o.attr['name']] = tmp[0]
        # Print out the model coverage
        for n,loc in og_locus.iteritems():
            print >>sys.stderr, '%s\t%d\t%s' % (n, loc.model_cov(), loc.category())

        # Parse user input
        z = raw_input('Action to take: ').strip()
        if z == '': continue
        inputcmd = z.strip().split(' ')
        if inputcmd[0] == 'REJECT':
            if len(inputcmd) == 1:
                # Only max will be kept
                st = sorted([loc for n,loc in og_locus.iteritems()], key=lambda x:x.model_cov(), reverse=True)[1:]
                loc_ids = [_.id for _ in st]
            elif len(inputcmd) == 2:
                loc_ids = inputcmd[1].split(',')
            else:
                assert False
            for loc_id in loc_ids:
                reject.add(og_locus[loc_id])
        elif inputcmd[0] == 'TANDEM':
            if len(inputcmd) == 1:
                assert len(og_locus)==2, 'More than 2 loci are present'
                tandem.append([loc for n,loc in og_locus.iteritems()])
            elif len(inputcmd) == 2:
                loc_ids = inputcmd[1].split('+')
                tandem.append([og_locus[loc_id] for loc_id in loc_ids])
            else:
                assert False
        elif inputcmd[0] == 'DIFF':
            n1,n2 = inputcmd[1].split('-')
            g1 = og_locus[n1]
            g2 = og_locus[n2]
            if g1.span()[0] < g2.span()[1]:
                g1.shorten(g2.span()[1]+20, g1.span()[1])
            elif g1.span()[1] < g2.span()[0]:
                g1.shorten(g1.span()[0], g2.span()[0]-20)
            else:
                print "no overlap!"
            print g1
        elif inputcmd[0] == 'IGNORE':
            continue
        else:
            assert False, 'Unknown command: "%s"' % inputcmd[0]

    # Remove rejected annotations
    for rloc in reject:
        if rloc in all_locs:
            all_locs.remove(rloc)
    
    # Create the tandem annotations
    for tgroup in tandem:
        tandem_loc = HERVLocus(id=tgroup[0].id)
        tandem_loc.internal = list(chain.from_iterable(loc.internal for loc in tgroup))
        if tandem_loc.strand() == '+':
            tandem_loc.internal.sort(key=lambda x:x.start)
        else:
            tandem_loc.internal.sort(key=lambda x:x.end, reverse=True)
    
        tandem_loc.find_ltr(args.ltr_files, 1000)
        tandem_loc.adjust_overlaps()
        tandem_loc.is_tandem = True
        all_locs.append(tandem_loc)
        # Remove from original
        for rloc in tgroup:
            all_locs.remove(rloc)
    
    print >>sys.stderr, "After overlap removal:"
    print >>sys.stderr, '\n'.join('%s%d' % (cat.ljust(20,' '),count) for cat,count in Counter(c.category() for c in all_locs).most_common())
    print >>sys.stderr, '%s%d' % ('Rejected'.ljust(20,' '), len(reject))
    if args.igv_preview and len(overlap_groups)>0: os.remove('tmp.gtf')
    
    ### Sort loci ########################################################################
    bychrom = defaultdict(list)
    for loc in all_locs:
        bychrom[loc.chrom()].append(loc)
    
    final_locs = []
    for chrom in utils.CHROMNAMES:
        if chrom in bychrom:
            for loc in sorted(bychrom[chrom], key=lambda x:x.span()[0]):
                final_locs.append(loc)
    
    for i,loc in enumerate(final_locs):
        loc.id = '%s_%04d' % (args.prefix, i+1)

    ### Rename loci according to cytoband #################################################
    # Create GTF with all_locs
    with open('tmp.gtf','w') as outh:
        for g in utils.sort_gtf(loc.span_gtf() for loc in final_locs):
            print >>outh, '\t'.join(g)    

    p1 = Popen('bedtools intersect -wo -a tmp.gtf -b ../other_sources/cytoband.gtf', shell=True, stdout=PIPE, stderr=PIPE)
    out,err = p1.communicate()
    os.remove('tmp.gtf')

    byband = defaultdict(list)
    for ll in out.strip('\n').split('\n'):
        f = ll.split('\t')
        g1 = GTFLine(f[:9])
        g2 = GTFLine(f[9:-1])
        band = '%s%s' % (g2.chrom.strip('chr'),g2.attr['gene_id'])
        byband[band].append(g1)

    namemap = {}
    for band,glist in byband.iteritems():
        if len(glist) == 1:
            namemap[glist[0].attr['name']] = '%s_%s' % (args.prefix, band)
        else:
            glist.sort(key=lambda x:x.start)
            for i,gl in enumerate(glist):
                namemap[gl.attr['name']] = '%s_%s%s' % (args.prefix, band, someletters[i])

    for loc in final_locs:
        loc.locus_name = namemap[loc.id]

    ### Create annotation files ##########################################################
    print >>sys.stderr, "Writing annotation files"
    with open('%s.gtf' % args.prefix,'w') as outh:
        liter = utils.sort_gtf(chain.from_iterable(loc.each_gtf() for loc in final_locs))
        print >>outh, '\n'.join('\t'.join(_) for _ in liter)
        # for loc in final_locs:
        #     print >>outh, '\n'.join('\t'.join(g) for g in loc.each_gtf())

    with open('%s_reject.gtf' % args.prefix,'w') as outh:
        liter = utils.sort_gtf(chain.from_iterable(loc.each_gtf() for loc in reject))
        print >>outh, '\n'.join('\t'.join(_) for _ in liter)        
        # for loc in reject:
        #     print >>outh, '\n'.join('\t'.join(g) for g in loc.each_gtf())
    
    with open('%s_span.gtf' % args.prefix,'w') as outh:
        for g in utils.sort_gtf(loc.span_gtf() for loc in final_locs):
            print >>outh, '\t'.join(g) 
    
    with open('%s_table.txt' % args.prefix,'w') as outh:
        print >>outh, '\t'.join(['locus_name','id','strand','chrom','start','end','strand','nfeats','width','model_cov','ltr5_model','int_model','ltr3_model'])
        for loc in final_locs:
            mgtf = GTFLine(loc.span_gtf())
            row = [loc.locus_name, loc.id, loc.category(),
                   mgtf.chrom, mgtf.start, mgtf.end, mgtf.strand,
                   mgtf.attr['nfeats'], loc.width(), loc.model_cov(),
                   loc.ltr_up_name(), loc.internal_name(), loc.ltr_down_name(),
                   ]
            print >>outh, '\t'.join(str(_) for _ in row)

    ### Extract sequences ################################################################
    if args.get_sequences:
        print >>sys.stderr, "Extracting sequences"
        genome_fasta = args.genome_fasta # '/Users/bendall/Projects/References/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
        genome = dict((s.id,s) for s in SeqIO.parse(genome_fasta,'fasta'))
        
        with open('%s.full.fasta' % args.prefix,'w') as outh:
            for loc in final_locs:
                gcoord = '%s:%d-%d(%s)' % (loc.chrom(), loc.span()[0], loc.span()[1], loc.strand())
                print >>outh, '>%s|%s|%s' % (loc.locus_name, loc.category(), gcoord)
                print >>outh, str(loc.entire_sequence(genome).seq)
        
        with open('%s.internal.fasta' % args.prefix,'w') as outh:
            for loc in final_locs:
                gcoord = '%s:%d-%d(%s)' % (loc.chrom(), min(p.start for p in loc.internal), max(p.end for p in loc.internal), loc.strand())
                print >>outh, '>%s_int|%s|%s|%s' % (loc.locus_name, loc.category(), gcoord, loc.format_print_clust())
                print >>outh, str(loc.internal_sequence(genome).seq)
        
        with open('%s.5ltr.fasta' % args.prefix,'w') as outh:
            for loc in final_locs:
                ltrseq = loc.ltr_up_sequence(genome)
                if ltrseq:
                    gcoord = '%s:%d-%d(%s)' % (loc.chrom(), min(p.start for p in loc.ltr_up), max(p.end for p in loc.ltr_up), loc.strand())
                    print >>outh, '>%s_5LTR|%s|%s' % (loc.locus_name, loc.ltr_up_name(), gcoord)
                    print >>outh, str(ltrseq.seq)
        
        with open('%s.3ltr.fasta' % args.prefix,'w') as outh:
            for loc in final_locs:
                ltrseq = loc.ltr_down_sequence(genome)
                if ltrseq:
                    gcoord = '%s:%d-%d(%s)' % (loc.chrom(), min(p.start for p in loc.ltr_down), max(p.end for p in loc.ltr_down), loc.strand())
                    print >>outh, '>%s_3LTR|%s|%s' % (loc.locus_name, loc.ltr_down_name(), gcoord)
                    print >>outh, str(ltrseq.seq)

    ### IGV snapshots ####################################################################
    if args.igv_snapshot:
            print >>sys.stderr, "Taking IGV snapshots"
            igv = IGV()
            igv.new()
            igv.genome('hg19')
            igv.load(os.path.join(os.getcwd(),'../other_sources/rmsk_LTR.hg19.gtf'))
            if os.path.isdir('tmp'):
                for compare_gtf in glob('tmp/*.gtf'):
                    igv.load(os.path.join(os.getcwd(), compare_gtf))
            
            igv.load(os.path.join(os.getcwd(),'%s.gtf' % args.prefix))
            igv.load(os.path.join(os.getcwd(),'%s_reject.gtf' % args.prefix))
            
            do_snapshots = True
            
            if do_snapshots:
                if not os.path.exists(os.path.join(os.getcwd(),'snapshots')):
                    os.mkdir(os.path.join(os.getcwd(),'snapshots'))
                if not os.path.exists(os.path.join(os.getcwd(),'reject')):
                    os.mkdir(os.path.join(os.getcwd(),'reject'))
                
                categories = ['prototype', 'oneside', 'internal']
                for cat in categories:
                    if not os.path.exists(os.path.join(os.getcwd(),'snapshots/%s' % cat)):
                        os.mkdir(os.path.join(os.getcwd(),'snapshots/%s' % cat))    
                    if not os.path.exists(os.path.join(os.getcwd(),'reject/%s' % cat)):
                        os.mkdir(os.path.join(os.getcwd(),'reject/%s' % cat))
                        
            for loc in final_locs:
                rc,lc = loc.span()
                locus_str = '%s:%d-%d' % (loc.chrom(), rc-5000, lc+5000)
                print >>sys.stderr, '%s\t%s\t%s' % (loc.locus_name, loc.category(), locus_str)
                igv.goto(locus_str)
                igv.expand()
                if do_snapshots:
                    igv.snapshotDirectory(os.path.join(os.getcwd(),'snapshots/%s' % loc.category().strip('*') ))
                    igv.snapshot(filename='%s.png' % loc.locus_name)
            
            for loc in reject:
                rc,lc = loc.span()
                locus_str = '%s:%d-%d' % (loc.chrom(), rc-5000, lc+5000)
                print >>sys.stderr, '%s\t%s\t%s' % (loc.id, loc.category(), locus_str)
                igv.goto(locus_str)
                igv.expand()
                if do_snapshots:
                    igv.snapshotDirectory(os.path.join(os.getcwd(),'reject/%s' % loc.category().strip('*') ))
                    igv.snapshot(filename='%s.png' % loc.id)    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Build the annotation')
    parser.add_argument('--longdist', type=int, default=10000, help='An extreme distance. Annotations with genomic distance greater than longdist are considered seperate loci.')
    parser.add_argument('--minpct', type=float, default=0.0, help='Minimum percentage of model coverage to be included as locus')
    parser.add_argument('--mincov', type=int, default=0, help='Minimum length of model coverage to be included as locus')
    parser.add_argument('--flank', type=int, default=1000, help='Distance to search for flanking LTRs')
    
    parser.add_argument('--igv_preview', action='store_true') 
    parser.add_argument('--igv_snapshot', action='store_true')
    parser.add_argument('--get_sequences', action='store_true')    
    parser.add_argument('--genome_fasta', default='/Users/bendall/Projects/References/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa')
    
    parser.add_argument('internal_file')
    parser.add_argument('ltr_files')
    parser.add_argument('prefix')    
    main(parser.parse_args())


"""




##########################################################################################
todelete = set()
with open('tmp.gtf','w') as outh:
    for g in utils.sort_gtf(clust.span_gtf() for clust in filter):
        print >>outh, '\t'.join(g)

p1 = Popen('bedtools intersect -wo -a tmp.gtf -b tmp.gtf', shell=True, stdout=PIPE, stderr=PIPE)
out,err = p1.communicate()
os.remove('tmp.gtf')

for ll in out.strip('\n').split('\n'):
    f = ll.split('\t')
    g1 = GTFLine(f[:9])
    g2 = GTFLine(f[9:-1])
    if g1.attr['name'] != g2.attr['name']:
        # Decide what to do
        if g1.start < g2.start < g1.end and g1.start < g2.end < g1.end:
            print >>sys.stderr, '%s contained within %s' % (g2.attr['name'],g1.attr['name'])
            todelete.add(g2.attr['name'])
        elif g2.start < g1.start < g2.end and g2.start < g1.end < g2.end:
            print >>sys.stderr, '%s contained within %s' % (g1.attr['name'],g2.attr['name'])
            todelete.add(g1.attr['name'])            
        elif g1.attr['category'].startswith('prototype') and not g2.attr['category'].startswith('prototype'):
            print >>sys.stderr, '%s not prototype, %s is prototype' % (g2.attr['name'],g1.attr['name'])
            todelete.add(g2.attr['name'])            
        elif g2.attr['category'].startswith('prototype') and not g1.attr['category'].startswith('prototype'):
            print >>sys.stderr, '%s not prototype, %s is prototype' % (g1.attr['name'],g2.attr['name'])
            todelete.add(g1.attr['name'])            
        else:
            print >>sys.stderr, 'Unusual case\n%s\n%s' % (g1,g2)

if len(todelete)>0:
    for clust in filter:
        if clust.id in todelete:
            reject.append(clust)
            filter.remove(clust)

if len(todelete)>0:
    print >>sys.stderr, "Filtered %d overlapping loci." % len(todelete)
    print >>sys.stderr, "After filtering:"
    print >>sys.stderr, '\n'.join('%s%d' % (cat.ljust(20,' '),count) for cat,count in Counter(c.category() for c in filter).most_common())

##########################################################################################
sys.exit()

### Sort loci ############################################################################
bychrom = defaultdict(list)
for loc in filter:
    bychrom[loc.chrom()].append(loc)

all_loci = []
for chrom in utils.CHROMNAMES:
    if chrom in bychrom:
        for loc in sorted(bychrom[chrom], key=lambda x:x.span()[0]):
            all_loci.append(loc)

for i,loc in enumerate(all_loci):
    loc.id = '%s_%04d' % (prefix, i+1)

### Rename loci according to cytoband ####################################################
with open('tmp.gtf','w') as outh:
    for loc in all_loci:
        print >>outh, '\t'.join(loc.span_gtf())

p1 = Popen('bedtools intersect -wo -a tmp.gtf -b ../other_sources/cytoband.gtf', shell=True, stdout=PIPE, stderr=PIPE)
out,err = p1.communicate()
os.remove('tmp.gtf')

byband = defaultdict(list)
for ll in out.strip('\n').split('\n'):
    f = ll.split('\t')
    g1 = GTFLine(f[:9])
    g2 = GTFLine(f[9:-1])
    locid = g1.attr['name']
    spos = g1.start
    band = '%s%s' % (g2.chrom.strip('chr'),g2.attr['gene_id'])
    byband[band].append((locid,spos))

from string import letters as someletters
namemap = {}
for band,llist in byband.iteritems():
    if len(llist) == 1:
        namemap[llist[0][0]] = '%s_%s' % (prefix,band)
    else:
        llist.sort(key=lambda x:x[1])
        for i,tup in enumerate(llist):
            namemap[tup[0]] = '%s_%s%s' % (prefix,band,someletters[i])

for loc in all_loci:
    loc.locus_name = namemap[loc.id]

print >>sys.stderr, "Writing annotation files"

with open('%s.gtf' % prefix,'w') as outh:
    for loc in all_loci:
        print >>outh, '\n'.join('\t'.join(g) for g in loc.each_gtf())

with open('%s_reject.gtf' % prefix,'w') as outh:
    for loc in reject:
        print >>outh, '\n'.join('\t'.join(g) for g in loc.each_gtf())

with open('%s_span.gtf' % prefix,'w') as outh:
    for loc in all_loci:
        print >>outh, '\t'.join(loc.span_gtf())

with open('%s_table.txt' % prefix,'w') as outh:
    print >>outh, '\t'.join(['locus_name','id','strand','chrom','start','end','strand','nfeats','width','model_cov','ltr5_model','int_model','ltr3_model'])
    for loc in all_loci:
        mgtf = GTFLine(loc.span_gtf())
        row = [loc.locus_name, loc.id, loc.category(),
               mgtf.chrom, mgtf.start, mgtf.end, mgtf.strand,
               mgtf.attr['nfeats'], loc.width(), loc.model_cov(),
               loc.ltr_up_name(), loc.internal_name(), loc.ltr_down_name(),
               ]
        print >>outh, '\t'.join(str(_) for _ in row)


### Extract sequences ####################################################################
print >>sys.stderr, "Extracting sequences"
genome_fasta = '/Users/bendall/Projects/References/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
genome = dict((s.id,s) for s in SeqIO.parse(genome_fasta,'fasta'))

with open('%s.full.fasta' % prefix,'w') as outh:
    for loc in all_loci:
        gcoord = '%s:%d-%d(%s)' % (loc.chrom(), loc.span()[0], loc.span()[1], loc.strand())
        print >>outh, '>%s|%s|%s' % (loc.locus_name, loc.category(), gcoord)
        print >>outh, str(loc.entire_sequence(genome).seq)

with open('%s.internal.fasta' % prefix,'w') as outh:
    for loc in all_loci:
        gcoord = '%s:%d-%d(%s)' % (loc.chrom(), min(p.start for p in loc.internal), max(p.end for p in loc.internal), loc.strand())
        print >>outh, '>%s_int|%s|%s|%s' % (loc.locus_name, loc.category(), gcoord, loc.format_print_clust())
        print >>outh, str(loc.internal_sequence(genome).seq)

with open('%s.5ltr.fasta' % prefix,'w') as outh:
    for loc in all_loci:
        ltrseq = loc.ltr_up_sequence(genome)
        if ltrseq:
            gcoord = '%s:%d-%d(%s)' % (loc.chrom(), min(p.start for p in loc.ltr_up), max(p.end for p in loc.ltr_up), loc.strand())
            print >>outh, '>%s_5LTR|%s|%s' % (loc.locus_name, loc.ltr_up_name(), gcoord)
            print >>outh, str(ltrseq.seq)

with open('%s.3ltr.fasta' % prefix,'w') as outh:
    for loc in all_loci:
        ltrseq = loc.ltr_down_sequence(genome)
        if ltrseq:
            gcoord = '%s:%d-%d(%s)' % (loc.chrom(), min(p.start for p in loc.ltr_down), max(p.end for p in loc.ltr_down), loc.strand())
            print >>outh, '>%s_3LTR|%s|%s' % (loc.locus_name, loc.ltr_down_name(), gcoord)
            print >>outh, str(ltrseq.seq)

# IGV
print >>sys.stderr, "Taking IGV snapshots"
from IGV import *

igv = IGV()
igv.new()
igv.genome('hg19')
igv.load(os.path.join(os.getcwd(),'../other_sources/rmsk_LTR.hg19.gtf'))
# igv.load(os.path.join(os.getcwd(),'../HERV9/tmp/concat.gtf'))
igv.load(os.path.join(os.getcwd(),'../HERV9/tmp/prototype.gtf'))
igv.load(os.path.join(os.getcwd(),'../HERV9/tmp/oneside.gtf'))
igv.load(os.path.join(os.getcwd(),'../HERV9/tmp/soloint.gtf'))
igv.load(os.path.join(os.getcwd(),'../HERV9/tmp/sololtr.gtf'))
igv.load(os.path.join(os.getcwd(),'../HERV9/tmp/unusual.gtf'))
igv.load(os.path.join(os.getcwd(),'%s.gtf' % prefix))
igv.load(os.path.join(os.getcwd(),'%s_reject.gtf' % prefix))

do_snapshots = True

if do_snapshots:
    if not os.path.exists(os.path.join(os.getcwd(),'snapshots')):
        os.mkdir(os.path.join(os.getcwd(),'snapshots'))
    if not os.path.exists(os.path.join(os.getcwd(),'reject')):
        os.mkdir(os.path.join(os.getcwd(),'reject'))
    
    categories = ['prototype', 'oneside', 'internal']
    for cat in categories:
        if not os.path.exists(os.path.join(os.getcwd(),'snapshots/%s' % cat)):
            os.mkdir(os.path.join(os.getcwd(),'snapshots/%s' % cat))    
        if not os.path.exists(os.path.join(os.getcwd(),'reject/%s' % cat)):
            os.mkdir(os.path.join(os.getcwd(),'reject/%s' % cat))
            
for loc in all_loci:
    rc,lc = loc.span()
    locus_str = '%s:%d-%d' % (loc.chrom(), rc-5000, lc+5000)
    print >>sys.stderr, '%s\t%s\t%s' % (loc.locus_name, loc.category(), locus_str)
    igv.goto(locus_str)
    igv.expand()
    # z = raw_input('Press enter to see next annotation')
    if do_snapshots:
        igv.snapshotDirectory(os.path.join(os.getcwd(),'snapshots/%s' % loc.category().strip('*') ))
        igv.snapshot(filename='%s.png' % loc.locus_name)

for loc in reject:
    rc,lc = loc.span()
    locus_str = '%s:%d-%d' % (loc.chrom(), rc-5000, lc+5000)
    print >>sys.stderr, '%s\t%s\t%s' % (loc.id, loc.category(), locus_str)
    igv.goto(locus_str)
    igv.expand()
    # z = raw_input('Press enter to see next annotation')
    if do_snapshots:
        igv.snapshotDirectory(os.path.join(os.getcwd(),'reject/%s' % loc.category().strip('*') ))
        igv.snapshot(filename='%s.png' % loc.id)




"""
