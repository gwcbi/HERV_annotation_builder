import sys
import re
import utils
from collections import defaultdict, Counter
from subprocess import Popen, PIPE

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
from glob import glob

GTFCOLS = ['chrom','source','feature','start','end','score','strand','frame']

def calculate_model_lengths(gtflines):
    ''' Calculate length of models from lines where repLeft == 0 '''
    ret = {}
    for g in gtflines:
        if g.attr['repLeft']==0:
            if g.attr['repName'] in ret:
                ret[g.attr['repName']] = max(ret[g.attr['repName']], g.attr['repEnd']) #assert ret[g.attr['repName']] == g.attr['repEnd']
            else:
                ret[g.attr['repName']] = g.attr['repEnd']
    return ret

def calculate_model_lengths2(gtflines):
    ''' Calculate length of models from all lines and take most common'''
    lengthlist = defaultdict(list)
    for g in gtflines:
        if g.strand == '+':
            lengthlist[g.attr['repName']].append(g.attr['repEnd'] - g.attr['repLeft'])
        else:
            lengthlist[g.attr['repName']].append(g.attr['repEnd'] - g.attr['repStart'])
    # for k,v in lengthlist.iteritems():
    #     print Counter(v).most_common()
    #     
    return dict((k,Counter(v).most_common()[0][0]) for k,v in lengthlist.iteritems())

def correct_model_coordinates(gtflines, modellen):
    ''' Fix incorrect model coordinates using model length '''
    for g in gtflines:
        if g.strand == '+':
            trueend = modellen[g.attr['repName']] + g.attr['repLeft']
            if trueend != g.attr['repEnd']:
                replen = g.attr['repEnd'] - g.attr['repStart']
                g.attr['repEnd'] = trueend
                g.attr['repStart'] = trueend - replen
        else:
            trueend = modellen[g.attr['repName']] + g.attr['repStart']
            if trueend != g.attr['repEnd']:
                replen = g.attr['repEnd'] - g.attr['repLeft']
                g.attr['repEnd'] = trueend
                g.attr['repLeft'] = trueend - replen

def flank_left(g,d):
    ret = g[:]
    ret[3] = '%d' % max(1, (int(g[3]) - d))
    ret[4] = '%d' % max(1, (int(g[3]) - 1))
    return ret

def flank_right(g,d):
    ret = g[:]
    ret[3] = '%d' % min(utils.CHROMSIZE[g[0]], (int(g[4]) + 1 ))
    ret[4] = '%d' % min(utils.CHROMSIZE[g[0]], (int(g[4]) + d ))
    return ret


class GTFLine:
    def __init__(self,row):
        for c,v in zip(GTFCOLS, row[:8]):
            try:
                setattr(self,c,int(v))
            except ValueError:
                setattr(self,c,v)
        _attrd = dict(re.findall('(\S+)\s+"([\s\S]+?)";',row[8]))
        self.attr = {}
        for k,v in _attrd.iteritems():
            try:
                self.attr[k] = int(v)
            except ValueError:
                self.attr[k] = v
    
    def fmt(self):
        _attrs = ' '.join('%s "%s";' % (k,v) for k,v in self.attr.iteritems())
        return [str(getattr(self,c)) for c in GTFCOLS] + [_attrs]
    
    def sequence(self, gdict):
        if self.strand == '+':
            return gdict[self.chrom][self.start:self.end]
        else:
            return gdict[self.chrom][self.start:self.end].reverse_complement()
    
    def __str__(self):
        return '\t'.join(self.fmt())

class HERVLocus:
    def __init__(self, id=None):
        self.id = id
        self.internal = []
        self.ltr_up   = []
        self.ltr_down = []
        self.ltr_int  = []
        self.locus_name = None 
        self.is_tandem = False       

    def chrom(self):
        return self.internal[0].chrom

    def strand(self):
        return self.internal[0].strand
    
    def category(self):
        if self.is_tandem: return 'tandem'
        if len(self.ltr_up)==0 and len(self.ltr_down)==0:
            ret = 'internal'
        elif len(self.ltr_up)>0 and len(self.ltr_down)>0:
            ret = 'prototype'
        elif len(self.ltr_up)>0 or len(self.ltr_down)>0:
            ret = 'oneside'
        else:
            ret = 'unknown'
        if len(self.ltr_int)>0:
            ret += '*'
        return ret

    def internal_name(self):
        return self.internal[0].attr['repName']
    
    def ltr_up_name(self):
        if len(self.ltr_up)==0: return ''
        longest = sorted(self.ltr_up, key=lambda x: x.end-x.start, reverse=True)[0]
        return longest.attr['repName']

    def ltr_down_name(self):
        if len(self.ltr_down)==0: return ''
        longest = sorted(self.ltr_down, key=lambda x: x.end-x.start, reverse=True)[0]
        return longest.attr['repName']
    
    def span(self):
        _allann = self.ltr_up + self.internal + self.ltr_int + self.ltr_down
        _spos = min(v.start for v in _allann)
        _epos = max(v.end for v in _allann)
        return  min(v.start for v in _allann),max(v.end for v in _allann)

    def width(self):
        l,r = self.span()
        return r - l
    
    def adjust_overlaps(self):
        _allann = self.ltr_up + self.internal + self.ltr_int + self.ltr_down
        if self.strand() == '+':
            _allann.sort(key=lambda x:x.start)
            for i,p1 in enumerate(_allann[1:]):
                p0 = _allann[i]
                gdist = p1.start - p0.end
                if gdist < 0:
                    p1.start = p1.start - gdist
                    p1.attr['repStart'] = p1.attr['repStart'] - gdist                   
        else:
            _allann.sort(key=lambda x:x.end, reverse=True)
            for i,p1 in enumerate(_allann[1:]):
                p0 = _allann[i]
                gdist = p0.start - p1.end
                if gdist < 0:
                    p1.end = p1.end + gdist 
                    p1.attr['repLeft'] = p1.attr['repLeft'] - gdist
    
    def adjust_genomic_overlaps(self):
        if self.strand() == '+':
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                gdist = p1.start - p0.end
                if gdist < 0:
                    p1.start = p1.start - gdist
                    p1.attr['repStart'] = p1.attr['repStart'] - gdist
        else:
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                gdist = p0.start - p1.end
                if gdist < 0:
                    p1.end = p1.end + gdist
                    p1.attr['repLeft'] = p1.attr['repLeft'] - gdist
    
    def adjust_ltr_overlaps(self):
        if self.strand() == '+':
            for i,p1 in enumerate(self.ltr_up[1:]):
                p0 = self.ltr_up[i]
                gdist = p1.start - p0.end
                if gdist < 0:
                    p1.start = p1.start - gdist
                    p1.attr['repStart'] = p1.attr['repStart'] - gdist
            for i,p1 in enumerate(self.ltr_down[1:]):
                p0 = self.ltr_down[i]
                gdist = p1.start - p0.end
                if gdist < 0:
                    p1.start = p1.start - gdist
                    p1.attr['repStart'] = p1.attr['repStart'] - gdist
        else:
            for i,p1 in enumerate(self.ltr_up[1:]):
                p0 = self.ltr_up[i]
                gdist = p0.start - p1.end
                if gdist < 0:
                    p1.end = p1.end + gdist
                    p1.attr['repLeft'] = p1.attr['repLeft'] - gdist
            for i,p1 in enumerate(self.ltr_down[1:]):
                p0 = self.ltr_down[i]
                gdist = p0.start - p1.end
                if gdist < 0:
                    p1.end = p1.end + gdist
                    p1.attr['repLeft'] = p1.attr['repLeft'] - gdist
    
    def shorten(self, spos, epos):
        _allann = self.ltr_up + self.internal + self.ltr_int + self.ltr_down
        for p0 in _allann:
            if p0.start < spos:
                diff = spos - p0.start + 1
                p0.start = p0.start + diff
                if p0.strand == '+':
                    p0.attr['repStart'] = p0.attr['repStart'] + diff
                else:
                    p0.attr['repEnd'] = p0.attr['repEnd'] - diff
            if p0.end > epos:
                diff = p0.end - epos - 1
                p0.end = p0.end - diff
                if p0.strand == '+':
                    p0.attr['repEnd'] = p0.attr['repEnd'] - diff
                else:
                    p0.attr['repLeft'] = p0.attr['repLeft'] + diff
    
    def format_print_clust(self):
        if self.strand() == '+':
            s = '[%s-%s]' % (self.internal[0].attr['repStart'], self.internal[0].attr['repEnd'])
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                gdist = p1.start - p0.end
                # mdist = p1.attr['repStart'] - p0.attr['repEnd']
                if gdist != 0:
                    s += '--%d--' % gdist
                s += '[%s-%s]' % (p1.attr['repStart'], p1.attr['repEnd'])
        else:
            s = '[%s-%s]' % (self.internal[0].attr['repLeft'], self.internal[0].attr['repEnd'])
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                gdist = p0.start - p1.end
                # mdist = p1.attr['repLeft'] - p0.attr['repEnd']
                if gdist != 0:
                    s += '--%d--' % gdist
                s += '[%s-%s]' % (p1.attr['repLeft'], p1.attr['repEnd'])
        return s
    
    def model_cov0(self):
        if self.strand() == '+':
            ret = self.internal[0].attr['repEnd'] - self.internal[0].attr['repStart']
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                ret += p1.attr['repEnd'] - p1.attr['repStart']
                if p1.attr['repStart'] < p0.attr['repEnd']:
                    ret += p1.attr['repStart'] - p0.attr['repEnd']
        else:
            ret = self.internal[0].attr['repEnd'] - self.internal[0].attr['repLeft']
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                ret += p1.attr['repEnd'] - p1.attr['repLeft']
                if p1.attr['repLeft'] < p0.attr['repEnd']:
                    ret += p1.attr['repLeft'] - p0.attr['repEnd']
        return ret

    def model_cov(self):
        ret = 0    
        if self.strand() == '+':
            for p0 in self.internal:
                ret += p0.attr['repEnd'] - p0.attr['repStart']
        else:
            for p0 in self.internal:
                ret += p0.attr['repEnd'] - p0.attr['repLeft']
        return ret

    
    '''
    def internal_sequence(self,gdict):
        if self.strand() == '+':
            iseq = self.internal[0].sequence(gdict).upper()
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                assert not p1.start < p0.end, "error: overlapping intervals:\n%s\n%s" % (p0,p1)
                if p1.start > p0.end:
                    iseq += gdict[p1.chrom][p0.end:p1.start].lower()
                iseq += p1.sequence(gdict).upper()
        else:
            iseq = self.internal[0].sequence(gdict).upper()
            for i,p1 in enumerate(self.internal[1:]):
                p0 = self.internal[i]
                assert not p0.start < p1.end, "error: overlapping intervals:\n%s\n%s" % (p0,p1)
                if p0.start > p1.end:
                    iseq += gdict[p1.chrom][p1.end:p0.start].reverse_complement().lower()
                iseq += p1.sequence(gdict).upper()
        return iseq
    '''
    def _get_seqlist(self, slist, gdict):
        if len(slist)==0: return Seq('',alphabet=SingleLetterAlphabet())
        if self.strand() == '+':
            iseq = slist[0].sequence(gdict).upper()
            for i,p1 in enumerate(slist[1:]):
                p0 = slist[i]
                assert not p1.start < p0.end, "error: overlapping intervals:\n%s\n%s" % (p0,p1)
                if p1.start > p0.end:
                    iseq += gdict[p1.chrom][p0.end:p1.start].lower()
                iseq += p1.sequence(gdict).upper()
        else:
            iseq = slist[0].sequence(gdict).upper()
            for i,p1 in enumerate(slist[1:]):
                p0 = slist[i]
                assert not p0.start < p1.end, "error: overlapping intervals:\n%s\n%s" % (p0,p1)
                if p0.start > p1.end:
                    iseq += gdict[p1.chrom][p1.end:p0.start].reverse_complement().lower()
                iseq += p1.sequence(gdict).upper()
        return iseq

    def internal_sequence(self,gdict):
        return self._get_seqlist(self.internal, gdict)

    def ltr_up_sequence(self, gdict):
        return self._get_seqlist(self.ltr_up, gdict)

    def ltr_down_sequence(self, gdict):
        return self._get_seqlist(self.ltr_down, gdict)

    def entire_sequence(self, gdict):
        l5seq = self.ltr_up_sequence(gdict)
        if l5seq:
            if self.strand() == '+':
                lc = max(p.end for p in self.ltr_up)
                rc = min(p.start for p in self.internal)
                if lc < rc:
                    l5seq += gdict[self.chrom()][lc:rc].lower()
            else:
                lc = max(p.end for p in self.internal)
                rc = min(p.start for p in self.ltr_up)
                if lc < rc:
                    l5seq += gdict[self.chrom()][lc:rc].reverse_complement().lower()

        l3seq = self.ltr_down_sequence(gdict)
        if l3seq:
            if self.strand() == '+':
                lc = max(p.end for p in self.internal)
                rc = min(p.start for p in self.ltr_down)
                if lc < rc:
                    l3seq = gdict[self.chrom()][lc:rc].lower() + l3seq
            else:
                lc = max(p.end for p in self.ltr_down)
                rc = min(p.start for p in self.internal)
                if lc < rc:
                    l3seq = gdict[self.chrom()][lc:rc].reverse_complement().lower() + l3seq

        return l5seq + self.internal_sequence(gdict) + l3seq
    
    def span_gtf(self):
        _allann = self.ltr_up + self.internal + self.ltr_int + self.ltr_down
        _spos = min(v.start for v in _allann)
        _epos = max(v.end for v in _allann)
        _attrs = 'name "%s"; category "%s"; nfeats "%d"; length "%d"; cov "%d";' % (self.id, self.category(), len(_allann), (_epos-_spos), self.model_cov())
        if self.locus_name is not None:
            _attrs += ' locus "%s";' % self.locus_name
        
        mergeline = [self.chrom(), 'merged', 'gene', str(_spos), str(_epos),
                     '.', self.strand(), '.', _attrs,
                    ]
        return mergeline

    def each_gtf(self):
        ret = []
        _allann = self.ltr_up + self.internal + self.ltr_int + self.ltr_down
        _allann.sort(key=lambda x:x.start)
        for g in _allann:
            gfmt = g.fmt()
            gfmt[1] = self.category()
            if self.locus_name is None:
                gfmt[8] = 'gene_id "%s"; transcript_id "%s"; %s' % (self.id, self.id, gfmt[8])
            else:
                gfmt[8] = 'gene_id "%s"; transcript_id "%s"; %s' % (self.locus_name, self.locus_name, gfmt[8])
                gfmt[8] += ' locus "%s";' % self.locus_name
            ret.append(gfmt)
        return ret
    
    def find_ltr(self, files, dist=1000):
        if self.category() != 'internal': return

        m = self.span_gtf()
        if self.strand()=='+':
            up_input   = '\t'.join(flank_left(m,dist))
            down_input = '\t'.join(flank_right(m,dist))
        else:
            up_input   = '\t'.join(flank_right(m,dist))
            down_input = '\t'.join(flank_left(m,dist))
        
        # Bedtools command for intersection
        intersect_cmd = 'bedtools intersect -wo -s -a -  -b %s' % files
        
        # Upstream
        p1 = Popen(intersect_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out,err = p1.communicate(input=up_input)
        if out != '':
            for l in out.strip('\n').split('\n'):
                f = l.split('\t')
                overlap = int(f[-1])
                self.ltr_up.append(GTFLine(f[-10:-1]))

        # Check overlap with internal
        if self.strand() == '+':
            self.ltr_up.sort(key=lambda x:x.start)
            for ltr in self.ltr_up:
                if ltr.end > min(p.start for p in self.internal):
                    ltr.end = min(p.start for p in self.internal)
        else:
            self.ltr_up.sort(key=lambda x:x.end, reverse=True)
            for ltr in self.ltr_up:
                if ltr.start < max(p.end for p in self.internal):
                    ltr.start = max(p.end for p in self.internal)

        # Downstream
        p1 = Popen(intersect_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out,err = p1.communicate(input=down_input)
        if out != '':
            for l in out.strip('\n').split('\n'):
                f = l.split('\t')
                overlap = int(f[-1])
                self.ltr_down.append(GTFLine(f[-10:-1]))

        # Check overlap with internal
        if self.strand() == '+':
            self.ltr_down.sort(key=lambda x:x.start)
            for ltr in self.ltr_down:
                if ltr.start < max(p.end for p in self.internal):
                    ltr.start = max(p.end for p in self.internal)
        else:
            self.ltr_down.sort(key=lambda x:x.end, reverse=True)
            for ltr in self.ltr_down:
                if ltr.end > min(p.start for p in self.internal): 
                    ltr.end = min(p.start for p in self.internal)
        
        # LTRs already found
        found_ltr = [g.attr['id'] for g in self.ltr_up] + [g.attr['id'] for g in self.ltr_down]

        # Internal
        p1 = Popen(intersect_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out,err = p1.communicate(input='\t'.join(m))
        if out != '':
            for l in out.strip('\n').split('\n'):
                f = l.split('\t')
                g = GTFLine(f[-10:-1])
                if g.attr['id'] not in found_ltr:
                    self.ltr_int.append(g)

        if self.strand == '+':
            self.ltr_int.sort(key=lambda x:x.start)
        else:
            self.ltr_int.sort(key=lambda x:x.end, reverse=True)
    
    def __str__(self):
        ret = '\t'.join(self.span_gtf())
        ret += '\n'
        ret += '\n'.join('\t'.join(_) for _ in self.each_gtf())
        return ret
