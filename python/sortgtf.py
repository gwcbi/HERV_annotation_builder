#! /usr/bin/env python

import utils

def main(parser):
    args = parser.parse_args()
    lines = utils.tab_line_gen(args.infile)
    for l in utils.sort_gtf(lines):
        print >>args.outfile, '\t'.join(l)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='sort GTF')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    main(parser)
