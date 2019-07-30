#!/usr/bin/env python
"""comp.py infile outfile

A script to complement fastq files
"""
from Bio import SeqIO
import gzip
import sys

def comp_writer(infile, outfile):
    """Write the complement of a fastq File

    """

    def _comp_map(rec):
        rec.seq = rec.seq.complement()
        #comp.id = rec.id
        #comp.name = rec.name
        #comp.description = rec.description
        return rec

    def _writer(ihandle, outhandle):
        recs = SeqIO.parse(ihandle, 'fastq')
        comprecs = map(_comp_map, recs)
        SeqIO.write(comprecs, ohandle,'fastq')


    if infile.endswith('.gz'):
        with gzip.open(infile, 'rt') as ihandle:
            with gzip.open(outfile, 'wt') as ohandle:
                _writer(ihandle, ohandle)
    else:
      with open(infile, 'r') as ihandle:
          with gzip.open(outfile, 'wt') as ohandle:
              _writer(ihandle, ohandle)

def main():
    comp_writer(infile=sys.argv[1], outfile=sys.argv[2])

if __name__ == '__main__':
    main()
