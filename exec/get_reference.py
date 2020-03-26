#!/usr/bin/python

"""
Given a chromosome and position (or range), seek to the correct place in
the named reference .fasta file and extract the correct sequence.
"""

import mmap
import os
import sys


class reference():

    def __init__(self, fasta_filename):
        self._fastafilename = ''
        self._fd = None
        self._index = {} # chr -> (first byte, length, linelength, offset_to_page)
        self._mmapped_chrome = {}
        self._is_gzipped = False
        self._fastafilename = fasta_filename
        self._read_fasta_index(fasta_filename + '.fai')


    def _read_fasta_index(self, filename):
        for line in open(filename):
            chrome, length, start, linelength, _ = line.split()
            length, start, linelength = int(length), int(start), int(linelength)
            aligned_start = int(start / mmap.PAGESIZE) * mmap.PAGESIZE
            offset = start - aligned_start
            self._index[chrome] = (start, length, linelength, offset)


    def get_sequence(self, chrome, start, end=None):
        """eg chrome='chr4', start and end should be integers.
            For backwards compatibility, we also accept start='12345'
            or start='12345-67890'. """
        if type(start) == str:
            if start.count('-'):
                # we've been given a range
                start, end = start.split('-')
            else:
                start, end = start, start
            start, end = int(start), int(end)
        if not end: end = start
        if end < start:
            tmp = start
            start = end
            end = tmp

        try:
            chr_start, chr_num_bases, chr_linelength, offset = self._index[chrome]
        except KeyError:
            print >>sys.stderr, ' '
            print >>sys.stderr, '**ERROR** incompatibility detected between input VCF and FASTA file, likely caused by no "chr" prefixes in %s' % self._fastafilename
            print >>sys.stderr, 'try setting "genome_build" arg in annotate_VCF() to "GRCh37'
            print >>sys.stderr, ' '
            print >>sys.stderr, 'missing contig in %s and/or index?' % self._fastafilename
            raise

        try:
            buf = self._mmapped_chrome[chrome]
        except KeyError:
            if not self._fd:
                self._fd = os.open(self._fastafilename, os.O_RDONLY)
            # have to align on PAGESIZE bytes, so round down to nearest multiple
            aligned_chr_start = int(chr_start / mmap.PAGESIZE) * mmap.PAGESIZE
            # account for newlines
            chr_length = chr_num_bases + (chr_num_bases / chr_linelength)
            buf = mmap.mmap(self._fd, chr_length+offset, flags=mmap.MAP_PRIVATE, prot=mmap.PROT_READ, offset=aligned_chr_start)
            self._mmapped_chrome[chrome] = buf

        # number of newlines between chr's start byte and our seq start byte
        num_newlines = int((start-1) / chr_linelength)
        # move to start of data. -1 since buf indexed from 0, but chrome:pos
        # indexed from 1
        # offset = byte where chromosome starts in the .fa file
        seq_to = offset + start-1 + num_newlines
        seq_length = end-(start-1) # +1 because we want to include final base
        try:
            buf.seek(seq_to, os.SEEK_SET)
        except ValueError, v: # customise error message
            v.message = 'get_reference error seeking .fasta; requested %s:%s-%s, tried seeking to %d' % (chrome, start, end, seq_to)
            raise v
        if self._is_gzipped:
            print "TODO!"
        else:
            # adjust for newlines in fasta file
            num_newlines = ((start-1) % chr_linelength + seq_length) / chr_linelength
            seq_length += int(num_newlines)
            seq = buf.read(seq_length)
        #buf.close()
        #os.close(fd)
        return seq.replace('\n', '').upper()

    def _test(self):
        for expected, params in (
            # hg18:
    #        ('T', ('chr1', '49')),
    #        ('ACCTATCAGCAGGATGTGGGTGGGAGCAGATTAGAGAATAAAAGCAGACTGC', ('chrY', '49-100')),
    #        ('CCTATCAGCAGGATGTGGGTGGGAGCAGATTAGAGAATAAAAGCAGACTGC', ('chrY', '50-100')),
    #        ('CTATCAGCAGGATGTGGGTGGGAGCAGATTAGAGAATAAAAGCAGACTGC', ('chrY', '51-100')),
            # hg19:
            ('TGAAGGCTTTTTTGGTTGCTTCTAAATTTGAGCTGTAAATATTCAGTGCATGTCTCCTGG', ('chr17', '60000-60059')),
            ('TGAAGGCTTTTTTGGTTGCTTCTAAATTTGAGCTGTAAATATTCAGTGCATGTCTCCTGGA', ('chr17', '60000-60060')),
            ('TGAAGGCTTTTTTGGTTGCTTCTAAATTTGAGCTGTAAATATTCAGTGCATGTCTCCTGGAG', ('chr17', '60000-60061')),
        ):
            value = get_sequence(*params)
            print params, expected==value
            print ' expected',expected,'got',value

if __name__ == '__main__':
    fa_path = '/data/public/reference_genomes/hg19/hg19lite_ucsc_chrM.fa'
    if len(sys.argv) == 1:
        print '%s [-f fastafile] [position [...]]' % sys.argv[0]
        print 'position is in "chr:pos1[-pos2]" format'
        print 'if no positions are given on command line, they are read from stdin'
        print '\neg: %s %s chr12:25,249,447-25,295,121' % (sys.argv[0], fa_path)
        sys.exit(1)
    positions = []
    if len(sys.argv) > 1:
        if sys.argv[1] == '-f':
            fa_path = sys.argv[2]
            positions = sys.argv[3:]
        else:
            positions = sys.argv[1:]
    if not positions:
        for line in sys.stdin:
            positions.append(line.rstrip())

    if fa_path.startswith('~/'):
        fa_path = os.environ['HOME'] + fa_path[1:]
    reference_obj = reference(fa_path)
    #_test() ; sys.exit(0)

    positions = map(lambda p: tuple(p.replace(',', '').split(':', 1)), positions)
    # positions is now a list of (chr, pos[-pos]) tuples
    
    print positions
    for pos in positions:
        print reference_obj.get_sequence(*pos)
