#!/usr/bin/python

"""take chr:pos points from stdin, and list if they are in target or not."""

import bisect
import subprocess
import sys
import os

class ranges:
    "this class loads and indexes a .bed file, allowing lookup by chr & pos"

    def __init__(self, targetfile):
        "loads the given .bed file."
        target_starts = {}
        target_ends = {}
        target_range_name = {} # CCDS/gene name etc

        # note - chrM will be filtered out if it isn't in the target file
        for i in range(1, 23) + ['X', 'Y', 'M']:
            c = "chr%s" % i
            i = str(i)
            target_starts[c] = []       # UCSC-style chr names
            target_ends[c] = []
            target_range_name[c] = []
            if i == 'M': i = 'MT'
            target_starts[i] = []       # NIH/GRC-style chr names for h37
            target_ends[i] = []
            target_range_name[i] = []
        self.target_filename = targetfile
        if targetfile.endswith('.gz'):
            import gzip
            IN = gzip.open(targetfile)
            targetfile = targetfile[:-3] # for gff/bed detection
        else:
            IN = open(targetfile)
        lastchr = None
        laststart = 0
        unknown_chr_warned = {}
        if targetfile.endswith('.gff') or targetfile.endswith('.gff3'):
            start0based = 0
        elif targetfile.endswith('.bed'):
            start0based = 1
        else:
            print >>sys.stderr, 'assuming .bed format with 0-based start coords.'
            start0based = 1
        for line in IN:
            if line.startswith('browser ') or line.startswith('track '):
                continue # .bed format header
            chrome, start, end, name = line.rstrip('\r\n').split('\t', 3)
            start = int(start) + start0based
            if chrome == lastchr and start < laststart:
                print >>sys.stderr, 'Error: target file is not sorted:', line
                sys.exit(1)
            try:
                target_starts[chrome].append(start)
                target_ends[chrome].append(int(end))
                target_range_name[chrome].append(name)
            except KeyError: # unknown chromosome
                if not unknown_chr_warned.has_key(chrome):
                    print >>sys.stderr, 'ignoring unmatched chromosome %s in %s' % (chrome, targetfile)
                    unknown_chr_warned[chrome] = True
            lastchr = chrome
            laststart = start
        IN.close()

        self.target_starts = target_starts
        self.target_ends = target_ends
        self.target_range_name = target_range_name


    def get_target_name(self, chrome, pos=None):
        """search through list of target positions on this chromosome.
            Either chromosome and position are given separately as 2 arguments,
            or provide a single string argument like "chr2:123456".
            returns the name of the matching target region, or None if not
            in target."""
        if pos==None:
            # assume we were given a single arg like "chr1:1234"
            chrome, pos = chrome.split(':')
        pos = int(pos)

        try:
            starts = self.target_starts[chrome]
        except KeyError:
            return None # unknown chromosome in target. eg chrUn_gl000229

        if not starts:
            # eg chrM is not in target list, so no target positions
            return None
        index = bisect.bisect_left(starts, pos)
        # the read's start position is before this target range[index]
        # special case 1:
        if index == len(starts):
            # the position is after the last range's start. Return name if it
            # is inside (ie before the end of this target range).
            if pos <= self.target_ends[chrome][index-1]:
                return self.target_range_name[chrome][index-1]
            return None
        # if here, then the position occurs before the start of this target
        # range. now check the previous target range
        # special case 2: we're in the first range
        if index == 0:
            if pos >= self.target_starts[chrome][0] and pos <= self.target_ends[chrome][0]:
                return self.target_range_name[chrome][index]
            return None
        if pos <= self.target_ends[chrome][index-1]:
            # read starts before the end of the previous target range
            return self.target_range_name[chrome][index-1]
        # see if we are at the very start of the current target range
        if pos == self.target_starts[chrome][index]:
            return self.target_range_name[chrome][index]

        # if here, the position is before the current range, but after the
        # previous range. IE between ranges:
        # search all the ranges before -- add by Mini
        while index > 0:
            if pos <= self.target_ends[chrome][index-1]:
                return self.target_range_name[chrome][index-1]
                break
            index -= 1

        return None

    def get_all_targets(self, chrome, pos):
        """search through list of target positions on this chromosome.
            Either chromosome and position are given separately as 2 arguments,
            or provide a single string argument like "chr2:123456".
            returns a list of matching target region names, or None if not
            in target."""
        matches = []
        try:
            starts = self.target_starts[chrome]
        except KeyError:
            return None # unknown chromosome in target. eg chrUn_gl000229

        if not starts:
            # eg chrM is not in target list, so no target positions
            return None
        # _left: "all e in a[:i] have e < x, and all e in a[i:] have e >= x"
        index = bisect.bisect_left(starts, pos)
        # the read's start position is before this target range[index]
        # special case 1:
        if index == len(starts):
            # the position is after the last range's start. Return name if it
            # is inside (ie before the end of this target range).
            if pos <= self.target_ends[chrome][index-1]:
                return [self.target_range_name[chrome][index-1]]
            return None
        # special case 2: we're in the first range
        if index == 0:
            if pos < starts[0]:
                return None # pos is before first range
            # if here, pos== start of range
            matches.append(self.target_range_name[chrome][index])
        else: # if we are at first position of a range
            if pos >= starts[index] and pos <= self.target_ends[chrome][index]:
                matches.append(self.target_range_name[chrome][index])
        # check current and all leftwards ranges
        index -= 1
        # can't assume all containing ranges will be next to each other,
        # because only sorted by start pos, not end pos. go through *all*
        # potentially overlapping ranges, for now.
        while index >= 0:
            if pos <= self.target_ends[chrome][index]:
                matches.append(self.target_range_name[chrome][index])
            index -= 1
        return matches
