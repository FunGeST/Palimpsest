#!/usr/bin/python

"""make_spectra_snvs.py [options] simple.file

creates a single mutation counts file for all the samples contained in the input file(s).

options:
--cachedir path     where intermediate files will be saved (default "$PWD/cache/catalogs")
--fasta path     the path to the reference fasta file
--genome | --exome  type of sequencing used for input files
--style dukenus|ludmil  TAB-separated or comma-separated output [ludmil]
--output path       default "mutation_counts-$REGION-$NUM.csv" (region is 
                    "genome" for --genome input, or "exome" for --exome input; 
                    number is the 96, 192, 1536 mutaion types)

Copyright (C) 2018 Steven G. Rozen and Mi Ni Huang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see
    <https://www.gnu.org/licenses/>

Contact: steverozen@gmail.com

"""

import sys, os
from getopt import getopt


script_dir = os.path.dirname(sys.argv[0])
annotation_dir = script_dir + '/annotation'
sys.path.append(annotation_dir)
import lookups

options, paths = getopt(sys.argv[1:], 'hc:', ('classes=', 'cachedir=', 'fasta=', 'output=', 'style=', 'filter-out=', 'dry-run', 'genome', 'exome', 'help'))
options = dict(options)
if '-h' in options or '--help' in options:
    print >>sys.stderr, __doc__
    sys.exit(0)
dry_run = '--dry-run' in options

input_type = None
if options.has_key('--genome'):
    input_type = 'genome'
elif options.has_key('--exome'):
    input_type = 'exome'
else:
    print >>sys.stderr, 'either --genome or --exome must be given.\n'
    print >>sys.stderr, __doc__
    sys.exit(1)
output_filename = options.get('--output')
style = options.get('--style', 'ludmil')
if not output_filename:
    if input_type == 'genome':
        output_filename = 'mutation_counts-genome-%s.csv'
    else:
        output_filename = 'mutation_counts-exome-%s.csv'

cachedir = options.get('--cachedir', 'cache/catalogs')
fasta_refdir = options.get('--fasta')

filters = {} # key => lambda_func
if options.get('--filter-out'):
    # check syntax
    # ENCODE => ENCODE_DAC or ENCODE_DUKE
    clauses = options['--filter-out'].split(',')
    for clause in clauses:
        clause = clause.strip() # whitespace
        k = 0
        while clause[k].isalpha(): k += 1
        key = clause[:k]
        k2 = k
        while not clause[k2].isalnum(): k2 += 1
        op = clause[k:k2].strip() # ignore whitespace
        value = clause[k2:]
        if op in ('<', '>', '<=', '>='):
            value = float(value)
            # put value as default arg and not in function body, because
            # body doesn't get evaluated until it's executed but default
            # args are evaluated when function is defined.
            # func returns True if we get a match (and should filter)
            if op == '<':
                func = lambda v, thresh=value: float(v) < thresh
            elif op == '>':
                func = lambda v, thresh=value: float(v) > thresh
            elif op == '<=':
                func = lambda v, thresh=value: float(v) <= thresh
            elif op == '>=':
                func = lambda v, thresh=value: float(v) >= thresh
        else:
            value_lc = value.lower()
            if value_lc == 'true':
                func = lambda v: True # any value
            elif value_lc == 'false':
                raise Exception('cannot currently test for absence of value')
            else: # assume only integer if testing for equality
                value = int(value)
                func = lambda v, thresh=value: int(v) == thresh
        if key == 'ENCODE': # treat as both
            for key in ('ENCODE_DAC', 'ENCODE_DUKE'):
                filters[key] = func
        else:
            filters[key] = func

if not os.path.isdir(cachedir):
    os.makedirs(cachedir)

lookups.init(fasta_refdir, cachedir)

########
features96, counts96 = None, None
# if genome, these will be filtered to in-transcript only
features192, counts192 = None, None
features1536, counts1536 = None, None # +-2 base context (1536)

for filename in paths:
    features_file96, counts_file96 = lookups.get_class_counts(filename, input_type, 96, variant_filters=filters)
    features_file1536, counts_file1536 = lookups.get_class_counts(filename, input_type, 96*16, variant_filters=filters)

    if input_type == 'genome':
        features_file192, counts_file192 = lookups.get_class_counts(filename, input_type, 192, filter_to="transcript", variant_filters=filters)
    else: # exome
        # assume input file is already restricted to exonic positions... TODO
        features_file192, counts_file192 = lookups.get_class_counts(filename, input_type, 192, variant_filters=filters)

    # merge in genome/exome counts
    if not features96: # first input file
        features96, counts96 = features_file96, counts_file96
    else:
        # need to merge samples from this input file with previous
        if features96 != features_file96: # TODO - just add any missing rows?
            print >>sys.stderr, 'mismatch between features:'
            print >>sys.stderr, features96
            print >>sys.stderr, features_file96
            sys.exit(1)
        for k,v in counts_file96.items():
            # sanity check
            if counts96.has_key(k):
                print >>sys.stderr, 'warning, overwriting existing sample', k
            counts96[k] = v

    if not features192: # first input file
        features192, counts192 = features_file192, counts_file192
    else:
        # need to merge samples from this input file with previous
        if features192 != features_file192: # TODO - just add any missing rows?
            print >>sys.stderr, 'mismatch between features:'
            print >>sys.stderr, features192
            print >>sys.stderr, features_file192
            sys.exit(1)
        for k,v in counts_file192.items():
            # sanity check
            if counts192.has_key(k):
                print >>sys.stderr, 'warning, overwriting existing sample', k
            counts192[k] = v

    if not features1536:
        features1536, counts1536 = features_file1536, counts_file1536
    else:
        if features1536 != features_file1536:
            print >>sys.stderr, 'mismatch between features:'
            print >>sys.stderr, features1536
            print >>sys.stderr, features_file1536
            sys.exit(1)
        for k,v in counts_file1536.items():
            if counts1536.has_key(k):
                print >>sys.stderr, 'warning, overwriting existing sample', k
            counts1536[k] = v

if output_filename.find('%s') > -1:
    filename96 = output_filename % '96'
    filename192 = output_filename % '192'
    filename1536 = output_filename % '1536'
else:
    filename96 = output_filename + '-96'
    filename192 = output_filename + '-192'
    filename1536 = output_filename + '-1536'

if style == 'dukenus':
    lookups.write_as_tsv(features96, counts96, filename96)
    lookups.write_as_tsv(features192, counts192, filename192)
    lookups.write_as_tsv(features1536, counts1536, filename1536)
else:
    lookups.write_as_csv(features96, counts96, filename96)
    lookups.write_as_csv(features192, counts192, filename192)
    lookups.write_as_csv(features1536, counts1536, filename1536)
