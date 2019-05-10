#!/usr/bin/python

"""make_spectra_indels.py [options] simple.file

creates a single 'indel' classification counts file for all the samples
contained in the input file.

options:
--cachedir path     where intermediate files will be saved (default "$PWD/cache/catalogs")
--fasta path     the path to the reference fasta file
--genome | --exome  type of sequencing used for input files
--style dukenus|ludmil  TAB-separated or comma-separated output [ludmil]
--output path       default "indel_counts.csv" 

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
import lookups_indels

options, paths = getopt(sys.argv[1:], 'hc:', ('classes=', 'cachedir=', 'fasta=', 'output=', 'style=', 'filter-out=', 'dry-run', 'genome', 'exome', 'help'))
options = dict(options)
if '-h' in options or '--help' in options or not paths:
    print >>sys.stderr, __doc__
    sys.exit(0)
dry_run = '--dry-run' in options

input_type = None
if options.has_key('--genome'):
    input_type = 'genome'
elif options.has_key('--exome'):
    input_type = 'exome'
else:
    #print >>sys.stderr, 'either --genome or --exome must be given.'
    #sys.exit(1)
    input_type = 'genome' # currently makes no difference for indels... TODO

output_filename = options.get('--output')
style = options.get('--style', 'ludmil')
if not output_filename:
    output_filename = 'indel_counts%s.csv'

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

lookups_indels.init(fasta_refdir, cachedir)

########
features_plus_strand, counts_plus_strand = None, None
## if genome, these will be filtered to in-transcript only
#features_sense_strand, counts_sense_strand = None, None

for filename in paths:
    file_features_plus, file_counts_plus = lookups_indels.get_indel_class_counts(filename, input_type, '+', variant_filters=filters)

    # don't worry about sense strand for indels, for now.
    #if input_type == 'genome':
    #    file_features_sense, file_counts_sense = lookups.get_indel_class_counts(filename, input_type, 'sense', filter_to="transcript", variant_filters=filters)
    #else: # exome
    #    # assume input file is already restricted to exonic positions... TODO
    #    file_features_sense, file_counts_sense = lookups.get_indel_class_counts(filename, input_type, 'sense', variant_filters=filters)

    # merge in genome/exome counts
    if not features_plus_strand: # first input file
        features_plus_strand, counts_plus_strand = file_features_plus, file_counts_plus
    else:
        # need to merge samples from this input file with previous
        if features_plus_strand != file_features_plus: # TODO - just add any missing rows?
            print >>sys.stderr, 'mismatch between features:'
            print >>sys.stderr, features_plus_strand
            print >>sys.stderr, file_features_plus
            sys.exit(1)
        for k,v in file_counts_plus.items():
            # sanity check
            if counts_plus_strand.has_key(k):
                print >>sys.stderr, 'warning, overwriting existing sample', k
            counts_plus_strand[k] = v
    """
    if not features_sense_strand: # first input file
        features_sense_strand, counts_sense_strand = file_features_sense, file_counts_sense
    else:
        # need to merge samples from this input file with previous
        if features_sense_strand != file_features_sense: # TODO - just add any missing rows?
            print >>sys.stderr, 'mismatch between features:'
            print >>sys.stderr, features_sense_strand
            print >>sys.stderr, file_features_sense
            sys.exit(1)
        for k,v in file_counts_sense.items():
            # sanity check
            if counts_sense_strand.has_key(k):
                print >>sys.stderr, 'warning, overwriting existing sample', k
            counts_sense_strand[k] = v
    """
## ignore strand for now
if output_filename.find('%s') > -1:
    filename_plus_strand = output_filename % ''
    #filename_plus_strand = output_filename % '+strand'
    #filename_sense_strand = output_filename % 'sense.strand'
else:
    filename_plus_strand = output_filename + '-indel'
    #filename_plus_strand = output_filename + '_+strand'
    #filename_sense_strand = output_filename + '_sense.strand'

if style == 'dukenus':
    lookups_indels.write_indels_as_tsv(features_plus_strand, counts_plus_strand, filename_plus_strand)
    #lookups.write_indels_as_tsv(features_sense_strand, counts_sense_strand, filename_sense_strand)
else:
    lookups_indels.write_indels_as_csv(features_plus_strand, counts_plus_strand, filename_plus_strand)
    #lookups.write_indels_as_csv(features_sense_strand, counts_sense_strand, filename_sense_strand)
