#!/usr/bin/python

"""make_spectra_dinucs.py [options] simple.file

creates a single 'dinucleotide' classification counts file for all the samples
contained in the input file.

options:
--output path       default "dinucs_counts.txt"

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
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

Contact: steverozen@gmail.com

"""

import getopt
import os
import sys
import numpy as np
import pandas as pd
from itertools import product
import functools
import gzip



args, filenames = getopt.getopt(sys.argv[1:], 'hc:', ('output=', 'help'))

args = dict(args)

if '-h' in args or '--help' in args or not filenames:
    # no args given or help requested
    print >>sys.stderr, __doc__
    sys.exit(0)

output_filename = args.get('--output')
if not output_filename:
    output_filename = 'dinucs_counts-%s.csv'

if output_filename.find('%s') > -1:
    output_filename = output_filename % 'no-strand'
else:
    output_filename = output_filename + '-dinuc'




pair = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
bases = ['A','C','G','T']


############# INDEX 144
## 16 dinucleotide type
combo = list(product(bases, repeat=2))
comboName = [''.join(x) for x in combo]

var_combo = list()
for dinc in combo:
    sub_var_combo = list(product([x for x in bases if x!=dinc[0]], [x for x in bases if x!=dinc[1]]))
    var_combo.append(sub_var_combo)

var_comboName = [''.join(x) for sublist in var_combo for x in sublist]

## index for dataframe
tuples = zip([x for x in comboName for i in range(9)], var_comboName)
index_144 = pd.MultiIndex.from_tuples(tuples, names=['Ref', 'Var'])


############# INDEX 78
## for strand analysis, there is only 4 dinucleotide type and their complements
ref_comb = ['AC', 'AT', 'CC', 'CG', 'CT', 'GC', 'TA', 'TC', 'TG', 'TT']
ref_comb_complement = [pair[x[1]] + pair[x[0]] for x in ref_comb]

var_comb_tuple = list()
for dinc in ref_comb:
    sub_var_comb = list(product([x for x in bases if x!=dinc[0]], [x for x in bases if x!=dinc[1]]))
    var_comb_tuple.append(sub_var_comb)

var_comb = [''.join(x) for sublist in var_comb_tuple for x in sublist]
var_comb_complement = [pair[x[1]] + pair[x[0]] for x in var_comb]

## index for dataframe
tuples = zip([x for x in ref_comb for i in range(9)], var_comb)
tuples_complement = zip([x for x in ref_comb_complement for i in range(9)], var_comb_complement)
## complementary dinucs are right next 
tuple_index = [item for x in zip(tuples, tuples_complement) for item in x]
## remove duplicates
unique = functools.reduce(lambda l, x: l+[x] if x not in l else l, tuple_index, [])
## remove symetrical ones
sym_list = [('AT', 'CG'), ('AT', 'GC'), ('AT', 'TA'), ('CG', 'AT'), ('CG', 'GC'), \
('CG', 'TA'), ('GC', 'AT'), ('GC', 'CG'), ('GC', 'TA'), ('TA', 'AT'), ('TA', 'CG'), ('TA', 'GC')]
tuple_index = [x for x in unique if x not in sym_list]

tuple_index = tuple_index + sym_list
index_78 = pd.MultiIndex.from_tuples(tuple_index, names=['Ref', 'Var'])


#### get counts for 144 types of DNPs
def get_dinucleotide_counts(filename):
    if filename.endswith('.gz'):
        handle = gzip.open(filename)
    else:
        handle = open(filename)    
        
    disease_idx = 0 # cancer type
    tumour_idx = 1
    vartype_idx = 4
    chrome_idx = 5
    pos_idx = 6
    refbase_idx = 8
    varbase_idx = 9

    line_before = ''
    chr_pos_before = tuple(['chr', 0])
    ref_before = ''
    var_before = ''
    
    df_counts = pd.DataFrame(index=index_144)
    
    for line in handle:
        line_value = line.strip().split('\t')
        cancer_type = line_value[disease_idx]
        uniq_id = line_value[tumour_idx]
        sample_id = (uniq_id)

        if sample_id not in df_counts.columns:
            ## All 144 cobinations, start with 0
            counts = [0]*144
            df_counts[sample_id] = pd.Series(counts, index=index_144)

        vartype = line_value[vartype_idx]
        # 'SNP' used by TCGA, 'single base substitution' used by ICGC
        # 'SNV' used in CCA lit. paper
        if vartype not in ('SNP', 'SNV', 'single base substitution'):
            continue # ignore indels for now

        chrome = str(line_value[chrome_idx])
        if chrome[0:3] != 'chr':
            chrome = 'chr' + chrome
        pos = int(float(line_value[pos_idx]))
        chr_pos = tuple([chrome, pos])
        ref = str(line_value[refbase_idx])
        var = str(line_value[varbase_idx])
       
        if chr_pos_before[0] == chr_pos[0] and chr_pos_before[1]+1 == chr_pos[1]: 
            idx_ref_dinc = combo.index((ref_before, ref))
            idx_var_dinc = var_combo[idx_ref_dinc].index((var_before, var))
            df_idx = idx_ref_dinc*9 + idx_var_dinc
            df_counts[sample_id][df_idx] += 1

        line_before = line_value
        chr_pos_before = chr_pos
        ref_before = ref
        var_before = var
    
    handle.close()
    
    return (df_counts)



#### get counts for 144 types of DNPs
df = pd.DataFrame(index=index_144)
for filename in filenames:
    df_counts = get_dinucleotide_counts(filename)
    df = pd.concat([df, df_counts], axis=1)

#### get counts for 78 folded types of DNPs
df = df.reindex(index_78)
df1 = df[:132]
df2 = df[132:]
df11 = df1[::2]
df12 = df1[1::2]
df12.index = df11.index
df1new = df11 + df12
df_out = df1new.append(df2)
#df_out = df_out.sort()
df_out = df_out.sort_index() ### later version of pandas


df_out.to_csv(output_filename)

exit(0)
