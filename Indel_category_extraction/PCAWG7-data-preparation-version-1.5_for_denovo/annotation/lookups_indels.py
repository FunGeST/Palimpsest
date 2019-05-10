"""
module for looking up reference .fasta files and gene .bed files.
"""

import json, os, sys
import subprocess

import get_reference
import in_target

_beddir = os.path.dirname(sys.argv[0])
_beddir = _beddir + '/annotation'
# _beddir = '/gpfs/nextgen1/pipeline/targets'
grch37_bed_file = _beddir + '/GENCODE_transcripts.v19.hs37.bed'
grch38_bed_file = _beddir + '/GENCODE_transcripts.v20.grch38.bed'

# 0=errors only, 1=errors+info, 2=verbose
debug_level = 0

#######

_cachedir = './cache'

_hg18_reference = None
_hg19_reference = None   # UCSC ref with leading "chr" prefix and older chrM
_grch37_reference = None # no leading "chr" prefix and newer MT
_grch38_reference = None


_complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
_de_IUB = {
    ('R','A'): 'G',
    ('R','G'): 'A',
    ('Y','C'): 'T',
    ('Y','T'): 'C',
    ('K','T'): 'G',
    ('K','G'): 'T',
    ('M','A'): 'C',
    ('M','C'): 'A',
    ('W','A'): 'T',
    ('W','T'): 'A',
    ('S','G'): 'C',
    ('S','C'): 'G',
}

def init(fasta_search_paths, cachedir=None):
    """load the .fasta.fai files for reference genomes found in the
        'fasta_search_paths' sequence. Entries in this sequence may be either
        paths to .fasta files, or paths to folders to be searched.

        Eg: init(['/gpfs/public/reference_genomes', '/home/ref/hg19.fasta'])
    """
    global _hg18_reference, _hg19_reference, _grch37_reference, _grch38_reference, _cachedir
    if cachedir: _cachedir = cachedir # use caller's preference
    if type(fasta_search_paths) == str: # if only single str was given
        fasta_search_paths = [fasta_search_paths]
    for entry in fasta_search_paths:
        if os.path.isdir(entry):
            fasta_search_paths += [entry + '/' + e for e in os.listdir(entry)]
            continue
        if not entry.endswith('.fasta') and not entry.endswith('.fa'):
            continue
        n = os.path.basename(entry).lower()
        if not _hg18_reference and n.count('hg18'):
            _hg18_reference = get_reference.reference(entry)
        elif not _hg19_reference and n.count('hg19'):
            _hg19_reference = get_reference.reference(entry)
        elif not _grch37_reference and (n.count('grch37') or n.count('hs37') or n.count('assembly19')):
            _grch37_reference = get_reference.reference(entry)
        elif not _grch38_reference and (n.count('h38') or n.count('hg38')):
            _grch38_reference = get_reference.reference(entry)
	# add demo reference for test run: chr22 in TCGA_LAML data
        elif not _grch37_reference and n.count('demo'):
            _grch37_reference = get_reference.reference(entry)


def get_reference_seq(ref_genome, chromosome, start, end, filepath=''):
    """gets the sequence from 'start' to 'end' (inclusive).
        'filepath' is only used in debugging output, if set."""
    sequence = None
    try:
        # '37' = TCGA, 'GRCh37' = ICGC
        if ref_genome in ('37', 'hs37d5', 'GRCh37', 'GRCh37-lite'):
            if chromosome == 'M': chromosome = 'MT' # match our hs37 ref
            sequence = _grch37_reference.get_sequence(chromosome, start, end)
        elif ref_genome in ('hg19'):
            # shouldn't happen, but some .maf file entries (eg from
            # bcgsc.ca) say hg19 but use GRCh37 chr names...
            if chromosome.startswith('GL'):
                which = int(chromosome[5:8])
                mapping = {191:1, 192:1, 193:4, 194:4, 195:7, 196:8, 197:8, 198:9, 199:9, 200:9, 201:9, 202:11, 203:17, 204:17, 205:17, 206:17, 207:18, 208:19, 209:19, 210:21}
                if mapping.has_key(which):
                    chromosome = 'chr%d_gl%s_random' % (mapping[which], chromosome[2:8])
                else:
                    # skip leading GL, truncate possible ".1" trailer
                    chromosome = 'chrUn_gl%s' % chromosome[2:8]
            if not chromosome.startswith('chr'):
                if chromosome == 'MT': chromosome = 'M'
                chromosome = 'chr' + chromosome
            sequence = _hg19_reference.get_sequence(chromosome, start, end)
        elif ref_genome in ('GRCh38', 'grch38'):
            # 1000genomes GRCh38 ref has leading 'chr' prefix, so does
            # UCSC hg38 ref.
            if not chromosome.startswith('chr'):
                chromosome = 'chr' + chromosome
            sequence = _grch38_reference.get_sequence(chromosome, start, end)
        elif ref_genome in ('36', 'hg18'):
            # we don't have GRCh36 - map to hg18 names instead
            if not chromosome.startswith('chr'):
                chromosome = 'chr' + chromosome
            if chromosome in ('chrM','chrMT'):
                return # silently ignore, not in our reference
            sequence = _hg18_reference.get_sequence(chromosome, start, end)
        else:
            print >>sys.stderr, (ref_genome, 'is not hg18/GRCh36 or hg19/GRCh37 ??')
    except ValueError: # invalid position for chromosome?
        raise # serious error, calling script needs to handle this
    except KeyError: # chromosome not present in ref?
        print >>sys.stderr, ('missing chr for context sequence at %s %s:%d-%d in %s' % (ref_genome, chromosome, start, end, filepath))
        # continue
    return sequence

_context_cache = {}
_context_cache_updated = False

def get_context(ref_genome, chromosome, pos, n=1, filepath=''):
    """get the +-n base reference sequence around this position.
        'filepath' is only used in debugging output, if set."""
    global _context_cache_updated
    try:
        # assumes we don't ask for a different 'n' !!
        return _context_cache[ref_genome][chromosome][pos]
    except KeyError: # this chromosome/position not currently in cache
        pass
    sequence = get_reference_seq(ref_genome, chromosome, pos-n, pos+n, filepath)
    if _context_cache.has_key(ref_genome):
        if not _context_cache[ref_genome].has_key(chromosome):
            _context_cache[ref_genome][chromosome] = {}
    else:
         _context_cache[ref_genome] = {chromosome: {}}
    _context_cache[ref_genome][chromosome][pos] = sequence
    _context_cache_updated = True
    return sequence


class list_i(list):
    "indexable class where failure returns -1 instead of raising ValueError."
    def __init__(self, l):
        self.extend(l)
    def get(self, v):
        try:
            return self.index(v)
        except ValueError:
            return -1

# list of target coordinates for transcripts, 1 dict for each reference
_transcripts_for_genome = {} # eg {'hg19': {'chr1':{...}, 'chr2':{...}}

def annotate_simple_indel(input_path, output_path):
    """calculate annotations (repeats, microhomology etc)
        and save info to file. The simple file format is 1-based for TCGA, need to standardize for the rest..."""

    in_proc = None
    if input_path.endswith('.gz'):
        in_proc = subprocess.Popen(('gzip', '-dc', input_path), stdout=subprocess.PIPE)
        in_handle = in_proc.stdout
    else:
        in_handle = open(input_path)

    # simple file has 'disease id cohort assembly snp/indel chr start end ref var annotations'
    disease_idx = 0 # cancer type
    tumour_idx = 1
    assembly_idx = 3
    vartype_idx = 4
    chrome_idx = 5
    pos_idx = 6
    endpos_idx = 7
    refbase_idx = 8
    varbase_idx = 9
    annotations_idx = 10
    strand_idx = -1

    blacklist = set() # ignore samples with invalid coordinates
    sample_ids_seen = set()
    # in case we have duplicate rows (eg multiple transcripts for a variant)
    patient_coordinates_seen = set()
    ignored_complex = 0
    ranges, last_genome = None, None # if checking coords are in a transcript
    out_proc = None
    if output_path.endswith('.gz'):
        out_handle_real = open(output_path + '-tmp', 'w')
        out_proc = subprocess.Popen(('gzip', '-c'), stdin=subprocess.PIPE, stdout=out_handle_real)
        out_handle = out_proc.stdin
    else:
        out_handle = open(output_path + '-tmp', 'w')
    for line in in_handle:
        fields = line.rstrip('\n').split('\t')
        cancer_type = fields[disease_idx]
        # for ICGC, this is a tissue sample ID, not really a patient ID
        uniq_id = fields[tumour_idx]
        sample_id = (uniq_id)
        sample_ids_seen.add(sample_id) # ignores duplicates
        if sample_id in blacklist: continue # ignore all further mutations
        vartype = fields[vartype_idx]
        if vartype not in ('INS', 'DEL', 'INDEL'):
            continue
        chromosome = fields[chrome_idx]
        # first affected base for del, last unaffected base 5' for ins
        start_pos = int(fields[pos_idx])
        # last affected base for del, first unaffected base 3' for ins
        end_pos = int(fields[endpos_idx])

        if (sample_id, chromosome, start_pos) in patient_coordinates_seen:
            if debug_level:
                print >>sys.stderr, 'ignoring duplicate indel in %s (%s) at %s:%d' % (sample_id, input_path, chromosome, start_pos)
            continue
        patient_coordinates_seen.add((sample_id, chromosome, start_pos))

        try:
            ref = fields[refbase_idx]
            var = fields[varbase_idx]
            if strand_idx > -1:
                strand = fields[strand_idx]
            else:
                strand = '' # we will need to calculate it
            ref_genome = fields[assembly_idx]
        except IndexError:
            print >>sys.stderr, ('in', path)
            print >>sys.stderr, ('malformed line?', line)
            #raise
            continue
        try:
            annotations = fields[annotations_idx]
        except IndexError: # some simple files don't have this column?
            annotations = ''

        # get transcript
        if ref_genome in ('hg19', 'hs37d5', '37', 'GRCh37', 'GRCh37-lite'):
            lookup_ref_genome = 'hs37'
            if ref_genome == 'hg19':
                if chromosome.startswith('chr'):
                    test_chromosome = chromosome[4:]
                else:
                    test_chromosome = chromosome
                if test_chromosome == 'M': test_chromosome = 'MT'
            else: # assume no leading 'chr' prefix
                test_chromosome = chromosome
        elif ref_genome in ('GRCh38', 'hs38'):
            lookup_ref_genome = 'grch38'
            # should all start with 'chr'
            if chromosome.startswith('chr'):
                test_chromosome = chromosome
            else:
                test_chromosome = 'chr' + chromosome
            # shouldn't happen, but does... :(
            if test_chromosome == 'chr23': test_chromosome = 'chrX'
            elif test_chromosome == 'chr24': test_chromosome = 'chrY'
        else: # we don't actually support any other references (yet)
            # TODO
            lookup_ref_genome = ref_genome
            test_chromosome = chromosome
        if last_genome != lookup_ref_genome:
            # need to get correct dict of coordinates for transcripts
            last_genome = lookup_ref_genome
            transcript_ranges = _transcripts_for_genome.get(lookup_ref_genome)
            if not transcript_ranges and lookup_ref_genome in ('hs37', 'grch38'):
                # not initialised; load ranges from a .bed file now
                if lookup_ref_genome == 'hs37':
                    bed_file = grch37_bed_file
                elif lookup_ref_genome == 'grch38':
                    bed_file = grch38_bed_file
                transcript_ranges = in_target.ranges(bed_file)
                _transcripts_for_genome[lookup_ref_genome] = transcript_ranges
        range_descr = transcript_ranges.get_target_name(test_chromosome, start_pos)
        # TODO - if not in a target, check the end position?
        if range_descr: # first found transcript
            # columns 4,5,6 of our .bed files
            # look for all transcripts, in case on both strands
            range_descrs = transcript_ranges.get_all_targets(test_chromosome, start_pos)
            try:
                # can't use zip directly because we give it a sequence,
                # not the contents of the sequence
                genes, transcripts, strands = apply(zip,[r.split('\t', 2) for r in range_descrs])
            except ValueError:
                print >>sys.stderr, lookup_ref_genome, `test_chromosome`, `start_pos`, `range_descr`
                print >>sys.stderr, range_descrs
                print >>sys.stderr, transcript_ranges.target_filename, len(transcript_ranges.target_starts[test_chromosome])
                raise
            # group by strand
            indices = range(len(strands))
            indices.sort(lambda a,b:cmp(strands[a],strands[b]))
            gene = ','.join(genes[i] for i in indices)
            transcript = ','.join(transcripts[i] for i in indices)
            strands = [strands[i] for i in indices]
            # collapse to either '-', '+', '+,-'
            if strands[0] == strands[-1]: # all on same strand
                strand = strands[0]
            else:
                strand = strands[0] + ',' + strands[-1]
        else:
            gene, transcript, strand = '', '', ''

        # get the context sequence around this variant, to look for indel classification:
        ## vartype * indel_size * indel_type * i

        if vartype == 'DEL':
            if var != '-': ## complex indels??
                print >>sys.stderr, 'ignore complex indel: %s, %s, %s, %s, %s|%s|%s' % (chromosome, start_pos, end_pos, ref, LHS_seq, var, RHS_seq)
                continue

            # note - some simple files don't include entire deleted sequence,
            # eg "GRCh37 INDEL 3 83299029 83299050 GAGAG<12>GAAGA  -"
            # signifies 12 bases aren't shown. get from reference now.
            if ref.find('<') > -1:
                try:
                    ref = get_reference_seq(ref_genome, chromosome, start_pos, end_pos-start_pos+1)
                except ValueError:
                    print >>sys.stderr, '%s (%s) now blacklisted due to invalid position %s %s:%s' % (sample_id, input_path, ref_genome, chromosome, start_pos)
                    blacklist.add(sample_id)
                    continue

            del_length = len(ref)
            ## indel_size : 1,2,3,4,5+
            if del_length <5: indel_size = str(del_length)
            else: indel_size = '5+' 

            try:    ## check for right and left side of 5 times long
                LHS_seq = get_reference_seq(ref_genome, test_chromosome, start_pos-5*del_length, start_pos-1)
                RHS_seq = get_reference_seq(ref_genome, test_chromosome, end_pos+1, end_pos+5*del_length)
            except ValueError:
                # blacklist this sample
                print >>sys.stderr, '%s (%s) now blacklisted due to invalid position %s %s:%s' % (sample_id, input_path, ref_genome, chromosome, start_pos)
                blacklist.add(sample_id)
                continue

            ## classification for 1bp del: indel_type is C or T; strand is T (-), N (+), U(untranscribed/both direction)
            if del_length == 1:
                if ref in ['A','G']: 
                    indel_type = _complement[ref]
                    # if strand == '+': indel_type = indel_type + '_T'
                    # elif strand == '-': indel_type = indel_type + '_N'
                    # else: indel_type = indel_type + '_U'
                else: 
                    indel_type = ref
                    # if strand == '+': indel_type = indel_type + '_N'
                    # elif strand == '-': indel_type = indel_type + '_T'
                    # else: indel_type = indel_type + '_U'

                ## get the repeat size i
                rhs_i = 0 # bases 3' of deletion that match the deletion
                lhs_i = 0 # bases 5' of deletion that match the deletion
                while True:
                    try:
                        if ref != RHS_seq[rhs_i]: break
                    except TypeError:
                        print >>sys.stderr, ref, RHS_seq, line
                        print >>sys.stderr, `ref_genome`, `test_chromosome`, `end_pos+1`, `end_pos+11`
                        raise
                    rhs_i += 1
                    if rhs_i == len(RHS_seq): break
                lhs_end = len(LHS_seq) - 1
                while True:
                    if ref != LHS_seq[lhs_end-lhs_i]: break
                    lhs_i += 1
                    if lhs_i == len(LHS_seq): break 
                # i = lhs_i + rhs_i 
                i = lhs_i + rhs_i + 1 

            ## check for indel_type: repeats
            elif ref == RHS_seq[0:del_length] or ref == LHS_seq[4*del_length:5*del_length]:
                indel_type = 'repeats'
                rhs_i = 0 # bases 3' of deletion that match the deletion
                lhs_i = 0 # bases 5' of deletion that match the deletion
                while True:
                    try:
                        if ref != RHS_seq[(0+rhs_i*del_length) : (1+rhs_i)*del_length]: break
                    except TypeError:
                        print >>sys.stderr, ref, RHS_seq, line
                        print >>sys.stderr, `ref_genome`, `test_chromosome`, `end_pos+1`, `end_pos+11`
                        raise
                    rhs_i += 1
                    if rhs_i == 5: break
                lhs_end = len(LHS_seq)
                while True:
                    if ref != LHS_seq[(len(LHS_seq)-(lhs_i+1)*del_length) : (len(LHS_seq)-lhs_i*del_length)]: break
                    lhs_i += 1
                    if lhs_i == 5: break
                # i = lhs_i + rhs_i # use sum or max? TODO
                i = lhs_i + rhs_i + 1

            ## check for indel_type: micro-homology
            else:
                indel_type = 'MH'    
                rhs_i = 0 # bases 3' of deletion that match beginning of deletion
                lhs_i = 0 # bases 5' of deletion that match end of deletion
                while True:
                    try:
                        if ref[rhs_i] != RHS_seq[rhs_i]: break
                    except TypeError:
                        print >>sys.stderr, ref, RHS_seq, line
                        print >>sys.stderr, `ref_genome`, `test_chromosome`, `end_pos+1`, `end_pos+11`
                        raise
                    rhs_i += 1
                    if rhs_i == len(ref): break
                ref_end = len(ref) - 1
                lhs_end = len(LHS_seq) - 1
                while True:
                    if ref[ref_end-lhs_i] != LHS_seq[lhs_end-lhs_i]: break
                    lhs_i += 1
                    if lhs_i > ref_end: break 
                i = max(lhs_i, rhs_i) # use sum or max? TODO
                if i == 0: 
                    indel_type = 'repeats'
                    i = 1

        
        else: # vartype == 'INS'
            if ref != '-': ## complex indels??
                print >>sys.stderr, 'ignore complex indel: %s, %s, %s, %s, %s|%s|%s' % (chromosome, start_pos, end_pos, ref, LHS_seq, var, RHS_seq)
                continue
            if len(var) == 0: ## lack insertion info
                print >>sys.stderr, 'Lack insertion info: %s, %s, %s, %s, %s|%s|%s' % (chromosome, start_pos, end_pos, ref, LHS_seq, var, RHS_seq)
                continue
            ins_length = len(var)
            ## indel_size : 1,2,3,4,5+
            if ins_length <5: indel_size = str(ins_length)
            else: indel_size = '5+'

            try:
                LHS_seq = get_reference_seq(ref_genome, test_chromosome, start_pos-5*ins_length, start_pos)
                RHS_seq = get_reference_seq(ref_genome, test_chromosome, end_pos, end_pos+5*ins_length-1)
            except ValueError:
                # blacklist this sample
                print >>sys.stderr, '%s (%s) now blacklisted due to invalid position %s %s:%s' % (sample_id, input_path, ref_genome, chromosome, start_pos)
                blacklist.add(sample_id)
                continue

            if ins_length == 1:
                if var in ['A','G']: 
                    indel_type = _complement[var]
                    #if strand == '+': indel_type = indel_type + '_T'
                    #elif strand == '-': indel_type = indel_type + '_N'
                    #else: indel_type = indel_type + '_U'
                else: 
                    indel_type = var
                    #if strand == '+': indel_type = indel_type + '_N'
                    #elif strand == '-': indel_type = indel_type + '_T'
                    #else: indel_type = indel_type + '_U'

                ## get the repeat size i
                rhs_i = 0 # bases 3' of deletion that match the deletion
                lhs_i = 0 # bases 5' of deletion that match the deletion
                while True:
                    try:
                        if var != RHS_seq[rhs_i]: break
                    except TypeError:
                        print >>sys.stderr, ref, RHS_seq, line
                        print >>sys.stderr, `ref_genome`, `test_chromosome`, `end_pos+1`, `end_pos+11`
                        raise
                    rhs_i += 1
                    if rhs_i == len(RHS_seq): break
                lhs_end = len(LHS_seq) - 1
                while True:
                    if var != LHS_seq[lhs_end-lhs_i]: break
                    lhs_i += 1
                    if lhs_i == len(LHS_seq): break 
                i = lhs_i + rhs_i 

            ## check for indel_type: repeats
            # elif var == RHS_seq[0:ins_length] or var == LHS_seq[4*ins_length:5*ins_length]:
            else:
                indel_type = 'repeats'
                rhs_i = 0 # bases 3' of deletion that match the deletion
                lhs_i = 0 # bases 5' of deletion that match the deletion
                while True:
                    try:
                        if var != RHS_seq[(0+rhs_i*ins_length) : (1+rhs_i)*ins_length]: break
                    except TypeError:
                        print >>sys.stderr, ref, RHS_seq, line
                        print >>sys.stderr, `ref_genome`, `test_chromosome`, `end_pos+1`, `end_pos+11`
                        raise
                    rhs_i += 1
                    if rhs_i == 5: break
                lhs_end = len(LHS_seq)
                while True:
                    if var != LHS_seq[(len(LHS_seq)-(lhs_i+1)*ins_length) : (len(LHS_seq)-lhs_i*ins_length)]: break
                    lhs_i += 1
                    if lhs_i == 5: break
                i = lhs_i + rhs_i # use sum or max? TODO

            ############### We don't care about micro-homology Insertions
            ## check for indel_type: micro-homology
            # else:
            #     indel_type = 'MH'
            #     rhs_i = 0 # bases 3' of insertion that match beginning of insertion
            #     lhs_i = 0 # bases 5' of insertion that match end of insertion
            #     while True:
            #         try:
            #             if var[rhs_i] != RHS_seq[rhs_i]: break
            #         except TypeError:
            #             print >>sys.stderr, ref, RHS_seq, line
            #             print >>sys.stderr, `ref_genome`, `test_chromosome`, `end_pos`, `end_pos+10`
            #             raise
            #         rhs_i += 1
            #         if rhs_i == len(var): break
            #     var_end = len(var) - 1
            #     lhs_end = len(LHS_seq) - 1
            #     while True:
            #         if var[var_end-lhs_i] != LHS_seq[lhs_end-lhs_i]: break
            #         lhs_i += 1
            #         if lhs_i > var_end: break
            #     i = max(lhs_i, rhs_i) # use sum or max? TODO
            #     if i == 0: indel_type = 'repeats'

        ## combine all features: vartype * indel_size * indel_type * i
        if vartype == 'DEL' and indel_type != 'MH':
            if i < 6: repeat_size = str(i)
            else: repeat_size = '6+'
        
        else:
            if i < 5: repeat_size = str(i)
            else: repeat_size = '5+'

        features = '%s,%s,%s,%s' % (vartype, indel_type, indel_size, repeat_size)
        
        if debug_level:
            ref_full_base_context = get_reference_seq(ref_genome, test_chromosome, start_pos-10, end_pos+10)
            print >>sys.stderr, '%s,%s,%s,%s,%s,%s,%s,%s,%s|%s|%s' % (features, chromosome, start_pos, end_pos, vartype, ref, var, ref_full_base_context, LHS_seq, ref, RHS_seq)

        out_handle.write(
            '\t'.join([sample_id, lookup_ref_genome, chromosome, str(start_pos), ref, var, features, transcript, strand, gene, annotations]) + '\n'
        )

    out_handle.close()
    if out_proc:
        out_proc.wait()
        out_handle_real.close()
    in_handle.close()
    if in_proc: in_proc.wait() # for gzip

    if ignored_complex and debug_level:
        print >>sys.stderr, 'ignored %d complex insertions with deletions from %s' % (ignored_complex, input_path)
    if blacklist:
        # filter out any mutations for blacklisted sample(s) that were written
        # before we blacklisted it
        in_proc, out_proc = None, None
        if output_path.endswith('.gz'):
            in_proc = subprocess.Popen(('gzip', '-dc', output_path + '-tmp'), stdout=subprocess.PIPE)
            in_handle = in_proc.stdout
            out_handle_real = open(output_path, 'w')
            out_proc = subprocess.Popen(('gzip', '-c'), stdin=subprocess.PIPE, stdout=out_handle_real)
            out_handle = out_proc.stdin
        else:
            in_handle = open(output_path + '-tmp')
            out_handle = open(output_path, 'w')
        for line in in_handle:
            sample_id = line[:line.find('\t')]
            if sample_id not in blacklist:
                out_handle.write(line)
        out_handle.close()
        in_handle.close()
        if in_proc:
            in_proc.wait()
            out_proc.wait()
        os.unlink(output_path + '-tmp')
    else:
        os.rename(output_path + '-tmp', output_path)


def _aggregate_counts_indel(annotated_cache_file, input_filename, aggr_strand='+', filter_to=None, variant_filters=None):
    """ aggr strand is either '+' or 'sense'... it should only be 'sense'
        if filter_to is also set (to 'transcript'). All callers(?) give
        indels relative to reference + strand (?).

        returns ([features], {sample: [counts]})."""

    # save to cache, or read existing cached file if present
    counts_cache_file = input_filename
    if counts_cache_file.startswith('./'):
        counts_cache_file = counts_cache_file[2:]
    counts_cache_file = _cachedir + '/counts/' + counts_cache_file.replace('/','%2F') + '-' + aggr_strand
    if filter_to: # 'exome' or 'transcript'
        counts_cache_file += '-' + filter_to

    observed_counts = {} # sample: {feature: count}
    # TODO - check if cache is outdated
    # Don't use cache if we are filtering variants by any criteria
    if os.path.exists(counts_cache_file) and not variant_filters:
        with open(counts_cache_file) as in_handle:
            cached = json.load(in_handle)
            feature_names, observed_counts = cached['keys'], cached['counts']
    else:
        # read in the annotated_cache file
        proc = None
        if annotated_cache_file.endswith('.gz'):
            proc = subprocess.Popen(('gzip', '-dc', annotated_cache_file), stdout=subprocess.PIPE)
            in_handle = proc.stdout
        else:
            in_handle = open(annotated_cache_file)
        ignored_annotations = set()
        for line in in_handle:
            # note - refbase and varbase will be complemented if - strand
            sample_id, ref_genome, chrome, pos, refseq, varseq, features, transcript, strand, gene, annotations = line.rstrip('\n').split('\t')
            # ignore intergenic var:
            if filter_to == 'transcript' and not transcript: continue
            # TODO! need to determine if in an exon!
            if filter_to == 'exome' and not transcript: continue

            if variant_filters and annotations:
                skip = False
                for ann in annotations.split(';'):
                    ann = ann.strip() # ignore wrapping whitespace
                    try:
                        key,value = ann.split('=')
                    except ValueError:
                        # some simple files from Ludmil have existing
                        # annotations like "SOMATIC" or "dbSNP v147". ignore
                        if ann not in ignored_annotations:
                            if ann != 'NA' and debug_level:
                                print >>sys.stderr, 'ignoring annotation "%s"' % ann
                            ignored_annotations.add(ann)
                        continue
                    # special case
                    if key == 'ncallers':
                        value, callers = value.split(',')
                    func = variant_filters.get(key)
                    if func and func(value): # if True, should ignore var
                        skip = True
                        break
                if skip: continue # ignore this variant
            # context is 5 bases
            if aggr_strand == 'sense':
                # TODO - flip tandem repeat base(s)
                if len(strand) != 1: # ie '-,+' instead of '+' or '-'
                    # ignore position transcribed on both strands
                    if debug_level:
                        print >>sys.stderr, 'ignoring %s %s %s:%s for 192 classes since transcribed on both strands' % (input_filename, ref_genome, chrome, pos)
                    continue
                pass # TODO
            for key in features.split(';'):
                try:
                    observed_counts[sample_id][key] = observed_counts[sample_id].get(key,0) + 1
                except KeyError:
                    observed_counts[sample_id] = {key: 1}
        in_handle.close()
        if proc: proc.wait()
        # change values from a hash to a list, and sort that each sample has
        # them in the same order.
        feature_names = set()
        for counts in observed_counts.values(): # for each sample's features
            feature_names.update(counts.keys()) # add names of features
        feature_names = sorted(feature_names) # set -> sorted list
        for sample_id in observed_counts.keys():
            observed_counts[sample_id] = [observed_counts[sample_id].get(k,0) for k in feature_names]

        # save to cache if we aren't filtering mutations
        if not variant_filters:
            tmpfilename = counts_cache_file + '-tmp'
            try:
                out_handle = open(tmpfilename, 'w')
            except:
                if not os.path.exists(_cachedir + '/counts'):
                    os.mkdir(_cachedir + '/counts')
                    out_handle = open(tmpfilename, 'w')
                else:
                    raise
            json.dump({'keys':feature_names, 'counts':observed_counts}, out_handle)
            out_handle.close()
            os.rename(tmpfilename, counts_cache_file)

    return feature_names, observed_counts

## generate reference indel catalogs:
# a = ["{:s},{:s}_{:s},{:s},{:s}".format(x, y, m, z, k) for x in ['DEL', 'INS'] for y in ['C', 'T'] 
#  for z in ['1'] for k in ['0','1','2','3','4','5+'] for m in ['T', 'N', 'U']]
a1 = ["{:s},{:s},{:s},{:s}".format(x, y, z, k) for x in ['DEL'] for y in ['C', 'T'] 
 for z in ['1'] for k in ['1','2','3','4','5','6+']]
a2 = ["{:s},{:s},{:s},{:s}".format(x, y, z, k) for x in ['INS'] for y in ['C', 'T'] 
 for z in ['1'] for k in ['0','1','2','3','4','5+']]
b1 = ["{:s},{:s},{:s},{:s}".format(x, y, z, k) for x in ['DEL'] for y in ['repeats'] 
 for z in ['2','3','4','5+'] for k in ['1','2','3','4','5','6+']]
b2 = ["{:s},{:s},{:s},{:s}".format(x, y, z, k) for x in ['INS'] for y in ['repeats'] 
 for z in ['2','3','4','5+'] for k in ['0','1','2','3','4','5+']]
c = ["{:s},{:s},{:s}".format(x, y, z) for x in ['DEL'] for y in ['MH'] 
 for z in ['2,1', '3,1', '3,2', '4,1', '4,2', '4,3', '5+,1', '5+,2', '5+,3', '5+,4', '5+,5+']]
ref_features = a1 + a2 + b1 + b2 + c

def write_indels_as_tsv(features, observed_counts, output_filename):
    """dumps feature counts to named output file."""
    sample_ids = sorted(observed_counts.keys())
    tmpfilename = output_filename + '-tmp'
    handle = open(tmpfilename, 'w')
    handle.write('Type\tSubtype\tIndel_size\tRepeat_MH_size\t' + '\t'.join(sample_ids) + '\n')

    for key in ref_features:
        feture_set = str.split(key, ",")
        handle.write("\t".join(feture_set))
        for patient in sample_ids:
            try:
                i = features.index(key)
            	handle.write('\t%d' % observed_counts[patient][i])
            except:
                handle.write('\t0')
        handle.write('\n')
    handle.close()
    os.rename(tmpfilename, output_filename)


def write_indels_as_csv(features, observed_counts, output_filename):
    """dumps feature counts to named output file."""
    sample_ids = sorted(observed_counts.keys())
    tmpfilename = output_filename + '-tmp'
    handle = open(tmpfilename, 'w')
    handle.write('Type,Subtype,Indel_size,Repeat_MH_size,' + ','.join(sample_ids) + '\n')
    
    for key in ref_features:
        feture_set = str.split(key, ",")
        handle.write(",".join(feture_set))
        for patient in sample_ids:
            try:
                i = features.index(key)
                handle.write(',%d' % observed_counts[patient][i])
            except:
                handle.write(',0')
        handle.write('\n')
    handle.close()
    os.rename(tmpfilename, output_filename)


def get_indel_class_counts(input_filename, genomic_region, strand='+', filter_to=None, variant_filters=None):
    """
        Given a .simple file containing indels, calculate various class counts.

        genome_region = "exome" or "genome"
        filter_to = "transcript", "exome", or None|"genome"

        Note - if genome_region=exome we assume that the input is already
            only in exonic regions, we don't check.
        returns a list of tuples (for Before/Ref/After/Var keys) and
        a dict of {sample: [counts]} where the counts are in the same order
        as the key's list.
    """

    # look for cached annotated file; add transcript,strand + context if needed
    annotated_cache_file = input_filename
    if annotated_cache_file.startswith('./'):
        annotated_cache_file = annotated_cache_file[2:]
    annotated_cache_file = _cachedir + '/annotated/' + annotated_cache_file.replace('/','%2F') + '-indels.gz'

    if not os.path.exists(annotated_cache_file):
        if not os.path.exists(_cachedir + '/annotated'):
            if not os.path.exists(_cachedir):
                os.mkdir(_cachedir)
            os.mkdir(_cachedir + '/annotated')
        try:
            if input_filename.find('.simple') > -1:
                annotate_simple_indel(input_filename, annotated_cache_file)
        except:
            print >>sys.stderr, ('failed while annotating', input_filename)
            raise

    return _aggregate_counts_indel(annotated_cache_file, input_filename, strand, filter_to, variant_filters)
