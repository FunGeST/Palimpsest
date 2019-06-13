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
# TODO - make these configurable instead of hard-coding
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
        elif not _grch38_reference and (n.count('h38') or n.count('hg38') or n.count('assembly19')):
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

def annotate_simple(input_path, output_path):
    """get annotation (+-2 base context) and save info to file."""

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
    refbase_idx = 8
    varbase_idx = 9
    annotations_idx = 10
    strand_idx = -1

    blacklist = set() # ignore samples with invalid coordinates
    sample_ids_seen = set()
    # in case we have duplicate rows (eg multiple transcripts for a variant)
    patient_coordinates_seen = set()
    ranges, last_genome = None, None # if checking coords are in a transcript
    out_handle = open(output_path + '-tmp', 'w')
    for line in in_handle:
        fields = line.rstrip('\n').split('\t')
        # print fields.    ## turn on to find lines in vcf that cause errors ############################################
        cancer_type = fields[disease_idx]
        # for ICGC, this is a tissue sample ID, not really a patient ID
        uniq_id = fields[tumour_idx]
        sample_id = (uniq_id)
        sample_ids_seen.add(sample_id) # ignores duplicates
        if sample_id in blacklist: continue # ignore all further mutations
        vartype = fields[vartype_idx]
        # 'SNP' used by TCGA, 'single base substitution' used by ICGC
        # 'SNV' used in CCA lit. paper
        if vartype not in ('SNP', 'SNV', 'single base substitution'):
            continue # ignore indels for now
        chromosome = fields[chrome_idx]
        #print fields[pos_idx] #prints fuck loads of stuff be careful
        pos = int(float(fields[pos_idx]))

        if (sample_id, chromosome, pos) in patient_coordinates_seen:
            if debug_level:
                print >>sys.stderr, 'ignoring duplicate mutation in %s (%s) at %s:%d' % (sample_id, input_path, chromosome, pos)
            continue
        patient_coordinates_seen.add((sample_id, chromosome, pos))

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
            # should all start with 'chr'?? 
            ## should remove 'chr' ....
            if chromosome.startswith('chr'):
                #test_chromosome = chromosome
                test_chromosome = chromosome[4:]
            else:
                #test_chromosome = 'chr' + chromosome
                test_chromosome = chromosome
            # shouldn't happen, but does... :(
            # if test_chromosome == 'chr23': test_chromosome = 'chrX'
            # elif test_chromosome == 'chr24': test_chromosome = 'chrY'
            if test_chromosome == 'chr23': test_chromosome = 'X'
            elif test_chromosome == 'chr24': test_chromosome = 'Y'
        else: # we don't actually support any other references (yet)
            # TODO
            lookup_ref_genome = ref_genome
            test_chromosome = chromosome
        if last_genome != lookup_ref_genome:
            # need to get correct dict of coordinates for transcripts
            last_genome = lookup_ref_genome
            ranges = _transcripts_for_genome.get(lookup_ref_genome)
            if not ranges and lookup_ref_genome in ('hs37', 'grch38'):
                # not initialised; load ranges from a .bed file now
                if lookup_ref_genome == 'hs37':
                    bed_file = grch37_bed_file
                elif lookup_ref_genome == 'grch38':
                    bed_file = grch38_bed_file
                ranges = in_target.ranges(bed_file)
                _transcripts_for_genome[lookup_ref_genome] = ranges
        range_descr = ranges.get_target_name(test_chromosome, pos)
        if range_descr: # first found transcript
            # columns 4,5,6 of our .bed files
            # look for all transcripts, in case on both strands
            range_descrs = ranges.get_all_targets(test_chromosome, pos)
            try:
                # can't use zip directly because we give it a sequence,
                # not the contents of the sequence
                genes, transcripts, strands = apply(zip,[r.split('\t', 2) for r in range_descrs])
            except ValueError:
                print >>sys.stderr, lookup_ref_genome, `test_chromosome`, `pos`, `range_descr`
                print >>sys.stderr, range_descrs
                print >>sys.stderr, ranges.target_filename, len(ranges.target_starts[test_chromosome])
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

        if len(var) > 1: # it contains 2 bases... remove ref base
            var = var.replace(ref, '')
        if var not in 'ACGT': # it's an ambiguity code
            try:
                var = _de_IUB[var,ref]
            except:
                print path
                print line
                raise

        # get the context sequence around this variant
        try:
            sequence = get_context(ref_genome, chromosome, pos, 2)
        except ValueError:
            # blacklist this sample
            print >>sys.stderr, '%s (%s) now blacklisted due to invalid position %s %s:%s' % (sample_id, input_path, ref_genome, chromosome, pos)
            blacklist.add(sample_id)
            continue

        if not sequence: # error already logged by get_context() function
            print >>sys.stderr, 'no context for %s %s:%s in %s' % (ref_genome, chromosome, pos, input_path)
            continue
        if 'N' in sequence:
            print >>sys.stderr, '%s (%s) now blacklisted due to N in %s %s:%s (%s)' % (sample_id, input_path, ref_genome, chromosome, pos, sequence)
            blacklist.add(sample_id)
            continue
        # TODO - check that the strand from the file actually matches the
        # strand that this gene is on. some maf files put everything relative
        # to + strand.
        if strand == '-':
            sequence = ''.join(_complement[b] for b in sequence[::-1])
        ## put ref and var consistent with the sequence context
        if ref == _complement[sequence[2]]:
            ref = _complement[ref]
            var = _complement[var]
        out_handle.write(
            '\t'.join([sample_id, lookup_ref_genome, chromosome, str(pos), ref, var, sequence, transcript, strand, gene, annotations]) + '\n'
        )

    out_handle.close()
    in_handle.close()
    if in_proc: in_proc.wait() # for gzip
    if blacklist:
        # filter out any mutations for blacklisted sample(s) that were written
        # before we blacklisted it
        in_handle = open(output_path + '-tmp')
        out_handle = open(output_path, 'w')
        for line in in_handle:
            sample_id = line[:line.find('\t')]
            if sample_id not in blacklist:
                out_handle.write(line)
        out_handle.close()
        in_handle.close()
        os.unlink(output_path + '-tmp')
    else:
        os.rename(output_path + '-tmp', output_path)


def _aggregate_counts(annotated_cache_file, input_filename, num_classes, filter_to=None, variant_filters=None):
    """returns ([keys], {sample: [counts]})."""

    # read in the annotated_cache file

    counts_cache_file = _cachedir + '/counts/' + input_filename.replace('/','%2F') + ('-%d' % num_classes)
    if filter_to: # 'exome' or 'transcript'
        counts_cache_file += '-' + filter_to

    observed_counts = {} # sample: {B,R,A,V: count}

    # TODO - check if cache is outdated
    # Don't use cache if we are filtering variants by any criteria
    if os.path.exists(counts_cache_file) and not variant_filters:
        with open(counts_cache_file) as in_handle:
            cached = json.load(in_handle)
            features, observed_counts = cached['keys'], cached['counts']
    else:
        in_handle = open(annotated_cache_file)
        ignored_annotations = set()
        for line in in_handle:
            # note - refbase and varbase will be complemented if - strand
            sample_id, ref_genome, chrome, pos, refbase, varbase, context, transcript, strand, gene, annotations = line.rstrip('\n').split('\t')
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
            if num_classes == 1536: # +- 2 bases
                before, after = context[:2], context[3:5]
                if refbase in 'AG': # fold to pyrimidine
                    before, after = _complement[after[1]]+_complement[after[0]], _complement[before[1]] + _complement[before[0]]
                    refbase, varbase = _complement[refbase], _complement[varbase]
            else: # +- 1 base
                if num_classes == 192 and len(strand) != 1:
                    # ignore position transcribed on both strands
                    if debug_level:
                        print >>sys.stderr, 'ignoring %s %s %s:%s for sense-strand counts since transcribed on both strands' % (input_filename, ref_genome, chrome, pos)
                    continue
                before, after = context[1], context[3]
                if num_classes == 96 and refbase in 'AG': # fold to pyrimidine
                    before, refbase, after, varbase = _complement[after], _complement[refbase], _complement[before], _complement[varbase]
            key = (before, refbase, after, varbase)
            try:
                observed_counts[sample_id][key] = observed_counts[sample_id].get(key,0) + 1
            except KeyError:
                observed_counts[sample_id] = {key: 1}
        in_handle.close()
        # change values from a hash to a list, and sort so that each sample has
        # them in the same order.
        # keys are B,R,A,V
        features = []
        bases = 'ACGT'
        if num_classes == 96 or num_classes == 1536: # not stranded
            refbases = 'CT'
        else:
            refbases = bases
        # sort by Ref, then Var, then Before, then After
        for ref in refbases:
            for var in bases:
                if ref == var: continue
                if num_classes == 1536:
                    contextbases = [b1+b2 for b1 in bases for b2 in bases]
                else:
                    contextbases = bases
                for before in contextbases:
                    for after in contextbases:
                        features.append((before,ref,after,var)) 
        #features.sort(lambda k1,k2: cmp(k1[1],k2[1]) or cmp(k1[3],k2[3]) or cmp(k1[0],k2[0]) or cmp(k1[2],k2[2]))
        for sample_id in observed_counts.keys():
            observed_counts[sample_id] = [observed_counts[sample_id].get(k,0) for k in features]

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
            json.dump({'keys':features, 'counts':observed_counts}, out_handle)
            out_handle.close()
            os.rename(tmpfilename, counts_cache_file)

    return features, observed_counts


def write_as_tsv(features, observed_counts, output_filename):
    """dumps feature counts to named output file."""
    sample_ids = sorted(observed_counts.keys())
    tmpfilename = output_filename + '-tmp'
    handle = open(tmpfilename, 'w')
    # TODO - if features is 5 columns instead of 4
    handle.write('Before\tRef\tAfter\tVar\t' + '\t'.join(sample_ids) + '\n')
    for i in range(len(features)):
        key = features[i]
        try:
            handle.write('%s\t%s\t%s\t%s' % tuple(key))
        except:
            print key
            raise
        for patient in sample_ids:
            # if we have a patient with no SNVs, then they won't
            # have a key created in the observed_counts dict.
            handle.write('\t%d' % observed_counts[patient][i])
        handle.write('\n')
    handle.close()
    os.rename(tmpfilename, output_filename)


def write_as_csv(features, observed_counts, output_filename):
    """dumps feature counts to named output file, formatted the same way
        as Ludmil's catalogs."""
    sample_ids = sorted(observed_counts.keys())
    tmpfilename = output_filename + '-tmp'
    handle = open(tmpfilename, 'w')
    num_features = len(features)
    if num_features <= 96:
        handle.write('Mutation type,Trinucleotide,' + ','.join(sample_ids) + '\n')
    elif num_features <= 192:
        handle.write('Strand,Mutation type,Trinucleotide,' + ','.join(sample_ids) + '\n')
    else: # 1536?
        handle.write('Mutation type,Pentanucleotide,' + ','.join(sample_ids) + '\n')
    
    for i in range(num_features):
        key = features[i]
        try:
            before,ref,after,var = key
            if num_features > 96 and num_features <= 192:
                if ref in 'CT':
                    strand = 'U'
                else:
                    before, after = _complement[after], _complement[before]
                    ref = _complement[ref]
                    var = _complement[var]
                    strand = 'T'
                handle.write('%s,%s>%s,%s%s%s' % (strand,ref,var,before,ref,after))
            else:
                handle.write('%s>%s,%s%s%s' % (ref,var,before,ref,after))
        except:
            print key
            raise
        for patient in sample_ids:
            # if we have a patient with no SNVs, then they won't
            # have a key created in the observed_counts dict.
            handle.write(',%d' % observed_counts[patient][i])
        handle.write('\n')
    handle.close()
    os.rename(tmpfilename, output_filename)


def get_class_counts(input_filename, genomic_region, num_classes=96, filter_to=None, variant_filters=None):
    """
        Given a .simple file containing SNVs, calculate mutation class counts.

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
        annotated_cache_file = annotated._cache_file[2:]
    annotated_cache_file = _cachedir + '/annotated/' + annotated_cache_file.replace('/','%2F')

    # ignore cached file if input file has been updated
    #if os.path.exists(counts_cache_file):
    #    cached_mtime = os.stat(counts_cache_file).st_mtime
    #    orig_mtime = os.stat(path).st_mtime
    #    if orig_mtime < cached_mtime:
    #        have_cached = True

    if not os.path.exists(annotated_cache_file):
        if not os.path.exists(_cachedir + '/annotated'):
            if not os.path.exists(_cachedir):
                os.mkdir(_cachedir)
            os.mkdir(_cachedir + '/annotated')
        try:
            if input_filename.find('.simple') > -1:
                annotate_simple(input_filename, annotated_cache_file)
        except:
            print >>sys.stderr, ('failed while annotating', input_filename)
            raise

    # for TCGA sample names: remap:
    # sample_id = _reformat_samplename(sample_id) # TODO
    return _aggregate_counts(annotated_cache_file, input_filename, num_classes, filter_to, variant_filters)

