#! /usr/bin/env python2
# -*- coding: utf-8 -*-
#
# This file is part of the Bacterial and Archaeal Genome Analyser
# Copyright (C) 2015-16 David Williams
# david.williams.at.liv.d-dub.org.uk
# License GPLv3+: GNU GPL version 3 or later
# This is free software: you are free to change and redistribute it
# There is NO WARRANTY, to the extent permitted by law
# 
# Work on this software was started at The University of Liverpool, UK 
# with funding from The Wellcome Trust (093306/Z/10) awarded to:
# Dr Steve Paterson (The University of Liverpool, UK)
# Dr Craig Winstanley (The University of Liverpool, UK)
# Dr Michael A Brockhurst (The University of York, UK)
#
'''
CallVariants module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains functions to infer single nucleotide polymorphisms, 
insertions and deletions from short reads mapped to a reference genome using 
the Genome Analysis Toolkit from the Broad Institute.
'''

# stdlib
from time import sleep as _sleep
from cStringIO import StringIO as _StringIO
from random import sample as _sample
from collections import OrderedDict as _OrderedDict
from collections import defaultdict as _defaultdict
from collections import Counter as _Counter
from glob import glob as _glob

from baga import _subprocess
from baga import _os
from baga import _multiprocessing
from baga import _cPickle
from baga import _gzip
from baga import _re
from baga import _tarfile
from baga import _array
from baga import _json
from baga import _time
from baga import _md5

# external Python modules
import pysam as _pysam
from Bio import SeqIO as _SeqIO
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord
from Bio.Data.CodonTable import TranslationError as _TranslationError

# package functions
from baga import decide_max_processes as _decide_max_processes
from baga import get_exe_path as _get_exe_path
from baga import report_time as _report_time
def main():
    pass

def parseVCF(path_to_VCF):
    '''returns (header, header_section_order, colnames, variants)'''
    header_section_order = _OrderedDict()
    header = _defaultdict(list)
    pattern = _re.compile('##([^=]+)=(.+)$')
    variants = []
    collect_variants = False
    for line in open(path_to_VCF):
        if collect_variants:
            variants += [line.rstrip()]
        elif line[:6] == '#CHROM':
            colnames = line.rstrip().split('\t')
            collect_variants = True
        else:
            section, value = _re.match(pattern, line.rstrip()).groups()
            header_section_order[section] = True
            header[section] += [line.rstrip()]

    header = dict(header)
    header_section_order = list(header_section_order)
    return(header, header_section_order, colnames, variants)

def dictify_vcf_header(header):
    '''convert VCF header to dict using re for quotes'''
    pattern1 = _re.compile(''',(?=(?:[^'"]|'[^']*'|"[^"]*")*$)''')
    pattern2 = _re.compile('''=(?=(?:[^'"]|'[^']*'|"[^"]*")*$)''')
    headerdict = {}
    for headersection,headerrows in header.items():
        for headerrow in headerrows:
            if headerrow[len('##{}='.format(headersection))] == '<':
                # split by commas ignoring quotes
                thisheaderdict = _re.split(pattern1, headerrow[len('##{}=<'.format(headersection)):-1])
                thisheaderdict = dict([_re.split(pattern2, i) for i in thisheaderdict])
                try:
                    headerdict[headersection] += [thisheaderdict]
                except KeyError:
                    headerdict[headersection] = [thisheaderdict]

    return(headerdict)

def sortVariantsKeepFilter(header, colnames, variantrows):
    '''
    Given a VCF file contents divided up by CallVariants.parseVCF(), return a 
    dict for tabulation of variants in which each is marked with appropriate 
    filters as described in the VCF header.
    Per sample filters stored in INFO columns must be described in header 
    ##INFO entries and have Descriptions starting "FILTER:".
    '''
    # identify chromosome for this VCF <== this will fail if running against > 1 contigs . . .
    pattern = _re.compile('##contig=<ID=([A-Za-z0-9_\.]+),length=([0-9]+)>')
    # for multiple chromosomes iterate here

    if 'contig' in header:
        genome_ids = []
        genome_lengths = []
        for contig in header['contig']:
            match = _re.match(pattern, contig)
            if match is None:
                raise Exception("Failed to identify which chromosome the variants were called on (couldn't find '##contig=' and/or recognise name)")
            
            genome_id, genome_length = match.groups()
            genome_lengths += [int(genome_length)]
            genome_ids += [genome_id]
        
        for genome_length, genome_id in zip(genome_lengths, genome_ids):
            print('Variants were called against {:,} bp genome: {}\n'.format(genome_length, genome_id))
    else:
        print('WARNING: no contig information found in header (reference free?)')


    parameter_names = colnames[:colnames.index('FORMAT')+1]
    parameter_names[0] = parameter_names[0].lstrip('#')
    sample_names = colnames[colnames.index('FORMAT')+1:]

    headerdict = dictify_vcf_header(header)

    FILTERfilters = set()
    for FILTER in headerdict['FILTER']:
        FILTERfilters.add(FILTER['ID'].strip('\'"'))

    INFOfilters = set()
    for INFO in headerdict['INFO']:
        if INFO['Description'].strip('"')[:len('FILTER:')] == 'FILTER:':
            INFOfilters.add(INFO['ID'].strip('\'"'))

    allfilters = FILTERfilters | INFOfilters

    if len(variantrows) == 0:
        variants = {}
        return(variants, allfilters)    

    cols = dict(zip(parameter_names, variantrows[0].rstrip().split('\t')[:len(parameter_names)]))

    e = 'Could not find FORMAT column in this VCF file. Probably no genotype data.'
    assert 'FORMAT' in cols, e
    e = 'Could not find GT under FORMAT column as per http://samtools.github.io/hts-specs/VCFv4.2.pdf. Probably no genotype data.'
    assert 'GT' in cols['FORMAT'], e

    sample_variant_info = cols['FORMAT'].split(':')

    # by sample chromosome position ref,query | filterlist

    # three level automatic dicts
    #                       sample               chromosome           position = [[ref,query], [filterlist]]
    variants = _defaultdict(lambda: _defaultdict(lambda: _defaultdict(dict)))

    for row in variantrows:
        # variant info
        cols = dict(zip(parameter_names, row.rstrip().split('\t')[:len(parameter_names)]))
        # collect per sample filters for this row
        INFO = dict([i.split('=') for i in cols['INFO'].split(';') if '=' in i])
        samples_filtered = {}
        for f in INFOfilters:
            if f in INFO:
                indexes = map(int, INFO[f].split(','))
                for i in indexes:
                    try:
                        samples_filtered[sample_names[i]].add(f)
                    except KeyError:
                        samples_filtered[sample_names[i]] = set([f])
        
        # collect sample wide filters
        thesefilterflags = set(cols['FILTER'].split(',')) - set(['PASS'])
        for sample_name in sample_names:
            try:
                samples_filtered[sample_name].update(thesefilterflags)
            except KeyError:
                samples_filtered[sample_name] = thesefilterflags
        
        # sample info
        sample_data = dict(zip(sample_names, row.rstrip().split('\t')[len(parameter_names):]))
        sample_data = dict([(s, dict(zip(sample_variant_info, d.split(':')))) for s,d in sample_data.items()])
        
        for sample,data in sample_data.items():
            if data['GT'] != '.':
                # there can be more than one variant per row
                # rare for SNPs but can happen for indels
                freqs = _Counter()
                alleles = cols['ALT'].split(',')
                allele_nums = map(int,data['GT'].split('/'))
                for allele_num in allele_nums:
                    if allele_num > 0:
                        freqs[alleles[allele_num-1]] += 1
                
                for ALT,freq in freqs.items():
                    if len(allele_nums) == 1:
                        # store ALT as single variant character
                        variants[sample][cols['CHROM']][int(cols['POS'])] = [[cols['REF'], ALT], samples_filtered[sample]]
                    else:
                        # store ALT as tuple of variant character, frequency, population_size
                        variants[sample][cols['CHROM']][int(cols['POS'])] = [[cols['REF'], (ALT,freq,len(allele_nums))], samples_filtered[sample]]
            elif data['GT'] == '.':
                # store ALT as '.' meaning insufficient information to call variant
                variants[sample][cols['CHROM']][int(cols['POS'])] = [[cols['REF'], '.'], samples_filtered[sample]]

    # convert nested defaultdicts to dicts
    variants = dict(variants)
    for k1,v1 in variants.items():
        variants[k1] = dict(v1)
        for k2,v2 in variants[k1].items():
            variants[k1][k2] = dict(v2)


    return(variants, allfilters)

def sortAmongBetweenReference(variants, sample_size):
    '''Separate "among sample" and "to reference" variants

    Takes a dictionary produced by .sortVariantsKeepFilter() as input, returns two 
    of the same shape in a single dictionary
    '''
    variant_freqs = _Counter()
    for sample, chromosomes in variants.items():
        # print('==> checking {}'.format(sample))
        for chromosome, positions in chromosomes.items():
            # iterate through variants by position
            for position, ((reference,query),filters) in sorted(positions.items()):
                variant_freqs[chromosome,position,query] += 1

    # sort by frequency
    among = set()
    to_ref = set()
    for info,freq in variant_freqs.items():
        if freq < sample_size:
            among.add(info)
        else:
            to_ref.add(info)

    print('among: {}; to reference: {}'.format(len(among),len(to_ref)))

    # divide up
    variants_among = {}
    variants_to_ref = {}
    for sample, chromosomes in variants.items():
        variants_among[sample] = {}
        variants_to_ref[sample] = {}
        for chromosome, positions in chromosomes.items():
            variants_among[sample][chromosome] = {}
            variants_to_ref[sample][chromosome] = {}
            # iterate through variants by position
            for position, ((reference,query),filters) in sorted(positions.items()):
                if (chromosome,position,query) in among:
                    variants_among[sample][chromosome][position] = ((reference,query),filters)
                elif (chromosome,position,query) in to_ref:
                    variants_to_ref[sample][chromosome][position] = ((reference,query),filters)
                else:
                    raise Exception('dividing variants between and among failed!')

    return({'among':variants_among, 'to_reference':variants_to_ref})

def to_by_position_filtered(variants, filters_applied, summarise = True):
    '''Sort variants by position and divide into non-filtered and filtered'''
    by_position = _defaultdict(_Counter)
    by_position_filtered = _defaultdict(_Counter)
    for sample, chromosomes in variants.items():
        for chromosome, positions in chromosomes.items():
            # iterate through variants by position
            for position, ((reference,query),filters) in sorted(positions.items()):
                if len(filters & set(filters_applied)) == 0:
                    # retain variants without any filters flagged (of those we are interested in)
                    by_position[chromosome][(position,reference,query,None)] += 1
                else:
                    for f in filters & set(filters_applied):
                        # also retain those with a filter flag, separately for each filter
                        by_position_filtered[chromosome][(position,reference,query,f)] += 1

    if summarise:
        for f1 in sorted(filters_applied):
            print('-- {} --'.format(f1))
            these = []
            if len(variants):
                for position,reference,query,f2 in sorted(by_position_filtered[chromosome]):
                    if f1 == f2:
                        these += ['{}: {} => {}, {}'.format(position,reference,query,f2)]
            
            print('Total: {}'.format(len(these)))
            print('\n'.join(these))

    return(by_position,by_position_filtered)

def reportCumulative(filter_order, reference_id, VCFs, VCFs_indels = False):
    '''
    Generate simple table of cumulative totals after filters applied to a VCF file

    This function parses VCF files already created by BAGA with various filters 
    applied and writes the total of each class of variant (SNP, indel) removed
    by each filter to a .csv file.
    '''
    # cumulative affect of filters

    ## build table column names
    colnames = ["Cumulative filters applied"]
    # dataset == reads_names
    for dataset in sorted(VCFs):
        for v in ['SNPs', 'InDels']:
            colnames += ["{} {}".format(v, dataset)]

    variant_type_order = ['SNPs', 'InDels']

    for v in variant_type_order:
        colnames += ["{} all samples".format(v)]



    # start with no filters
    filters_applied_ordered = [()]
    # then add some
    for filtername in filter_order:
        filters_applied_ordered += [short2fullnames[filtername]]

    collect_baga_filters = [f for f in filter_order if 'GATK' not in f]
    from glob import glob as _glob
    # find actual VCFs . . find the VCFs with suffixes that match the requested filters
    VCFs_use = {}
    for dataset,varianttypes in sorted(VCFs.items()):
        VCFs_use[dataset] = {}
        for varianttype,fname in varianttypes.items():
            VCFs_use[dataset][varianttype] = []
            if not isinstance(fname, list):
                # --calleach --calljoint
                filenames = [fname]
            else:
                # --callsingles
                filenames = fname
            for filename in filenames:
                bits = filename.split(_os.path.extsep)
                pattern = _os.path.extsep.join(bits[:-1]) + '*' + bits[-1]
                for checkthis in _glob(pattern):
                    filter_present = []
                    for fltr in collect_baga_filters:
                        if 'F_'+fltr in checkthis:
                            filter_present += [fltr]
                    
                    if set(filter_present) >= set(collect_baga_filters):
                        # OK if additional filters included in a VCF
                        VCFs_use[dataset][varianttype] += [checkthis]

    ### need to know (i) how many samples per dataset which may span VCF files or may not . . .

    # build table

    cumulative_filters = set()

    variant_groups = ('all', 'among', 'to_reference')
    rows = {group_name:list() for group_name in variant_groups}

    for filters in filters_applied_ordered:
        cumulative_filters.update(filters)
        
        this_row = {}
        for group_name in variant_groups:
            try:
                this_row[group_name] = [filter_names[filters]]
            except KeyError:
                this_row[group_name] = ['None']
        
        # must be saved as sets to only count variants once each
        totals_by_type = {}
        for group_name in variant_groups:
            totals_by_type[group_name] = {}
            for varianttype in variant_type_order:
                totals_by_type[group_name][varianttype] = set()
        
        for dataset,varianttypes in sorted(VCFs_use.items()):
            print('dataset: {}'.format(dataset))
            for varianttype in variant_type_order:
                for filename in varianttypes[varianttype]:
                    header, header_section_order, these_colnames, variantrows = parseVCF(filename)
                    variants, allfilters = sortVariantsKeepFilter(header, these_colnames, variantrows)
                    # divide variants into those among sample only, those between sample
                    # and reference
                    variants_divided = sortAmongBetweenReference(variants, sample_size = len(these_colnames[9:]))
                    variants_divided['all'] = variants
                    # cumulative filters applied here
                    for group_name in variant_groups:
                        by_position, by_position_filtered = to_by_position_filtered(
                                variants_divided[group_name], cumulative_filters)
                        
                        # reference_id is the chromosome ID
                        print('{} {}'.format(len(by_position[reference_id]), varianttype))
                        this_row[group_name] += [len(by_position[reference_id])]
                        totals_by_type[group_name][varianttype].update([info[0] for info in by_position[reference_id]])
        
        # add totals for variant class in columns corresponding to variant_type_order
        for group_name in variant_groups:
            this_row[group_name] += [len(totals_by_type[group_name][varianttype]) for varianttype in variant_type_order]
            rows[group_name] += [this_row[group_name]]


    for group_name in variant_groups:
        # just the totals by variant class
        outfilename = 'Table_of_cumulative_variant_totals_as_filters_applied_{}.csv'.format(group_name)
        print('Printing to {}'.format(outfilename))
        with open(outfilename, 'w') as fout:
            fout.write(','.join(['"{}"'.format(c) for c in colnames])+'\n')
            for row in rows[group_name]:
                fout.write(','.join(['"{}"'.format(row[0])]+[str(c) for c in row[1:]])+'\n')

def reportLists(include_filters, reference_id, VCFs, VCFs_indels = False):
    '''
    Generate simple table of variants and the filters applied from a VCF file

    This function parses VCF files already created by BAGA with various filters 
    applied and writes each variant with names of filters if masked to a .csv
    file.
    '''

    ## Frequencies with filters applied divided by group that variants were called in
    ## Single position might appear more than once if affected by a filter in only
    ## some samples. 

    ## arbitrary group e.g. within versus between clades?
    ## VCFs are currently divided as called so this is not easily implemented
    ## would need to specify a second list of groups

    ## build table column names
    colnames = ["Position", "Reference", "Variant", "Frequency", "Sample Group", "Filter"]

    # dataset == reads_name == Sample Group
    # also include combined frequencies

    variant_type_order = ['SNPs', 'InDels']

    # GATK filters added upstream so are in all BAGA produced VCFs if GATK used
    collect_baga_filters = [f for f in include_filters if 'GATK' not in f]
    from glob import glob as _glob
    # find actual VCFs . . find the VCFs with suffixes that match the requested filters
    VCFs_use = {}
    filters_per_VCF = {}
    checked = []
    for dataset,varianttypes in sorted(VCFs.items()):
        VCFs_use[dataset] = {}
        for varianttype,filename in varianttypes.items():
            bits = filename.split(_os.path.extsep)
            pattern = _os.path.extsep.join(bits[:-1]) + '*' + bits[-1]
            print(pattern)
            for checkthis in _glob(pattern):
                filter_present = []
                for fltr in include_filters:
                    if 'F_'+fltr in checkthis:
                        filter_present += [fltr]
                
                checked += [checkthis]
                if set(filter_present) >= set(collect_baga_filters):
                    # OK if additional filters included in a VCF
                    VCFs_use[dataset][varianttype] = checkthis
                    # retain requested filters present in this VCF
                    filters_per_VCF[checkthis] = set(filter_present) & set(collect_baga_filters)
                    print('Selecting: {}'.format(checkthis))

    checked.sort()
    assert len(filters_per_VCF) > 0, 'No suitable VCFs located for filters: {}. '\
            'These were considered based on content of the CallVariants.CallerGATK '\
            'object:\n{}'.format(', '.join(collect_baga_filters), '\n'.join(checked))

    # build table

    # expand multi-part filters
    include_filters2 = [a for b in [short2fullnames[f] for f in include_filters] for a in b]

    variant_groups = ('all', 'among', 'to_reference')
    rows = {group_name:list() for group_name in variant_groups}

    for reads_name,varianttypes in sorted(VCFs_use.items()):
        print('=> Dataset: {}'.format(dataset))
        for varianttype in variant_type_order:
            print('==> Variant class: {}'.format(varianttype))
            filename = varianttypes[varianttype]
            header, header_section_order, these_colnames, variantrows = parseVCF(filename)
            variants, allfilters = sortVariantsKeepFilter(header, these_colnames, variantrows)
            # divide variants into those among sample only, those between sample
            # and reference
            variants_divided = sortAmongBetweenReference(variants, sample_size = len(these_colnames[9:]))
            variants_divided['all'] = variants
            # filters applied here
            for group_name in variant_groups:
                print('===> Filtered in variant group: {}'.format(group_name))
                by_position, by_position_filtered = to_by_position_filtered(
                        variants_divided[group_name], include_filters2)
                
                ### '"Position","Reference","Variant","Frequency","Sample Group","Filter"' <======
                
                # reference_id is the chromosome ID
                for info,frequency in sorted(by_position[reference_id].items()):
                    position,ref_char,query,this_filter = info
                    rows[group_name] += ['{},"{}","{}",{},"{}","retained"'.format(
                            position,ref_char,query,frequency,reads_name)]
                for filter_of_interest in include_filters2:
                    for info,frequency in sorted(by_position_filtered[reference_id].items()):
                        position,ref_char,query,this_filter = info
                        if filter_of_interest == this_filter:
                            rows[group_name] += ['{},"{}","{}",{},"{}","{}"'.format(
                                    position,ref_char,query,frequency,reads_name,this_filter)]

    for group_name in variant_groups:
        # just the totals by variant class
        outfilename = 'Table_of_variants_with_filters_if_masked_{}.csv'.format(group_name)
        print('Printing to {}'.format(outfilename))
        with open(outfilename, 'w') as fout:
            fout.write(','.join(['"{}"'.format(c) for c in colnames])+'\n')
            fout.write('\n'.join(rows[group_name])+'\n')

# hard coded information for known filter types
known_filters = {}
# per sample filters have custom INFO added to allow for 
# sample-specific per-site filtering

known_filters['rearrangements'] = {}
# if dict of ranges e.g. for rearrangements regions extended
# by no or low read depth mapping, need a list here with each
# INFO entry
known_filters['rearrangements']['string'] = [
'##INFO=<ID=rearrangements1,Number=.,Type=Integer,Description="FILTER: Within a region affected by rearrangements between reference genome and sample. Sample-specific. INFO field is list of base-0 indexes for failed samples in column order">',
'##INFO=<ID=rearrangements2,Number=.,Type=Integer,Description="FILTER: Adjacent to region affected by rearrangements between reference genome and sample and with >50% without reads mapped, i.e., absent in query or possibly insufficient mapping quality. Sample-specific. INFO field is list of base-0 indexes for failed samples in column order">']
known_filters['rearrangements']['per_sample'] = True

# reference genome specific filters need a conventional FILTER entry so all variants at a position are excluded
known_filters['genome_repeats'] = {}
known_filters['genome_repeats']['string'] = [
'##FILTER=<ID=genome_repeats,Description="Within a long repeat unit in the reference genome">'
]
known_filters['genome_repeats']['per_sample'] = False


# non-BAGA filters (e.g., from GATK) that are worth knowing about
other_filters = {'GATK':('LowQual', 'standard_hard_filter')}

# determine order of filters to include

# GATK's LowQual and standard_hard_filter
short2fullnames = dict(other_filters.items())
short2fullnames['genome_repeats'] = ('genome_repeats',)
short2fullnames['rearrangements'] = ('rearrangements1', 'rearrangements2')
# names for table
filter_names = {
            ('LowQual', 'standard_hard_filter'): "GATK standard 'hard' filter", 
            ('genome_repeats',): 'baga reference genome repeats', 
            ('rearrangements1', 'rearrangements2'): 'baga genome rearrangements'
}

class CallerGATK:
    '''
    Wrapper around Broad Institute's Genome Analysis Tool Kit for variant calling

    Requires a collection of short read datasets that have been aligned to the 
    same genome sequence in Sequence Alignment/Map (SAM) format for variant 
    calling using Genome Analysis Tool Kit (GATK).
    '''
    def __init__(self, alignments = False, baga = False):
        '''
        Initialise with:
        a baga.AlignReads.SAMs object
        or
        path to a saved baga.CallVariants.CallerGATK object
        '''
        assert alignments or baga, 'Instantiate with alignments or the path to a '\
                'previously saved CallerGATK'
        assert not (alignments and baga), 'Instantiate with alignments OR the path '\
                'to a previously saved CallerGATK!'
        
        if alignments:
            try:
                self.ready_BAMs = alignments.ready_BAMs
            except AttributeError:
                print('ERROR: baga.CallVariants.CallerGATK needs a baga.AlignReads.SAMs '\
                        'object with a "ready_BAMs" attribute. This can be obtained with '\
                        'the "IndelRealignGATK()" method of the AlignReads module.')
            try:
                self.genome_sequence = alignments.genome_sequence
                self.genome_id = alignments.genome_id
            except AttributeError:
                print('ERROR: baga.CallVariants.CallerGATK needs a baga.AlignReads.SAMs '\
                        'object with a "genome_sequence" attribute. This can be '\
                        'obtained by running all methods in the AlignReads module.')
        elif baga:
            with _tarfile.open(baga, "r:gz") as tar:
                for member in tar:
                    contents = _StringIO(tar.extractfile(member).read())
                    try:
                        # either json serialised conventional objects
                        contents = _json.loads(contents.getvalue())
                    except ValueError:
                        #print('json failed: {}'.format(member.name))
                        # or longer python array.array objects
                        contents = _array('c', contents.getvalue())
                    
                    setattr(self, member.name, contents)


    def saveLocal(self, name):
        '''
        Save processed object info to a local compressed archive of json strings.
        'name' should exclude extension: .baga will be added
        '''
        fileout = 'baga.CallVariants.CallerGATK-%s.baga' % name
        with _tarfile.open(fileout, "w:gz") as tar:
            print('Writing to {} . . . '.format(fileout))
            for att_name, att in self.__dict__.items():
                if isinstance(att, _array):
                    io = _StringIO(att.tostring())
                    io.seek(0, _os.SEEK_END)
                    length = io.tell()
                    io.seek(0)
                    thisone = _tarfile.TarInfo(name = att_name)
                    thisone.size = length
                    tar.addfile(tarinfo = thisone, fileobj = io)
                else:
                    # try saving everything else here by jsoning
                    try:
                        io = _StringIO()
                        _json.dump(att, io)
                        io.seek(0, _os.SEEK_END)
                        length = io.tell()
                        io.seek(0)
                        thisone = _tarfile.TarInfo(name = att_name)
                        thisone.size = length
                        tar.addfile(tarinfo = thisone, fileobj = io)
                    except TypeError:
                        # ignore non-jsonable things like functions
                        # include unicodes, strings, lists etc etc
                        #print('omitting {}'.format(att_name))
                        pass

    def CallVCFsGATK(self, 
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            local_variants_path = ['variants'],
            use_java = 'java',
            force = False,
            mem_num_gigs = 8, 
            max_cpus = -1,
            arguments = False):
        '''
        Part of GATK "Best Practices" for DNA sequencing variant calling
        https://www.broadinstitute.org/gatk/guide/best-practices/?bpm=DNAseq
        
        this method follows the single sample variant calling described here:
        https://www.broadinstitute.org/gatk/guide/article?id=2803
        not the joint genotyping (see .CallgVCFsGATK() and .GenotypeGVCFsGATK())
        
        An output after this method is a list of str to the
        paths of the VCF for each sample. Because there is a 
        recalibration step, the list of str is added to a list (of lists) here:
        self.path_to_unfiltered_VCF
        
        If this is a list of lists (not str), downstream steps can infer whether a 
        separate genotyping (not joint) analysis is being run.
        
        max_cpus for this GATK module is "cpu threads per data thread"
        '''

        print(self.genome_id)

        genome_fna = 'genome_sequences/%s.fna' % self.genome_id
        if not _os.path.exists(genome_fna):
            _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                    genome_fna, 
                    'fasta')

        jar = _os.path.sep.join(jar)
        local_variants_path = _os.path.sep.join(local_variants_path)
        if not _os.path.exists(local_variants_path):
            _os.makedirs(local_variants_path)

        local_variants_path_genome = _os.path.sep.join([
                                        local_variants_path,
                                        self.genome_id])

        if not _os.path.exists(local_variants_path_genome):
            _os.makedirs(local_variants_path_genome)

        max_processes = _decide_max_processes( max_cpus )

        start_time = _time.time()
        paths_to_raw_VCFs = []
        exe = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar]
        # call the last set of ready BAMs added
        for cnum,BAM in enumerate(self.ready_BAMs[-1]):
            VCF_out = BAM[:-4] + '_unfiltered.VCF'
            VCF_out = _os.path.sep.join([local_variants_path_genome, VCF_out.split(_os.path.sep)[-1]])
            if not _os.path.exists(VCF_out) or force:
                # cmd should be built as 'option':[argument list] dictionary
                # with None as values for flag options
                cmd = {}
                # '-T', 'HaplotypeCaller', '-R', genome_fna, '-I', BAM, #'-L', '20', 
                cmd['-T'] = ['HaplotypeCaller']
                cmd['-R'] = [genome_fna]
                cmd['-I'] = [BAM]
                # '--genotyping_mode', 'DISCOVERY',
                cmd['--genotyping_mode'] = ['DISCOVERY']
                #'--sample_ploidy', '1',
                cmd['--sample_ploidy'] = ['1']
                #'--heterozygosity', '0.0001',         # this is less expected diversity than the default 0.001
                cmd['--heterozygosity'] = ['0.0001']
                #'--indel_heterozygosity', '0.00001',  # this is less expected diversity than the default 0.0001
                cmd['--indel_heterozygosity'] = ['0.00001']
                #'--emitRefConfidence', 'GVCF',        # not for single sample calling
                #'--variant_index_type', 'LINEAR',
                cmd['--variant_index_type'] = ['LINEAR']
                #'--variant_index_parameter', '128000',
                cmd['--variant_index_parameter'] = ['128000']
                #'-nct',  str(max_processes),
                cmd['-nct'] = [str(max_processes)]
                #'-stand_emit_conf', '10', 
                cmd['-stand_emit_conf'] = ['10']
                #'-stand_call_conf', '20', 
                cmd['-stand_call_conf'] = ['20']
                #'-o', VCF_out]
                cmd['-o'] = [VCF_out]
                
                if arguments:
                    # overwrite defaults with direct arguments
                    # (e.g. via -A/--arguments cli)
                    from baga import parse_new_arguments
                    cmd = parse_new_arguments(arguments, cmd)
                
                # make commands into a list suitable for subprocess
                cmds = []
                for opt,arg in cmd.items():
                    cmds += [opt]
                    if arg is not None:
                        cmds += arg
                
                print('Called: %s' % (' '.join(map(str, exe + cmds))))
                _subprocess.call(exe + cmds)
                
            else:
                print('Found:')
                print(VCF_out)
                print('use "force = True" to overwrite')
            
            paths_to_raw_VCFs += [VCF_out]
            
            # report durations, time left etc
            _report_time(start_time, cnum, len(self.ready_BAMs[-1]))

        # add to a list because this is done twice
        if hasattr(self, 'path_to_unfiltered_VCF'):
            self.path_to_unfiltered_VCF += [paths_to_raw_VCFs]
        else:
            self.path_to_unfiltered_VCF = [paths_to_raw_VCFs]

    def CallgVCFsGATK(self, 
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            local_variants_path = ['variants'],
            use_java = 'java',
            force = False,
            mem_num_gigs = 8, 
            max_cpus = -1,
            arguments = False):
        '''
        Part of GATK "Best Practices" for DNA sequencing variant calling
        https://www.broadinstitute.org/gatk/guide/best-practices/?bpm=DNAseq
        
        this method follows the joint genotyping calling described here:
        https://www.broadinstitute.org/gatk/guide/article?id=3893
        the method .GenotypeGVCFsGATK() should be called after this.
        
        An output after this method and .GenotypeGVCFsGATK() is a str to the
        path of the combined joint genotyped VCF. Because there is a 
        recalibration step, the str is added to a list here:
        self.path_to_unfiltered_VCF
        
        If this is a list of str (not lists), downstream steps can infer whether a 
        joint genotyping (not separate) analysis is being run.
        
        max_cpus for this GATK module is "cpu threads per data thread"
        '''

        print(self.genome_id)

        genome_fna = 'genome_sequences/%s.fna' % self.genome_id
        if not _os.path.exists(genome_fna):
            _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                    genome_fna, 
                    'fasta')

        jar = _os.path.sep.join(jar)
        local_variants_path = _os.path.sep.join(local_variants_path)
        if not _os.path.exists(local_variants_path):
            _os.makedirs(local_variants_path)

        local_variants_path_genome = _os.path.sep.join([
                                        local_variants_path,
                                        self.genome_id])

        if not _os.path.exists(local_variants_path_genome):
            _os.makedirs(local_variants_path_genome)

        max_processes = _decide_max_processes( max_cpus )

        start_time = _time.time()
        paths_to_raw_gVCFs = []
        exe = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar]
        # call the last set of ready BAMs added
        for cnum,BAM in enumerate(self.ready_BAMs[-1]):
            VCF_out = BAM[:-4] + '_unfiltered.gVCF'
            VCF_out = _os.path.sep.join([local_variants_path_genome, VCF_out.split(_os.path.sep)[-1]])
            if not _os.path.exists(VCF_out) or force:
                # cmd should be built as 'option':[argument list] dictionary
                # with None as values for flag options
                cmd = {}
                # '-T', 'HaplotypeCaller', '-R', genome_fna, '-I', BAM, #'-L', '20', 
                cmd['-T'] = ['HaplotypeCaller']
                cmd['-R'] = [genome_fna]
                cmd['-I'] = [BAM]
                # '--genotyping_mode', 'DISCOVERY',
                cmd['--genotyping_mode'] = ['DISCOVERY']
                #'--sample_ploidy', '1',
                cmd['--sample_ploidy'] = ['1']
                #'--heterozygosity', '0.0001',         # this is less expected diversity than the default 0.001
                cmd['--heterozygosity'] = ['0.0001']
                #'--indel_heterozygosity', '0.00001',  # this is less expected diversity than the default 0.0001
                cmd['--indel_heterozygosity'] = ['0.00001']
                #'--emitRefConfidence', 'GVCF',        # make vcfs appropriate for doing GenotypeGVCFs after <========
                cmd['--emitRefConfidence'] = ['GVCF']  # make vcfs appropriate for doing GenotypeGVCFs after
                #'--variant_index_type', 'LINEAR',
                cmd['--variant_index_type'] = ['LINEAR']
                #'--variant_index_parameter', '128000',
                cmd['--variant_index_parameter'] = ['128000']
                #'-nct',  str(max_processes),
                cmd['-nct'] = [str(max_processes)]
                #'-stand_emit_conf', '10', 
                cmd['-stand_emit_conf'] = ['10']
                #'-stand_call_conf', '20', 
                cmd['-stand_call_conf'] = ['20']
                #'-o', VCF_out]
                cmd['-o'] = [VCF_out]
                
                if arguments:
                    # overwrite defaults with direct arguments
                    # (e.g. via -A/--arguments cli)
                    from baga import parse_new_arguments
                    cmd = parse_new_arguments(arguments, cmd)
                
                # make commands into a list suitable for subprocess
                cmds = []
                for opt,arg in cmd.items():
                    cmds += [opt]
                    if arg is not None:
                        cmds += arg
                
                print('Called: %s' % (' '.join(map(str, exe + cmds))))
                _subprocess.call(exe + cmds)
                
            else:
                print('Found:')
                print(VCF_out)
                print('use "force = True" to overwrite')
            
            paths_to_raw_gVCFs += [VCF_out]
            
            # report durations, time left etc
            _report_time(start_time, cnum, len(self.ready_BAMs[-1]))

        # add to a list because this is done twice
        if hasattr(self, 'paths_to_raw_gVCFs'):
            self.paths_to_raw_gVCFs += [paths_to_raw_gVCFs]
        else:
            self.paths_to_raw_gVCFs = [paths_to_raw_gVCFs]

    def GenotypeGVCFsGATK(self, data_group_name,
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            local_variants_path = ['variants'],
            use_java = 'java',
            force = False,
            mem_num_gigs = 8,
            arguments = False):
        
        jar = _os.path.sep.join(jar)
        local_variants_path = _os.path.sep.join(local_variants_path)
        local_variants_path_genome = _os.path.sep.join([
                                        local_variants_path,
                                        self.genome_id])

        genome_fna = 'genome_sequences/%s.fna' % self.genome_id
        if not _os.path.exists(genome_fna):
            _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                    genome_fna, 
                    'fasta')

        e1 = 'Could not find "paths_to_raw_gVCFs" attribute. \
        Before starting performing joint GATK analysis, variants must be called. \
        Please run:\n\
        CallgVCFsGATK()\n\
        method on this SAMs instance.'

        assert hasattr(self, 'paths_to_raw_gVCFs'), 'Could not find "paths_to_raw_gVCFs" '\
                'attribute. Before starting performing joint GATK analysis, variants must '\
                'be called. Please run:\nCallgVCFsGATK()\nmethod on this SAMs instance.'

        # use the last VCFs called
        open('variants.list', 'w').write('\n'.join(self.paths_to_raw_gVCFs[-1]))

        # this method can be called prior or post base score recalibration
        # so give output a number corresponding to how many times variants called <===== need to never let total to go over two . . . <== fail with a warning/error?
        use_name = '{}_{}_samples_unfiltered.vcf'.format(data_group_name, len(self.paths_to_raw_gVCFs))
        VCF_out = _os.path.sep.join([local_variants_path_genome, use_name])

        exe = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar]

        # cmd should be built as 'option':[argument list] dictionary
        # with None as values for flag options
        cmd = {}
        cmd['-T'] = ['GenotypeGVCFs']
        cmd['-R'] = [genome_fna]
        cmd['--heterozygosity'] = ['0.0001']        # 650 total indels prior in all samples i.e., population
        cmd['--indel_heterozygosity'] = ['0.00001'] # 65 total indels prior in all samples i.e., population
        cmd['-stand_emit_conf'] = ['10']
        cmd['-stand_call_conf'] = ['20']
        cmd['-V'] = ['variants.list']
        cmd['-o'] = [VCF_out]

        if arguments:
            # overwrite defaults with direct arguments
            # (e.g. via -A/--arguments cli)
            from baga import parse_new_arguments
            cmd = parse_new_arguments(arguments, cmd)

        # make commands into a list suitable for subprocess
        cmds = []
        for opt,arg in cmd.items():
            cmds += [opt]
            if arg is not None:
                cmds += arg

        print('Called: %s' % (' '.join(map(str, exe + cmds))))
        _subprocess.call(exe + cmds)


        # add to a list because this is done twice
        if hasattr(self, 'path_to_unfiltered_VCF'):
            self.path_to_unfiltered_VCF += [VCF_out]
        else:
            self.path_to_unfiltered_VCF = [VCF_out]

    def hardfilterSNPsGATK(self, 
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            use_java = 'java',
            force = False):
        
        jar = _os.path.sep.join(jar)
        genome_fna = 'genome_sequences/%s.fna' % self.genome_id
        if not _os.path.exists(genome_fna):
            _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                    genome_fna, 
                    'fasta')

        e1 = 'Could not find "path_to_unfiltered_VCF" attribute. \
        Before filtering, joint calling of variants is necessary. \
        Please run:\n\
        CallgVCFsGATK and GenotypeGVCFsGATK() or just CallVCFsGATK\n\
        method on this SAMs instance.'

        assert hasattr(self, 'path_to_unfiltered_VCF'), e1

        ## filtering must be done differently if
        ## a single joint called VCF is present
        ## or a set of VCFs are present

        if isinstance(self.path_to_unfiltered_VCF[-1],str) or \
                isinstance(self.path_to_unfiltered_VCF[-1],unicode):
            # single item to joint called VCF
            # extract the SNPs
            raw_SNPs = self.path_to_unfiltered_VCF[-1][:-4] + '_SNPs.vcf'
            cmd = [use_java, '-jar', jar,
                '-T', 'SelectVariants',
                '-R', genome_fna,
                '-V', self.path_to_unfiltered_VCF[-1],
                #'-L', '20',
                '-selectType', 'SNP',
                '-o', raw_SNPs]
            print(' '.join(cmd))
            _subprocess.call(cmd)
            # filter the SNPs
            hf_SNPs = (self.path_to_unfiltered_VCF[-1][:-4] + '_SNPs.vcf').replace('unfiltered','hardfiltered')
            cmd = [use_java, '-jar', jar,
                '-T', 'VariantFiltration',
                '-R', genome_fna,
                '-V', raw_SNPs,
                '--filterExpression', 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0',
                '--filterName', 'standard_hard_filter',
                '-o', hf_SNPs]
            print(' '.join(cmd))
            _subprocess.call(cmd)
        else:
            # must be a list
            hf_SNPs = []
            for unfilteredVCF in self.path_to_unfiltered_VCF:
                # single item to joint called VCF
                # extract the SNPs
                raw_SNPs = unfilteredVCF[:-4] + '_SNPs.vcf'
                cmd = [use_java, '-jar', jar,
                    '-T', 'SelectVariants',
                    '-R', genome_fna,
                    '-V', unfilteredVCF,
                    #'-L', '20',
                    '-selectType', 'SNP',
                    '-o', raw_SNPs]
                print(' '.join(cmd))
                _subprocess.call(cmd)
                # filter the SNPs
                this_hf_SNPs = (unfilteredVCF[:-4] + '_SNPs.vcf').replace('unfiltered','hardfiltered')
                cmd = [use_java, '-jar', jar,
                    '-T', 'VariantFiltration',
                    '-R', genome_fna,
                    '-V', raw_SNPs,
                    '--filterExpression', 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0',
                    '--filterName', 'standard_hard_filter',
                    '-o', this_hf_SNPs]
                hf_SNPs += [this_hf_SNPs]
                print(' '.join(cmd))
                _subprocess.call(cmd)

        # add to a list because this is done twice
        if hasattr(self, 'path_to_hardfiltered_SNPs'):
            self.path_to_hardfiltered_SNPs += [hf_SNPs]
        else:
            self.path_to_hardfiltered_SNPs = [hf_SNPs]

    def hardfilterINDELsGATK(self, 
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            use_java = 'java',
            force = False):
        jar = _os.path.sep.join(jar)
        genome_fna = 'genome_sequences/%s.fna' % self.genome_id
        if not _os.path.exists(genome_fna):
            _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                    genome_fna, 
                    'fasta')

        e1 = 'Could not find "path_to_unfiltered_VCF" attribute. \
        Before filtering, joint calling of variants is necessary. \
        Please run:\n\
        CallgVCFsGATK and GenotypeGVCFsGATK() or just CallVCFsGATK\n\
        method on this SAMs instance.'

        assert hasattr(self, 'path_to_unfiltered_VCF'), e1

        if isinstance(self.path_to_unfiltered_VCF[-1],str) or \
                isinstance(self.path_to_unfiltered_VCF[-1],unicode):
            # single item to joint called VCF
            # extract the INDELs
            raw_INDELs = self.path_to_unfiltered_VCF[-1][:-4] + '_INDELs.vcf'
            cmd = [use_java, '-jar', jar,
                '-T', 'SelectVariants',
                '-R', genome_fna,
                '-V', self.path_to_unfiltered_VCF[-1],
                #'-L', '20',
                '-selectType', 'INDEL',
                '-o', raw_INDELs]
            print(' '.join(cmd))
            _subprocess.call(cmd)
            # filter the INDELs
            hf_INDELs = (self.path_to_unfiltered_VCF[-1][:-4] + '_INDELs.vcf').replace('unfiltered','hardfiltered')
            cmd = [use_java, '-jar', jar,
                '-T', 'VariantFiltration',
                '-R', genome_fna,
                '-V', raw_INDELs,
                '--filterExpression', 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0',
                '--filterName', 'standard_indel_hard_filter',
                '-o', hf_INDELs]
            print(' '.join(cmd))
            _subprocess.call(cmd)
        else:
            # must be a list
            hf_INDELs = []
            for unfilteredVCF in self.path_to_unfiltered_VCF:
                # single item to joint called VCF
                # extract the INDELs
                raw_INDELs = unfilteredVCF[:-4] + '_INDELs.vcf'
                cmd = [use_java, '-jar', jar,
                    '-T', 'SelectVariants',
                    '-R', genome_fna,
                    '-V', unfilteredVCF,
                    #'-L', '20',
                    '-selectType', 'INDEL',
                    '-o', raw_INDELs]
                print(' '.join(cmd))
                _subprocess.call(cmd)
                # filter the INDELs
                this_hf_INDELs = (unfilteredVCF[:-4] + '_INDELs.vcf').replace('unfiltered','hardfiltered')
                cmd = [use_java, '-jar', jar,
                    '-T', 'VariantFiltration',
                    '-R', genome_fna,
                    '-V', raw_INDELs,
                '--filterExpression', 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0',
                '--filterName', 'standard_indel_hard_filter',
                    '-o', this_hf_INDELs]
                hf_INDELs += [this_hf_INDELs]
                print(' '.join(cmd))
                _subprocess.call(cmd)


        # add to a list because this is done twice
        if hasattr(self, 'path_to_hardfiltered_INDELs'):
            self.path_to_hardfiltered_INDELs += [hf_INDELs]
        else:
            self.path_to_hardfiltered_INDELs = [hf_INDELs]

    def recalibBaseScoresGATK(self, 
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            samtools_exe = ['external_programs', 'samtools', 'samtools'],
            local_variants_path = ['variants'],
            use_java = 'java',
            force = False,
            mem_num_gigs = 8, 
            max_cpus = -1):
        '''
        https://www.broadinstitute.org/gatk/guide/best-practices/?bpm=DNAseq
        max_cpus for this GATK module is "cpu threads per data thread"
        '''
        jar = _os.path.sep.join(jar)
        samtools_exe = _os.path.sep.join(samtools_exe)
        genome_fna = 'genome_sequences/%s.fna' % self.genome_id
        if not _os.path.exists(genome_fna):
            _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                    genome_fna, 
                    'fasta')

        local_variants_path = _os.path.sep.join(local_variants_path)
        if not _os.path.exists(local_variants_path):
            _os.makedirs(local_variants_path)

        local_variants_path_genome = _os.path.sep.join([
                                        local_variants_path,
                                        self.genome_id])

        if not _os.path.exists(local_variants_path_genome):
            _os.makedirs(local_variants_path_genome)

        paths_to_recalibrated_BAMs = []

        max_processes = _decide_max_processes( max_cpus )

        start_time = _time.time()

        for cnum,BAM in enumerate(self.ready_BAMs[-1]):
            table_out_pre = BAM[:-4] + '_baserecal_pre.table'
            if isinstance(self.path_to_hardfiltered_SNPs[-1],str) or \
                    isinstance(self.path_to_hardfiltered_SNPs[-1],unicode):
                # joint calling was used
                knownSitesVCF = self.path_to_hardfiltered_SNPs[-1]
            else:
                # per sample single calling was used
                knownSitesVCF = self.path_to_hardfiltered_SNPs[-1][cnum]
            
            if not _os.path.exists(table_out_pre) or force:
                cmd = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar,
                    '-T', 'BaseRecalibrator',
                    '-R', genome_fna,
                    '-I', BAM,
                    #'-L', '20',
                    '-nct',  str(max_processes),
                    '-knownSites', knownSitesVCF,
                    #'--validation_strictness', 'LENIENT',
                    '-o', table_out_pre]
                
                print('Called: %s' % (' '.join(map(str, cmd))))
                _subprocess.call(cmd)
            else:
                print('Found:')
                print(table_out_pre)
                print('use "force = True" to overwrite')
            
            table_out_post = BAM[:-4] + '_baserecal_post.table'
            if not _os.path.exists(table_out_post) or force:
                cmd = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar,
                    '-T', 'BaseRecalibrator',
                    '-R', genome_fna,
                    '-I', BAM,
                    #'-L', '20',
                    '-nct',  str(max_processes),
                    '-knownSites', knownSitesVCF,
                    '-BQSR', table_out_pre,
                    '-o', table_out_post]
                
                print('Called: %s' % (' '.join(map(str, cmd))))
                _subprocess.call(cmd)
                
            else:
                print('Found:')
                print(table_out_post)
                print('use "force = True" to overwrite')
            
            BAM_out = BAM[:-4] + '_baserecal.bam'
            if not _os.path.exists(BAM_out) or force:
                cmd = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar,
                    '-T', 'PrintReads',
                    '-R', genome_fna,
                    '-I', BAM,
                    '-nct',  str(max_processes),
                    '-BQSR', table_out_post,
                    '-o', BAM_out]
                
                print('Called: %s' % (' '.join(map(str, cmd))))
                _subprocess.call(cmd)
                
            else:
                print('Found:')
                print(BAM_out)
                print('use "force = True" to overwrite')
            
            cmd = [samtools_exe, 'index', BAM_out]        
            _subprocess.call(cmd)
            paths_to_recalibrated_BAMs += [BAM_out]
            
            # report durations, time left etc
            _report_time(start_time, cnum, len(self.ready_BAMs[-1]))

        # the last list of BAMs in ready_BAMs is input for CallgVCFsGATK
        # both IndelRealignGATK and recalibBaseScoresGATK put here
        self.ready_BAMs += [paths_to_recalibrated_BAMs]

class CallerDiscoSNP:
    '''
    Wrapper around DiscoSNP++ for SNP and InDel calling from short reads.

    Requires a collection of short read datasets. Alignment and reference-free.
    Can produce a VCF file if a reference genome provided.
    '''
    def __init__(self, reads = False, 
                       genome = False, 
                       path_to_baga = False):
        '''
        Initialise with:
        a list of one or more baga.PrepareReads.Reads object and optionally,
        a baga.CollectData.Genome object.
        
        OR
        
        a path to baga.CallVariants.CallerDiscoSNP object (like this one) that 
        was previously saved.
        '''
        assert reads or path_to_baga, 'Instantiate with reads or a previously saved CallerDiscoSNP'
        assert not (reads and path_to_baga), 'Instantiate with reads or a previously saved CallerDiscoSNP!'
        
        if reads:
            reads_paths = []
            for read_group in reads:
                reads_paths += [{}]
                for attribute_name in dir(read_group):
                    if 'read_files' in attribute_name:
                        reads_paths[-1][attribute_name] = getattr(read_group, attribute_name)
            
            self.reads = reads_paths
            if genome:
                self.genome_sequence = genome.sequence
                self.genome_id = genome.id
        elif path_to_baga:
            with _tarfile.open(path_to_baga, "r:gz") as tar:
                for member in tar:
                    contents = _StringIO(tar.extractfile(member).read())
                    try:
                        # either json serialised conventional objects
                        contents = _json.loads(contents.getvalue())
                    except ValueError:
                        #print('json failed: {}'.format(member.name))
                        # or longer python array.array objects
                        contents = _array('c', contents.getvalue())
                    
                    setattr(self, member.name, contents)


    def saveLocal(self, name):
        '''
        Save processed object info to a local compressed archive of json strings.
        'name' should exclude extension: .baga will be added
        '''
        fileout = 'baga.CallVariants.CallerDiscoSNP-%s.baga' % name
        with _tarfile.open(fileout, "w:gz") as tar:
            print('Writing to {} . . . '.format(fileout))
            for att_name, att in self.__dict__.items():
                if isinstance(att, _array):
                    io = _StringIO(att.tostring())
                    io.seek(0, _os.SEEK_END)
                    length = io.tell()
                    io.seek(0)
                    thisone = _tarfile.TarInfo(name = att_name)
                    thisone.size = length
                    tar.addfile(tarinfo = thisone, fileobj = io)
                else:
                    # try saving everything else here by jsoning
                    try:
                        io = _StringIO()
                        _json.dump(att, io)
                        io.seek(0, _os.SEEK_END)
                        length = io.tell()
                        io.seek(0)
                        thisone = _tarfile.TarInfo(name = att_name)
                        thisone.size = length
                        tar.addfile(tarinfo = thisone, fileobj = io)
                    except TypeError:
                        # ignore non-jsonable things like functions
                        # include unicodes, strings, lists etc etc
                        #print('omitting {}'.format(att_name))
                        pass

    def call(self, 
            local_variants_path = ['variants_DiscoSNP'],
            path_to_bwa = False,
            use_existing_graph = False,
            add_prefix = False,
            force = False,
            max_cpus = -1,
            arguments = False):
        '''
        Call SNPs and InDels using kissnp2 and kissreads2 of DiscoSNP++
        '''

        path_to_exe = _get_exe_path('discosnp')

        local_variants_path = _os.path.sep.join(local_variants_path)
        if not _os.path.exists(local_variants_path):
            _os.makedirs(local_variants_path)

        if hasattr(self, 'genome_id'):
            use_genome_name = self.genome_id
        else:
            use_genome_name = 'no_reference'

        local_variants_path_genome = _os.path.sep.join([
                                        local_variants_path,
                                        use_genome_name])

        if not _os.path.exists(local_variants_path_genome):
            _os.makedirs(local_variants_path_genome)

        max_processes = _decide_max_processes( max_cpus )

        # select read files to use

        use_reads = {}
        for n,reads_group in enumerate(self.reads):
            if 'trimmed_read_files' in reads_group:
                for pairname,pair in reads_group['trimmed_read_files'].items():
                    use_reads[pairname] = pair
            elif 'adaptorcut_read_files' in reads_group:
                print('WARNING: "trimmed_read_files" attribute not found in '\
                        'baga.PreparedReads.Reads object {}'.format(n+1))
                for pairname,pair in reads_group['adaptorcut_read_files'].items():
                    use_reads[pairname] = pair
                print('Using read files that have not been trimmed by position score')
            else:
                print('WARNING: "trimmed_read_files" nor "adaptorcut_read_files" attribute '\
                        'not found in baga.PreparedReads.Reads object {}'.format(n+1))
                assert 'read_files' in reads_group, 'Try recreating '\
                        'baga.PreparedReads.Reads . . . could not find '\
                        'any read files.'
                for pairname,pair in reads_group['read_files'].items():
                    use_reads[pairname] = pair
                print('Using read files that have not been trimmed by position score nor had '\
                        'potential library preparation artifacts removed ("adaptor" sequences '\
                        'etc) {}'.format(n+1))


        pair_file_names = []
        for sample_name, pair in use_reads.items():
            this_file_name = _os.path.sep.join([
                    local_variants_path_genome, 
                    'readspair_paths_for_{}.txt'.format(sample_name)])
            pair_file_names += [this_file_name]
            with open(this_file_name, 'w') as fout:
                fout.write('{}\n{}\n'.format(
                        _os.path.abspath(pair[1]),
                        _os.path.abspath(pair[2])))

        this_file_name = _os.path.sep.join([local_variants_path_genome, 
                'readpairs_for_DiscoSNP.txt'])

        with open(this_file_name, 'w') as fout:
            for pair_file_name in pair_file_names:
                fout.write('{}\n'.format(_os.path.abspath(pair_file_name)))

        try:
            # previous failed runs can leave this file which causes problems
            _os.unlink(this_file_name+'_removemeplease')
        except OSError:
            pass

        # cmd should be built as 'option':[argument list] dictionary
        # with None as values for flag options
        cmd = {}
        cmd['-r'] = [this_file_name]
        cmd['-T'] = None

        if hasattr(self, 'genome_sequence') and hasattr(self, 'genome_id'):
            print('Will map DiscoSNP++ variants to {} ({:,} bp)'.format(
                    self.genome_id, len(self.genome_sequence)))
            if not path_to_bwa:
                path_to_bwa = _get_exe_path('bwa')
            if _os.path.isfile(path_to_bwa):
                # VCF maker only wants path to bwa exe, not the exe itself
                path_to_bwa = path_to_bwa[:-4]
            genome_fna = 'genome_sequences/%s.fna' % self.genome_id
            if not _os.path.exists(genome_fna):
                _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                        genome_fna, 
                        'fasta')
            
            cmd['-G'] = [genome_fna]
            cmd['-B'] = [path_to_bwa]

        if use_existing_graph:
            print('Will attempt to use existing graph: be sure your input reads match the graph!')
            # discoRes_k_31_c_auto.h5 for default settings
            cmd['-g'] = None

        if add_prefix:
            cmd['-p'] = [add_prefix+'__DiscoSNP']

        if arguments:
            # overwrite defaults with direct arguments
            # (e.g. via -A/--arguments cli)
            from baga import parse_new_arguments
            cmd = parse_new_arguments(arguments, cmd)

        # make commands into a list suitable for subprocess
        cmds = []
        for opt,arg in cmd.items():
            cmds += [opt]
            if arg is not None:
                cmds += arg

        print(' '.join([path_to_exe] + cmds))
        _subprocess.call([path_to_exe] + cmds)


class Summariser:
    '''
    Summarise variants in a VCF in various ways.

    Parse VCF files and create .csv files for viewing in spreadsheets
    '''
    def __init__(self, VCF_paths = False, path_to_caller = False, genomes = False):
        '''
        A CallVariants.Summariser object can be instantiated with:
            - a list of paths to VCF files
        or
            - path to a CallVariants.CallerGATK object
        or
            - path to a CallVariants.CallerDiscoSNP object
        '''
        
        assert bool(VCF_paths) ^ bool(path_to_caller), 'cannot instantiate with both '\
                'paths to VCFs and a Caller object'
        
        if VCF_paths:
            for VCF in VCF_paths:
                assert _os.path.exists(VCF), 'Could not find {}.\nPlease ensure all '\
                        'files exist'.format(VCF)
            self.VCF_paths = VCF_paths
        elif path_to_caller:
            raise NotImplementedError('only direct paths to VCFs is currently implemented')
        
        # allow genomes to be provided as single object for 1 replicon
        # or list or tuple for one or more
        # or False
        if isinstance(genomes, list) or isinstance(genomes, tuple):
            self.genomes = {genome.id:genome for genome in genomes}
        elif genomes:
            self.genomes = {genomes.id:genomes}
        else:
            self.genomes = genomes
            
            
        


    def collect_variants(self):
        '''
        Collect variants from VCF files for 
        '''

        all_variants = {}
        all_headerdicts = {}

        for VCF in self.VCF_paths: #break
            header, header_section_order, these_colnames, variantrows = parseVCF(VCF)
            if len(variantrows) == 0:
                print('WARNING: no variants found in {}'.format(VCF))
                continue
            headerdict = dictify_vcf_header(header)
            all_headerdicts[VCF] = headerdict
            #print(headerdict)
            variants, allfilters = sortVariantsKeepFilter(header, these_colnames, variantrows)
            #print(variants)
            for sample,chromosomes in variants.items():
                if sample not in all_variants:
                    all_variants[sample] = {}
                for chromosome,positions in chromosomes.items():
                    if chromosome not in all_variants[sample]:
                        all_variants[sample][chromosome] = {}
                    for pos1,info in positions.items():
                        if pos1 in all_variants[sample][chromosome]:
                            if all_variants[sample][chromosome][pos1] != info:
                                print('WARNING: VCFs in conflict for {} in chromosome '\
                                        '{} at position {}. You could try processing '\
                                        'VCFs one at a time and compare tables produced. '\
                                        'A variant at this position was omitted.'\
                                        ''.format(sample,chromosome,pos1))
                        else:
                            all_variants[sample][chromosome][pos1] = info

        return(all_variants, all_headerdicts)


    def simple(self, nodata_char = '?'):
        '''
        List all variants with rows corresponding to those in the VCF file(s) in a .csv file
        '''

        all_variants, all_headerdicts = self.collect_variants()

        if len(all_variants) == 0:
            print('WARNING: no variants found the supplied VCF files!')
            return()

        # ['ALT', 'FILTER', 'FORMAT', 'GATKCommandLine', 'INFO', 'contig']
        genome_IDs = set()
        for VCF,headerdict in all_headerdicts.items():
            try:
                genome_IDs.update([contig['ID'] for contig in headerdict['contig']])
            except KeyError:
                print('WARNING: no contig information found in header of {}. Reference free?'.format(VCF))

        if len(genome_IDs):
            print('Chromosomes found: {}'.format(', '.join(genome_IDs)))
        else:
            print('No chromsome information found in: {}'.format(', '.join(sorted(all_headerdicts))))


        # assert len(genome_IDs) == 1, "only one reference sequence at a time is "\
                # "implemented for simple summary. {} found: {}".format(len(genome_IDs), 
                # ', '.join(genome_IDs))

        if self.genomes:
            annotate_with_these = set(self.genomes) & genome_IDs
            requested_not_found = set(self.genomes) - genome_IDs
            found_not_requested = genome_IDs - set(self.genomes)
            if len(requested_not_found) == 1:
                print('WARNING: {} was requested for use with annotations, but not found in VCFs'\
                        ''.format(sorted(requested_not_found)[0]))
            elif len(requested_not_found) > 1:
                print('WARNING: {} were requested for use with annotations, but not found in VCFs'\
                        ''.format(', '.join(sorted(requested_not_found))))
            if len(found_not_requested) == 1:
                print('WARNING: {} was found in VCFs, but not provided for use with annotations'\
                        ''.format(sorted(found_not_requested)[0]))
            elif len(found_not_requested) > 1:
                print('WARNING: {} were found in VCFs, but not provided for use with annotations'\
                        ''.format(', '.join(sorted(found_not_requested))))
            if len(genome_IDs) == 0:
                annotate_with_these = set(self.genomes)
                print('WARNING: Assuming VCFs were called reference free because no chromosome '\
                        'information found in VCFs. Cannot check provided chromosomes match '\
                        'variants VCF: please double check annotations!')
            if len(annotate_with_these):
                print('Using {} for annotations'.format(', '.join(sorted(annotate_with_these))))


        # collect by chromsome,position
        by_position = _defaultdict(dict)
        by_position_nodata = _defaultdict(dict)
        by_position_freqs = _defaultdict(dict)
        annotations = {}
        for sample,chromosomes in all_variants.items():
            for chromosome,positions in chromosomes.items():
                for pos1,((r,q),filters) in positions.items():
                    ## this was to keep sortVariantsKeepFilter()
                    ## compatible with other uses that assume single clone, always 1/1
                    if isinstance(q, tuple):
                        ALT,freq,total = q
                        use_q = ALT
                    else:
                        freq,total = 1,1
                        use_q = q
                    
                    if use_q == '.':
                        # . in VCF is insufficient data to call:
                        # the data is missing (unknown)
                        # use desired "no data" character
                        by_position_nodata[chromosome,pos1][sample] = nodata_char
                        # no need to attempt codon prediction or
                        # add a row as a real variant would below
                        continue
                    
                    try:
                        by_position[chromosome,pos1][(r,use_q)][sample] = filters
                        by_position_freqs[chromosome,pos1][(r,use_q)][sample] = freq,total
                    except KeyError:
                        by_position[chromosome,pos1][(r,use_q)] = {}
                        by_position[chromosome,pos1][(r,use_q)][sample] = filters
                        by_position_freqs[chromosome,pos1][(r,use_q)] = {}
                        by_position_freqs[chromosome,pos1][(r,use_q)][sample] = freq,total
                    
                    try:
                        ORFs_info = [(ORF_id,s,e,st,gene_name) for ORF_id,(s,e,st,gene_name) in \
                                sorted(self.genomes[chromosome].ORF_ranges.items()) if s < pos1 < e]
                        if len(ORFs_info) > 1:
                            print('WARNING: variant in more than one ORF (overlapping). '\
                                    'Detailed annotations not yet implemented (marked as "multi")')
                            multi_ORF_IDs = ','.join([ORF_id for ORF_id,s,e,st,gene_name in ORFs_info])
                            multi_ORF_names = ','.join([gene_name for ORF_id,s,e,st,gene_name in ORFs_info])
                            annotations[(chromosome,pos1,r,use_q)] = ("multi", "multi", 
                                    "multi", "multi", "multi", multi_ORF_IDs, multi_ORF_names, "multi")
                        elif len(ORFs_info) == 1:
                            ORF_id,s,e,strand,gene_name = ORFs_info[0]
                            ORF_seq = self.genomes[chromosome].sequence[s:e].tostring()
                            ORF0 = pos1 - 1 - s
                            ORF0_codon_start = ORF0 - ORF0 % 3
                            ref_codon = ORF_seq[ORF0_codon_start:ORF0_codon_start+3]
                            frame1 = ORF0 - ORF0_codon_start + 1
                            if len(r) == len(use_q) == 1:
                                # substitution
                                var_codon = list(ref_codon)
                                var_codon[frame1-1] = use_q
                                var_codon = ''.join(var_codon)
                                assert self.genomes[chromosome].sequence[pos1-1] == ref_codon[frame1-1]
                                if st == 1:
                                    var_AA = str(_Seq(var_codon).translate())
                                    ref_AA = str(_Seq(ref_codon).translate())
                                else:
                                    var_AA = str(_Seq(var_codon).reverse_complement().translate())
                                    ref_AA = str(_Seq(ref_codon).reverse_complement().translate())
                                    var_codon = str(_Seq(var_codon).reverse_complement())
                                    ref_codon = str(_Seq(ref_codon).reverse_complement())
                                    frame1 = 3 - (frame1-1)
                            else:
                                # indel
                                var_codon = '-'
                                var_AA = '-'
                                if st == -1:
                                    frame1 = 3 - (frame1-1)
                                    ref_AA = str(_Seq(ref_codon).reverse_complement().translate())
                                    ref_codon = str(_Seq(ref_codon).reverse_complement())
                                else:
                                    ref_AA = str(_Seq(ref_codon).translate())
                            
                            annotations[(chromosome,pos1,r,use_q)] = (ref_codon, var_codon, 
                                    ref_AA, var_AA, frame1, ORF_id, gene_name, strand)
                    except (KeyError, TypeError):
                        pass
                    except _TranslationError:
                        print('WARNING: problem encountered translating a DNA sequence '\
                                'causing a variant to be omitted from the summary. '\
                                'This may be caused by unexpected characters including "*" '\
                                'sometimes introduced by variant callers.')
                        print('Problem variant was: "{}" at position {} in seq. ID "{}", '\
                                'sample "{}"'.format(q, pos1, chromosome, sample))
                        pass

        # allow KeyError that defaultdict does not raise
        by_position_nodata = dict(by_position_nodata)

        all_chromosomes = sorted(set([a for b in all_variants.values() for a in b]))

        filenameout = 'Simple_summary_for_{}_and_{}_others__{}.csv'.format(sorted(all_variants)[0], 
                len(all_variants)-1, '_'.join(all_chromosomes))

        sample_order = sorted(all_variants)

        def quote(v):
            if str(v).replace('-','').replace('.','').isdigit():
                return(str(v))
            else:
                return('"'+v+'"')

        colnames = ['Chromosome','Position (bp)', 'Reference', \
                'Variant', 'Gene ID', 'Gene name', 'Strand', 'Reference codon', 'Variant codon', \
                'Reference AA', 'Variant AA', 'Codon position'] + sample_order + \
                ['Frequency', 'Frequency (filtered)']

        with open(filenameout, 'w') as fout:
            print('Writing to {}'.format(filenameout))
            fout.write(','.join(map(quote,colnames))+'\n')
            for (chromosome,pos1),variants in sorted(by_position.items()):
                for (r,q),samples in variants.items():
                    try:
                        (ref_codon, var_codon, ref_AA, var_AA, frame1, ORF_id, \
                                gene_name, strand) = annotations[(chromosome,pos1,r,q)]
                    except KeyError:
                        (ref_codon, var_codon, ref_AA, var_AA, frame1, ORF_id, \
                                gene_name, strand) = '', '', '', '', -1, '', '', 0
                    
                    this_row = [chromosome,pos1,r,q] + [ORF_id, gene_name, strand, ref_codon, \
                            var_codon, ref_AA, var_AA, frame1]
                    freq_all = 0
                    freq_filtered = 0
                    for sample in sample_order:
                        try:
                            this_row += [by_position_nodata[chromosome,pos1][sample]]
                        except KeyError:
                            # not no data ambiguity so get freq
                            try:
                                freq,total = by_position_freqs[chromosome,pos1][(r,q)][sample]
                                if freq == total == 1:
                                    val = 1
                                else:
                                    val = float(freq) / float(total)
                                
                                filters = samples[sample]
                                # what does the {'.'} here mean for filters?
                                if filters == set(['.']) or \
                                        filters == set(['PASS']) or \
                                        len(filters) == 0:
                                    this_row += [val]
                                    freq_filtered += 1
                                else:
                                    this_row += ['+'.join(filters)]
                                freq_all += 1
                            except KeyError:
                                this_row += [0]
                    
                    this_row += [freq_all, freq_filtered]
                    fout.write(','.join(map(quote,this_row))+'\n')

class Checker:
    '''
    Check variants in a VCF against regional de novo assemblies.

    Extract short reads from a BAM aligned around each variant, de novo 
    assemble using SPAdes, then optimal global align contigs back to those 
    regions and check for variants.
    '''
    def __init__(self, VCF_paths, BAM_paths, genome):
        '''
        A CallVariants.Checker object must be instantiated with:
            - a list of paths to VCF files
            - a CollectData genome object
        '''
        
        e = 'Could not find %s.\nPlease ensure all files exist'
        for VCF in VCF_paths:
            assert _os.path.exists(VCF), e % VCF
        
        for BAM in BAM_paths:
            assert _os.path.exists(VCF), e % BAM
        
        self.VCF_paths = VCF_paths
        self.BAM_paths = BAM_paths
        self.genome = genome

    def saveLocal(self, name):
        '''
        Save processed object info to a local compressed archive of json strings.
        'name' should exclude extension: .baga will be added
        '''
        fileout = 'baga.CallVariants.Checker-%s.baga' % name
        with _tarfile.open(fileout, "w:gz") as tar:
            print('Writing to {} . . . '.format(fileout))
            for att_name, att in self.__dict__.items():
                if isinstance(att, _array):
                    io = _StringIO(att.tostring())
                    io.seek(0, _os.SEEK_END)
                    length = io.tell()
                    io.seek(0)
                    thisone = _tarfile.TarInfo(name = att_name)
                    thisone.size = length
                    tar.addfile(tarinfo = thisone, fileobj = io)
                else:
                    # try saving everything else here by jsoning
                    try:
                        io = _StringIO()
                        _json.dump(att, io)
                        io.seek(0, _os.SEEK_END)
                        length = io.tell()
                        io.seek(0)
                        thisone = _tarfile.TarInfo(name = att_name)
                        thisone.size = length
                        tar.addfile(tarinfo = thisone, fileobj = io)
                    except TypeError:
                        # ignore non-jsonable things like functions
                        # include unicodes, strings, lists etc etc
                        #print('omitting {}'.format(att_name))
                        pass

    def collect_variants(self):
        '''
        Collect variants from VCF files for 
        '''

        all_variants = {}

        for VCF in self.VCF_paths: #break
            header, header_section_order, these_colnames, variantrows = parseVCF(VCF)
            headerdict = dictify_vcf_header(header)
            #print(headerdict)
            variants, allfilters = sortVariantsKeepFilter(header, these_colnames, variantrows)
            #print(variants)
            all_variants = dict(all_variants.items() + variants.items())

        return(all_variants)
            

    def doCheck(self, num_padding = 2000, max_memory = False, force = False):
        '''
        De novo assemble variant regions and compare with variant calls
        '''
        # VCFs and genome provided at instantiation
        # do some checks on provided genome and VCFs
        genome_ids = {}
        genome_lengths = {}
        pattern = _re.compile('##contig=<ID=([A-Za-z0-9\._]+),length=([0-9]+)>')
        for VCF in self.VCF_paths: #break
            for line in open(VCF):
                
                if line[:9] == '##contig=':
                    try:
                        genome_id, genome_length = _re.match(pattern, line).groups()
                    except AttributeError:
                        print('Failed to parse genome information from {}'.format(line))
                    genome_lengths[int(genome_length)] = VCF
                    genome_ids[genome_id] = VCF
                    #print(genome_id, genome_length)
                    identified = True
                    break
                
            e = "Failed to identify which chromosome the variants in {} were called on (couldn't find '##contig=')".format(VCF)
            assert identified, e

        e = 'Differing reference genome among provided VCFs? {}'.format(genome_ids.items())
        assert len(genome_ids) == 1, e

        e = 'Genome provided ({}) does not match genome ID in provided VCFs: {}'.format(self.genome.id, genome_ids.keys()[0])
        assert self.genome.id == genome_ids.keys()[0], e
        print('Variants were called against {:,} bp genome: {}\n'.format(int(genome_length), genome_id))


        all_variants = self.collect_variants()

        BAMs_by_ids = {}
        for BAM in self.BAM_paths:
            header = _pysam.Samfile(BAM, 'rb').header
            BAMs_by_ids[(header['RG'][0]['ID'],header['SQ'][0]['SN'])] = BAM
            assert header['SQ'][0]['SN'] == self.genome.id, 'mismatch between '\
                    'BAM genome ({}) and genome used by BAGA ({})'.format(
                    header['SQ'][0]['SN'], self.genome.id)

        use_samples = set(all_variants) & set([sample for sample,chromosome in BAMs_by_ids])

        print('Found {} samples common to supplied VCFs and BAMs:\n{}'.format(
                len(use_samples),', '.join(use_samples)))


        from baga import Structure

        try:
            _os.mkdir('variant_checks')
        except OSError:
            pass

        path_to_variant_checks = ['variant_checks', self.genome.id]
        # inconsistant requirements for path (list or str) below
        path_to_variant_checks_str = _os.path.sep.join(path_to_variant_checks)

        try:
            _os.mkdir(path_to_variant_checks_str)
        except OSError:
            pass

        import baga

        denovo_info = {}
        all_contigs = {}
        for sample in use_samples:
            denovo_info[sample] = {}
            all_contigs[sample] = {}
            chromosomes = all_variants[sample]
            for chromosome,variants in chromosomes.items():
                collector = Structure.Collector(BAMs_by_ids[(sample, chromosome)])
                single_assembly = False
                R1 = 'variant_checks/{}/{}__{}_unmapped_R1.fastq'.format(chromosome,sample,chromosome)
                R2 = 'variant_checks/{}/{}__{}_unmapped_R2.fastq'.format(chromosome,sample,chromosome)
                RS = 'variant_checks/{}/{}__{}_unmapped_S.fastq'.format(chromosome,sample,chromosome)
                if _os.path.exists(R1) and _os.path.exists(R2) and \
                        _os.path.getsize(R1) > 0 and _os.path.getsize(R2) > 0 \
                        and not force:
                    print('Found unmapped reads at {} and {}\nUse --force/-F to overwrite'.format(R1,R2))
                    r1_out_path_um, r2_out_path_um, rS_out_path_um = R1, R2, RS
                else:
                    print('Extracting poorly and unaligned reads for sample {}'.format(sample))
                    collector.getUnmapped()
                    r1_out_path_um, r2_out_path_um, rS_out_path_um = \
                            collector.writeUnmapped(path_to_variant_checks_str)
                import AssembleReads
                if max_memory:
                    use_mem_gigs = max_memory
                else:
                    # round down available GBs
                    use_mem_gigs = int(baga.get_available_memory())
                    # unless to zero!
                    if use_mem_gigs == 0:
                        use_mem_gigs = 1
                
                reads_path_unmapped = {}
                output_folder_um = '_'.join(
                        r1_out_path_um.split('_')[:-1]).split(_os.path.sep)[-1]
                reads_path_unmapped[output_folder_um] = (r1_out_path_um, 
                        r2_out_path_um, rS_out_path_um)
                path_to_bad_unmapped_contigs = _os.path.sep.join([
                        path_to_variant_checks_str, output_folder_um, 
                        'contigs.fasta'])
                if _os.path.exists(path_to_bad_unmapped_contigs) and \
                        _os.path.getsize(path_to_bad_unmapped_contigs) > 0 and \
                        not force:
                    print('Found assembly at {}\nUse --force/-F to overwrite. '\
                            'Skipping . . .'.format(path_to_bad_unmapped_contigs))
                else:
                    if not force:
                        # if args.force specified, don't need to say anything
                        # either way, do assembly
                        print('Nothing found at {}. Doing assembly.'.format(
                                path_to_bad_unmapped_contigs))
                    
                    reads = AssembleReads.DeNovo(paths_to_reads = reads_path_unmapped)
                    reads.SPAdes(output_folder = path_to_variant_checks, 
                            mem_num_gigs = use_mem_gigs, only_assembler = True, 
                            careful = False)
                # assemble read from each region with poorly/unmapped
                reads_paths = {}
                # make a second dict of reads for assembly, all values for unmapped reads
                # that need to be included in each assembly
                # first, collect reads, make fastqs, recording in a dict
                reads_path_unmapped = {}
                assemblies_by_variant = {}
                
                # join regions with multiple variants to avoid redundacy in de novo assemblies
                join_dist = 10000
                position_regions = []
                allpositions = sorted(variants)
                this_region = [allpositions[0]]
                for pn,pos1 in enumerate(allpositions[:-1]):
                    if pos1 + join_dist > allpositions[pn+1]:
                        this_region += [allpositions[pn+1]]
                    else:
                        position_regions += [this_region]
                        this_region = [allpositions[pn+1]]
                
                position_regions += [this_region]
                regions_for_de_novo = []
                for these_positions in position_regions:
                    regions_for_de_novo += [(these_positions[0]-500,these_positions[-1]+500)]
                    if regions_for_de_novo[-1][0] < 0:
                        regions_for_de_novo[-1][0] = 0
                    if regions_for_de_novo[-1][1] > len(self.genome.sequence):
                        ## this bit not compatible with multiple chromosomes
                        regions_for_de_novo[-1][1] = len(self.genome.sequence)
                
                # ensure no very long contigs else PW alignment will not be possible
                # longer regions more likely for samples more divergent from reference
                # use ajoining regions with small overlaps instead
                regions_for_de_novo2 = []
                # this seems to be maximum length before seq-align seg-faults
                # about 16GB memory required for these long alignments
                max_len = 20000
                for s,e in regions_for_de_novo:
                    if e - s > max_len:
                        num_pieces = ((e-s)//float(max_len)+1)
                        newlen = int((e-s) / num_pieces)
                        for i in range(int(num_pieces)):
                            # make internal join overlap incase variant close to a join
                            pre_pad = 0
                            end_pad = 0
                            if i > 0:
                                pre_pad = 100
                            if i < num_pieces:
                                end_pad = 100
                            new_s,new_e = s+(newlen*i)-pre_pad,s+(newlen*(i+1)+end_pad)
                            regions_for_de_novo2 += [[new_s,new_e]]
                    else:
                        regions_for_de_novo2 += [[s,e]]
                
                # print('allpositions',len(allpositions))
                # print('regions_for_de_novo',len(regions_for_de_novo))
                # print('regions_for_de_novo2',len(regions_for_de_novo2))
                num_padding = 0
                #for pos1,info in sorted(variants.items()):
                    # print('Extracting reads aligned near variant at {} in sample {}'\
                            # ''.format(pos1, sample))
                    # collector.makeCollection(pos1 - 1, pos1, num_padding)
                for s,e in regions_for_de_novo2:
                    print('Extracting reads aligned over variant from {} to {} in sample {}'\
                            ''.format(s, e, sample))
                    collector.makeCollection(s, e, 0)
                    r1_out_path, r2_out_path, rS_out_path = collector.writeCollection(
                            path_to_variant_checks_str)
                    if not r1_out_path:
                        # if no reads found, False returned
                        print('WARNING: No reads found in region +/- {} positions around '\
                                'variant at {}'.format(num_padding, pos1))
                        continue
                    # put assembly in folder with same name as read files
                    output_folder = '_'.join(r1_out_path.split('_')[:-1]).split(
                            _os.path.sep)[-1]
                    print(output_folder)
                    path_to_contigs = _os.path.sep.join([path_to_variant_checks_str,
                                                        output_folder, 
                                                        'contigs.fasta'])
                    # collect all contigs paths to be assembled and aligning
                    # for summarising below
                    assemblies_by_variant[s, e] = path_to_contigs
                    if _os.path.exists(path_to_contigs) and \
                            _os.path.getsize(path_to_contigs) > 0 and \
                            not force:
                        print('Found assembly at {}\nUse --force/-F to overwrite. '\
                                'Skipping . . .'.format(path_to_contigs))
                    else:
                        # if omitted from this "reads_paths" dict, no assembly done
                        # but aligning still done if in "assemblies_by_variant"
                        reads_paths[output_folder] = (r1_out_path, r2_out_path, 
                                rS_out_path)
                        reads_path_unmapped[output_folder] = (r1_out_path_um, 
                                r2_out_path_um, rS_out_path_um)
                
                # second, run each assembly in single call to
                # AssembleReads.DeNovo.SPAdes
                print('Assemble reads for each region around variants')
                reads = AssembleReads.DeNovo(paths_to_reads = reads_paths, 
                        paths_to_reads2 = reads_path_unmapped)
                reads.SPAdes(output_folder = path_to_variant_checks, 
                        mem_num_gigs = use_mem_gigs, single_assembly = single_assembly, 
                        only_assembler = True, careful = False)
                
                all_contigs[sample][chromosome] = assemblies_by_variant
                # a dict of paths to contigs per region
                aligner = Structure.Aligner(self.genome)
                unmappedfasta = _os.path.sep.join([path_to_variant_checks_str, 
                                                  output_folder_um, 
                                                  'contigs.fasta'])
                if _os.path.exists(unmappedfasta) and _os.path.getsize(unmappedfasta) > 0:
                    # provide dict of range tuples
                    aligner.alignRegions(assemblies_by_variant, num_padding, 
                            path_to_omit_sequences = unmappedfasta, 
                            single_assembly = single_assembly, min_region_length = 0)
                    denovo_info[sample][chromosome] = aligner.reportAlignments()
                else:
                    print('WARNING: no assembled unmapped and poorly mapped reads found at:\n{}'.format(unmappedfasta))
                    try:
                        r1_size = _os.path.getsize(r1_out_path_um)
                        r2_size = _os.path.getsize(r2_out_path_um)
                        print('but reads, {} ({:,} bytes) and {} ({:,} bytes), exist . . check SPAdes assembly log in {}'.format(
                                                                r1_out_path_um,
                                                                r1_size,
                                                                r2_out_path_um,
                                                                r2_size,
                                                                unmappedfasta.replace('contigs.fasta','')))
                    except IOError:
                        print('WARNING: could not find unmapped and poorly '\
                        'aligned reads at:\n{}\n{}\nthis is unexpected but '\
                        'conceivable (if ALL reads really did map to reference!).'.format(
                                r1_out_path_um,r2_out_path_um))
                    print('proceeding with alignment of assembled putatively '\
                            'rearranged regions to reference nonetheless')
                    aligner.alignRegions(assemblies_by_variant, num_padding,
                            single_assembly = single_assembly, min_region_length = 0, 
                            force = False)
                    denovo_info[sample][chromosome] = aligner.reportAlignments()


        num_padding = 0

        # collect all the de novo assembly called variants and compare
        table_outname = 'variant_checks/{}/Table_of_variants_with_de_novo_comparison__{}_and_{}_others.csv'\
                ''.format(chromosome,sorted(use_samples)[0],len(use_samples)-1)
        colnames = ['Sample','Chromosome','Position (bp base-1)','Reference','Variant','Status','Filters','Alignment']
        with open(table_outname,'w') as fout:
            fout.write(','.join(['"{}"'.format(a) for a in colnames])+'\n')
            zeropadding = len(str(len(self.genome.sequence)))
            agree = {}
            missing = {}
            disagree = {}
            for sample in use_samples:
                chromosomes = all_variants[sample]
                agree[sample] = {}
                missing[sample] = {}
                disagree[sample] = {}
                for chromosome,variants in chromosomes.items():
                    these_agree = {}
                    these_missing = {}
                    these_disagree = {}
                    # first collect all variants among do novo alignments
                    all_de_novo_vars = {}
                    vars2alignment = {}
                    for ref_region_id,alnd_contigs in denovo_info[sample][chromosome]['variants_by_contig'].items():
                        alignment_name = 'variant_checks/{1}/{0}__{1}_{2}_multi_alnd.fna'.format(
                                sample,
                                chromosome,
                                ref_region_id[4:])
                        for contig_id,pos0s in alnd_contigs.items():
                            for pos0,(ref_char,alnd_char) in pos0s.items():
                                vars2alignment[pos0+1] = alignment_name
                                try:
                                    all_de_novo_vars[pos0+1].add((ref_char,alnd_char))
                                except KeyError:
                                    all_de_novo_vars[pos0+1] = set([(ref_char,alnd_char)])
                    
                    # then compare to variants and write result to table
                    for pos1,((ref,query),filters) in sorted(variants.items()):
                        if pos1 in all_de_novo_vars:
                            if len(all_de_novo_vars[pos1]) == 1:
                                ref_char,alnd_char = all_de_novo_vars[pos1].pop()
                                #assert ref_char == ref, 'reference character mismatch: {},{} {},{}'.format(ref,query,ref_char,alnd_char)
                                if ref_char != ref:
                                    print('WARNING: reference character mismatch: {},{} {},{}'.format(ref,query,ref_char,alnd_char))
                                    print(pos1,vars2alignment[pos1])
                                if alnd_char == query:
                                    status = "corroborated"
                                    these_agree[pos1] = (ref,query),filters
                                else:
                                    status = "negative"
                                    these_disagree[pos1] = (ref,query),filters
                            else:
                                status = "ambiguous"
                                these_disagree[pos1] = (ref,query),filters
                        else:
                            status = "absent"
                            these_missing[pos1] = (ref,query),filters
                        
                        try:
                            alignment = vars2alignment[pos1]
                        except KeyError:
                            alignment = 'none'
                        
                        fout.write('"{}","{}",{},"{}","{}","{}","{}","{}"\n'.format(
                                sample,
                                chromosome,
                                pos1,
                                ref,
                                query,
                                status,
                                '+'.join(filters),
                                alignment
                                ))
                    
                    agree[sample][chromosome] = these_agree
                    disagree[sample][chromosome] = these_disagree
                    missing[sample][chromosome] = these_missing

        print('Wrote variants and whether thay were found in their respective de '\
                'novo assembled contigs in:\n{}'.format(table_outname))

        # #print('summary')
        for sample,chromosomes in agree.items():
            for chromosome,these_agree in chromosomes.items():
                print('{} aligned against {}: {} SNPs called also in de novo '\
                        'assemblies; {} not.'.format(sample, chromosome, 
                        len(these_agree), len(missing[sample][chromosome])))

        self.agree = agree
        self.disagree = disagree
        self.missing = missing

        # generate name from genome + used samples in sorted order
        # should be unique to this combination without needing a specific --read_group X --genome Y combination
        # currently arbitrary VCFs and BAMs can be provided
        # does not allow use of combined chromosomes
        # could be standardised baga-wide?
        hasher = _md5()
        hasher.update(str([self.genome.id]+sorted(use_samples)))
        unique_name = hasher.hexdigest()

        self.saveLocal(unique_name)



class Filter:
    '''
    Methods to remove variant calls from VCFs according to position specific 
    filters inferred using the Structure and Repeats modules of the Bacterial 
    and Archaeal Genome Analyser.
    '''

    def __init__(self, VCF_paths, genome):  #, reads):
        '''
        A CallVariants.Filter object must be instantiated with:
            - a list of paths to VCF files
            - a CollectData genome object
        '''
        
        e = 'Could not find %s.\nPlease ensure all files exist'
        for VCF in VCF_paths:
            assert _os.path.exists(VCF), e % VCF
        
        self.VCF_paths = VCF_paths
        
        self.ORF_ranges = genome.ORF_ranges
        
        # http://www.1000genomes.org/node/101
        # FILTER filter: PASS if this position has passed all filters, i.e. a call is made at this position.
        # Otherwise, if the site has not passed all filters, a semicolon-separated list of codes for filters
        # that fail. e.g. "q10;s50" might indicate that at this site the quality is below 10 and the number
        # of samples with data is below 50% of the total number of samples. "0" is reserved and should not
        # be used as a filter String.



    def apply_filter_by_ranges(self, variant_rows, ranges, filter_id, sample_index = False):
        '''
        return variant VCF rows changing:
            FILTER cell if sample_index = False (applies to all samples)
            adding an entry to INFO describing which sample filtered if sample_index provided

        sample_index is base-0 index of sample column order
        '''
        new_rows = []
        filtered = {}
        for row in variant_rows:
            cells = row.split('\t')
            pos1 = int(cells[1])
            filter_status = cells[6]
            info_cell = cells[7]
            fail = False
            for s, e in ranges:
                if s < pos1 <= e:
                    # this variant is within region to reject
                    fail = True
                    break
            
            if fail: #break
                filtered[pos1] = cells[3:5]
                if sample_index is not False:
                    # applies just to this sample
                    infos = info_cell.split(';')
                    # ignore entries without an "="
                    infos = dict([i.split('=') for i in infos if '=' in i])
                    if filter_id in infos:
                        # add this sample index to existing of previously filtered samples at this position
                        new_filtered_index_list = map(int,infos[filter_id].split(',')) + [sample_index]
                        # could assert this sample not already filtered at this position . . .
                        # for now silently overwrite (!)
                        infos[filter_id] = ','.join(map(str,sorted(set(new_filtered_index_list))))
                    else:
                        # make new entry for this index
                        infos[filter_id] = str(sample_index)
                    
                    info_cell = ';'.join(['='.join([k,v]) for k,v in sorted(infos.items())])
                    cells[7] = info_cell
                else:
                    # applies to all: change FILTER
                    if filter_id in filter_status:
                        print('WARNING: {} already has {}'.format(cells[0], filter_id))
                    elif filter_status == 'PASS':
                        filter_status = filter_id
                    else:
                        filter_status = ';'.join(filter_status.split(';') + [filter_id])
                    
                    cells[6] = filter_status
                    
                #print(filter_id,pos1,s,e,cells[6],info_cell.split(';')[-1])
                
                new_rows += ['\t'.join(cells)]
            else:
                new_rows += [row]

        return(new_rows, filtered)

    def markVariants(self, filters_to_apply):
        '''
        Given details of one or more filter to apply, write new VCFs with variants marked
        '''
        all_filtered = {}
        for VCF_path in self.VCF_paths:
            all_filtered[VCF_path] = {}
            use_VCF_path = VCF_path
            for filter_id, filter_info in filters_to_apply.items(): #break
                # new VCF with suffix saved for each filter applied
                print('Applying {} filter to variants in:\n{}'.format(filter_id, use_VCF_path))
                header, header_section_order, colnames, variants = parseVCF(use_VCF_path)
                print('{} variant positions found'.format(len(variants)))
                #print('variants',variants)
                sample_order = colnames[9:]
                these_filtered = {}
                # filter_id = 'rearrangements'
                # filter_info = filters_to_apply[filter_id]
                if filter_info['per_sample']:
                    for sample,ranges in sorted(filter_info['ranges'].items()): #break
                        # either vcf per sample or multi-sample vcfs
                        # try all samples in all vcfs to handle either scenario
                        if sample in sample_order:
                            these_filtered[sample] = {}
                            infos_to_add = set()
                            if isinstance(ranges, dict):
                                # sometimes 'extended' versions of filters e.g., no or few reads adjacent to disrupted regions
                                # add as filtername1, filtername2 etc
                                for n,(filter_variant,these_ranges) in enumerate(sorted(ranges.items())): #break
                                    # variants list of rows gets iteratively added to per sample
                                    variants, filtered = self.apply_filter_by_ranges(variants, 
                                                                                these_ranges, 
                                                                                filter_id+str(n+1), 
                                                                                sample_index = sample_order.index(sample))
                                    
                                    # manually check that changes were applied
                                    #[t.split('\t')[:8] for t in variants if int(t.split('\t')[1]) in filtered]
                                    infos_to_add.add(filter_info['string'][n])
                                    these_filtered[sample][filter_id+str(n+1)] = filtered
                                    
                            else:
                                # one set of filter ranges, per sample
                                variants, filtered = self.apply_filter_by_ranges(variants, 
                                                                            ranges, 
                                                                            filter_id, 
                                                                            sample_index = sample_order.index(sample))
                                infos_to_add.add(filter_info['string'])
                                these_filtered[sample] = filtered
                    
                    # add filter info as INFO
                    header['INFO'] += list(infos_to_add)
                else:
                    # just a single reference-genome specific filter to be applied to all samples via the FILTER property
                    variants, filtered = self.apply_filter_by_ranges(variants, filter_info['ranges'], filter_id)
                    # record filtered positions per sample even though determined by reference genome
                    these_filtered = dict([(sample,filtered) for sample in sample_order])
                    header['FILTER'] += filter_info['string']
                
                all_filtered[VCF_path][filter_id] = these_filtered
                
                newname = _os.path.extsep.join((use_VCF_path.split(_os.path.extsep)[:-1]))
                newname = _os.path.extsep.join([newname + '__F_' + filter_id, use_VCF_path.split(_os.path.extsep)[-1]])
                print('Writing all variants with filter information to:\n{}'.format(newname))
                new_VCF_content = []
                for header_section in header_section_order:
                    new_VCF_content += ['\n'.join(header[header_section])]
                
                new_VCF_content += ['\t'.join(colnames)]
                new_VCF_content += ['\n'.join(variants)]
                new_VCF_content = '\n'.join(new_VCF_content)
                
                open(newname, 'w').write(new_VCF_content)
                use_VCF_path = newname

        self.all_filtered = all_filtered


    def reportFiltered(self, to_csv = True):
        '''
        Generate a "comma separated values" (.csv) text file for loading into a 
        spreadsheet or importing into a document that summarises variants in one or 
        more VCF files. If to_csv = False, parsed variants are only stored as an 
        attribute for further analysis.
        '''
        #self.known_filters['genome_repeats']['per_sample']
        for VCF, filters in sorted(self.all_filtered.items()):
            print('Organising filtered variants from:\n{}'.format(VCF))
            per_sample_per_position_info = {}
            for this_filter, info in filters.items():
                for sample, info2 in info.items():
                    if sample not in per_sample_per_position_info:
                        per_sample_per_position_info[sample] = {}
                    # if info2 dict contains dicts, a subfilter was applied
                    # else a single main filter which needs a single iteration
                    try:
                        a_value = info2.values()[0]
                    except IndexError:
                        # could be empty though: nothing to report
                        continue
                    
                    if isinstance(a_value, dict):
                        for sub_filter, positions in info2.items():
                            for position, (ref_char_state, sample_char_state) in positions.items():
                                if position in per_sample_per_position_info[sample]:
                                    per_sample_per_position_info[sample][position] += [(ref_char_state, sample_char_state, sub_filter)]
                                else:
                                    per_sample_per_position_info[sample][position] = [(ref_char_state, sample_char_state, sub_filter)]
                    else:
                        for position, (ref_char_state, sample_char_state) in info2.items():
                            if position in per_sample_per_position_info[sample]:
                                per_sample_per_position_info[sample][position] += [(ref_char_state, sample_char_state, this_filter)]
                            else:
                                per_sample_per_position_info[sample][position] = [(ref_char_state, sample_char_state, this_filter)]
            
            # save for general use
            self.per_sample_per_position_info = per_sample_per_position_info
            
            if to_csv:
                outfilename = VCF[:-3] + 'filtered.csv'
                print('Writing list of filtered variants to:\n{}'.format(outfilename))
                with open(outfilename, 'w') as fout:
                    colnames = ['sample', 'position', 'CDS', 'reference', 'variant', 'filter']
                    fout.write(','.join(['"'+c+'"' for c in colnames])+'\n')
                    for sample, positions in sorted(per_sample_per_position_info.items()):
                        for position, info in sorted(positions.items()):
                            ORFs = [ORF for ORF,(s, e, d, name) in self.ORF_ranges.items() if s < position <= e]
                            ORFnames = []
                            for ORF in ORFs:
                                if len(name):
                                    ORFnames += ['{} ({})'.format(ORF,name)]
                                else:
                                    ORFnames += ['{}'.format(ORF)]
                            
                            if len(ORFnames) > 0:
                                ORF = '"'+','.join(ORFnames)+'"'
                            else:
                                ORF = '""'
                            
                            for (ref_char_state, sample_char_state, this_filter) in info:
                                row_cells = ['"'+sample+'"', str(position), ORF, '"'+ref_char_state+'"', '"'+sample_char_state+'"', '"'+this_filter+'"']
                                fout.write(','.join(row_cells)+'\n')



    def doFiltering(self, filters):
        '''
        filter dict keys must be in known_filters dict below
        values are a list or dict of ranges
        '''
        # do some checks on provided genome and VCFs
        genome_ids = {}
        genome_lengths = {}
        pattern = _re.compile('##contig=<ID=([A-Za-z0-9\._]+),length=([0-9]+)>')
        for VCF in self.VCF_paths: #break
            for line in open(VCF):
                
                if line[:9] == '##contig=':
                    try:
                        genome_id, genome_length = _re.match(pattern, line).groups()
                    except AttributeError:
                        print('Failed to parse genome information from {}'.format(line))
                    genome_lengths[int(genome_length)] = VCF
                    genome_ids[genome_id] = VCF
                    #print(genome_id, genome_length)
                    identified = True
                    break
                
            e = "Failed to identify which chromosome the variants in {} were called on (couldn't find '##contig=')".format(VCF)
            assert identified, e

        e = 'Differing reference genome among provided VCFs? {}'.format(genome_ids.items())
        assert len(genome_ids) == 1, e

        ## is this only part that uses genome (i.e. not necessary)
        # e = 'Genome provided ({}) does not match genome ID in provided VCFs: {}'.format(self.genome_genbank_record.id, genome_ids.keys()[0])
        # assert self.genome_genbank_record.id == genome_ids.keys()[0], e
        # print('Variants were called against {:,} bp genome: {}\n'.format(int(genome_length), genome_id))


        filters_to_apply = {}
        for this_filter, ranges in filters.items():
            try:
                filters_to_apply[this_filter] = known_filters[this_filter]
                filters_to_apply[this_filter]['ranges'] = ranges
            except KeyError:
                print('Unknown filter type: {}. Choose from {}'.format(this_filter, ', '.join(known_filters) ))

        self.markVariants(filters_to_apply)

        self.reportFiltered()


class Linkage:
    '''Methods to measure co-incidence of alleles on the same reads or fragments.

    Currently baga.CollectData.Genome only supports single chromosome genomes
    This test should be run separately for each reference genome (chromosome) against which 
    reads are mapped against.

    Alleles on the same reads (and therefore same chromosomes) called at 
    polymorphisms in a sample of pooled genomic DNA. Infrequent co-incidence of 
    variants on the same read in nearby polymorphisms implies variants occuring 
    in different genomes in the sample (separate lineages) and has been described 
    as a "multidiverse" signature in:

    Lieberman, T. D., Flett, K. B., Yelin, I., Martin, T. R., McAdam, A. J., Priebe, 
    G. P. & Kishony, R. 
    Genetic variation of a bacterial pathogen within individuals with cystic 
    fibrosis provides a record of selective pressures.
    Nature Genetics, 2013, 46, 82-87
    '''
    def __init__(self, vcf_paths = False, 
                       alignment_paths = False, 
                       genome = False, 
                       baga = False):
        """
        Instatiate a baga.CallVariants.Linkage object.
        
        Requires either:
        - a list of baga genomes (current implementation means: chromosomes)
        (genomes = )
        - a list of BAM file paths (alignment_paths = )
        - a list to VCF file paths (VCF_paths = )
        or
        - a path to a saved baga.CallVariants.Linkage object (baga =)
        """
        
        e = 'Instantiate with an alignment object and paths to VCF files or ' + \
        'a previously saved CallerGATK'
        assert ((genome and alignment_paths and vcf_paths) and not baga) or \
               (not (genome and alignment_paths and vcf_paths) and baga), e
        
        if alignment_paths:
            e = 'Could not find file: "%s".\nPlease ensure all files exist'
            for VCF in vcf_paths:
                assert _os.path.exists(VCF), e % VCF
            
            self.VCF_paths = vcf_paths
            
            for BAM in alignment_paths:
                assert _os.path.exists(BAM), e % BAM
            
            self.alignment_paths = alignment_paths
            
            self.genome = genome
        
        elif baga:
            with _tarfile.open(baga, "r:gz") as tar:
                for member in tar:
                    contents = _StringIO(tar.extractfile(member).read())
                    try:
                        # either json serialised conventional objects
                        contents = _json.loads(contents.getvalue())
                    except ValueError:
                        #print('json failed: {}'.format(member.name))
                        # or longer python array.array objects
                        contents = _array('c', contents.getvalue())
                    
                    setattr(self, member.name, contents)

    def saveLocal(self, name):
        '''
        Save processed baga object info to a local compressed pickle file.
        
        'name' should exclude extension: .baga will be added
        '''
        fileout = 'baga.CallVariants.Linkage-%s.baga' % name
        with _tarfile.open(fileout, "w:gz") as tar:
            print('Writing to {} . . . '.format(fileout))
            for att_name, att in self.__dict__.items():
                if isinstance(att, _array):
                    io = _StringIO(att.tostring())
                    io.seek(0, _os.SEEK_END)
                    length = io.tell()
                    io.seek(0)
                    thisone = _tarfile.TarInfo(name = att_name)
                    thisone.size = length
                    tar.addfile(tarinfo = thisone, fileobj = io)
                else:
                    # try saving everything else here by jsoning
                    try:
                        io = _StringIO()
                        _json.dump(att, io)
                        io.seek(0, _os.SEEK_END)
                        length = io.tell()
                        io.seek(0)
                        thisone = _tarfile.TarInfo(name = att_name)
                        thisone.size = length
                        tar.addfile(tarinfo = thisone, fileobj = io)
                    except TypeError:
                        # ignore non-jsonable things like functions
                        # include unicodes, strings, lists etc etc
                        #print('omitting {}'.format(att_name))
                        pass

    def parsePooledVCF(self, minGQ = 0):
        '''extract variant information from a VCF with ploidy > 1 e.g. pooled'''

        def do_type(v):
          try:
            return int(v)
          except ValueError:
            try:
              return float(v)
            except ValueError:
              return v

        pooled_variants = {}

        for VCF in self.VCF_paths:
            header, header_section_order, colnames, variants = parseVCF(VCF)
            headerdict = dictify_vcf_header(header)
            these_variants = {}
            for variantline in variants:
                chrm, pos, x, ref_char, alt_chars, s, filter_status, info1, info2keys, info2values = variantline.rstrip().split('\t')
                if filter_status == 'PASS':
                    # parse VCF line
                    info1_extra = set([a for a in info1.split(';') if '=' not in a])
                    info1 = dict([a.split('=') for a in info1.split(';') if '=' in a])
                    del info1['set']
                    info2 = dict(zip(info2keys.split(':'),info2values.split(':')))
                    allinfo = {}
                    for k,v in info1.items()+info2.items():
                        allinfo[k] = do_type(v)
                    # first in list is ref, others are alts
                    allinfo['reference'] = ref_char
                    allinfo['variants'] = alt_chars.split(',')
                    allinfo['extra'] = set(info1_extra)
                    allinfo['GT'] = tuple(map(int,allinfo['GT'].split('/')))
                    ### set stringency here: Phred-scaled confidence for GT 
                    ### <10 is <90%, <20 is <99%, >20 is generally desirable
                    ### however GT isn't well suited to the artificially high
                    ### ploidy samples containing e.g. >10 samples because
                    ### of the greater demand on precision versus presence/absence
                    ### i.e. a low GQ for a 20/40 SNP doesn't imply the SNP call
                    ### itself is a false positive
                    if allinfo['GQ'] >= minGQ:
                        try:
                            these_variants[chrm][int(pos)] = allinfo
                        except KeyError:
                            these_variants[chrm] = {}
                            these_variants[chrm][int(pos)] = allinfo
            
            pooled_variants[VCF] = these_variants
                
        self.pooled_variants = pooled_variants



    def collectAdjacentPolymorphisms(self, dist = 1000):
        '''collect polymorphsims within a specific distance on chromosome'''

        clusters = {}
        for VCF,chromosomes in sorted(self.pooled_variants.items()):
            clusters[VCF] = {}
            for chromosome,variants in chromosomes.items():
                positions = sorted(variants)
                done = set()
                these_clusters = []
                for n,p1 in enumerate(positions):
                    this_cluster = []
                    for p2 in positions[(n+1):]:
                        if p2 - p1 < dist:
                            this_cluster += [p2]
                        else:
                            if len(this_cluster) > 0:
                                this_cluster += [p1]
                                this_cluster = sorted(set(this_cluster) - done)
                                if len(this_cluster) > 0:
                                    if len(this_cluster) == 1:
                                        this_cluster += [these_clusters[-1][-1]]
                        
                                    these_clusters += [sorted(this_cluster)]
                                    done.update(this_cluster)
                            break
                
                clusters[VCF][chromosome] = these_clusters

        self.clusters = clusters

    def check_loci_in_read(self, check_loci_pos1, refseq, refread0_to_varread0, chrm1_to_refread0, r):
        '''Check coincidence of alleles on single reads from pooled gDNA samples'''
        alleles_per_loci = {}
        for n,(locus_pos1,alleles) in enumerate(check_loci_pos1):
            ref_char = alleles[0]
            var_chars = alleles[1:]
            
            ## get this aligned region of ref chromosome
            # (equivalent to read without any variant positions)
            # start position first
            refread_start0 = chrm1_to_refread0[locus_pos1]
            # end of each piece is either next variant or the end of the sequence
            if n == len(check_loci_pos1) - 1:
                # this is last variant: select to end of alignment (of a read with no variants)
                refread_end0 = len(refseq) - 1
            else:
                # next variant does not align in a read if it is spanned by a deletion
                # at this or another variant <== not tested yet
                # ==> if long del absent: no problem; if del present, r.reference_end
                # will be extended because it is defined by the length of alignment, not
                # read length
                refread_end0 = chrm1_to_refread0[check_loci_pos1[n+1][0]]
                len(chrm1_to_refread0)
            
            this_piece_ref = str(refseq[refread_start0:refread_end0])
            
            ## get this aligned region of read using provided alignment
            # (containing none, some or all variants in an unknown combination of alleles
            # for comparison with refread)
            # start of a slice is at pos0
            # If var read has a deletion at start of this segment, those positions will
            # be without key (refread0 index) in refread0_to_varread0 because no
            # homologous sequence in refread0.
            # To collect appropriate, variant-containing, segment from r.query_sequence
            # as a slice: refread0-1 => varread0; then varread0+1 is correct slice start.
            # If var read has an insertion at start of this segment, those positions
            # will be without value (varread0 index) in refread0_to_varread0 because no
            # homologous sequence in varread0.
            # to collect appropriate, variant-containing, segment from r.query_sequence
            # as a slice: refread0-1 => varread0; then varread0+1 is correct slice start.
            
            if refread_start0 == 0:
                # for SNPs at first position
                varread_start0 = refread0_to_varread0[refread_start0]
            else:
                try:
                    # this should fail on indels without: -1 in ref index, +1 in read index
                    varread_start0 = refread0_to_varread0[refread_start0-1]+1
                except KeyError:
                    alleles_per_loci[locus_pos1] = 'noisy segment: unexpected read alignment'
                    continue
            
            try:
                # Fail e.g. 1) deletion up to end -1 of read causes
                # refread0_to_varread0[refread_end0-1]+1 to fail, just omitting -1 +1
                # fixes this.
                # Not encountered: if next variant is a deletion, the pos0 ending slice 
                # will be at a position not aligned between ref and read, so -1 +1 also
                # needed (did offsetting indels cause this?).
                varread_end0 = refread0_to_varread0[refread_end0-1]+1
            except KeyError:
                varread_end0 = refread0_to_varread0[refread_end0]
            
            # Other failures outside of this should be due to noisy reads and can be
            # ignored (i.e., store read as noisy).
            this_piece_read = r.query_sequence[varread_start0:varread_end0]
            this_piece_read_qualities = r.query_qualities[varread_start0:varread_end0]
            # A) does ref segment == (naive) read segment?
            variant_found = False
            if this_piece_ref == this_piece_read:
                # no variants here, record as so, then continue to next segement
                if all([len(char) == 1 for char in var_chars]):
                    alleles_per_loci[locus_pos1] = 'no mutation'
                else:
                    # revert position to pre-indel as in VCFs for recording
                    alleles_per_loci[locus_pos1-1] = 'no mutation'
                continue
                
            else:
                # This bit applies the actual changes: which variant (allele) makes ref
                # segment == var segment?
                # For reads of a population, need to try each allele to see if any match
                allele_found = False
                for var_char in var_chars:
                    # Test the region affected by latest varinat only (independent from
                    # prior variants).
                    # This work for insertions i.e. replace '' at start of ref segment
                    # does the insertion at the beginning.
                    this_piece_ref_with_var = this_piece_ref.replace(ref_char, var_char, 1)
                    # Is mutated refseq like the read? (aligned bit of chromosome plus
                    # leading unaligned bit after mutation i.e., variant present)
                    # initial compare lengths.
                    if len(this_piece_ref_with_var) != len(this_piece_read):
                        # This variant not here, check next variant (allele)
                        # Assume because this allele is absent but could also be because
                        # of apparent mutation elsewhere in segment (noise because wasn't
                        # called across all the reads in the 'pileup')
                        continue
                    else:
                        # given equal length after 'mutating' ref, after omitting low
                        # confidence positions, does ref segment with variant applied
                        # equal read segment?
                        this_piece_ref_with_var_qualpass = []
                        this_piece_read_qual_pass = []
                        for i,(char,qual) in enumerate(zip(this_piece_read, this_piece_read_qualities)):
                            if qual >= 20:
                                this_piece_read_qual_pass += [char]
                                this_piece_ref_with_var_qualpass += [this_piece_ref_with_var[i]]
                            else:
                                # these will be invisible in comparison
                                this_piece_read_qual_pass += ['-']
                                this_piece_ref_with_var_qualpass += ['-']
                        
                        this_piece_read_qual_pass = ''.join(this_piece_read_qual_pass)
                        this_piece_ref_with_var_qualpass = ''.join(this_piece_ref_with_var_qualpass)
                        if this_piece_read_qual_pass == this_piece_ref_with_var_qualpass:
                            # this is the variant here, record locus and allele of this variant
                            allele_found = True
                            if all([len(char) == 1 for char in var_chars]):
                                alleles_per_loci[locus_pos1] = var_char
                            else:
                                # revert position to pre-indel as in VCFs for recording
                                alleles_per_loci[locus_pos1-1] = var_char
                            # pieces_this_read += [this_piece_read]
                            # on to next segment
                            break
                        else:
                            # this variant not here, check next variant (allele)
                            # assume because this allele is absent but could also be
                            # because of apparent mutation elsewhere in segment (noise
                            # because wasn't called across all the reads in the 'pileup')
                            continue
                
                # This bit is currently assuming noise in an earlier segment is only a
                # problem for that segment: not later segements.
                # Previously, whole read abandoned.
                if not allele_found:
                    if all([len(char) == 1 for char in var_chars]):
                        alleles_per_loci[locus_pos1] = 'noisy segment: undetermined'
                    else:
                        # revert position to pre-indel as in VCFs for recording
                        alleles_per_loci[locus_pos1-1] = 'noisy segment: undetermined'

        return(alleles_per_loci)
    def check_within_frags(self, spanning_frags, corrected_indel_alleles, these_reads, these_reads_bypos0_indel_offsets, minMappingQuality = 60):
        '''Check coincidence of alleles on paired reads (fragements) from pooled gDNA samples'''
        alleles_per_loci_per_frag = {}
        for fragID in spanning_frags:    # break
            r1 = these_reads[fragID, True]
            r2 = these_reads[fragID, False]
            # check alignment quality is adequate
            if r1.mapq < minMappingQuality or r2.mapq < minMappingQuality:
                alleles_per_loci_per_frag[fragID] = 'low quality alignment'
                continue
            
            # 1) b) trim out read-aligned region of chromosome for attempting to apply all
            # reported variants and recording which are present.
            refseq1 = str(self.genome.sequence[r1.reference_start:r1.reference_end].tostring())
            refseq2 = str(self.genome.sequence[r2.reference_start:r2.reference_end].tostring())
            # 2) a) convert base-1 chromosome index of variant to base-0 read as aligned to
            # chromosome forward strand ==> a position mapping dictionary chrm1_to_refread0.
            chrm1_to_refread0_r1 = dict(zip(range(r1.reference_start+1,r1.reference_end+1),
                                            range(r1.reference_end-r1.reference_start)))
            refread0_to_chrm1_r1 = dict([(v,k) for k,v in chrm1_to_refread0_r1.items()])
            varread0_to_chrm1_r1 = dict(
                [(read0,chrm0+1) for read0,chrm0 in r1.aligned_pairs if \
                None not in (read0,chrm0)])
            chrm1_to_varread0_r1 = dict(
                [(chrm0+1,read0) for read0,chrm0 in r1.aligned_pairs if \
                None not in (read0,chrm0)])
            varread0_to_refread0_r1 = dict(
                [(read0,chrm1_to_refread0_r1[chrm0+1]) for read0,chrm0 in r1.aligned_pairs if \
                None not in (read0,chrm0)])
            refread0_to_varread0_r1 = dict(
                [(chrm1_to_refread0_r1[chrm0+1],read0) for read0,chrm0 in r1.aligned_pairs if \
                None not in (read0,chrm0)])
            
            chrm1_to_refread0_r2 = dict(zip(range(r2.reference_start+1,r2.reference_end+1),
                                            range(r2.reference_end-r2.reference_start)))
            refread0_to_chrm1_r2 = dict([(v,k) for k,v in chrm1_to_refread0_r2.items()])
            varread0_to_chrm1_r2 = dict(
                [(read0,chrm0+1) for read0,chrm0 in r2.aligned_pairs if \
                None not in (read0,chrm0)])
            chrm1_to_varread0_r2 = dict(
                [(chrm0+1,read0) for read0,chrm0 in r2.aligned_pairs if \
                None not in (read0,chrm0)])
            varread0_to_refread0_r2 = dict(
                [(read0,chrm1_to_refread0_r2[chrm0+1]) for read0,chrm0 in r2.aligned_pairs if \
                None not in (read0,chrm0)])
            refread0_to_varread0_r2 = dict(
                [(chrm1_to_refread0_r2[chrm0+1],read0) for read0,chrm0 in r2.aligned_pairs if \
                None not in (read0,chrm0)])
            
            # 2) b) collect polymorphic loci with alleles potentially in this read
            # can be multiple alleles per locus
            # confirm at least 2 polymorphic loci are spanned by this read
            
            these_pos0 = [pos0 for pos0,rIDs in these_reads_bypos0_indel_offsets.items() if \
                                ((fragID, True) in rIDs) or ((fragID, False) in rIDs)]
            if len(these_pos0) < 2:
                print(
                '*** problem with this read: spans {} variants ***'.format(
                                                                len(these_pos0))
                )
            
            # Collect the adjusted strings where indels occurred (omit the preceeding,
            # identical character) for this read.
            check_loci_pos1 = sorted(
                [(pos1,these_allele_strings) for pos1,these_allele_strings in \
                corrected_indel_alleles.items() if \
                pos1-1 in these_pos0])
            
            # 5) apply SNPs and indels to aligned region of chromosome and compare with
            # each query read sequence:
            # Compile list of read segments the first of which is always like reference
            # and starts list.
            # Then check each subsequent for variants at start.
            # Slicing at ORF0 position includes each ORF0 position at beginning of each
            # piece.
            # pieces_orig = [str(refseq[:chrm1_to_refread0[check_loci_pos1[0][0]]])]
            # pieces_this_read = [str(refseq[:chrm1_to_refread0[check_loci_pos1[0][0]]])]
            
            # some noisy reads might still cause exceptions but need to double check which
            # exceptions are due to noise/alignment ambiguity . . 
            
            check_loci_pos1_use = [(pos1,alleles) for pos1,alleles in check_loci_pos1 if \
                                    pos1 in chrm1_to_refread0_r1]
            alleles_per_loci_r1 = self.check_loci_in_read(check_loci_pos1_use, 
                                                          refseq1, 
                                                          refread0_to_varread0_r1, 
                                                          chrm1_to_refread0_r1, 
                                                          r1)
            
            check_loci_pos1_use = [(pos1,alleles) for pos1,alleles in check_loci_pos1 if \
                                    pos1 in chrm1_to_refread0_r2]
            alleles_per_loci_r2 = self.check_loci_in_read(check_loci_pos1_use, 
                                                          refseq2, 
                                                          refread0_to_varread0_r2, 
                                                          chrm1_to_refread0_r2, 
                                                          r2)
            
            # There shouldn't now be any 'empty' results: all polymorphic loci should
            # report be either 'no mutation', 'noisy' or the allele.
            alleles_per_loci_per_frag[fragID] = dict(alleles_per_loci_r1.items() + \
                                                     alleles_per_loci_r2.items())

        return(alleles_per_loci_per_frag)
    def check_within_reads(self, spanning_reads, corrected_indel_alleles, these_reads, these_reads_bypos0_indel_offsets, minMappingQuality = 60):
        '''
        Check coincidence of alleles on single reads from pooled gDNA samples

        Reads may be from paired-end fragments
        '''
        alleles_per_loci_per_read = {}
        for fragID, is_read1 in spanning_reads:
            r = these_reads[fragID, is_read1]
            # Check alignment quality is adequate
            if r.mapq < minMappingQuality:
                alleles_per_loci_per_read[(fragID, is_read1)] = 'low quality alignment'
                continue
            
            # 1) b) trim out read-aligned region of chromosome for attempting to apply
            # all reported variants and recording which are present.
            refseq = str(self.genome.sequence[r.reference_start:r.reference_end].tostring())
            # 2) a) convert base-1 chromosome index of variant to base-0 read as
            # aligned to chromosome forward strand ==> a position mapping dictionary
            # chrm1_to_refread0.
            #
            # reference_start = 0-based leftmost coordinate == left end of slice
            # reference_end = aligned reference position of the read on the reference genome
            # reference_end = aend which points to one past the last aligned
            # residue == right end of a (pos0) slice.
            # r.reference_start,r.reference_end form a python slice as prepared by
            # pySAM, not an inclusive list of indexes.
            # range(5) == [0,1,2,3,4]
            # range(range5[0]+1,range5[-1]+2) == [1,2,3,4,5]
            # range(1,4)
            # range(3) == [0,1,2,3,4][0:3] == [0,1,2]
            #
            # r.qstart not necessary here because this is 'ref read'; r.qstart is
            # accounted for by using r.aligned_pairs below (I think... non-0 r.qstart not tested)
            chrm1_to_refread0 = dict(zip(range(r.reference_start+1,r.reference_end+1),
                                         range(r.reference_end-r.reference_start)))
            refread0_to_chrm1 = dict([(v,k) for k,v in chrm1_to_refread0.items()])
            # Using the pysam-supplied .aligned_pairs
            # (22, 4015415),
            # (23, None),         <== omit from dicts
            # (24, 4015416),
            varread0_to_chrm1 = dict(
                    [(read0,chrm0+1) for read0,chrm0 in r.aligned_pairs if \
                    None not in (read0,chrm0)])
            chrm1_to_varread0 = dict(
                    [(chrm0+1,read0) for read0,chrm0 in r.aligned_pairs if \
                    None not in (read0,chrm0)])
            varread0_to_refread0 = dict(
                    [(read0,chrm1_to_refread0[chrm0+1]) for read0,chrm0 in r.aligned_pairs if \
                    None not in (read0,chrm0)])
            refread0_to_varread0 = dict(
                    [(chrm1_to_refread0[chrm0+1],read0) for read0,chrm0 in r.aligned_pairs if \
                    None not in (read0,chrm0)])
            
            # 2) b) collect polymorphic loci with alleles potentially in this read.
            # Can be multiple alleles per locus.
            # Confirm at least 2 polymorphic loci are spanned by this read.
            these_pos0 = [pos0 for pos0,rIDs in these_reads_bypos0_indel_offsets.items() if (fragID, is_read1) in rIDs]
            ### reads always span a maximum of two variants even when total to check is >>2
            #print('len(these_pos0) == {}; len(these_reads_bypos0_indel_offsets) == {}'.format(len(these_pos0), len(these_reads_bypos0_indel_offsets)))
            ### more than two screw up check_loci_in_read() but never seems to happen

            if len(these_pos0) < 2:
                print(
                '*** problem with this read: spans {} variants ***'.format(
                                                                len(these_pos0))
                )
            
            # Collect the adjusted strings where indels occurred (omit the preceeding,
            # identical character) for this read.
            check_loci_pos1 = sorted(
                [(pos1,these_allele_strings) for pos1,these_allele_strings in \
                corrected_indel_alleles.items() if \
                pos1-1 in these_pos0])
            
            #print('len(check_loci_pos1) == {}'.format(len(check_loci_pos1)))
            
            # 5) apply SNPs and indels to aligned region of chromosome and compare with
            # each query read sequence:
            # Compile list of read segments the first of which is always like reference
            # and starts list.
            # Then check each subsequent for variants at start.
            # Slicing at ORF0 position includes each ORF0 position at beginning of each piece
            # pieces_orig = [str(refseq[:chrm1_to_refread0[check_loci_pos1[0][0]]])]
            # pieces_this_read = [str(refseq[:chrm1_to_refread0[check_loci_pos1[0][0]]])]
            
            # some noisy reads might still cause exceptions but need to double check which
            # exceptions are due to noise/alignment ambiguity . . 
            alleles_per_loci = self.check_loci_in_read(check_loci_pos1, 
                                                       refseq, 
                                                       refread0_to_varread0, 
                                                       chrm1_to_refread0, 
                                                       r)
            # there shouldn't now be any 'empty' results: all polymorphic loci should
            # report be either 'no mutation', 'noisy' or the allele.
            alleles_per_loci_per_read[fragID, is_read1] = alleles_per_loci

        return(alleles_per_loci_per_read)

    def checkAlignments(self):
        '''
        parse BAM files checking for variants on same read or read pair (sequenced fragment)


        These notes are 'draft' :-) 

        |-----template_length,query_length---|
        -R1---\
        *============*======
                     \---R2-   if r.is_read2, query_length; pos


        if r.is_read1:  r.template_length + r.query_length == r.pnext - r.reference_start

        really same as r.is_read1?
        if r.is_read2:  r.template_length + r.query_length == r.pnext - r.reference_start

        insert end position
        r.reference_start + r.template_length

        length of this read
        r.query_length

        length of other read
        r.template_length - (r.pnext - r.reference_start)


        Track each variant in a read from start to finish building a 1-to-1 mapping of
        positions. This makes checking presence or absence variants much easier.

        0) with potentially linked variants (pos,(r,q)) and many read/read pairs in hand.

        1) a) collect single reads and fragments to be analysed based on spanning of at
        least two clustered variants. Then iterate through reads:
        1) b) trim out read-aligned region of chromosome for attempting to apply all
        reported variants and recording which are present.

        2) a) convert base-1 chromosome index of variant to base-0 read as aligned to
        chromosome forward strand ==> a position mapping dictionary chrm1_to_refread0.
        2) b) collect variants with alleles potentially in this read.

        3) adjust method (position) of indels as reported: position is where indel
        happened and one string is length zero.

        4) make another mapping dict refread0_to_chrm1: chrm1_to_refread0 and 
        refread0_to_chrm1 need to be updated as refread0 gets deletions (as varread0 
        gets insertions)

        5) apply SNPs and indels to query ORF sequence:
          thisORFseq_variant, 
          refread0_2_varread0, 
          varread0_2_refread0 = applyVariantsGetMappings(thisORFseq, 
                                                         chrm1_to_refread0, 
                                                         check_variants)

        Variants reported relative to ref read can be converted to position in variant 
        read and vice versa.
        '''

        # NB:
        # pysam is pos0
        # VCFs are pos1

        # really inserts because read pairs share ID
        num_reads_by_pos1_indel_offsets = {}

        # by VCF, by chromosome, by cluster, by [reads and/or fragment]
        polymorphism_linkage_bypop = {}

        # clusters[VCF][chromosome]

        # assumes lists of VCFs and BAMs correspond by name: not checked
        VCF2BAM = dict(zip(sorted(self.VCF_paths),sorted(self.alignment_paths)))

        for VCF,chromosomes in self.clusters.items():
            polymorphism_linkage_bypop[VCF] = {}
            for chromosome,near_vars in chromosomes.items():
                # Currently baga.CollectData.Genome only supports single chromosome
                # genomes.
                # This test should be run separately for each reference genome
                # (chromosome) mapped against
                if chromosome != self.genome.id:
                    polymorphism_linkage_bypop[VCF][chromosome] = {}
                    print(
                    "Skipping chromosome {} present in BAM because "
                    "it doesn't match supplied genome: {}".format(
                        chromosome, self.genome.id))
                    continue
                
                # 0) with potentially linked variants (pos,(r,q)) and many read/read
                # pairs in hand . . .
                print('Collecting variants for {} in {}'.format(chromosome, VCF))
                ## lists to tuples for use as keys and +1 for indels
                # # for use as key to check results at the end
                near_vars_tuples = []
                # for use in analysing sets of reads                
                near_vars_indel_offsets = []
                # for collecting all the required reads
                to_get_reads_indel_offsets = []
                corrected_indel_alleles = {}
                # for getting distance beyond a potential deletion that read needs to span
                # (not fully implemented yet?)
                indel_max_allele_lengths = {'del':{},'ins':{}}
                for these_vars in near_vars:
                    these_vars_indel_offsets = []
                    for pos1 in these_vars:
                        these_allele_strings = self.pooled_variants[VCF][chromosome][pos1]['variants']
                        longest_length = max(map(len,these_allele_strings))
                        if longest_length > 1:
                            to_get_reads_indel_offsets += [pos1+1]
                            these_vars_indel_offsets += [pos1+1]
                            corrected_indel_alleles[pos1+1] = [v[1:] for v in these_allele_strings]
                            
                            if len(these_allele_strings[0]) > 1:
                                indel_max_allele_lengths['del'][pos1+1] = len(
                                                                these_allele_strings[0]
                                                                            )-1
                            
                            if len(these_allele_strings[0]) < longest_length:
                                indel_max_allele_lengths['ins'][pos1+1] = longest_length-1
                            
                        else:
                            to_get_reads_indel_offsets += [pos1]
                            these_vars_indel_offsets += [pos1]
                            corrected_indel_alleles[pos1] = these_allele_strings
                    
                    near_vars_indel_offsets += [these_vars_indel_offsets]
                    near_vars_tuples += [tuple(these_vars)]
                
                # get the reads out of the BAM file
                reads = _pysam.Samfile(VCF2BAM[VCF])
                # all variants positions
                these_reads_bypos0_indel_offsets = {}
                these_reads = {}
                # Fetch all reads spanning all variant positions for this sample
                # => some variant positions will span reads from chromosomes with
                # deletions (which should be indicated in called variants).
                # -> i)  indels will be treated as occurring at +1 their reported
                # positions which must be accounted for when collecting reads.
                # -> ii) deletions will cause alignment length increase, insertions
                # alignments to ref will decrease in length.
                # => sometimes a fragment, identifiable by r.qname, appears twice in the
                # over the region of interest as r.is_read1 and also as r.is_read1.
                # -> therefore reads must be stored using fragment name (query_name) and
                # is_read1 (yes or no, if no its read 2).
                
                for pos1 in to_get_reads_indel_offsets:
                    # a pos0 slice from a pos1 index
                    reads_iter = reads.fetch( chromosome, pos1-1, pos1)
                    these_reads_bypos0_indel_offsets[pos1-1] = {}
                    try:
                        read_end_offset = indel_max_allele_lengths['del'][pos1]
                    except KeyError:
                        read_end_offset = 0
                    
                    for r in reads_iter:
                        # filter reads that don't span potential deletions. Need all of
                        # a potential deletion present within aligned region to replace
                        # in ref read.
                        # Rd:        ===========ddddd==
                        # Rf:    -------------V--------    <== reads must not only span
                        #                                      V but also beyond end of 
                        #                                      ddddd like Rd here
                        if pos1 + read_end_offset <= r.reference_end:
                            # pos0 end slice == pos1 end inclusive and must be present in
                            # read alignment.
                            these_reads_bypos0_indel_offsets[pos1-1][(r.query_name,r.is_read1)] = r
                            these_reads[(r.query_name,r.is_read1)] = r
                    
                    num_reads_by_pos1_indel_offsets[pos1] = len(
                                                    these_reads_bypos0_indel_offsets[pos1-1])
                
                linkages = {}
                for these_positions,these_positions_indel_offsets in zip(near_vars_tuples,
                                                                         near_vars_indel_offsets):
                    # 1) a) collect single reads and fragments to be analysed based on
                    # spanning of at least two clustered variants.
                    
                    # if len(set(these_positions_indel_offsets) & \
                           # set(indel_max_allele_lengths['ins'])):
                        # print(
                        # 'insertion at %s' % ','.join(map(str,
                                                         # set(these_positions_indel_offsets) & \
                                                         # set(indel_max_allele_lengths['ins'])
                                                         # )))
                    # elif len(set(these_positions_indel_offsets) & \
                             # set(indel_max_allele_lengths['del'])):
                        # print(
                        # 'deletion at %s' % ','.join(map(str,
                                                        # set(these_positions_indel_offsets) & \
                                                        # set(indel_max_allele_lengths['del'])
                                                        # )))
                    
                    # This assumes all variants are recorded in near_vars - additional low
                    # confidence/noise will render reads unusable . . .
                    # Not if e.g. SNP prevents recognition of mutated ref segment but
                    # also != without the mutation.
                    # Collect read pairs spanning at least two variants (potential linkage
                    # info; either in a single read or both each end of a fragment)
                    
                    # sort locus and read info per fragment
                    read_at_locus_per_frag = _defaultdict(list)
                    for pos1 in these_positions_indel_offsets:
                        for fragID, is_read1 in these_reads_bypos0_indel_offsets[pos1-1]:
                            read_at_locus_per_frag[fragID] += [(pos1,is_read1)]
                    
                    read_at_locus_per_frag = dict(read_at_locus_per_frag)
                    
                    # collect all reads spanning two or more loci (read 1 or 2)
                    spanning_reads = []
                    spanning_frags = []
                    for fragment, reads in read_at_locus_per_frag.items():
                        read1s = set()
                        read2s = set()
                        for locus,isread1 in reads:
                            if isread1:
                                read1s.add(locus)
                            else:
                                read2s.add(locus)
                        
                        if len(read1s) > 1:
                            spanning_reads += [(fragment,True)]
                        
                        if len(read2s) > 1:
                            spanning_reads += [(fragment,False)]
                        
                        # Some of these fragments will include reads with two loci
                        # themselves; 
                        # Will need merging after analyses avoiding double counts of
                        # linkage etc. Need to ensure each read spans a _different_
                        # polymorphism if collecting a fragment.
                        if len(read1s - read2s) >= 1 and len(read2s - read1s) >= 1:
                            spanning_frags += [fragment]
                    
                    spanning_reads = set(spanning_reads)
                    spanning_frags = set(spanning_frags)
                    
                    print('{}:\n\tpositions: {} ({} bp):\n'
                          '\treads spanning >=2 pos {};\n'
                          '\tread pairs spanning >=2 pos {};\n'
                          '\twith >=2 in a read {}'.format(
                                VCF,
                                ','.join(map(str,these_positions)),
                                these_positions[-1]-these_positions[0],
                                len(spanning_reads),
                                len(spanning_frags),
                                len(set([a for a,b in spanning_reads]) & spanning_frags)))
                    
                    ## check for linkage of alleles spanned by single reads
                    alleles_per_loci_per_read = self.check_within_reads(
                                                        spanning_reads, 
                                                        corrected_indel_alleles, 
                                                        these_reads, 
                                                        these_reads_bypos0_indel_offsets)
                    
                    ## check for linkage of alleles spanned by 2 reads of same fragments
                    alleles_per_loci_per_frag = self.check_within_frags(
                                                        spanning_frags,
                                                        corrected_indel_alleles,
                                                        these_reads,
                                                        these_reads_bypos0_indel_offsets)
                    
                    ## check whether fragment-level counts and read-level counts include same variant
                    ## merge results . . . not implemented yet <==============
                    overlap = set(
                        [a[0] for a in alleles_per_loci_per_read if isinstance(a,tuple)]) & \
                        set(alleles_per_loci_per_frag)
                    overlap2 = set(
                        [a for a in alleles_per_loci_per_read if \
                            isinstance(a,tuple) and \
                            a[0] in alleles_per_loci_per_frag])
                    if len(overlap) > 0:
                        print(
                        '*** overlap between within read and fragment linkage reports ***'
                        )
                        print('{}\n{}'.format(len(overlap),len(overlap2)))
                        #print('{}\n{}\n{}'.format(overlap,overlap2,alleles_per_loci_per_frag))
                    
                    ## at this stage reads are saved with keys as just the fragment name (ID)
                    ## if on a single read, or as tuple with is_read1 if on a fragment
                    ## ==> should these have same type of key here? add is_read1 to single reads?
                    linkages[these_positions] = dict(
                                                  alleles_per_loci_per_read.items() + \
                                                  alleles_per_loci_per_frag.items())
                    
                    
                polymorphism_linkage_bypop[VCF][chromosome] = linkages

        self.polymorphism_linkage_bypop = polymorphism_linkage_bypop


    def compare_allele_freqs(self, reads, VCF, chromosome):
        '''Count allele frequencies with selected reads for tabulation'''
        some_freqs = []
        num_reads = 0
        for oligo,alleles in reads.items():
            if isinstance(alleles,str):
                # read not mapped with high enough quality (MQ = 60)
                # ==> do not include in total read count
                continue
            
            if 'noisy segment: undetermined' in alleles:
                # alignment at one or other position has low position score or alignment is ambiguous
                # ==> do not include in total read count
                continue
            
            if isinstance(oligo,tuple):
                fragID = oligo[0]
            else:
                fragID = oligo
            
            num_reads += 1
            for pos1,allele in alleles.items():
                # correct for differences in indel position reporting ... 4015416 includes insertions but is not at pos1+1 ... expected all indels to be pos1+1
                if pos1 in self.pooled_variants[VCF][chromosome]:
                    use_pos1 = pos1
                else:
                    use_pos1 = pos1 - 1
                
                # correct how indels are reported for comparisons to VCF info
                these_alleles = self.pooled_variants[VCF][chromosome][use_pos1]['variants']
                longest_length = max(map(len,these_alleles))
                if longest_length > 1:
                    these_alleles = [a[1:] for a in these_alleles]
                            
                if allele in these_alleles:
                    some_freqs += [use_pos1]

        some_freqs = _Counter(some_freqs)
        for pos1,b in some_freqs.items():
            freq_at_this_locus_in_reads = int(round(40*(b/float(num_reads))))
            freq_at_this_locus_in_calls = len(
                                           filter(
                                            lambda x: x != 0, 
                                            self.pooled_variants[VCF][chromosome][pos1]['GT']
                                          ))
            print('At {}, counted in reads: {}/40, called {}/40'.format(
                                                           pos1,
                                                           freq_at_this_locus_in_reads,
                                                           freq_at_this_locus_in_calls))

        return(some_freqs)

    def tabulateResults(self, min_spanning_read_depth = 20):
        '''prepare tables containing linkage information.

        One table compares frequencies of each allele included (which checks accuracy of 
        counting method). The column headers are:
        "Chromosome position"
        "Observed Frequencies"
        "Called Frequencies"

        The other table contains biological info i.e., which bit of chromosome, its function, 
        groups of linked variants (only pairs in this case so row per pair). The column 
        headers are:
        "Variant pair positions"
        "Total reads spanning pair"
        "Reads variant at both positions"
        "Reads variant at first position only"
        "Reads variant at second position only"
        "Annotations"
        '''


        no_variants = set(('no mutation','noisy segment: undetermined'))

        tables_freqcheck = {}
        tables_linkage = {}
        linkages_undetermined = {}

        some_unaccounted = {}
        for VCF,chromosomes in self.polymorphism_linkage_bypop.items():
            for chromosome, positions in chromosomes.items():
                # Currently baga.CollectData.Genome only supports single chromosome
                # genomes.
                # This test should be run separately for each reference genome
                # (chromosome) mapped against
                if chromosome != self.genome.id:
                    # not blanking for now - not sure why/if this was necessary
                    # self.polymorphism_linkage_bypop[VCF][chromosome] = {}
                    print(
                    "Skipping chromosome {} present in BAM because "
                    "it doesn't match supplied genome: {}. This analysis "
                    "should be run separately for each reference genome "
                    "(chromosome) mapped against".format(chromosome, self.genome.id))
                    continue
                
                table_freqcheck = []
                table_linkage = []
                linkages_undetermined[VCF] = {}
                for these_positions,reads in positions.items():
                    # some summaries
                    # (tests for > 2 variants close enough to be tested by e.g. 500 bp
                    # insert fragments)
                    no_polymorphisms = set()
                    one_polymorphism = set()
                    linked_polymorphism = {}
                    low_quality_alignment = set()
                    too_noisy_to_test_linkage = set()
                    # "oligo" key is each str of fragemnt ID if only one read checked
                    # or tuple of is, is_read1 if both reads in a pair were checked
                    for oligo,alleles in reads.items():
                        if alleles == 'low quality alignment':
                            low_quality_alignment.add(oligo)
                        elif len([v for v in alleles.values() if \
                                            'noisy' not in v]) < 2:
                            too_noisy_to_test_linkage.add(oligo)
                        elif all([v == 'no mutation' for v in alleles.values() if \
                                            'noisy' not in v]):
                            no_polymorphisms.add(oligo)
                        elif len([v for v in alleles.values() if \
                                            v != 'no mutation' and \
                                            'noisy' not in v]) == 1:
                            one_polymorphism.add(oligo)
                        else:
                            linked_polymorphism[oligo] = dict(
                                        [(pos1,v) for pos1,v in alleles.items() if \
                                                v != 'no mutation' and \
                                                'noisy' not in v]
                                                         )
                    
                    print('Number of positions: {}'.format(len(these_positions)))
                    if len(linked_polymorphism) > 0 and len(these_positions) > 2:
                        print(max(map(len,linked_polymorphism.values())))
                        if max(map(len,linked_polymorphism.values())) > 2:
                            print('Three way linkage: not implemented')
                    
                    print('Pop: {}'.format(VCF))
                    print('\treads: {}'.format(len(reads)))
                    print('\tpositions: {}'.format(','.join(map(str,these_positions))))
                    print('\tlinked: {}'.format(len(linked_polymorphism)))
                    print('\tsingles: {}'.format(len(one_polymorphism)))
                    print('\tno mutations: {}'.format(len(no_polymorphisms)))
                    print('\t(not mapped: {}; too noisy: {})'.format(
                                    len(low_quality_alignment), 
                                    len(too_noisy_to_test_linkage)))
                    
                    #### linkage table ####
                    # divide up read findings by pairs of loci ==> rows
                    
                    ## variants at this stage are not always just pairs, some of these frags have >2 <== WHY?
                    ## i.e. 'by_pair' is a wrong assumption
                    ## would actually be better to scale up to arbitrary number of adjacent variants
                    
                    by_cluster = _defaultdict(list)
                    for loci in reads.values():
                        if isinstance(loci,dict):
                            variant_positions = tuple(sorted(loci))
                            muts = map(loci.get,sorted(loci))
                            by_cluster[variant_positions] += [muts]
                        else:
                            # could count unusable reads here also . . .
                            pass
                    
                    ## minimum read depth limit applied here:
                    by_cluster = dict([(cluster,combos) for cluster,combos in by_cluster.items() if \
                                            len(combos) >= min_spanning_read_depth])
                    
                    for cluster,combos in by_cluster.items():
                        # get ORF annotation(s) for this pair
                        # Not implemented: would need baga.CollectData.Genome to retain more
                        # information.
                        ## could optionally provide a gbk or DL via accession for locus info?
                        # these_features = {}
                        # for p in pair:
                            # for ORF,(s,e) in ORFslices['PAO1'].items():
                                # if s < p <= e:
                                    # these_features[p] = [ORF]
                                    # break
                                
                                # these_features[p] = [-1]
                        
                        # annos_for_loci_affected = set(collect_ORF_anno(these_features))
                        
                        # annos = []
                        # for anno in annos_for_loci_affected:
                            # if len(anno) == 3:
                                # this_annotation = '%s: %s (%s)' % tuple(anno)
                            # elif len(anno) == 2:
                                # # gene name not available
                                # this_annotation = '%s: %s' % tuple(anno)
                            # else:
                                # this_annotation = anno
                            
                            # annos += [this_annotation]
                        
                        # if len(set(annos)) == 1:
                            # annotation = annos[0]
                        # else:
                            # annotation = '; '.join(annos)
                        
                        annotation = 'not available'
                        
                        ABlinked = 0
                        Aonly = 0
                        Bonly = 0
                        wild_type = 0
                        # 1 or both variants noisy so linkage cannot be determined
                        linkage_undetermined = 0
                        
                        #### this bit needs updating to handle more than just pairs . . .
                        if len(cluster) > 2:
                            e = 'More than two alleles per read or fragment not implemented! '\
                            'Raise an issue at github.com/daveuu/baga if you need this feature. '\
                            'Polymorphisms in your data are sufficiently close for >2 to be spanned '\
                            'by single reads/fragments.'
                            raise NotImplementedError(e)
                        
                        for a,b in combos:
                            if a not in no_variants and b not in no_variants:
                                ABlinked += 1
                            elif 'noisy segment: undetermined' in (a,b):
                                linkage_undetermined += 1
                            elif a not in no_variants:
                                Aonly += 1
                            elif b not in no_variants:
                                Bonly += 1
                            else:
                                wild_type += 0
                        
                        #print(ABlinked, Aonly, Bonly, wild_type, linkage_undetermined)
                        table_linkage += [[ # "Variant pair pos"
                                            ', '.join(map(str,sorted(cluster))), 
                                            #  "Total infm reads spanning"
                                            ABlinked + Aonly + Bonly + wild_type, 
                                            # "variant at both"
                                            ABlinked, 
                                            # "variant at first"
                                            Aonly, 
                                            # "variant at second"
                                            Bonly, 
                                            # "Annotations"
                                            annotation]]
                        linkages_undetermined[VCF][tuple(sorted(cluster))] = linkage_undetermined
                    
                    #### frequency check table ####
                    # How can mismatches be explained?
                    # Reads with low alignment scores (<60), ambiguous indel alignments?
                    # "Chromosome position", "Observed Frequencies", "Called Frequencies"
                    some_freqs = self.compare_allele_freqs(reads, VCF, chromosome)
                    for pos1,b in some_freqs.items():
                        # correct for unchecked reads from total here?
                        # i.e. -len(low_quality_alignment)-len(too_noisy_to_test_linkage)
                        # correct for uncheckable reads from total read count
                        freq_at_this_locus_in_reads = int(round(40*(b/float(len(reads) - \
                                                          len(low_quality_alignment) - \
                                                          len(too_noisy_to_test_linkage)))))
                        freq_at_this_locus_in_calls = len(filter(lambda x: x != 0, 
                                                                 pooled_variants[VCF]['NC_002516.2'][pos1]['GT']))
                        table_freqcheck += [[ # "Chromosome position"
                                              pos1, 
                                              # "Observed Frequencies"
                                              freq_at_this_locus_in_reads,
                                              # "Called Frequencies"
                                              freq_at_this_locus_in_calls]]
                    
                    if len(reads) != len(linked_polymorphism) + \
                                     len(one_polymorphism) + \
                                     len(no_polymorphisms) + \
                                     len(low_quality_alignment) + \
                                     len(too_noisy_to_test_linkage):
                        print('SOME READS UNACCOUNTED FOR')
                        some_unaccounted[VCF,these_positions] = sorted(set(reads) - (linked_polymorphism | \
                                                                                     one_polymorphism | \
                                                                                     linked_polymorphism | \
                                                                                     low_quality_alignment | \
                                                                                     too_noisy_to_test_linkage))
                tables_freqcheck[VCF] = table_freqcheck
                tables_linkage[VCF] = table_linkage

        self.tables_freqcheck = tables_freqcheck
        self.tables_linkage = tables_linkage

    def writeTables(self, freq_table_name = 'freq_table.csv', 
                          linkage_table_name = 'linkage_table.csv'):
        '''make csv tables containing linkage information'''
        # frequencies of each allele
        headers_freqcheck = ["Sample",
                             "Chromosome position", 
                             "Observed Frequencies", 
                             "Called Frequencies"]

        # biological info i.e., which bit of chromosome, its function,
        # groups of linked variants (only pairs in this case so row per pair)
        headers_linkage = [ "Sample",
                            "Variant pair positions", 
                            "Total informative reads spanning pair",
                            "Reads variant at both positions", 
                            "Reads variant at first position only", 
                            "Reads variant at second position only", 
                            "Annotations"]

        with open(freq_table_name, 'w') as fout:
            print('Writing table of observed frequencies (reads counted) with frequencies '
                  'provided in the VCFs to {}'.format(freq_table_name))
            fout.write(','.join(headers_freqcheck)+'\n')
            for VCF, clusters in self.tables_freqcheck.items():
                sample = VCF.split(_os.path.sep)[-1]
                for pos, observed_freqs, called_freqs in clusters:
                    #print(sample, pos, observed_freqs, called_freqs)
                    row = '"{}","{}",{},{},{},{},"{}"\n'.format(sample, pos, 
                                                                        observed_freqs, 
                                                                        called_freqs)
                    fout.write(row)

        with open(linkage_table_name, 'w') as fout:
            print('Writing table of linkage status of nearby polymorphisms to '
                  '{}'.format(linkage_table_name))
            fout.write(','.join(headers_linkage)+'\n')
            for VCF, clusters in self.tables_linkage.items():
                sample = VCF.split(_os.path.sep)[-1]
                for pos, total_reads, reads_w_both, reads_w_A, reads_w_B, annotation in clusters:
                    #print(sample, pos, total_reads, reads_w_both, reads_w_A, reads_w_B, annotation)
                    if sum([total_reads, reads_w_both, reads_w_A, reads_w_B]) == 0:
                        print('Reads spanning polymorphisms at {} too noisy for inference ({})'.format(pos, sample))
                        continue
                    row = '"{}","{}",{},{},{},{},"{}"\n'.format(sample, pos, total_reads, 
                                                                             reads_w_both, 
                                                                             reads_w_A, 
                                                                             reads_w_B, 
                                                                             annotation)
                    fout.write(row)

    def doLinkageCheck(self, dist = 1000):
        '''Call various methods to perform linkage testing'''


        self.parsePooledVCF()
        print('Collecting nearby variants')
        self.collectAdjacentPolymorphisms(dist = dist)

        for BAM in self.alignment_paths:
            indexfile = _os.path.extsep.join([BAM,'bai'])
            if not(_os.path.exists(indexfile) and _os.path.getsize(indexfile) > 0):
                print('indexing {}'.format(BAM))
                _pysam.index(BAM)

        print('Checking read-reference alignments')
        self.checkAlignments()

        self.tabulateResults()

        self.writeTables()

if __name__ == '__main__':
    main()
