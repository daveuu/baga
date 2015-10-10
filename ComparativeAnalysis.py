#! /usr/bin/env python2
# -*- coding: utf-8 -*-@language python
#
# This file is part of the Bacterial and Archaeal Genome Analyser
# Copyright (C) 2015 David Williams
# david.williams.at.liv.d-dub.org.uk
# License GPLv3+: GNU GPL version 3 or later
# This is free software: you are free to change and redistribute it
# There is NO WARRANTY, to the extent permitted by law
# 
# Work on this software was started at The University of Liverpool, UK 
# with funding from The Wellcome Trust (093306/Z/10) awarded to:
# Dr Steve Paterson (The University of Liverpool, UK)
# Dr Craig Winstanley (The University of Liverpool, UK)
# Dr Michael A Brockhurst (University of York, UK)
#
'''
ComparativeAnalysis module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains functions for phylogenetic reconstruction and inference 
of homologous recombination from multiple sequence analysis.
'''
# stdlib
from baga import _os
from baga import _sys
from baga import _cPickle
from baga import _gzip
from baga import _subprocess
from baga import _time
from baga import _re
from baga import _array
from baga import _tarfile
from baga import _StringIO
import baga

from random import sample as _sample
from array import array as _array
import re as _re
from collections import defaultdict as _defaultdict
from collections import Counter as _Counter


# external Python modules <= could make these per Class? import at instantiation time?
from Bio import SeqIO as _SeqIO
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord
from Bio.Align import MultipleSeqAlignment as _MultipleSeqAlignment
from Bio import AlignIO as _AlignIO

from baga import get_exe_path as _get_exe_path
from baga import report_time as _report_time
from baga import CallVariants

# for non-stdlib modules that are only required by certain Classes
# issue warnings here if not found

try:
    import pysam as _pysam
except ImportError:
    print('WARNING: the MultipleSequenceAlignment Class requires PySAM but it could not be imported')
    print('Try "{} Dependencies --get pysam" to install locally'.format(_sys.argv[0]))


try:
    import dendropy as _dendropy
except ImportError:
    print('WARNING: the Plotter Class requires the DendroPy version 4 module but DendroPy could not be imported')
    print('Try "{} Dependencies --get dendropy" to install locally'.format(_sys.argv[0]))

try:
    import svgwrite as _svgwrite
except ImportError:
    print('WARNING: the Plotter Class requires the svgwrite module but it could not be imported')
    print('Try "{} Dependencies --get svgwrite" to install locally'.format(_sys.argv[0]))


def main():
    pass

def makeRanges(disjoint_consecs):
  ranges = []
  s = disjoint_consecs[0]
  for n,i in enumerate(disjoint_consecs):
    if n == 0 or n == len(disjoint_consecs) + 1:
      continue
    if i - 1 != disjoint_consecs[n-1]:
      ranges += [(s,disjoint_consecs[n-1]+1)]
      s = disjoint_consecs[n]
  
  ranges += [(s,disjoint_consecs[-1]+1)]
  return(tuple(ranges))

class MultipleSequenceAlignment:
    '''
    The MultipleSequenceAlignment class of the ComparativeAnalyses module contains 
    methods to construct genome-scale multiple sequence alignments from a 
    collection of VCF files or a single VCF file containing variants for multiple 
    samples. Read coverage information from BAM files is used to correctly 
    introduce gaps in the alignment, usually interpretted as missing data.
    '''
    # could add de novo MSA from sequences using e.g. muscle and Entrez downloads?

    def __init__(self, paths_to_VCFs):
        '''
        A MultipleSequenceAlignment Builder object must be instantiated with:
            - list of path(s) to VCF file(s)
        
        '''
        
        for VCF in paths_to_VCFs:
            try:
                f = open(VCF)
            except IOError:
                e = 'Could not access {}.\nPlease ensure all files exist and are accessible'.format(VCF)
        
        self.paths_to_VCFs = paths_to_VCFs


    def old__init__(self, path_to_VCFs, path_to_INDEL_VCFs = False):
        '''
        A MultipleSequenceAlignment Builder object must be instantiated with:
            
            - path(s) to VCF file(s)
        '''
        
        e = 'Please supply 1 VCF path for mixed variants or 2 VCF paths for separate SNPs and InDels.\n{} paths supplied'.format(len(path_to_VCFs))
        assert len(path_to_VCFs) < 3, e
        
        e = 'Could not find %s.\nPlease ensure all files exist'
        for VCF in path_to_VCFs:
            assert _os.path.exists(VCF), e % VCF
        
        if path_to_INDEL_VCFs:
            for VCF in path_to_INDEL_VCFs:
                assert _os.path.exists(VCF), e % VCF
        
        if path_to_INDEL_VCFs:
            self.path_to_SNPs_VCFs = path_to_VCFs
            self.path_to_InDels_VCFs = path_to_INDEL_VCFs
        else:
            self.mixed_VCFs = path_to_VCFs


    def collectVariants(self, 
                        samples_to_include = [], 
                        samples_to_exclude = [],
                        filters = ['rearrangements','genome_repeats','LowQual','standard_hard_filter'],
                        force_inclusion_of_invariants = False, 
                        show_totals = True):
        
        '''
        Given list of VCFs, parse them obeying specified filters and return an optional subset of samples.
        filters is a list of filters to exclude indicated in either INFO or FILTER column.
        e.g. filters = ['genome_repeats', rearrangments]
        filters must be described in CallVarinats.known_filters
        '''

        # some filters like rearrangements are expanded into two sub filters: rearrangements1 and 2
        obeyfilters_INFO = set()
        obeyfilters_FILTER = set()
        pattern = _re.compile('ID=([0-9A-Za-z_]+),')
        for this_filter in filters:
            if this_filter in CallVariants.known_filters:
                for vcf_header_string in CallVariants.known_filters[this_filter]['string']:
                    filter_name = _re.findall(pattern, vcf_header_string)[0]
                    if vcf_header_string[:6] == '##INFO':
                        obeyfilters_INFO.add(filter_name)
                    else:
                        obeyfilters_FILTER.add(filter_name)
            else:
                print('Unknown filter type: {}. Choose from {}'.format(this_filter, ', '.join(CallVariants.known_filters) ))

        genome_ids = {}
        genome_lengths = {}

        SNPs = _defaultdict(dict)
        InDels = _defaultdict(dict)
        filtered_log = []
        for VCF_path in self.paths_to_VCFs:
            header, header_section_order, colnames, these_variants = CallVariants.parseVCF(VCF_path)
            headerdict = CallVariants.dictify_vcf_header(header)
            FILTERfilters = set()
            if 'FILTER' in headerdict:
                for FILTER in headerdict['FILTER']:
                    FILTERfilters.add(FILTER['ID'].strip('\'"'))
            
            INFOfilters = set()
            if 'INFO' in headerdict:
                for INFO in headerdict['INFO']:
                    if INFO['Description'].strip('"')[:len('FILTER:')] == 'FILTER:':
                        INFOfilters.add(INFO['ID'].strip('\'"'))
            
            ### hardwired for single chromosome (contig)
            genome_ids[headerdict['contig'][0]['ID']] = VCF_path
            genome_lengths[headerdict['contig'][0]['length']] = VCF_path
            
            sample_names = colnames[9:]
            for line in these_variants:
                bits = line.split('\t')
                chromosome, pos1, ID, ref, query, qual, FILTER, INFO, FORMAT = bits[:9]
                # (polymorphisms separated with commas)
                query_char_states = query.split(',')
                FILTER = set(FILTER.split(';'))
                if len(obeyfilters_FILTER & FILTER):
                    # at least some filters present
                    filtered_log += ['Omitted variant {} to {} at {} from {} because of "{}" filter'.format(
                                        ref,
                                        query,
                                        pos1,
                                        ','.join(sample_names),
                                        ','.join(FILTER))]
                    continue
                
                # collect per sample filters for this row
                sample_INFOs = bits[9:]
                GTindex = FORMAT.split(':').index('GT')
                INFO = dict([i.split('=') for i in INFO.split(';') if '=' in i])
                #print(INFO,sample_INFOs)
                samples_filtered = {}
                for f in obeyfilters_INFO:
                    if f in INFO:
                        indexes = map(int, INFO[f].split(','))
                        for i in indexes:
                            try:
                                samples_filtered[sample_names[i]].add(f)
                            except KeyError:
                                samples_filtered[sample_names[i]] = set([f])
                
                for s,info in zip(sample_names,sample_INFOs):
                    try:
                        filtered_log += ['Omitted variant {} to {} at {} from {} because of "{}" filter'.format(
                                        ref,
                                        query,
                                        pos1,
                                        s,
                                        ','.join(samples_filtered[s]))]
                    except KeyError:
                        # not filtered so record variant
                        GTstate = info.split(':')[GTindex]
                        assert '/' not in GTstate, 'Pooled data with allele frequencies not implemented for MSAs ({})'.format(
                                                                    VCF_path)
                        if GTstate == '.':
                            try:
                                if int(INFO['DP']) != 0:
                                    print('WARNING: absent genotype but some reads present: {}'.format(INFO['DP']))
                            except KeyError:
                                # if not DP, no depth, no alignment: indel
                                query = '-'
                        else:
                            # comma separated variants, but 0 in GT is not variant
                            # so 1 should index first query, 0 causes exclusion of all variants
                            query = query_char_states[int(GTstate)-1]
                        
                        if GTstate != '0':
                            if len(ref) == len(query) == 1 and query != '-':
                                SNPs[s][int(pos1)] = (ref, query)
                            else:
                                InDels[s][int(pos1)] = (ref, query)


        self.genome_id = genome_ids.keys()[0]
        self.genome_length = genome_lengths.keys()[0]


        # one day log this
        # print('\n'.join(filtered_log))


        e = 'Differing reference genome among provided VCFs? {}'.format(genome_ids.items())
        assert len(genome_lengths) == 1, e

        if len(genome_ids) > 1:
            print('WARNING: differing genome IDs among VCFs: {}, but lengths equal ({:,}) so assuming alternative IDs for same genome.'.format(
                    ', '.join(genome_ids), int(genome_lengths.keys()[0])))


        allsamples_found = set(SNPs) | set(InDels)
        if samples_to_include:
            # retain only requested samples
            assert len(set(samples_to_include) & allsamples_found), 'none of requested samples ({}) found in VCFs ({})'.format(
                                                                                        ','.join(samples_to_include),
                                                                                        ','.join(allsamples_found))
            SNPs = dict([(a,b) for a,b in SNPs.items() if a in samples_to_include])
            InDels = dict([(a,b) for a,b in InDels.items() if a in samples_to_include])

        # ensure samples with no variants are included in variant dicts (and MSA)
        if force_inclusion_of_invariants:
            if samples_to_include:
                use_samples = samples_to_include
            else:
                use_samples = allsamples_found
            
            for sample in use_samples:
                try:
                    SNPs[sample] = SNPs[sample]
                except KeyError:
                    print('WARNING: {} requested but no SNPs found in VCF . . including as invariant because force_inclusion_of_invariants == True'.format(sample))
                    SNPs[sample] = {}
                try:
                    InDels[sample] = InDels[sample]
                except KeyError:
                    print('WARNING: {} requested but no InDels found in VCF . . including as invariant because force_inclusion_of_invariants == True'.format(sample))
                    InDels[sample] = {}

        if samples_to_exclude:
            # exclude samples if requested
            initial_samples_SNPs = set(SNPs)
            SNPs = dict([(a,b) for a,b in SNPs.items() if a not in samples_to_exclude])
            InDels = dict([(a,b) for a,b in InDels.items() if a not in samples_to_exclude])
            excluded_from_SNPs = initial_samples_SNPs - set(SNPs)
            if len(excluded_from_SNPs):
                print('Excluded from SNPs:\n{}'.format('\n'.format(excluded_from_SNPs)))
                if len(samples_to_exclude) > len(excluded_from_SNPs):
                    print('WARNING: Requested to exclude {}, but not found in VCFs'.format(
                                ', '.join(set(samples_to_exclude) - excluded_from_SNPs)))
            else:
                print('WARNING: Requested to exclude {}, but none found in VCFs'.format(
                                ', '.join(samples_to_exclude)))



        if show_totals:
            print('-SNPs-\n{}\n'.format('\n'.join(
                ['{}: {}'.format(a,len(b)) for a,b in sorted(SNPs.items())]
            )))
            print('-InDels-\n{}\n'.format('\n'.join(
                ['{}: {}'.format(a,len(b)) for a,b in sorted(InDels.items())]
            )))


        self.SNPs = SNPs
        self.InDels = InDels

    def old_parseVCFs(self, samples_to_include = [], 
                        filters_include = ['rearrangements1','rearrangements2'], 
                        filters_exclude = [], 
                        force_inclusion_of_invariants = False, 
                        show_totals = True):
        '''
        Parse VCFs.
        filters_include is a list of filter names described in INFO column.
        filters_exclude is a list of filters in the FILTER column to ignore,
        e.g. filters_exclude = ['genome_repeats']
        '''
        
        # identify chromosome for this VCF <== this will fail if running against > 1 contigs . . .
        identified = False
        genome_ids = {}
        genome_lengths = {}
        pattern = _re.compile('##contig=<ID=([A-Za-z0-9]+\.[0-9]+),length=([0-9]+)>')
        for VCF in self.paths_to_VCFs:
            for line in open(VCF):
                if line[:9] == '##contig=':
                    genome_id, genome_length = _re.match(pattern, line).groups()
                    genome_lengths[int(genome_length)] = VCF
                    genome_ids[genome_id] = VCF
                    #print(genome_id, genome_length)
                    identified = True
                    break
        
        e = 'Differing reference genome among provided VCFs? {}'.format(genome_ids.items())
        assert len(genome_ids) == 1, e
        
        e = "Failed to identify which chromosome the variants in {} were called on (couldn't find '##contig=')".format(VCFs['SNPs'])
        assert identified, e
        
        self.genome_id = genome_ids.keys()[0]
        self.genome_length = genome_lengths.keys()[0]
        
        print('Variants were called against {:,} bp genome: {}\n'.format(self.genome_length, self.genome_id))
        
        variants = {'SNPs': {}, 'InDels': {}}
        
        all_missing_samples = set()
        all_found_samples = set()
        for SNPs_VCF, InDels_VCF in zip(VCFs['SNPs'], VCFs['InDels']):
            
            # get sample names and order
            for line in open(SNPs_VCF):
                if '#CHROM' in line:
                    #print(line)
                    parameter_names = line.rstrip().split('\t')[1:9]
                    sample_names = line.rstrip().split('\t')[9:]
                    break
            
            # some sanity checks
            e = 'Requested to include {} samples, but only {} found in {}'.format(
                                                        len(samples_to_include), 
                                                        len(sample_names), 
                                                        SNPs_VCF)
            assert len(samples_to_include) <= len(sample_names), e
            
            if len(samples_to_include) > 0:
                
                # check we found them all
                samples_here = set(samples_to_include) & set(sample_names)
                all_found_samples.update(samples_here)
                #print('all_found_samples',all_found_samples)
                #print('samples_here',samples_here)
                all_missing_samples.difference_update(samples_here)
                #print('all_missing_samples',all_missing_samples)
                missing_samples = sorted(set(samples_to_include) - set(sample_names) - all_found_samples)
                #print('missing_samples',missing_samples)
                all_missing_samples.update(missing_samples)
                #print('all_missing_samples',all_missing_samples)
                e = 'Following samples not found in VCF: {}'.format(','.join(all_missing_samples))
                assert len(all_missing_samples) == 0, e
                
                print(missing_samples, SNPs_VCF)
                
                if len(samples_here) == 0:
                    continue
                
                samples_found_to_include = sorted(set(samples_to_include) & set(sample_names))
                
                samples_use_indexes = []
                for sample in samples_found_to_include:
                    samples_use_indexes += [sample_names.index(sample)]
                
                
            else:
                
                samples_use_indexes = range(len(sample_names))
                samples_found_to_include = sample_names
            
            
            for var_type, VCF_path in zip(['SNPs', 'InDels'], [SNPs_VCF, InDels_VCF]):
                these_variants = _defaultdict(dict)
                for line in open(VCF_path):
                    
                    if line[0] == '#':
                        continue
                    
                    cols = line.rstrip().split('\t')
                    #print(cols)
                    
                    # remove filters to ignore from FILTER column
                    use_FILTER = set(cols[6].split(',')) - set(filters_exclude)
                    if len(use_FILTER) == 0:
                        use_FILTER = 'PASS'
                        omitted = ', '.join(set(cols[6].split(',')) & set(filters_exclude))
                        print('Filters: ignoring "{}" at {} in {}'.format(omitted, cols[1], VCF_path))
                    else:
                        use_FILTER = ', '.join(use_FILTER)
                    
                    
                    if use_FILTER == 'PASS':
                        
                        # check any sample specific filters
                        INFO = dict([i.split('=') for i in cols[7].split(';')])
                        omit_samples = {}
                        for this_filter in filters_include:
                            try:
                                omit_samples_indexes = map(int, INFO[this_filter].split(','))
                            except KeyError:
                                continue
                            
                            if omit_samples_indexes:
                                for o in omit_samples_indexes:
                                    if sample_names[o] in omit_samples:
                                        omit_samples[sample_names[o]] += [this_filter]
                                    else:
                                        omit_samples[sample_names[o]] = [this_filter]
                        
                        # collect variant characteristics
                        pos1 = int(cols[1])
                        reference_char_state = cols[3]
                        # (polymorphisms separated with commas)
                        query_char_states = cols[4].split(',')
                        
                        # collect variant status for each sample
                        info_names = cols[8].split(':')
                        info = cols[9:]
                        for i,sample in zip(samples_use_indexes, samples_found_to_include):
                            if sample in omit_samples:
                                print('Omitting {} for {} because of filter(s) {}'.format(pos1, sample, ', '.join(omit_samples[sample])))
                                continue
                            
                            this_info = {}
                            for name, value in zip(info_names, info[i].split(':')):
                                this_info[name] = value
                            
                            # collect
                            if this_info['GT'] == '.':
                                try:
                                    if int(this_info['DP']) != 0:
                                        print('WARNING: absent genotype but some reads present: {}'.format(this_info['DP']))
                                except KeyError:
                                    # if not DP, no depth, no alignment: indel
                                    these_variants[sample][pos1] = (reference_char_state, '-')
                                
                            elif int(this_info['GT']) != 0:
                                these_variants[sample][pos1] = (reference_char_state, query_char_states[int(this_info['GT'])-1])
                    else:
                        print('Omitting {} because of filter(s): {}, for all samples in:\n{}'.format(cols[1], use_FILTER, VCF_path))
                
                these_variants = dict(these_variants)
                
                # ensure samples with no variants are included in the alignment if requested
                if len(samples_to_include) > 0 and force_inclusion_of_invariants:
                    for sample in samples_to_include:
                        try:
                            these_variants[sample] = these_variants[sample]
                        except KeyError:
                            print('WARNING: {} requested but not found in VCF . . including as invariant because force_inclusion_of_invariants == True'.format('sample'))
                            these_variants[sample] = {}
                
                variants[var_type] = dict(variants[var_type].items() + these_variants.items())
        
        if len(samples_to_include):
            e = 'Following samples not found in any VCFs: {}'.format(','.join(missing_samples))
            assert len(missing_samples) == 0, e
        
        if show_totals:
            print('-SNPs-\n{}\n'.format('\n'.join(
                ['{}: {}'.format(a,len(b)) for a,b in sorted(variants['SNPs'].items())]
            )))
            print('-InDels-\n{}\n'.format('\n'.join(
                ['{}: {}'.format(a,len(b)) for a,b in sorted(variants['InDels'].items())]
            )))
        
        self.SNPs = variants['SNPs']
        self.InDels = variants['InDels']




    def getCoverageRanges(self, path_to_BAMs = '', 
                                BAM_suffix = '', 
                                resolution = 10, 
                                min_mapping_quality = 30, 
                                load = True, 
                                save = True):
        '''
    collect per position along genome, aligned read coverage depth that pass a 
    quality standard and return ranges without reads.
    Defaults to trying load a previous set of scanned missing regions if possible.
        '''
        
        if len(path_to_BAMs):
            BAMs = {}
            reference = {}
            if isinstance(path_to_BAMs, str):
                checkBAMs = []
                for BAM in _os.listdir(path_to_BAMs):
                    if (len(BAM_suffix) == 0 and BAM[-3:] in ('bam','BAM')) or BAM[(-1*len(BAM_suffix)):] == BAM_suffix:
                        checkBAMs += [path_to_BAMs + _os.path.sep + BAM]
            else:
                checkBAMs = path_to_BAMs
            
            for BAM in checkBAMs:
                thisBAM = _pysam.Samfile(BAM)
                BAMs[thisBAM.header['RG'][0]['ID']] = thisBAM
                reference[thisBAM.header['SQ'][0]['SN']] = thisBAM.header['RG'][0]['ID']
            
            missing = set(self.SNPs) - set(BAMs)
            e = 'Samples with variants but no BAM file in {}: {}'.format(path_to_BAMs, ', '.join(missing))
            assert len(missing) == 0, e
            
            assert len(set(reference.values())), 'BAMs are relative to different references: {}'.format(', '.join(set(reference.values())))
            
            # only use BAMs for which we have SNPs
            BAMs = dict([(k,v) for k,v in BAMs.items() if k in self.SNPs])
        
        # start by collecting coverages: 
        # either done previously with Structure module
        # or by scanning here and saving
        missing_regions = {}
        start_time = _time.time()
        for snum,(sample,BAM) in enumerate(BAMs.items()):
            found_depths = False
            depths = {}
            if load:
                use_name = BAM.filename[:-4] + '_depths.baga'
                # this is a hack to save time: previous depths not found in _baserecal version
                filein = use_name.replace('_baserecal','')
                try:
                    with _tarfile.open(filein, "r:gz") as tar:
                        for member in tar:
                            contents = _StringIO(tar.extractfile(member).read())
                            depths[member.name] = _array('i', contents.getvalue())
                    
                    print("Found and loaded previously scanned coverage depths from {} which saves time!".format(filein))
                    found_depths = True
                except IOError:
                    print("Couldn't find previously scanned coverages at {}".format(filein))
            
            if not found_depths:
                print('Scanning coverages from BAM file . . .')
                # also collect proper pair depths to maintain compatibility with Structure module
                num_ref_positions = BAM.header['SQ'][0]['LN']
                depth_total = _array('i', (0,) * (int(num_ref_positions / resolution) + 1))
                depth_proper = _array('i', (0,) * (int(num_ref_positions / resolution) + 1))
                pl_iter = BAM.pileup( BAM.references[0] )
                for x in pl_iter:
                    # Coordinates in pysam are always 0-based.
                    # SAM text files use 1-based coordinates.
                    if x.pos % resolution == 0:
                        read_count_proper_pair = 0
                        read_count_non_dup = 0
                        for r in x.pileups:
                            if \
                              not r.alignment.is_qcfail and \
                              not r.alignment.is_duplicate and \
                              not r.alignment.is_secondary and \
                              r.alignment.mapping_quality >= min_mapping_quality:
                                
                                read_count_non_dup += 1
                                
                                if r.alignment.is_proper_pair:
                                    
                                    read_count_proper_pair += 1
                        
                        pos_in_arrays = x.pos / resolution
                        depth_total[pos_in_arrays] = read_count_non_dup
                        depth_proper[pos_in_arrays] = read_count_proper_pair
                
                depths['depth_total'] = depth_total
                depths['depth_proper'] = depth_proper
            
            if save and not found_depths:
                use_name = BAM.filename[:-4] + '_depths.baga'
                # this is a hack to save time: previous depths not found in _baserecal version
                fileout = use_name.replace('_baserecal','')
                try:
                    print("Saving scanned coverage depths at {} to save time if reanalysing".format(fileout))
                    def add_array(obj, name):
                        io = _StringIO(obj.tostring())
                        io.seek(0, _os.SEEK_END)
                        length = io.tell()
                        io.seek(0)
                        info = _tarfile.TarInfo(name = name)
                        info.size = length
                        tar.addfile(tarinfo = info, fileobj = io)
                    
                    with _tarfile.open(fileout, "w:gz") as tar:
                        print('Writing to {} . . . '.format(fileout))
                        add_array(depths['depth_total'], 'depth_total')
                        add_array(depths['depth_proper'], 'depth_proper')
                        
                except IOError:
                    print("Attempt to save scanned coverage depths at {} failed . . .".format(fileout))
            
            these_missing_regions = []
            collecting = False
            for pos,depth in enumerate(depths['depth_total']):
                if depth == 0:
                    if not collecting:
                        # zero depth, not collecting so start
                        collecting = True
                        these_missing_regions += [pos * resolution]
                else:
                    if collecting:
                        # reads and collecting so stop
                        collecting = False
                        these_missing_regions += [pos * resolution]
                        s, e = these_missing_regions[-2:]
                        print('{} bp at {}-{}'.format(e - s, s, e))
            
            these_missing_regions = [these_missing_regions[n:n+2] for n in range(0,len(these_missing_regions),2)]
            print('{} bp in {} gaps missing from {}'.format(sum([(e-s) for s,e in these_missing_regions]), len(these_missing_regions), sample))
            missing_regions[sample] = these_missing_regions
            # report durations, time left etc
            _report_time(start_time, snum, len(BAMs))
        
        self.missing_regions = missing_regions

    def writeMSA(self,  MSA_filename = 'multiple_alignment', 
                        strict_core = False, 
                        include_invariants = False, 
                        genome = None, 
                        missing_char = '?'):
        '''
        write SNPs to an alignment . . . 
        optionally including invariant sites and, 
        optionally restricting to core sites
        indels and missing pieces of genome are treated the same.
        '''
        
        if include_invariants:
            # do some checks for information needed for entire genome
            e = 'If including all positions in alignment, original genome must be provided'
            assert genome != None, e
            
            e = 'If including all positions in alignment, BAM files from which VCFs generated must be scanned for missing regions using getCoverageRanges() method'
            assert hasattr(self, 'missing_regions'), e
            
            ref_genome_seq = genome.sequence
            
            gaps = _Counter()
            for sample,ranges in self.missing_regions.items():
                for (s,e) in ranges:
                    for pos0 in xrange(s,e):
                        gaps[pos0] += 1
            
            #disjoint_consecs = [a/10 for a in sorted(gaps)]
            # turn gap positions into ranges
            gap_ranges = makeRanges(sorted(gaps))
            #gap_ranges = [(s*10,e*10) for s,e in gap_ranges]
            
        else:
            # prepare appropriate functions and an iterator for only variable positions
            all_SNP_reference = {}
            for s,info in self.SNPs.items():
                for pos1,(r,q) in info.items():
                    if len(r) == len(q) == 1:
                        all_SNP_reference[pos1] = r
            
            all_SNP_reference_sorted = sorted(all_SNP_reference)
            
            # get missing info 
            gaps = _Counter()
            for s,info in self.SNPs.items():
                for pos1,(r,q) in info.items():
                    if q == '-':
                        gaps[pos1] += 1
        
        if strict_core:
            print('Excluding {:,} bp from strict core'.format(len(gaps)))
        
        # attempt to add the variant for each variant column; if no variant, add the wildtype nucleotide
        alignment_arrays = {}
        # make a dict of variable column index mapping to reference genome position
        # needed to find e.g. variants shared by recombination
        variable_positions_pos1 = {}
        start_time = _time.time()
        excluded_sites = set()
        for snum,(sample,these_variants) in enumerate(self.SNPs.items()):
            dropped = []
            print('Building aligned sequence for {} ({} of {})'.format(sample, snum, len(self.SNPs)))
            this_sequence = _array('c')
            if include_invariants:
                variant_or_missing = sorted(set(gaps.keys() + these_variants.keys()))
                #variants_not_gaps = set(variant_or_missing) - set(gaps)
                start = 0
                for pos1 in variant_or_missing:
                    #print(s,pos1)
                    # add reference sequence up to this position
                    # seq = 'ccc-ddd-sss'
                    # seq[:4-1]
                    # seq[4:8-1]
                    # seq[8:]
                    if pos1 != start + 1:
                        # if not in consecutive missing pieces,
                        # add reference characters since end of last gap
                        this_sequence.extend(ref_genome_seq[start:pos1 - 1])
                    
                    start = pos1
                    # decide whether to add a variant character or skip position because not in core
                    if any([strt < pos1 <= end for strt,end in gap_ranges]):
                        # this position missing somewhere
                        if strict_core or gaps[pos1] == len(self.SNPs):
                            # skip a character because missing in at least one sample (non-core) and strict core requested
                            # or skip a character because missing in all the samples
                            dropped += ['Dropped {}, character {} because missing in {} samples'.format(pos1, ref_genome_seq[pos1 - 1], gaps[pos1])]
                            excluded_sites.add(pos1 - 1)
                        else:
                            # present in at least one sample
                            if any([strt < pos1 <= end for strt,end in self.missing_regions[sample]]):
                                # but missing here
                                this_sequence.append(missing_char)
                            else:
                                # not missing here
                                try:
                                    # add a variant character if possible
                                    this_sequence.append(these_variants[pos1][1])
                                    variable_positions_pos1[len(this_sequence)] = pos1
                                except KeyError:
                                    # add the wild-type character
                                    this_sequence.append(ref_genome_seq[pos1 - 1])
                    else:
                        # in core: either variant or not
                        try:
                            # add a variant character if possible
                            this_sequence.append(these_variants[pos1][1])
                            variable_positions_pos1[len(this_sequence)] = pos1
                        except KeyError:
                            this_sequence.append(ref_genome_seq[pos1 - 1])
                
                # add the remaining characters
                this_sequence.extend(ref_genome_seq[start:])
                
                _report_time(start_time, snum, len(self.SNPs))
                
            else:
                for pos1 in all_SNP_reference_sorted:
                    if strict_core:
                        # determine inclusion based on core strictness
                        if pos1 in gaps:
                            dropped += ['Dropped {}, character {} because missing in {} samples'.format(pos1, all_SNP_reference[pos1], gaps[pos1])]
                            excluded_sites.add(pos1)
                            continue
                    
                    try:
                        this_sequence.append(these_variants[pos1][1])
                        variable_positions_pos1[len(this_sequence)] = pos1
                    except KeyError:
                        this_sequence.append(all_SNP_reference[pos1])
            
            
            alignment_arrays[sample] = this_sequence
        
        # print('\n'.join(dropped))
        
        print('{} variable positions found (columns in multiple sequence alignment)'.format(len(variable_positions_pos1)))
        
        print("Saving alignment column to reference genome variant position mapping.")
        baga.bagasave(variable_positions_pos1, 'baga.ComparativeAnalysis.MSA.{}_dict2ref'.format(MSA_filename))
        
        print("{} sites excluded".format(len(excluded_sites)))
        if not strict_core:
            print("Excluded sites contained data for reference genome only")
        
        if include_invariants:
            if strict_core:
                alignment_arrays[self.genome_id] = _array('c', [char for pos0,char in enumerate(ref_genome_seq) if pos0 not in excluded_sites])
            else:
                # excluded sites if absent in all samples
                alignment_arrays[self.genome_id] = _array('c', [char for pos0,char in enumerate(ref_genome_seq) if pos0 not in excluded_sites])
        else:
            if strict_core:
                alignment_arrays[self.genome_id] = _array('c', [char for pos1,char in sorted(all_SNP_reference.items()) if pos1 not in excluded_sites])
            else:
                alignment_arrays[self.genome_id] = _array('c', [char for pos1,char in sorted(all_SNP_reference.items())])
        
        print('{} total columns in alignment'.format(len(alignment_arrays[self.genome_id])))
        
        # could add InDels but would only really contribute to missing data
        # unless end-user chooses to encode indels as characters?
        
        MSA_filename = MSA_filename.replace('.fna','').replace('.fasta','')
        forMSA = []
        for sample, sequence_array in sorted(alignment_arrays.items()):
            # without explicitly setting description = '', BioPython adds "<unknown description>" to sequence file
            # for each record. Some sequence file parsers include "<unknown description>" with id which can confuse
            # things e.g., ClonalFrameML where tree tip label will not match a sequence
            forMSA += [_SeqRecord( id = sample, name = '', description = '', seq = _Seq(''.join(sequence_array)) )]
            if len(sample) > 10:
                print('WARNING: {} may get truncated to ten characters ({}) in {}.phy because of Phylip specs'.format(
                                        sample, sample[:10], MSA_filename))
        
        MSA = _MultipleSeqAlignment(forMSA)
        
        print('Writing multiple nucleotide sequence alignment to Fasta file {}.fna'.format(MSA_filename))
        _SeqIO.write(MSA,'{}.fna'.format(MSA_filename), 'fasta')
        print('Writing multiple nucleotide sequence alignment to Phylip file {}.phy'.format(MSA_filename))
        _SeqIO.write(MSA,'{}.phy'.format(MSA_filename), 'phylip')
        print('Also writing a compressed Phylip file to {}.phy.gz.'.format(MSA_filename))
        print('This interleaved format compresses low diversity alignments very well and is a convenient means of archiving large alignments')
        _SeqIO.write(MSA, _gzip.open('{}.phy.gz'.format(MSA_filename), 'wb'), 'phylip')

class Phylogenetics:
    '''
    The Phylogenetics class of the ComparativeAnalyses module contains methods 
    to infer phylogenetic relationships from multiple sequence alignments using 
    programs like PhyML and for inferring recombination using ClonalFrameML. A
    multiple sequence alignment is required for input, such as one generated by 
    the Builder class of the MultipleAlignments module.
    '''

    def __init__(self, path_to_multiple_alignment, tree = False, path_to_tree = False):
        '''
        A ComparativeAnalyses Phylogenetics object must be instantiated with:
            
            - path to a multiple sequence alignment (MSA) in fasta or phylip format
            - optionally path to a Newick formated tree file corresponding to the MSA
        '''
        
        e = 'Could not find %s.\nPlease ensure file exist'
        assert _os.path.exists(path_to_multiple_alignment), e % path_to_multiple_alignment
        
        self.path_to_MSA = path_to_multiple_alignment
        
        if tree or path_to_tree:
            # check dendropy is around
            dependencies_module = 'The Dependencies module can install this package locally.'
            
            try:
                import dendropy as _dendropy
            except ImportError:
                # also check version
                print('Could not import the required DendroPy 4 package. {}'.format(dependencies_module))
            
            dendropy_major_version = int(_dendropy.__version__.split('.')[0])
            
            e = 'Please install version 4 of DendroPy (Current version is {}). {}'.format(
                                    dendropy_major_version,
                                    dependencies_module
                                    )
            assert dendropy_major_version == 4, e
        
        if tree:
            e = '"tree" must be a dendropy Tree, not a {}'.format(type(tree))
            assert type(tree) is _dendropy.datamodel.treemodel.Tree, e
            if path_to_tree:
                print('Ignoring path to tree because also supplied with a tree')
            
            self.tree = tree
            
        elif path_to_tree:
            
            e = 'Could not find %s.\nPlease ensure file exist'
            assert _os.path.exists(path_to_tree), e % path_to_tree
            
            
            self.tree = _dendropy.Tree.get_from_path(path_to_tree, 'newick')
        

    def estimate_phylogeny_PhyML(self, path_to_exe = False, num_bootstraps = 0, collect_previous = True):
        '''Infer a phylogeny using phyml and collect parameter estimates e.g. kappa (Tv/Ts ratio)'''
        try:
            with open(self.path_to_MSA + '_phyml_tree') as f:
                tree_exists = True
        except IOError:
            tree_exists = False

        try:
            with open(self.path_to_MSA + '_phyml_stats') as f:
                stats_exists = True
        except IOError:
            stats_exists = False

        if tree_exists and stats_exists and collect_previous:
            print('Collecting previously inferred phylogeny and stats from {}'.format(self.path_to_MSA + '_phyml_tree'))
            self.path_to_PhyML_tree = self.path_to_MSA + '_phyml_tree'
            self.collectPhyMLstats()
        elif collect_previous:
            print('Request to collect previous tree but . . .')
            if not tree_exists:
                print('PhyML tree not found at {}'.format(self.path_to_MSA + '_phyml_tree'))
            if not tree_exists:
                print('PhyML tree stats not found at {}'.format(self.path_to_MSA + '_phyml_stats'))
            
        else:
            if not path_to_exe:
                path_to_exe = _get_exe_path('phyml')
            
            firstline = open(self.path_to_MSA).next()
            if firstline[0] == '>':
                print('WARNING: {} looks like a FASTA file but PhyML requires Phylip sequence files')
            
            cmd = [path_to_exe]
            cmd += ['-i', self.path_to_MSA, '-o', 'tlr', '-s', 'BEST', '-t', 'e', '-d', 'nt', '-f', 'm', '-v', '0', '-b', str(num_bootstraps)]
            print('Called: {}'.format(' '.join(cmd)))
            _subprocess.call(cmd)
            self.path_to_PhyML_tree = self.path_to_MSA + '_phyml_tree'
            self.collectPhyMLstats()
            # initiate a new filename to save to
            self.savetreepath = self.path_to_PhyML_tree
    def collectPhyMLstats(self, path_to_stats_file = False):
        if path_to_stats_file:
            stats = path_to_stats_file
        else:
            stats = self.path_to_MSA + '_phyml_stats'

        print('Parsing {}'.format(stats))
        relative_rates = []
        nucleotide_frequencies = {}
        try:
            for line in open(stats):
                if 'Log-likelihood' in line:
                    Log_likelihood = float(line.split(' ')[-1])
                elif 'Unconstrained likelihood' in line:
                    Unconstrained_likelihood = float(line.split(' ')[-1])
                elif 'Parsimony' in line:
                    Parsimony = int(line.split(' ')[-1])
                elif 'Tree size' in line:
                    Tree_size = float(line.split(' ')[-1])
                elif 'Gamma shape parameter' in line:
                    Gamma_shape_parameter = float(line.split(' ')[-1])
                elif 'Relative rate in class' in line:
                    relative_rates += [float(line.split(' ')[-3].lstrip('\t'))]
                elif 'Transition/transversion ratio' in line:
                    Ts_Tv = float(line.split(' ')[-1])
                else:
                    m = _re.search('f\(([ACGT])\)= ([0-9]+\.[0-9]+)', line)
                    if m is not None:
                        nuc,freq = m.groups()
                        nucleotide_frequencies[nuc] = float(freq)
        except IOError:
            print("It seems {} couldn't be opened . . .".format(stats))
            print("An alternative can be specified with path_to_stats_file")

        try:
            self.relative_rates = relative_rates
            self.nucleotide_frequencies = nucleotide_frequencies
            self.Log_likelihood = Log_likelihood
            self.Unconstrained_likelihood = Unconstrained_likelihood
            self.Parsimony = Parsimony
            self.Tree_size = Tree_size
            self.Gamma_shape_parameter = Gamma_shape_parameter
            self.Ts_Tv = Ts_Tv
        except UnboundLocalError:
            _sys.exit('Failed to parse phyML stats file: {}'.format(stats))


    def load_tree(self, path_to_tree = ''):
        '''
        If attribute 'phyml_tree_path' is present, load tree as dendropy object
        Or load from path if provided
        '''


        if len(path_to_tree) == 0:
            e = '"path_to_PhyML_tree" attribute not found. Run estimate_phylogeny_PhyML().'
            assert hasattr(self, 'path_to_PhyML_tree'), e
            use_path = self.path_to_PhyML_tree
        else:
            e = 'Could not find %s.\nPlease ensure file exist'
            assert _os.path.exists(path_to_tree), e % path_to_tree
            use_path = path_to_tree

        dependencies_module = 'The Dependencies module can install this package locally.'

        try:
            import dendropy as _dendropy
        except ImportError:
            # also check version
            print('Could not import the required DendroPy 4 package. {}'.format(dependencies_module))

        dendropy_major_version = int(_dendropy.__version__.split('.')[0])

        e = 'Please install version 4 of DendroPy (Current version is {}). {}'.format(
                                dendropy_major_version,
                                dependencies_module
                                )
        assert dendropy_major_version == 4, e

        self.tree = _dendropy.Tree.get_from_path(use_path, 'newick')

        self.savetreepath = use_path

    def write_tree(self, schema = "newick"):
        '''
        If attribute 'tree' is present, write the tree as a newick string.
        '''


        e = 'Need a tree attribute! Try .load_tree() or .estimate_phylogeny_PhyML().'
        assert hasattr(self, 'tree'), e

        dependencies_module = 'The Dependencies module can install this package locally.'

        try:
            import dendropy as _dendropy
        except ImportError:
            # also check version
            print('Could not import the required DendroPy 4 package. {}'.format(dependencies_module))

        dendropy_major_version = int(_dendropy.__version__.split('.')[0])

        e = 'Please install version 4 of DendroPy (Current version is {}). {}'.format(
                                dendropy_major_version,
                                dependencies_module
                                )
        assert dendropy_major_version == 4, e

        self.tree.write(path = self.savetreepath, schema = schema)
        print("Wrote tree as '{}' to {}".format(schema,self.savetreepath))

    def reroot_to_outgroup(self, label_list, suffix = '_rooted'):
        '''
        Given a list of outgroup labels for rooting (often list of one), reroot phylogeny.
        '''
        out_group = set(self.tree.taxon_namespace.get_taxa(label_list, first_match_only = True))
        in_group = set(self.tree.taxon_namespace) - out_group
        in_group_member = self.tree.find_node_with_taxon(lambda t: t == list(in_group)[0])


        # this creates a root node which remains with length zero
        # even after .reroot_at_edge() below which should "suppress_unifurcations=True"
        # DendroPy 4 bug?
        # self.tree.reroot_at_edge(   in_group_member.edge, 
                                    # update_bipartitions=False)

        # side effect of rerooting at a node with a taxon is to drop that taxon from tree

        self.tree.reroot_at_node(   in_group_member.parent_node, 
                                    update_bipartitions=False)

        root_node = self.tree.mrca(taxa = out_group)

        self.tree.reroot_at_edge(   root_node.edge, 
                                    length1 = root_node.edge_length / 2, 
                                    length2 = root_node.edge_length / 2, 
                                    update_bipartitions=True)

        # insert suffix appropriate for operation
        bits = self.savetreepath.split(_os.path.extsep)
        bits[-2] = bits[-2] + suffix
        self.savetreepath = _os.path.extsep.join(bits)

    def infer_recombination(self, path_to_tree = False, path_to_exe = False, bootstraps = False, output_suffix = ''):
        '''
        Infer a recombaintion using ClonalFrameML
        Either a tree object acquired at instantiation or, an external tree is used 
        if a path is provided.
        -ignore_incomplete_sites       true   Ignore sites with any ambiguous bases.
        -use_incompatible_sites        true   Use homoplasious and multiallelic sites to correct branch lengths.
        '''

        e = 'Please run Phylogenetics.estimate_phylogeny_PhyML() first or supply estimate of transitions/transversions as "Ts_Tv attribute"'
        assert hasattr(self, 'Ts_Tv'), e

        # get path to phylogeny
        e = 'Either a path to a Newick tree file, a tree object provided at instantiation or a tree from Phylogenetics.estimate_phylogeny_PhyML() is required.'
        assert path_to_tree or hasattr(self, 'tree'), e

        if path_to_tree and hasattr(self, 'tree'):
            print('Both tree object available and path to tree file provided: using external tree in file')
        elif hasattr(self, 'tree') and not path_to_tree:
            outputname = _os.path.extsep.join(self.path_to_MSA.split(_os.path.extsep)[:-1])
            path_to_tree = outputname + '_for_ClonalFrameML_baga.tree'
            
            if _os.path.exists(path_to_tree):
                print('WARNING: overwriting {} for ClonalFrameML analysis'.format(path_to_tree))
            
            self.tree.write_to_path(path_to_tree,'newick')


        # get path to executable
        if not path_to_exe:
            path_to_exe = _get_exe_path('clonalframeml')

        # get path to multiple sequence alignment
        outputname = _os.path.extsep.join(self.path_to_MSA.split(_os.path.extsep)[:-1])
        # check if fasta version exists
        make_fasta = False
        if _os.path.isfile(_os.path.extsep.join([outputname,'fna'])):
            fasta = _os.path.extsep.join([outputname,'fna'])
        elif _os.path.isfile(_os.path.extsep.join([outputname,'fasta'])):
            fasta = _os.path.extsep.join([outputname,'fasta'])
        elif _os.path.isfile(_os.path.extsep.join([outputname,'fa'])):
            fasta = _os.path.extsep.join([outputname,'fa'])
        else:
            print("Cannot find fasta")
            make_fasta = True
            fasta = _os.path.extsep.join([outputname,'fna'])

        # check tree tips match fasta

        if not make_fasta:
            inMSA = sorted([rec.id.replace('_',' ') for rec in _SeqIO.parse(fasta, 'fasta')])
            inTree = sorted([n.taxon.label.replace('_',' ') for n in self.tree.leaf_nodes()])
            
            if inMSA != inTree:
                print("fasta seems out of date")
                make_fasta = True

        if make_fasta:
            fasta_out = open(fasta, 'w')
            print('Writing fasta file for ClonalFrameML . . .')
            n = 0
            total_tips = len(self.tree.leaf_nodes())
            for rec in _SeqIO.parse(self.path_to_MSA, 'phylip'):
                n += 1
                rec.name = ''
                rec.description = ''
                # else BioPython pollutes the fasta with <unknown description>
                _SeqIO.write(rec, fasta_out, 'fasta')
                print('{}: {} of {}'.format(rec.id, n, total_tips))
            
            fasta_out.close()


        cmd = [path_to_exe]
        cmd += [path_to_tree, fasta, outputname + output_suffix, '-kappa', str(self.Ts_Tv)]
        if bootstraps:
            e = 'please supply an integer for number of bootstraps'
            assert isinstance(bootstraps, int), e
            
            cmd += ['-emsim', str(bootstraps)]

        cmd += ['-ignore_incomplete_sites', 'true', '-use_incompatible_sites', 'true']

        print('Called: {}'.format(' '.join(cmd)))
        with open(outputname + '.log', "w") as out:
            _subprocess.call(cmd, stdout = out)






    def summarise_recombined_variants(self, output_suffix = '', reference_id = ''):
        '''
        Using dictionary generated when MSA constructed, translate ClonalFrameML 
        reported ordinates back to reference genome. A 1:1 relationship may not
        exist if only core genome aligned.
        output_suffix should match that used when calling ClonalFrame via .infer_recombination
        '''

        outputname = _os.path.extsep.join(self.path_to_MSA.split(_os.path.extsep)[:-1]) + output_suffix
        #outputname = _os.path.extsep.join(phylo_analyser.path_to_MSA.split(_os.path.extsep)[:-1]) + output_suffix

        # ClonalFrameML output
        #Both_all_core_positions_Fgr_Fr_rooted.cfml.pdf
        #Both_all_core_positions_Fgr_Fr_rooted.em.txt
        #Both_all_core_positions_Fgr_Fr_rooted.importation_status.txt
        #Both_all_core_positions_Fgr_Fr_rooted.labelled_tree.newick
        #Both_all_core_positions_Fgr_Fr_rooted.ML_sequence.fasta
        #Both_all_core_positions_Fgr_Fr_rooted.position_cross_reference.txt

        ClonalFrameML_tree = outputname + '.labelled_tree.newick'
        ClonalFrameML_regions = outputname + '.importation_status.txt'

        try:
            with open(ClonalFrameML_tree) as f:
                tree_exists = True
        except OSError:
            tree_exists = False

        try:
            with open(ClonalFrameML_regions) as f:
                importation_status_exists = True
        except OSError:
            importation_status_exists = False

        e = "Couldn't find {} nor {}. Please run .infer_recombination() first (else ensure output_suffix matches)".format(
                                                                outputname + '.labelled_tree.newick',
                                                                outputname + '.importation_status.txt'
                                                                )
        assert tree_exists and importation_status_exists, e

        # parse ClonalFrameML's "import" ranges:
        # pos1 inclusive ranges of MSA columns
        # tree node names: internal or terminal
        bits = [l.rstrip().split('\t') for l in open(ClonalFrameML_regions)]
        homoplasies = dict([(tuple(map(int,bit[1:])),[]) for bit in bits[1:]])
        affected = dict([(bit[0].replace('_',' '),[]) for bit in bits[1:]])
        for bit in bits[1:]:
            # newick underscores are converted to spaces
            nodename = bit[0].replace('_',' ')
            homoplasy = tuple(map(int,bit[1:]))
            homoplasies[homoplasy] += [nodename]
            affected[nodename] += [homoplasy]

        total_region = sum([(e-s) for s,e in homoplasies])
        print('Found {} regions totalling {} bp affecting {} samples or sample ancestors.'.format(
                            len(homoplasies), 
                            total_region,
                            len(affected)
                            ))

        # load MSA and the dict that maps MSA columns with reference genome positions
        MSA_filename = self.path_to_MSA.replace('.phy','')
        #MSA_filename = phylo_analyser.path_to_MSA.replace('.phy','')
        variable_positions_pos1 = baga.bagaload('baga.ComparativeAnalysis.MSA.{}_dict2ref'.format(MSA_filename))

        print("Loading multiple sequence alignment: {}".format(self.path_to_MSA))
        msa = _AlignIO.read(self.path_to_MSA, 'phylip')
        #msa = _AlignIO.read(phylo_analyser.path_to_MSA, 'phylip')

        # make a dict to look up MSA sequences by name
        label2MSAindex = dict([(s.id,n) for n,s in enumerate(msa)])

        # load the tree modified and labelled by ClonalFrameML
        tree = _dendropy.Tree.get_from_path(ClonalFrameML_tree, 'newick')

        # convert internal "NODE" labels to groups of terminal taxa labels 
        # (collected from appropriate clades)
        homoplasies_to_samples = {}
        for (s,e),nodes in homoplasies.items():
            homoplasies_to_samples[s,e] = []
            for node in nodes:
                this_node = tree.find_node(lambda x: x.label == node)
                if this_node:
                    homoplasies_to_samples[s,e] += [n.taxon.label for n in this_node.leaf_nodes()]
                else:
                    this_node = tree.find_node_with_taxon_label(node)
                    homoplasies_to_samples[s,e] += [this_node.taxon.label]

        # from ranges of recombinant chromosome described in terms of MSA regions and tree nodes
        # return SNPs with position in reference genome
        SNPs_by_homoplasies = {}
        for (s,e),samples in homoplasies_to_samples.items():
            SNPs_by_homoplasies[s,e] = {}
            ref_seq = str(msa[label2MSAindex[reference_id],s-1:e].seq)
            #print(s,e,variable_positions_pos1[s],variable_positions_pos1[e])
            for sample in samples: #break
                sample_seq = str(msa[label2MSAindex[sample],s-1:e].seq)
                for n,ch in enumerate(sample_seq):
                    if ch != '-' and ref_seq[n] != ch:
                        try:
                            SNPs_by_homoplasies[s,e][variable_positions_pos1[s+n], ref_seq[n], ch] += [sample]
                        except KeyError:
                            try:
                                print(variable_positions_pos1[s+n])
                                print(ref_seq[n])
                                print(ch)
                                SNPs_by_homoplasies[s,e][variable_positions_pos1[s+n], ref_seq[n], ch] = [sample]
                            except KeyError:
                                # not a SNP?
                                pass
                
                if len(SNPs_by_homoplasies[s,e]) == 0:
                    print(sample_seq)

        self.SNPs_by_homoplasies = SNPs_by_homoplasies

        # s,e = 881132, 881138 ####
        # [(a,b) for a,b in sorted(variable_positions_pos1.items()) if s - 10000 < a < e + 10000]
        # for sample,index in label2MSAindex.items():
            # A = str(msa[index,s-1:e].seq)
            # B = str(msa[label2MSAindex[reference_id],s-1:e].seq)
            # print(sample,all([a==b for a,b in zip(A,B) if a != '-']))

    def find_homoplasies(self, column_index = False):

        # collect variants from MSA
        # need to decide how to describe variants:
            # relative to an outgroup?
            # ancestral state?
            # arbitrary tip?
            # find outgroup?
            # midpoint_root?
            # longest terminal edge?
            # actually, do need rooted tree, and a focal 'in group' or clade in which to define homoplasies.

        assert hasattr(self, 'tree'), 'a tree attribute on which to find homoplasies is required (as a DendroPy Tree)'

        # load MSA
        print("Loading multiple sequence alignment: {}".format(self.path_to_MSA))
        try:
            msa = _AlignIO.read(self.path_to_MSA, 'phylip')
        except ValueError:
            assert self.path_to_MSA[-3:].lower() != 'phy', 'problem opening your Phylip file at {}'.format(self.path_to_MSA)
            try:
                msa = _AlignIO.read(self.path_to_MSA, 'fasta')
            except Exception:
                _sys.exit('There seems to be a problem opening your alignment, assumed to be a FASTA file: {}'.format(self.path_to_MSA))

        if column_index:
            # if requested, load mapping of alignment columns to chromosome positions (and therefore annotations)
            MSA_filename = self.path_to_MSA.replace('.phy','')
            #MSA_filename = phylo_analyser.path_to_MSA.replace('.phy','')
            variable_positions_pos1 = baga.bagaload('baga.ComparativeAnalysis.MSA.{}_dict2ref'.format(MSA_filename))


        # make a dict to look up MSA sequences by name
        label2MSAindex = {}
        for n,s in enumerate(msa):
            label2MSAindex[self.tree.taxon_namespace.get_taxon(label = s.id.replace('_',' '))] = n

        e = 'multiple sequence alignment and tree do not contain the same taxa'
        try:
            assert self.tree.taxon_namespace.taxa_bitmask(taxa = label2MSAindex) == self.tree.seed_node.tree_leafset_bitmask, e
        except KeyError:
            _sys.exit(e)

        for taxon,i in label2MSAindex.items():
            print(taxon, len(msa[i].seq))


        # collect deepest monophyletic clades that a variant is present in: more than one is a homoplasy
        # currently requires an outgroup of one to provide an effective ancestral state against which to compare

        # get ancestral state from outgroup
        nodes = [node for node in self.tree.seed_node.child_nodes() if node.taxon is not None]
        ## ancestral state could be provided as an argument: SeqRecord same length as MSA etc
        ## then user can either us outgroup, reference, reconstruct etc.
        assert len(nodes) == 1, 'require a single taxon outgroup to provide ancestral state'

        ancestral_state = str(msa[label2MSAindex[nodes[0].taxon]].seq)

        # use this to define whether a variant is exists or not
        print(ancestral_state)

        derived_variants_i = {}
        #for n in t.postorder_internal_node_iter():  # internal will only do pairs . . also interested in homoplasies among single pools
        for node in self.tree.postorder_node_iter():
            taxa = set([k.taxon for k in node.leaf_nodes()])
            if taxa == set(label2MSAindex):
                # not from root
                continue
            
            #### collect the positions common to the samples in these_samples
            
            these_variants = []
            for taxon in taxa:
                chars = str(msa[label2MSAindex[taxon]].seq)
                for n,(r,q) in enumerate(zip(ancestral_state,chars)):
                    if r != q:
                        these_variants += [n]
            
            these_variants = _Counter(these_variants)
            
            #print(these_variants)
            
            # these_variants = []
            # for sample in these_samples:
                # for v in in_ORF_variants_chrm_pos1[sample]:
                    # these_variants += [v]
                # for v in non_ORF_variants_chrm_pos1[sample]:
                    # these_variants += [v]
            
            # these_variants = Counter(these_variants)
            for v,c in these_variants.items():
                if c == len(taxa):
                    if v not in derived_variants_i:
                        derived_variants_i[v] = [taxa]
                    else:
                        keep_these = []
                        for g in derived_variants_i[v]:
                            if not g.issubset(taxa):
                                keep_these += [g]
                        
                        derived_variants_i[v] = keep_these + [taxa]

        # len(derived_variants_i)                                            # 1926

        print(len([(a,len(b)) for a,b in derived_variants_i.items() if len(b) > 1]))    # 23 homoplasies (losses or recombinations) across haplotypes (isolates)!

        for a,b in sorted([(a,len(b)) for a,b in derived_variants_i.items() if len(b) > 1]):
            print(a,b)

        return(derived_variants_i)

        # sorted([(a,b) for a,b in derived_variants_i.items() if len(b) > 1])



class Plotter:
    '''
    Plotter class of the ComparativeAnalysis module contains methods to plot 
    phylogenies as SVG files from .
    '''

    def __init__(self,  plot_output_path, 
                        path_to_tree = False, 
                        tree = False, 
                        width_cm = 30, height_cm = 30, 
                        viewbox_width_px = 1800, viewbox_height_px = 1800,
                        plot_width_prop = 0.8, plot_height_prop = 0.8, 
                        white_canvas = True):
        '''
        Plot a phylogeny.
        '''
        
        dependencies_module = 'The Dependencies module can install this package locally.'
        try:
            import dendropy as _dendropy
        except ImportError:
            # also check version
            print('Could not import the required DendroPy 4 package. {}'.format(dependencies_module))
        
        dendropy_major_version = int(_dendropy.__version__.split('.')[0])
        
        e = 'Please install version 4 of DendroPy (Current version is {}). {}'.format(
                                dendropy_major_version,
                                dependencies_module
                                )
        assert dendropy_major_version == 4, e
        
        try:
            import svgwrite as _svgwrite
        except ImportError:
            print('Could not import the required svgwrite package. {}'.format(dependencies_module))
        
        e = 'please supply a DendroPy tree object or a path to a newick formatted tree to load'
        assert path_to_tree or tree, e
        
        if path_to_tree:
            e = 'Could not find %s.\nPlease ensure file exist'
            assert _os.path.exists(path_to_tree), e % path_to_tree
            
            self.tree = _dendropy.Tree.get_from_path(path_to_tree, 'newick')
            if tree:
                print('WARNING: loading tree from path, but tree also supplied. Ignoring latter.')
        elif tree:
            assert type(tree) == _dendropy.datamodel.treemodel.Tree, '"tree" must be a DendroPy tree, not {}'.format(type(tree))
            self.tree = tree
        
        self.tree.update_bipartitions()
        
        if plot_output_path[-4:] not in ('.svg', '.SVG'):
            plot_output_path = plot_output_path + '.svg'
        
        dwg = _svgwrite.Drawing(plot_output_path, width='%scm' % width_cm, height='%scm' % height_cm,
                                profile='full', debug=True)
        
        dwg.viewbox(width = viewbox_width_px, height = viewbox_height_px)
        
        if white_canvas:
            dwg.add(_svgwrite.shapes.Rect(insert=(0, 0), size=(viewbox_width_px, viewbox_height_px), fill = _svgwrite.rgb(100, 100, 100, '%')))
        
        self.viewbox_width_px = viewbox_width_px
        self.viewbox_height_px = viewbox_height_px
        self.plot_width_prop = plot_width_prop
        self.plot_height_prop = plot_height_prop
        self.dwg = dwg
        self.rgb = _svgwrite.rgb
        self.contrasting_colours = ['#cc2727',
     '#ffdab0',
     '#4f6646',
     '#296fd9',
     '#8b1fa6',
     '#7f1818',
     '#ffba30',
     '#70ff70',
     '#184280',
     '#2d1633',
     '#a64949',
     '#a6791f',
     '#9ee6a8',
     '#5f90d9',
     '#f570ff',
     '#4c2222',
     '#4c380f',
     '#6a9970',
     '#96b0d9',
     '#af7bb3',
     '#e69e9e',
     '#d9b05f',
     '#136629',
     '#6a7d99',
     '#d929cd',
     '#ff4c30',
     '#735d32',
     '#43995a',
     '#465366',
     '#801879',
     '#330f0a',
     '#99896a',
     '#27cc7f',
     '#232a33',
     '#d95fc0',
     '#ff8370',
     '#4d4535',
     '#b0ffda',
     '#2958d9',
     '#4d2244',
     '#733b32',
     '#7f6b18',
     '#18805d',
     '#1d3e99',
     '#ff30ba',
     '#73544f',
     '#ffe270',
     '#0f4d38',
     '#0a1533',
     '#731654',
     '#4d3835',
     '#ffefb0',
     '#6a9989',
     '#182680',
     '#ffb0e5',
     '#cc5327',
     '#333023',
     '#354d45',
     '#0f174d',
     '#734f67',
     '#7f3418',
     '#fff130',
     '#30ffd6',
     '#7083ff',
     '#d92987',
     '#4c1f0f',
     '#bfb524',
     '#0a332b',
     '#22274d',
     '#994371',
     '#ff9670',
     '#4c480f',
     '#5fd9c0',
     '#b0baff',
     '#b37b99',
     '#995a43',
     '#666446',
     '#233330',
     '#6a7099',
     '#7f1842',
     '#d9a796',
     '#606613',
     '#188079',
     '#35384d',
     '#330a1a',
     '#99766a',
     '#939943',
     '#29cdd9',
     '#2929d9',
     '#ff70a9',
     '#ff8330',
     '#bbbf84',
     '#0f484d',
     '#2e1d99',
     '#4d2233',
     '#bf6224',
     '#8ba61f',
     '#9ee1e6',
     '#9670ff',
     '#593e49',
     '#663413',
     '#e2ff70',
     '#6a9699',
     '#6249a6',
     '#33232a',
     '#ffa970',
     '#7e8c61',
     '#29b6d9',
     '#443273',
     '#ff3068',
     '#996643',
     '#454d35',
     '#186b80',
     '#a796d9',
     '#a62043',
     '#4c3322',
     '#a0d95f',
     '#0a2b33',
     '#8729d9',
     '#4c0f1f',
     '#d9b096',
     '#557332',
     '#466066',
     '#2c2333',
     '#ffb0c5',
     '#735d4f',
     '#daffb0',
     '#299ed9',
     '#5d1880',
     '#ff304c',
     '#ff9f30',
     '#83ff30',
     '#185d80',
     '#380f4d',
     '#66131e',
     '#bf7724',
     '#274d0f',
     '#96c2d9',
     '#b05fd9',
     '#d95f70',
     '#734716',
     '#1a330a',
     '#6a8999',
     '#e5b0ff',
     '#7f3842',
     '#33200a',
     '#2a3323',
     '#35454d',
     '#7e618c',
     '#996a70',
     '#ffbc70',
     '#4ebf24',
     '#184f80',
     '#503e59',
     '#997143',
     '#348018',
     '#0f304d',
     '#d630ff']

    def reroot_to_outgroup(self, label_list):
        '''Given a list of outgroup labels for rooting (often list of one)'''

        # first, test rooting
        already_rooted_as_requested = False
        root_node = self.tree.mrca(taxa = [n.taxon for n in self.tree.leaf_nodes()])
        basal_groups = root_node.child_nodes()
        for basal_group in basal_groups:
            if set([n.taxon.label for n in basal_group.leaf_nodes()]) == set(label_list):
                already_rooted_as_requested = True
                break

        if already_rooted_as_requested and len(basal_groups) == 2:
            if not self.tree.is_rooted:
                self.tree.is_rooted = True
                self.tree.update_bipartitions()
        else:
            # although outgroup is basal, there is a multi-furcation: reroot
            out_group = set(self.tree.taxon_namespace.get_taxa(label_list, first_match_only = True))
            in_group = set(self.tree.taxon_namespace) - out_group
            in_group_member = self.tree.find_node_with_taxon(lambda t: t == list(in_group)[0])
            
            # this creates a root node which remains with length zero
            # even after .reroot_at_edge() below which should "suppress_unifurcations=True"
            # DendroPy 4 bug?
            # self.tree.reroot_at_edge(   in_group_member.edge, 
                                        # update_bipartitions=False)
            
            # side effect of rerooting at a node with a taxon is to drop that taxon from tree
            
            # finding some MRCAs make the code below fail . . . bug?
            #if not (self.tree.mrca(taxa = in_group).parent_node == self.tree.mrca(taxa = out_group).parent_node == self.tree.mrca(taxa = set(self.tree.taxon_namespace))):
                # else already rooted as requested
                
            self.tree.reroot_at_node(   in_group_member.parent_node, 
                                        update_bipartitions=True)
            
            root_node = self.tree.mrca(taxa = out_group)
            
            self.tree.reroot_at_edge(   root_node.edge, 
                                        length1 = root_node.edge_length / 4, 
                                        length2 = root_node.edge_length / 4 * 3, 
                                        update_bipartitions=True)

    def getTipOrder(self, ladderize = True):
        '''Use DendroPy to get a compatible tip order . . defaults to ladderised'''

        if ladderize:
            self.tree.ladderize(ascending=True)

        tree_string = self.tree.as_string( schema = 'newick',
                                                    suppress_internal_taxon_labels = True,
                                                    suppress_internal_node_labels = True,
                                                    suppress_rooting = True,
                                                    suppress_edge_lengths = True,
                                                    unquoted_underscores = True,
                                                    preserve_spaces = True).rstrip()

        # this may mess up underscores . . I think Newick converts spaces to underscores?
        self.tiporder = _re.findall("[^\(^\)^,^;^']+", tree_string)


    def getVGTs(self):
        '''Given a rooted Dendropy tree extract Vertical Genetic Transfers (VGTs) and put in a level ordered list'''
        # start at root node
        thisnode = self.tree.mrca(taxa=[n.taxon for n in self.tree.leaf_nodes()])
        #thisnode = phylo_plotter.tree.mrca(taxa=[n.taxon for n in phylo_plotter.tree.leaf_nodes()])
        tip_labels = [n.taxon.label for n in self.tree.leaf_nodes()]
        # tip_labels = [n.taxon.label for n in phylo_plotter.tree.leaf_nodes()]
        # stop when no children and no more edges left to follow
        arechildren = True
        # make lists of edges to return to of age (y), level (for calculating x), pair (for matching VGTs as bifurcations in x calcs)
        VGTstarts = [];VGTends = [];VGTlabels = [];VGTdepth = [];VGTsisters = [];VGTdescs = [];VGTancsisters = []
        # keep info of current position in tree. Currentpair is unique numerical ID = could be >2 == 'sisters'
        # Age is cumulative distance along edges from root and gets 'restored' when returning to a 'todo' node
        # Depth is number of edges from root and is used in x-displacement calculations and storing VGTs
        currentage = 0;currentdepth = 0;currentsisters = 0;ancsisters = None
        # store current info for each node still to do (to return to) so that 'current info' can be set accordingly
        nodestodoAge = [];nodestodoDepth = [];nodestodoSisters = [];nodestodoLSisters = [];nodestodo = []
        while arechildren:
            children = thisnode.child_nodes()                   # get children (2) => can be more if an edge is collapsed
            VGTancsisters.append(ancsisters)
            # if no children, it's a terminal node
            if len(children) == 0:
                VGTstarts.append(currentage)                    # record the start depth (y)
                assert thisnode.edge_length is not None, 'len collected: {}; len thisnode leaves: {}, {}'.format(len(VGTstarts), len(thisnode.leaf_nodes()), thisnode.taxon.label)
                currentage += thisnode.edge_length              # add edge length (edge to this node)
                VGTends.append(currentage)                      # record end depth
                VGTlabels.append(thisnode.taxon.label)          # record label (in this case a tip label)
                VGTdepth.append(currentdepth)                   # record the depth/level
                VGTsisters.append(currentsisters)                    # record pair id
                VGTdescs.append([thisnode.taxon.label])
                if len(nodestodo) == 0:
                    arechildren = False
                else:
                    # just did this node, so remove from 'to do'
                    thisnode = nodestodo.pop(0)
                    currentage = nodestodoAge.pop(0)
                    currentdepth = nodestodoDepth.pop(0)
                    currentsisters = nodestodoSisters.pop(0)
                    ancsisters = nodestodoLSisters.pop(0)
                
            else:
                # else it's internal: record 'locating' info and increment it for next VGT (edge)
                VGTdepth.append(currentdepth)
                VGTsisters.append(currentsisters)
                currentdepth += 1
                VGTstarts.append(currentage)
                if thisnode.edge_length is not None:
                    currentage += thisnode.edge_length
                    
                VGTends.append(currentage)
                if thisnode.label == None:
                    VGTlabels.append(None)                      # if present
                else:
                    #VGTlabels.append(float(thisnode.label))     # label will be a support value
                    # for combined support values pre-formatted as a string
                    VGTlabels.append(thisnode.label)     # label will be a support value
                    
                VGTdescs.append([a.taxon.label for a in thisnode.leaf_nodes()])
                # prepare for next VGT: set 'this node' to first child in list
                thisnode = children[0]
                # if bifurcation, 1 other node, but could be more: all share age and depth AND sister ID if multifurcation
                for c in children[1:]:
                    nodestodo.append(c)
                    nodestodoAge.append(currentage)
                    nodestodoDepth.append(currentdepth)
                    nodestodoSisters.append(max(VGTsisters)+1)     # the id of the other bifurcation will be that of the next
                    nodestodoLSisters.append(currentsisters)
                    
                ancsisters = currentsisters
                currentsisters = max(VGTsisters)+1                   # ids are unique so make it bigger than the biggest already existing

        # sort into levels
        alllevels = sorted(set(VGTdepth))
        VGTbylevel = []
        toadd = []              # if a VGT ends in a tip at one level (VGTlabels entry is in tip_labels), include it in all subsequent levels (include in the toadd set which gets added at each subsequent level)
        for l in alllevels:
            VGTbylevel.append([])
            VGTbylevel[l].extend(toadd)
            for n,s in enumerate(VGTstarts):
                if VGTdepth[n] == l:
                    VGTbylevel[l].append(dict({'s':VGTstarts[n],'e':VGTends[n],'l':VGTlabels[n],'d':VGTdescs[n],'S':VGTsisters[n],'a':VGTancsisters[n],'b':VGTdepth[n]}))
                    if VGTlabels[n] in tip_labels:
                        toadd.append(dict({'s':VGTstarts[n],'e':VGTends[n],'l':VGTlabels[n],'d':VGTdescs[n],'S':VGTsisters[n],'a':VGTancsisters[n],'b':VGTdepth[n]}))

        self.VGTbylevel = VGTbylevel


    def assignX_tip_distances(self, inclusionlevel = 1):
        '''
        with a given tip order, assign x's to the tips: using effective tips at a given depth
        only edges born in level >= inclusionlevel are considered for HGTs and have inner circle space ('inclusion level') 
        '''
        # this is a much simplified version of the circular plot
        # start by allocating tips an x position around inner across plot
        # all tips initial level because will span levels where not bifurcating
        # save VGTs separately to a dict indexing by an ID as they have initial x ordinates calculated
        # levels 0 is root: levels are adressed here in reverse. Descendent tips are in the highest level.
        ntips = float(len(self.VGTbylevel[-1]))
        tipwidth = 1/(ntips)
        ## trying to index 'done' using tuple of descendents i.e., split mask
        donex = {}      # this informs when spanning a level if no bifurcation and inheriting an x displacement from previous: indexed by label!
        VGTid = 0
        VGTs = {}
        # the highest level x order is determined by tiporder which in turn determines the x values for subsequent lower levels
        for n,tip in enumerate(self.tiporder):
            # find the VGT by its tip label
            thisVGT = [a for a in self.VGTbylevel[-1] if a['l']==tip][0]
            # assign even x spacing to each tip
            thisVGT['x'] = tipwidth*n
            # save x value for propogating down self.VGTbylevel where there are no bifurcations and therefore no change in x (using descendents as key)
            donex[frozenset(thisVGT['d'])] = tipwidth*n
            # store VGT with its x value under a unique ID key (incrementing integer . . . not required in this version: only value list returned)
            VGTs[VGTid] = thisVGT
            VGTid += 1

        ## go to each next level and assign x to each VGT as mid point of the VGTs
        nodelinks = []
        br = False
        for l in reversed(range(len(self.VGTbylevel)-1)):
            for VGT in self.VGTbylevel[l]:
                # a VGT may span multiple levels if fewer than max. splits along this tip to root path
                if frozenset(VGT['d']) in donex:
                    VGT['x'] = donex[frozenset(VGT['d'])]
                else:
                    # examine previous level to find descendents and X values
                    Xs_from_descendents = [thisVGT['x'] for thisVGT in self.VGTbylevel[l+1] if len(set(thisVGT['d']) & set(VGT['d'])) > 0]
                    # save this X as midpoint of previous
                    VGT['x'] = sum([max(Xs_from_descendents),min(Xs_from_descendents)])/2
                    # record which are done for X for reuse if necessary (if no bifurcations)
                    donex[frozenset(VGT['d'])] = VGT['x']
                    nodelinks.append(dict({'age':VGT['e'],'x':[min(Xs_from_descendents),max(Xs_from_descendents)]}))
                    VGTs[VGTid] = VGT
                    VGTid += 1

        self.VGTs = VGTs.values()
        self.nodelinks = nodelinks


    def addVGT( self,
                VGT, 
                bottom_left, 
                top_right, 
                max_VGT_height, 
                max_VGT_width, 
                stroke_width, 
                label_size, 
                add_to_group = False,
                colour = (10, 10, 16, '%'), 
                plotinnerlabel = False, 
                plottiplabel = False, 
                tip_label_offset = 5,
                plot_extra_line = False, 
                plotlabel_height = -1, 
                ends = 'round', 
                use_fontfamily = 'sans-serif'):
        '''
        given the bottom-left, top-right positions of tree plot area, off-set (position on x), 
        and parent and child depth (y positions) plots an edge in a phylogram

        '''

        # labels and edges could be further grouped separately (not yet implemented)
        if not add_to_group:
            add_to_this = self.dwg
        else:
            add_to_this = add_to_group

        # if it is a VGT leading to a tip, make the end square to avoid adding length with the round tip
        # with image viewbox dimensions, tree plot box areas, calculate xnorm and ynorm
        h = bottom_left[1] - top_right[1]
        w = top_right[0] - bottom_left[0]
        y_plotoffset = top_right[1]
        x_plotoffset = bottom_left[0]
        ynorm = h/float(max_VGT_height)
        xnorm = w/float(max_VGT_width)
        x = VGT['x']*xnorm + x_plotoffset

        s = h - VGT['s']*ynorm + y_plotoffset
        e = h - VGT['e']*ynorm + y_plotoffset

        # for terminal edge tip labels
        if (plottiplabel or plot_extra_line):
            textanc = 'start'      # if in to out
            if plotlabel_height == -1:
                y = e - tip_label_offset # in cm? px?
            else:
                y = h - plotlabel_height*ynorm - tip_label_offset + y_plotoffset
                if plot_extra_line:
                    add_to_this.add(self.dwg.line((x, e - 5), (x, y), 
                                                stroke = self.rgb(80, 80, 80, '%'), 
                                                style = 'stroke-linecap:butt', 
                                                stroke_width = '%spt' % stroke_width))
            
            if plottiplabel:
                # why doesn't this plot internal labels?
                tiplabel = self.dwg.text('', x = [x], y = [y], 
                                        fill='black', 
                                        font_family = use_fontfamily, 
                                        text_anchor = textanc, 
                                        font_size = '%spt' % label_size)
                tiplabel.add(self.dwg.tspan(VGT['l'], baseline_shift='-30%'))
                tiplabel.rotate(-90, center=(x,y))
                add_to_this.add(tiplabel)

        # for internal and group (family) labels
        if plotinnerlabel and VGT['l'] is not None:
            x_l, y_l = x+0.01, (e+s)/2
            add_to_this.add(self.dwg.text( VGT['l'], x=[x_l], y=[y_l], 
                                        fill = 'red', 
                                        font_family = use_fontfamily, 
                                        text_anchor='middle', 
                                        font_size = '10pt'))

        # some things have to be as CSS style strings for InkScape
        add_to_this.add(self.dwg.line(  (x, s), (x, e), 
                                stroke = self.rgb(*colour), 
                                style = 'stroke-linecap:%s' % ends, 
                                stroke_width = '%spt' % stroke_width))

    def addNodeLinks(   self, 
                        link, 
                        bottom_left, 
                        top_right, 
                        max_VGT_height, 
                        max_VGT_width, 
                        stroke_width, 
                        add_to_group = False):
        '''
        given the central coordinate, start and end off-set (position on x), 
        depth/age, link two nodes
        '''

        # labels and edges could be further grouped separately (not yet implemented)
        if not add_to_group:
            add_to_this = self.dwg
        else:
            add_to_this = add_to_group

        # the order of points given to the path command influences how it is drawn
        h = bottom_left[1] - top_right[1]
        w = top_right[0] - bottom_left[0]
        y_plotoffset = top_right[1]
        x_plotoffset = bottom_left[0]
        ynorm = h/float(max_VGT_height)
        xnorm = w/float(max_VGT_width)

        age = h - link['age']*ynorm + y_plotoffset

        xnode1 = link['x'][0]*xnorm + x_plotoffset
        xnode2 = link['x'][1]*xnorm + x_plotoffset

        if (xnode1-xnode2) > 0:
            xnode1,xnode2 = xnode2,xnode1
          
        add_to_this.add(self.dwg.line( (xnode1, age), 
                                    (xnode2, age), 
                                    stroke = self.rgb(10, 10, 16, '%'), 
                                    style = 'stroke-linecap:round', 
                                    stroke_width = '%spt' % stroke_width))

    def addHGT( self,
                HGT_coords, 
                bottom_left, 
                top_right, 
                max_VGT_height, 
                max_VGT_width, 
                position_on_edge = (1,1),
                stroke_width = 6, 
                total_tips = None,
                add_to_group = False,
                colour = (10, 10, 16, '%')):
        '''
        given the bottom-left, top-right positions of tree plot area, off-set (position 
        on x), and parent and child depth (y positions) plots linking lines between 
        homoplasies inferred to be caused by homologous recombination, or circles 
        if given one at a time.
        NB self.rgb() handling to hex string is not handled internally.
        '''

        # labels and edges could be further grouped separately (not yet implemented)
        if not add_to_group:
            add_to_this = self.dwg
        else:
            add_to_this = add_to_group

        # if it is a VGT leading to a tip, make the end square to avoid adding length with the round tip
        # with image viewbox dimensions, tree plot box areas, calculate xnorm and ynorm
        h = bottom_left[1] - top_right[1]
        w = top_right[0] - bottom_left[0]
        y_plotoffset = top_right[1]
        x_plotoffset = bottom_left[0]
        ynorm = h/float(max_VGT_height)
        xnorm = w/float(max_VGT_width)

        # determine radius of circles - this is repeated in addHGTkey()
        # 4% of page height
        radius_page = h * 0.04

        if total_tips:
            # radius relative to tips
            radius_tips = (1.0 / total_tips) * 0.65 * xnorm
            # select smallest (else circles get too big for few tips)
            radius = min(radius_page, radius_tips)
        else:
            radius = radius_page


        # translate initial coords (x: 0-1, y: phylogeny units)
        # update y offset for multiple homoplasies on one edge
        (thisHnum, totalHs) = position_on_edge
        HGT_coords_norm = []
        for x,y in HGT_coords:
            x_norm = x * xnorm + x_plotoffset
            y_norm = h - y * ynorm + y_plotoffset
            # update position on edge if other homoplases here
            if totalHs > 1:
                # start
                y_norm -= radius * (totalHs - 1)
                # update according to position
                y_norm += radius * 2 * thisHnum
            
            HGT_coords_norm += [[x_norm, y_norm]]


        if len(HGT_coords_norm) == 1:
            # obtained from outside of dataset: draw circle.
            # Radius should ideally be a bit less than 2 * distance between edges for which total needed.
            
            add_to_this.add(self.dwg.circle((HGT_coords_norm[0][0], HGT_coords_norm[0][1]), radius, 
                                        stroke = colour, 
                                        fill = colour, 
                                        fill_opacity = 0.6,
                                        stroke_width = '%spt' % (stroke_width * 0.2)))
        else:
            # shared within dataset: draw path linking.
            # currently not used: unclear when a line crosses many terminal edges
            commands = ['M {} {}'.format(*HGT_coords_norm[0])]
            
            for x,y in HGT_coords_norm[1:-1]:
                commands += ['L {} {}'.format(x,y)]
            
            commands += ['L {} {} z'.format(*HGT_coords_norm[-1])]
            
            add_to_this.add(self.dwg.path(d=' '.join(commands),
                                        stroke = colour, 
                                        stroke_opacity = 0.2,
                                        fill = 'none', 
                                        style = 'stroke-linecap:round', 
                                        stroke_width = '%spt' % stroke_width))


    def addHGTkey( self,
                homoplasies, 
                bottom_left, 
                top_right, 
                max_VGT_height, 
                max_VGT_width, 
                label_size, 
                stroke_width = 3, 
                radius = None, 
                total_tips = 10,
                use_fontfamily = 'sans-serif',
                add_to_group = False):
        '''
        given the bottom-left, top-right positions of tree plot area, off-set (position 
        on x), and parent and child depth (y positions) plots linking lines between 
        homoplasies inferred to be caused by homologous recombination
        NB self.rgb() handling to hex string is not handled internally.
        '''

        # labels and edges could be further grouped separately (not yet implemented)
        if not add_to_group:
            add_to_this = self.dwg
        else:
            add_to_this = add_to_group

        # with image viewbox dimensions, tree plot box areas, calculate xnorm and ynorm
        h = bottom_left[1] - top_right[1]
        w = top_right[0] - bottom_left[0]
        y_plotoffset = top_right[1]
        x_plotoffset = bottom_left[0]
        ynorm = h/float(max_VGT_height)
        xnorm = w/float(max_VGT_width)

        # determine radius of circles by total tips is radius not provided, defaults to 10 tips

        # for the key, fix this value
        radius_page = h * 0.012 # 0.04

        # if total_tips:
            # # radius relative to tips
            # radius_tips = (1.0 / total_tips) * 0.65 * xnorm
            # # select smallest (else circles get too big for few tips)
            # radius = min(radius_page, radius_tips)
        # else:
            # radius = radius_page

        radius = radius_page


        textanc = 'start'
        lineheight = radius * 2.5
        linestart = len(homoplasies) * lineheight
        thisline_y = linestart
        for n,homoplasy in enumerate(sorted(homoplasies)):
            add_to_this.add(self.dwg.circle((0, thisline_y), radius, 
                                        stroke = self.contrasting_colours[n], 
                                        fill = self.contrasting_colours[n], 
                                        fill_opacity = 0.6,
                                        stroke_width = '%spt' % (stroke_width * 0.2)))
            
            tiplabel = self.dwg.text('', 
                                    x = [radius * 2], 
                                    y = [thisline_y], 
                                    fill='black', 
                                    font_family = use_fontfamily, 
                                    text_anchor = textanc, 
                                    font_size = '%spt' % label_size)
            
            tiplabel.add(self.dwg.tspan('{:,} - {:,} bp'.format(*homoplasy), baseline_shift='-30%'))
            #tiplabel.rotate(-90, center=(x,y))
            add_to_this.add(tiplabel)
            thisline_y += lineheight


    def addScaleBar( self,
                bottom_left, 
                top_right, 
                max_VGT_height, 
                max_VGT_width, 
                prop_total_height = 0.1,
                stroke_width = 4, 
                label_size = 15, 
                genome_length = False,
                add_to_group = False,
                ends = 'butt', 
                use_fontfamily = 'sans-serif'):
        '''
        given the bottom-left, top-right positions of tree plot area, off-set (position on x), 
        and parent and child depth (y positions) plots an edge in a phylogram

        '''

        # labels and edges could be further grouped separately (not yet implemented)
        if not add_to_group:
            add_to_this = self.dwg
        else:
            add_to_this = add_to_group

        # if it is a VGT leading to a tip, make the end square to avoid adding length with the round tip
        # with image viewbox dimensions, tree plot box areas, calculate xnorm and ynorm
        h = bottom_left[1] - top_right[1]
        w = top_right[0] - bottom_left[0]
        y_plotoffset = top_right[1]
        x_plotoffset = bottom_left[0]
        ynorm = h/float(max_VGT_height)
        xnorm = w/float(max_VGT_width)

        # find unit length for a scale bar e.g. 10% height of page
        num_bar_plot_units = prop_total_height * h
        # ynorm is plot units per phylo units
        num_bar_phylo_units = num_bar_plot_units / ynorm
        # default unless genome length supplied
        bar_phylo_units = 'substitutions per site'

        if genome_length:
            num_bar_phylo_units *= genome_length
            bar_phylo_units = 'substitutions'
            # round to nearest whole number of substitutions
            use_num_bar_phylo_units = int(round(num_bar_phylo_units))
            rel_diff = use_num_bar_phylo_units / num_bar_phylo_units
            use_num_bar_plot_units = num_bar_plot_units * rel_diff
        else:
            if num_bar_phylo_units > 10:
                use_num_bar_phylo_units = int(round(num_bar_phylo_units, 0))
            elif num_bar_phylo_units > 1:
                use_num_bar_phylo_units = round(num_bar_phylo_units, 2)
            else:
                use_num_bar_phylo_units = '{:.2e}'.format(num_bar_phylo_units)
            
            use_num_bar_plot_units = num_bar_plot_units


        label_offset = h * 0.04

        add_to_this.add(self.dwg.line((0, 0), (use_num_bar_plot_units, 0), 
                                    stroke = 'black', 
                                    style = 'stroke-linecap:butt', 
                                    stroke_width = '%spt' % stroke_width))

        textanc = 'middle'

        scalelabel = self.dwg.text('', x = [use_num_bar_plot_units/2], y = [label_offset], 
                                fill='black', 
                                font_family = use_fontfamily, 
                                text_anchor = textanc, 
                                font_size = '%spt' % label_size)

        scalelabel.add(self.dwg.tspan('{} {}'.format(use_num_bar_phylo_units, bar_phylo_units), baseline_shift='-30%'))
        #scalelabel.rotate(-90, center=(x,y))
        add_to_this.add(scalelabel)

    def SVGtoPNG(self, outfilename, dpi = 300, format = 'png'):
        # https://download.gnome.org/sources/librsvg/2.40/librsvg-2.40.9.tar.xz
        # librsvg-2.40.9.sha256sum
        if outfilename[-3:] not in (format, format.upper()):
            outfilename = _os.path.extsep.join([outfilename,format])
        
        _subprocess.call(['rsvg-convert', self.dwg.filename, '--dpi-x', str(dpi), '--dpi-y', str(dpi), '--format', format, '--output', outfilename])

    def doPlot(self,    outgroup_label_list = [], 
                        stroke_width = 3, 
                        label_size = 10, 
                        plotinnerlabels = True,
                        plottiplabels = True,
                        tip_label_offset = 15,
                        plot_extra_lines = False,
                        direction = 'right',
                        use_names = None,
                        plot_transfers = False,
                        scale_bar = True,
                        genome_length = False,
                        raster_out = 'png'):
        '''
        Plot a phylogeny based on the information in an instantiated 
        ComparativeAnalyses.Plotter object. SVG dimensions and plot area 
        are determined at instantiation time (not here). Details of 
        plotting style can be provided as arguments.
        '''


        if len(outgroup_label_list):
            print('rerooting to outgroup {}'.format(', '.join(outgroup_label_list)))
            # replace underscores with spaces per Newick standard
            outgroup_label_list = [a.replace('_',' ') for a in outgroup_label_list]
            missing_outgroup_labels = []
            for o in outgroup_label_list:
                if self.tree.find_node_with_taxon_label(o) is None:
                    missing_outgroup_labels += [o]
            
            e = 'Requested outgroup member(s): {} not found in tree'.format(', '.join(missing_outgroup_labels))
            assert len(missing_outgroup_labels) == 0, e
            
            self.reroot_to_outgroup(outgroup_label_list)
        else:
            print('rerooting to midpoint')
            self.tree.reroot_at_midpoint(update_bipartitions=True)


        if plot_transfers:
            ## ClonalFrameML seems to set a minimum edge length of 1e-7
            ## even when PhyML provides edges of 0 and resolution to 1e-8
            ## seems to be a loss of information
            # only applies to internals
            # self.tree.collapse_unweighted_edges(threshold=1e-07, update_bipartitions=True)
            edges = [e for e in self.tree.postorder_edge_iter() if e.length is not None]
            shorties = [e for e in edges if e.length == 1e-7]
            if min([e.length for e in edges]) == 1e-7:
                print('Warning: rounding {} edge with lengths 1e-7 down to 0. These were probably introduced by ClonalFrameML'.format(len(shorties)))
                for e in shorties:
                    e.length = 0

        self.getTipOrder()

        # get start and finish info in y-axis for each edge
        self.getVGTs()


        # get horizontal displacement (for vertical tree)
        # this is relative positions between 0 and 1 to be converted to absolute positions below for plotting
        self.assignX_tip_distances()
        # padthese = [] is probably for giving internal nodes their own space

        # get dimenions of plotting area
        plotstart_x = (self.viewbox_width_px - (self.viewbox_width_px * self.plot_width_prop)) / 2
        plotend_x = plotstart_x + (self.viewbox_width_px * self.plot_width_prop)
        plottop_y = (self.viewbox_height_px - (self.viewbox_height_px * self.plot_height_prop)) / 2
        plotbottom_y = plottop_y + (self.viewbox_height_px * self.plot_height_prop)
        bottom_left = (plotstart_x, plotbottom_y)
        top_right = (plotend_x, plottop_y)
        max_VGT_height = max([a['e'] for a in self.VGTs])
        # if using lines to aligned tip labels
        # max_VGT_height += additional_line_dist
        max_VGT_width = max([a['x'] for a in self.VGTs])

        # add everything to a group make transformations easy
        main_group = self.dwg.add(self.dwg.g(id = 'phylogeny'))

        # for aligned labels:
        #plotlabel_height = max_VGT_height

        # for labels at tips
        plotlabel_height = -1

        tip_labels = set([n.taxon.label for n in self.tree.leaf_nodes()])

        #print(type(plot_transfers))

        if plot_transfers:
            # collect VGT by node name
            node_2_VGT = {}
            for n,VGT in enumerate(self.VGTs):
                node_2_VGT[VGT['l']] = VGT
            if isinstance(plot_transfers, str):
                # path to raw ClonalFrame output provided
                # get those affected and their relationships
                bits = [l.rstrip().split('\t') for l in open(plot_transfers)]
                homoplasies = dict([(tuple(map(int,bit[1:])),[]) for bit in bits[1:]])
                affected = dict([(bit[0].replace('_',' '),[]) for bit in bits[1:]])
                for bit in bits[1:]:
                    # newick underscores are converted to spaces
                    nodename = bit[0].replace('_',' ')
                    homoplasy = tuple(map(int,bit[1:]))
                    thisVGT = node_2_VGT[nodename]
                    homoplasies[homoplasy] += [(
                                     # name of node
                                     nodename, 
                                     # raw position coordinates
                                     [thisVGT['x'],thisVGT['s'] + (thisVGT['e'] - thisVGT['s']) * 0.5], 
                                     # number of other homoplasies here and which this one is
                                     # to be updated below once all homoplasies collected
                                     (0,1))]
                    affected[nodename] += [homoplasy]
            
            elif isinstance(plot_transfers, dict):
                # processed ClonalFrameML output provided
                total_region = 0
                affected = {}
                # this is a misnomer! Recombinant haplotypes is better . . .
                homoplasies = {}
                SNPs_per_homoplasy = {}
                for (MSA_s,MSA_e),SNPs in sorted(plot_transfers.items()):
                    these_samples = set()
                    for samples in SNPs.values():
                        these_samples.update(samples)
                    
                    assert these_samples == set(samples)
                    if len(these_samples) == 1:
                        nodename = these_samples.pop()
                        print(nodename)
                    else:
                        nodename = self.tree.mrca(taxon_labels = samples).label
                        print(nodename)
                    
                    thisVGT = node_2_VGT[nodename]
                    
                    # this time record in terms of reference genome positions
                    SNP_ref_positions = [SNP[0] for SNP in SNPs]
                    homoplasy = min(SNP_ref_positions),max(SNP_ref_positions)
                    SNPs_here = SNPs
                    if homoplasy not in homoplasies:
                        homoplasies[homoplasy] = []
                        SNPs_per_homoplasy[homoplasy] = SNPs
                    
                    homoplasies[homoplasy] += [(
                                     # name of node
                                     nodename, 
                                     # raw position coordinates
                                     [thisVGT['x'],thisVGT['s'] + (thisVGT['e'] - thisVGT['s']) * 0.5], 
                                     # number of other homoplasies here and which this one is
                                     # to be updated below once all homoplasies collected
                                     (0,1))]
                    try:
                        affected[nodename] += [homoplasy]
                    except KeyError:
                        affected[nodename] = [homoplasy]
            
            #print(homoplasies)
            total_region = sum([(e-s)+1 for s,e in homoplasies])
            print('Found {} regions totalling {} bp affecting {} samples or sample ancestors.'.format(
                                len(homoplasies), 
                                total_region,
                                len(affected)
                                ))
            
            # update "number of other homoplasies here and which this one is"
            for nodename,thesehomoplasies in affected.items():
                for n,homoplasy in enumerate(thesehomoplasies):
                    for n2, ((thisnodename, [x, y], [thisHnum, totalHs])) in enumerate(homoplasies[homoplasy]):
                        if nodename == thisnodename:
                            thisHnum, totalHs = n, len(thesehomoplasies)
                            homoplasies[homoplasy][n2] = ((thisnodename, [x, y], [thisHnum, totalHs]))
            
            for n, (homoplasy, infos) in enumerate(sorted(homoplasies.items())):
                # homoplasies share colours but are plotted by unkinked circles
                # would need to apply self.rgb() to make hex string before here
                this_colour = self.contrasting_colours[n]
                for (thisnodename, [x, y], [thisHnum, totalHs]) in infos:
                    # need thisHnum, totalHs to apply correct offset for
                    # multiple homoplasies on one edge
                    HGT_coords = [[x, y]]
                    self.addHGT( HGT_coords, 
                                 bottom_left, 
                                 top_right, 
                                 max_VGT_height, 
                                 max_VGT_width, 
                                 position_on_edge = (thisHnum, totalHs),
                                 stroke_width = 6, 
                                 total_tips = len(self.tree.leaf_nodes()),
                                 add_to_group = main_group,
                                 colour = this_colour)
            
            # add homoplasy key
            HGT_key_group = self.dwg.add(self.dwg.g(id = 'HGT_key_group'))
            self.addHGTkey(homoplasies, 
                        bottom_left, 
                        top_right, 
                        max_VGT_height, 
                        max_VGT_width, 
                        label_size = 20,
                        stroke_width = 6 * 0.2,
                        #radius = 10, 
                        total_tips = len(self.tree.leaf_nodes()),
                        use_fontfamily = 'sans-serif',
                        add_to_group = HGT_key_group)
            
            # move it bottom right
            HGT_key_group.translate(plotend_x * 0.88, plotbottom_y * 0.70)


        if use_names:
            d = _re.compile('[\t ]+')
            try:
                name_translate = dict([_re.split(d, l.rstrip(), maxsplit = 1) for l in open(use_names) if len(l) > 2])
            except IOError:
                _sys.exit('Could not access the requested --use_names file: {}'.format(use_name))
            
            # remove underscores for compatibility with DendroPy's Newick handling i.e., with Newick
            name_translate = dict([(a.replace('_',' '),b.replace('_',' ')) for a,b in name_translate.items()])
            
            for t in tip_labels:
                if t not in name_translate:
                    print('WARNING: cannot find tip label {} in supplied tip name translation file: {}'.format(
                                                        t, use_names))
            
            # current style is to retain original label (reads accession) in parentheses
            new_tip_labels = set()
            for t in tip_labels:
                for VGT in self.VGTs:
                    if VGT['l'] == t:
                        try:
                            VGT['l'] = '{} ({})'.format(name_translate[t], t)
                            new_tip_labels.add('{} ({})'.format(name_translate[t], t))
                        except KeyError:
                            # no translation supplied for this tip
                            new_tip_labels.add(t)
            
            tip_labels = new_tip_labels

        if scale_bar:
            
            scalebar_group = self.dwg.add(self.dwg.g(id = 'scale_bar'))
            
            self.addScaleBar(bottom_left, 
                        top_right, 
                        max_VGT_height, 
                        max_VGT_width, 
                        prop_total_height = 0.1,
                        stroke_width = 4, 
                        label_size = 20, 
                        genome_length = genome_length,
                        add_to_group = scalebar_group,
                        ends = 'butt', 
                        use_fontfamily = 'sans-serif')
            
            # move it bottom left
            scalebar_group.translate(plotend_x * 0.15, plotbottom_y * 0.98)



        for VGT in self.VGTs:
            if VGT['l'] in tip_labels:
                # its a tip
                self.addVGT(VGT, 
                            bottom_left, 
                            top_right, 
                            max_VGT_height, 
                            max_VGT_width, 
                            stroke_width, 
                            label_size, 
                            add_to_group = main_group,
                            plotinnerlabel = False, 
                            plottiplabel = plottiplabels, 
                            tip_label_offset = tip_label_offset,
                            plot_extra_line = plot_extra_lines, 
                            plotlabel_height = plotlabel_height, 
                            ends = 'butt')
            else:
                self.addVGT(VGT,
                            bottom_left, 
                            top_right, 
                            max_VGT_height, 
                            max_VGT_width, 
                            stroke_width, 
                            label_size, 
                            add_to_group = main_group,
                            plotinnerlabel = plotinnerlabels, 
                            plottiplabel = False, 
                            plot_extra_line = False, 
                            plotlabel_height = plotlabel_height, 
                            ends = 'round')

        for link in self.nodelinks:
            self.addNodeLinks(  link, 
                                bottom_left, 
                                top_right, 
                                max_VGT_height, 
                                max_VGT_width, 
                                stroke_width,
                                add_to_group = main_group)



        # do any transformations
        if direction != 'right':
            print('Direction: "{}", not implemented! Will plot to the right.'.format(direction))

        angle = 90
        main_group.rotate(  angle, 
                            center = (
                                (plotstart_x + plotend_x) / 2, 
                                (plottop_y + plotbottom_y) / 2
                            )
                          )


        self.dwg.save()

        if raster_out:
            self.SVGtoPNG(self.dwg.filename[:-4], format = raster_out)


if __name__ == '__main__':
    main()
