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

# external Python modules
from Bio import SeqIO as _SeqIO
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord

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
    dict for tabulation of variants sorted by filters described in the VCF 
    header.
    Per sample filters stored in INFO columns must be described in header 
    ##INFO entries and have Descriptions starting "FILTER:".
    '''
    # identify chromosome for this VCF <== this will fail if running against > 1 contigs . . .
    try:
        pattern = _re.compile('##contig=<ID=([A-Za-z0-9]+\.[0-9]+),length=([0-9]+)>')
        # for multiple chromosomes iterate here
        genome_id, genome_length = _re.match(pattern, header['contig'][0]).groups()
        genome_length = int(genome_length)
    except KeyError:
        print("Failed to identify which chromosome the variants were called on (couldn't find '##contig=')")

    print('Variants were called against {:,} bp genome: {}\n'.format(genome_length, genome_id))

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
                variants[sample][cols['CHROM']][int(cols['POS'])] = [[cols['REF'], cols['ALT'].split(',')[int(data['GT'])-1]], samples_filtered[sample]]

    # convert nested defaultdicts to dicts
    variants = dict(variants)
    for k1,v1 in variants.items():
        variants[k1] = dict(v1)
        for k2,v2 in variants[k1].items():
            variants[k1][k2] = dict(v2)

    allfilters = FILTERfilters | INFOfilters

    return(variants, allfilters)

def reportCumulative(filter_order, reference_id, VCFs, VCFs_indels = False):
    '''Generate simple table of cumulative effects of filters applied to a VCF file'''
    # cumulative affect of filters

    ## build table column names
    colnames = ["Cumulative filters applied"]
    for dataset in sorted(VCFs):
        for v in ['SNPs', 'InDels']:
            colnames += ["{} {}".format(v, dataset)]

    variant_type_order = ['SNPs', 'InDels']

    for v in variant_type_order:
        colnames += ["{} all samples".format(v)]


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

    # start with no filters
    filters_applied_ordered = [()]
    # then add some
    for filtername in filter_order:
        filters_applied_ordered += [short2fullnames[filtername]]

    collect_baga_filters = [f for f in filter_order if 'GATK' not in f]
    from glob import glob as _glob
    # find actual VCFs . . find the VCFs with suffices that match the requested filters
    VCFs_use = {}
    for dataset,varianttypes in sorted(VCFs.items()):
        VCFs_use[dataset] = {}
        for varianttype,filename in varianttypes.items():
            bits = filename.split(_os.path.extsep)
            pattern = _os.path.extsep.join(bits[:-1]) + '*' + bits[-1]
            for checkthis in _glob(pattern):
                filter_present = []
                for fltr in collect_baga_filters:
                    if 'F_'+fltr in checkthis:
                        filter_present += [fltr]
                
                if set(filter_present) >= set(collect_baga_filters):
                    # OK if additional filters included in a VCF
                    VCFs_use[dataset][varianttype] = checkthis


    # build table

    cumulative_filters = set()

    rows = []
    for filters in filters_applied_ordered:
        cumulative_filters.update(filters)
        try:
            this_row = [filter_names[filters]]
        except KeyError:
            this_row = ['None']
        
        # must be saved as sets to only count variants once each
        #totals_by_type = dict(zip(variant_type_order, [set()] * len(variant_type_order)))
        totals_by_type = {}
        for varianttype in variant_type_order:
            totals_by_type[varianttype] = set()
        
        for dataset,varianttypes in sorted(VCFs_use.items()):
            print(dataset)
            for varianttype in variant_type_order:
                filename = varianttypes[varianttype]
                header, header_section_order, these_colnames, variantrows = parseVCF(filename)
                variants, allfilters = sortVariantsKeepFilter(header, these_colnames, variantrows)
                by_position, by_position_filtered = to_by_position_filtered(variants, cumulative_filters)
                print(varianttype, len(by_position[reference_id]))
                this_row += [len(by_position[reference_id])]
                totals_by_type[varianttype].update([info[0] for info in by_position[reference_id]])
                print(totals_by_type[varianttype])
        
        this_row += [len(totals_by_type[varianttype]) for varianttype in variant_type_order]
        
        rows += [this_row]

    outfilename = 'Table_of_cumulative_variant_totals_with_filters.csv'
    print('Printing to {}'.format(outfilename))
    with open(outfilename, 'w') as fout:
        fout.write(','.join(['"{}"'.format(c) for c in colnames])+'\n')
        for row in rows:
            fout.write(','.join(['"{}"'.format(row[0])]+[str(c) for c in row[1:]])+'\n')


def to_by_position_filtered(variants, filters_applied, summarise = True):
    by_position = _defaultdict(_Counter)
    by_position_filtered = _defaultdict(_Counter)
    for sample, chromosomes in variants.items():
        for chromosome, positions in chromosomes.items():
            # iterate through variants by position
            for position, ((reference,query),filters) in sorted(positions.items()):
                if len(filters & filters_applied) == 0:
                    # retain variants without any filters flagged (of those we are interested in)
                    by_position[chromosome][(position,reference,query,None)] += 1
                else:
                    for f in filters & filters_applied:
                        # also retain those with a filter flag, seperately for each filter
                        by_position_filtered[chromosome][(position,reference,query,f)] += 1

    if summarise:
        for f1 in sorted(filters_applied):
            print('-- {} --'.format(f1))
            these = []
            for position,reference,query,f2 in sorted(by_position_filtered[chromosome]):
                if f1 == f2:
                    these += ['{}: {} => {}, {}'.format(position,reference,query,f2)]
            
            print('Total: {}'.format(len(these)))
            print('\n'.join(these))

    return(by_position,by_position_filtered)

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
class Caller:
    '''
    A collection of short read datasets that have been aligned to the same genome 
    sequence in Sequence Alignment/Map (SAM) format for variant calling using GATK.
    '''
    def __init__(self, alignments = False, baga = False):
        '''
        Initialise with:
        a baga.AlignReads.SAMs object
        or
        a saved baga.CallVariants.Caller object
        '''
        assert alignments or baga, 'Instantiate with alignments or a previously saved Caller'
        assert not (alignments and baga), 'Instantiate with alignments OR a previously saved Caller!'
        
        if alignments:
            try:
                self.ready_BAMs = alignments.ready_BAMs
            except AttributeError:
                print('''
        ERROR: baga.CallVariants.Caller needs a baga.AlignReads.SAMs object 
        with a "ready_BAMs" attribute. This can be obtained with the 
        "IndelRealignGATK()" method of the AlignReads module.
        ''')
            
            try:
                self.genome_sequence = alignments.genome_sequence
                self.genome_id = alignments.genome_id
            except AttributeError:
                print('''
        ERROR: baga.CallVariants.Caller needs a baga.AlignReads.SAMs object 
        with a "genome_sequence" attribute. This can be obtained by running all methods in the
        AlignReads module.
        ''')
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
        Save processed SAM file info to a local compressed pickle file.
        'name' should exclude extension: .baga will be added
        '''
        fileout = 'baga.CallVariants.Caller-%s.baga' % name
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

    def CallgVCFsGATK(self, 
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            local_variants_path = ['variants'],
            use_java = 'java',
            force = False,
            mem_num_gigs = 8, 
            max_cpus = -1):
        '''
        https://www.broadinstitute.org/gatk/guide/best-practices/?bpm=DNAseq
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
        # call the last set of ready BAMs added
        for cnum,BAM in enumerate(self.ready_BAMs[-1]):
            VCF_out = BAM[:-4] + '_unfiltered.gVCF'
            VCF_out = _os.path.sep.join([local_variants_path_genome, VCF_out.split(_os.path.sep)[-1]])
            if not _os.path.exists(VCF_out) or force:
                cmd = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar,
                '-T', 'HaplotypeCaller', '-R', genome_fna, '-I', BAM, #'-L', '20', 
                '--genotyping_mode', 'DISCOVERY',
                '--sample_ploidy', '1',
                '--heterozygosity', '0.0001',         # 650 total SNPs prior in all samples i.e., population <== make this an argument
                '--indel_heterozygosity', '0.00001',  # 65 total indels prior in all samples i.e., population
                '--emitRefConfidence', 'GVCF',        # make vcfs appropriate for doing GenotypeGVCFs after
                '--variant_index_type', 'LINEAR',
                '--variant_index_parameter', '128000',
                '-nct',  str(max_processes),
                '-stand_emit_conf', '10', 
                '-stand_call_conf', '20', 
                '-o', VCF_out]
                print('Called: %s' % (' '.join(map(str, cmd))))
                _subprocess.call(cmd)
                
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
            mem_num_gigs = 8):
        
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

        assert hasattr(self, 'paths_to_raw_gVCFs'), e1

        # use the last VCFs called
        open('variants.list', 'w').write('\n'.join(self.paths_to_raw_gVCFs[-1]))

        #rm /scratch/dw_temp/*_t_rawg.vcf

        # this method can be called prior or post base score recalibration
        # so give output a number corresponding to how many times variants called
        use_name = '{}_{}_samples_unfiltered.vcf'.format(data_group_name, len(self.paths_to_raw_gVCFs))
        VCF_out = _os.path.sep.join([local_variants_path_genome, use_name])

        cmd = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar, 
            '-T', 'GenotypeGVCFs', '-R', genome_fna,
            '--heterozygosity', '0.0001',         # 650 total indels prior in all samples i.e., population
            '--indel_heterozygosity', '0.00001',  # 65 total indels prior in all samples i.e., population
            '-stand_emit_conf', '10', 
            '-stand_call_conf', '20'] + \
            ['-V', 'variants.list'] + \
            ['-o', VCF_out]

        print('Called: %s' % (' '.join(cmd)))
        _subprocess.call(cmd)

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
        GenotypeGVCFsGATK()\n\
        method on this SAMs instance.'

        assert hasattr(self, 'path_to_unfiltered_VCF'), e1

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
            '--filterExpression', '"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"',
            '--filterName', '"standard_hard_filter"',
            '-o', hf_SNPs]

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
        GenotypeGVCFsGATK()\n\
        method on this SAMs instance.'

        assert hasattr(self, 'path_to_unfiltered_VCF'), e1

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
            '--filterExpression', '"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"',
            '--filterName', '"standard_indel_hard_filter"',
            '-o', hf_INDELs]

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
            if not _os.path.exists(table_out_pre) or force:
                cmd = [use_java, '-Xmx%sg' % mem_num_gigs, '-jar', jar,
                    '-T', 'BaseRecalibrator',
                    '-R', genome_fna,
                    '-I', BAM,
                    #'-L', '20',
                    '-nct',  str(max_processes),
                    '-knownSites', self.path_to_hardfiltered_SNPs[-1],
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
                    '-knownSites', self.path_to_hardfiltered_SNPs[-1],
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


if __name__ == '__main__':
    main()
