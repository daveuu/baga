#! /usr/bin/env python2
# -*- coding: utf-8 -*-
#
# This file is part of the Bacterial and Archaeal Genome Analyser
# Copyright (C) 2015-2016 David Williams
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
Structure module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains functions to detect structural rearrangements between a 
reference genome sequence and a query genome for which paired end reads for 
chromosome fragements are available. The intended use is to mark the regions 
likely affected by structural rearrangements for exclusion from a short read 
mapping experiment. Structurally altered regions, for example where a prophage 
has integrated into one genome but not the other, violate the assumption of the 
short read mapping method that reference and query sequences share 1-to-1 
orthology or homologous replacement.
'''

# stdlib
from baga import _os
from baga import _pickle
from baga import _gzip
from baga import _tarfile
from baga import _StringIO
from baga import _json
from baga import _subprocess
from baga import _sys
from baga import _logging

from glob import glob as _glob
from array import array as _array
import operator as _operator
import time as _time
import string as _string
from collections import Counter as _Counter

# external Python modules
import pysam as _pysam
from Bio import SeqIO as _SeqIO
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord

from baga import report_time as _report_time
from baga import get_exe_path as _get_exe_path
from baga import get_available_memory as _get_available_memory

from baga import MetaSample as _MetaSample
from baga import PROGRESS
from baga import PY3 as _PY3

from baga import load as _load
from baga import save as _save

if _PY3:
    _zip = zip
else:
    from itertools import izip as _zip

def main():
    pass


def moving_stats(values, window = 500, step = 1, resolution = 10):
    '''
    calculate:
        mean, 
        (commented out:
        sample variance (/(n-1)), 
        standard deviation)
    
    over window of specified width moving by specified step size. Input can have 
    a lower than 1:1 resolution e.g., 1 in 10 positions with resolution = 10
    '''
    mean = _array('f')
    #variance = _array('f')
    #st_dev = _array('f')
    for i in range(0, len(values) * resolution - window, resolution)[::step]:  #break
        # collect nt range allowing for resolution of ratios
        these = values[(i / resolution) : ((i + window) / resolution)]
        # exclude zero ratios from calculation of mean: 
        # associated with areas of low quality read alignments
        # also present at edges of large deletions so prevent them lowering window mean prematurely
        # require minimum half window length 
        these = [t for t in these if t > 0]
        if len(these) <= window / 2.0 / resolution:
            # when length of non-zero window is less than half that specified,
            # make up with zeros to decrease mean moderately close to edge of 
            # reference chromosome region with reads mapped
            these += [0] * int(window / 2.0 / resolution - len(these))
        
        # if len(these) >= 3:
            # these_mean = 0
        # else:
        these_mean = sum(these) / float(len(these))
        #these_variance = sum([(t - these_mean)**2 for t in these]) / (float(len(these) - 1))
        #these_st_dev = these_variance**0.5
        mean.append(these_mean)
        #variance.append(these_variance)
        #st_dev.append(these_st_dev)
    
    #return(mean, variance, st_dev)
    return(mean)


def loadCheckerInfo(filein):
    checker_info = {}
    with _tarfile.open(filein, "r:gz") as tar:
        for member in tar:
            #print(member.name)
            contents = _StringIO(tar.extractfile(member).read())
            try:
                # either json serialised conventional objects
                contents = _json.loads(contents.getvalue())
                
            except ValueError:
                # or longer python array.array objects
                if member.name == 'sequence':
                    contents = _array('c', contents.getvalue())
                elif 'ratio' in member.name:
                    contents = _array('f', contents.getvalue())
                else:
                    contents = _array('i', contents.getvalue())
            
            checker_info[member.name] = contents
    
    return(checker_info)

def checkStructure(BAMs, 
                   min_mapping_quality = 5, 
                   depth_resolution = 10, 
                   smoothed_resolution = 10, 
                   smoothed_window = False,
                   ratio_threshold = 0.15,
                   genome_name = False,
                   use_existing_bam_indexes = False,
                   force = False,
                   # if part of a pipeline with existing logger:
                   task_name = False,
                   console_verbosity_lvl = False,
                   log_folder = False):
    '''
    check for structural rearrangements . . .
    smoothed_window defaults to half of estimated fragment length
    '''

    # this is a module level function calling other objects from module

    module_name = __name__
    # pass log_folder on to objects.__init__()
    # also console_verbosity_lvl but not sure when used <== missing .set_verbosity bug? <==
    # analysis_path needed?

    if task_name:
        # if this has been called via the CLI it is part of a task
        logger = _logging.getLogger(task_name)
        # need to add adaptor after getting (handlers etc are retained from configureLogger)
        # in this namespace logger is the adaptor but not beyond
        logger = _logging.LoggerAdapter(logger, {'task': task_name})
    else:
        # else conform to conventional Python logging and use module name
        logger = _logging.getLogger(module_name)


    # instantiate Structure Checkers
    checkers = {}
    for BAM in BAMs:
        # ensure are indexed
        indexfile = _os.path.extsep.join([BAM,'bai'])
        if not use_existing_bam_indexes or \
                not(_os.path.exists(indexfile) and \
                _os.path.getsize(indexfile) > 0):
            # always re-index in case data has changed
            # unless explicitly told not to
            logger.info('indexing {}'.format(BAM))
            fail = False
            try:
                logger.debug('Using pySAM')
                _pysam.index(BAM)
            except TypeError as e:
                fail = True
                logger.debug("pysam didn't like something about this bam: {}"\
                        "".format(e))
            except IOError as e:
                fail = True
                logger.debug("pysam didn't like something about this bam: {}"\
                        "".format(e))
            if fail:
                path_to_exe = _get_exe_path('samtools')
                logger.debug("Trying '{} index {}' instead."\
                        "".format(path_to_exe,BAM))
                try:
                    logger.debug('Attempting to use SAMTools to index: {}'\
                            ''.format(BAM))
                    proc = _subprocess.Popen([path_to_exe,'index',BAM], 
                            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
                    stdout,stderr = proc.communicate()
                except OSError as e:
                    raise OSError("pysam and samtools didn't like something "\
                            "about this bam: {}".format(e))
                if 'fail to open' in stderr:
                    raise OSError("pysam and samtools didn't like something "\
                            "about this bam: {}".format(stderr))
        
        inherit_from = False
        if not force:
            # because not force (--force via CLI) can check if previous checker 
            # exists and that can be loaded from in which coverage depths etc.
            # are recorded
            # ===> this functionality might be duplicated in baga_cli.py
            reads = _pysam.Samfile(BAM, 'rb')
            # baga only compatible with one reads group per bam currently
            reads_name = reads.header['RG'][0]['ID']
            sample_name = '{}__{}'.format(reads_name, genome_name)
            if _os.path.exists('baga.Structure.Checker-{}.baga'.format(sample_name)):
                inherit_from = 'self'
        
        if inherit_from == 'self':
            logger.info('Will reloading from previous analysis: "{}"'.format(sample_name))
            this_checker = Checker(sample_name = sample_name, 
                    task_name = task_name, 
                    console_verbosity_lvl = console_verbosity_lvl, 
                    log_folder = log_folder, inherit_from = 'self')
        else:
            logger.info('Starting new analysis of "{}" with "{}"'\
                    ''.format(BAM, genome_name))
            this_checker = Checker(path_to_bam = BAM, genome_name = genome_name, 
                    task_name = task_name, 
                    console_verbosity_lvl = console_verbosity_lvl, 
                    log_folder = log_folder)
        checkers[this_checker.reads_name] = this_checker

    # mapped read pair inclusion:
    # not is_duplicate
    # not is_secondary
    # not is_qcfail

    start_time = _time.time()
    for cnum,(sample,checker) in enumerate(sorted(checkers.items())):
        logger.info('Collecting coverage for {} reads aligned to {}'\
                ''.format(sample, ', '.join(sorted(checker.genome_lengths))))
        # depth collection resolution impacts downstream resolutions
        # i.e., 1 in 10 here plus 1 in 10 smoothed ratio == 1 in 100 considered
        ## but why is this hardwired, but smoothed ratio provided with 10 as argument?
        ## what if these resolutions differ?
        ## check how getSmoothedRatios() uses smoothed ration carefully . . .
        depth_resolution = 10
        checker.getCoverageDepths(min_mapping_quality = min_mapping_quality, 
                resolution = depth_resolution, force = force)
        
        logger.info('Calculating mean insert size for {} reads aligned to {}'\
                ''.format(sample, ', '.join(sorted(checker.genome_lengths))))
        checker.getMeanInsertSize(force = force)
        if checker.mean_insert_size == -1:
            logger.info('Skipping: no aligned pairs found')
            continue
        logger.info('mean insert size == {:.1f}'.format(checker.mean_insert_size))
        logger.info('Collecting non-proper : proper pair ratios for {} reads aligned to {}'\
                ''.format(sample, ', '.join(sorted(checker.genome_lengths))))
        checker.getProperRatios(force = force)
        logger.info('Calculating smoothed ratios for {} reads aligned to {}'\
                ''.format(sample, ', '.join(sorted(checker.genome_lengths))))
        
        if not smoothed_window:
            # window for smoothing here is half estimated mean fragment length
            use_smoothed_window = int(round(checker.mean_insert_size)/2.0)
        else:
            use_smoothed_window = smoothed_window
        
        checker.getSmoothedRatios(window = use_smoothed_window, 
                resolution = smoothed_resolution)
        
        ## apparently not used?
        # omit zero depths and zero/infinite ratios for mean (~ omit NAs)
        # use_smoothed_ratios = {}
        # for name,info in checker.smoothed_ratios
            # use_smoothed_ratios = [m for m in info if m != 0]
            # mean_smoothed_ratios = sum(use_smoothed_ratios) / len(use_smoothed_ratios)
        
        # save for plotting
        checker.smoothed_resolution = smoothed_resolution
        checker.smoothed_window = use_smoothed_window
        checker.threshold = ratio_threshold
        
        #offset = int(round(depth_resolution * smoothed_resolution / 2.0))
        
        # also used in plotting
        offset = int(round(checker.smoothed_window / 2.0))
        
        filter_name = 'rearrangements'
        test = '>='
        checker.scan_threshold( checker.smoothed_ratios, 
                                threshold = checker.threshold, 
                                offset = offset,
                                test = test, 
                                filter_name = filter_name,
                                resolution = smoothed_resolution, 
                                buffer_distance = int(round(checker.mean_insert_size)))
        
        for chrm_name,suspect_regions in checker.suspect_regions[filter_name].items():
            t = len(suspect_regions)
            s = sum([(e - s) for s,e in suspect_regions])
            logger.info('For {} in sequence {} at threshold {} {}, '\
                    'excludes {} regions spanning {:,} basepairs'\
                    ''.format(sample, chrm_name, test, checker.threshold, t, s))
        
        # extend regions affected by sequence translocations if adjacent regions are 
        # without aligned reads or anomolous read alignment i.e., all proper-pairs
        filter_name = 'rearrangements_extended'
        filter_to_extend = 'rearrangements'
        # this was set to less than 0.5 because certain regions with read mapping gaps were 
        # not collected that were between two close disrupted regions
        checker.extend_regions(filter_name, filter_to_extend, 
                resolution = smoothed_resolution, threshold = 0.15)
        for chrm_name,suspect_regions in checker.suspect_regions[filter_name].items():
            t = len(suspect_regions)
            s = sum([(e - s) for s,e in suspect_regions])
            logger.info('For {} in sequence {} between high non-proper proportions, '\
                    'excludes {} regions spanning {:,} basepairs'.format(sample, 
                    chrm_name, t, s))
        
        checker.saveLocal()
        
        # report durations, time left etc <===
        # ==== this needs updating to work with logging
        _report_time(start_time, cnum, len(checkers))

    return(checkers)

class Checker(_MetaSample):
    '''
    The Checker class of the Structure module contains a method to check for 
    structural rearrangements between a reference genome sequence and a query 
    genome from which a set of paired end short reads are provided. Structural 
    rearrangements are detected by changes in the proportion of mapped reads 
    that were assigned to proper pairs by Burrows-Wheeler Aligner (BWA).
    '''

    def __init__(self, path_to_bam = False, genome_name = False, 
            sample_name = False, inherit_from = False, **kwargs):
        '''
        A structure checker object must be instantiated with either:
           
           sample_name of form <read group name>__<genome name>
           inferit_from = 'self'
        
        for restoring a previous analysis, or to start a new one from an existing 
        sorted BAM file: 
           
           path_to_bam
           genome_name
        '''
        module_name = __name__
        
        # argument combination options must be exclusive or for
        # A: start new analysis
        # b: load previous analysis
        option_A = ((path_to_bam and genome_name) and \
                not(sample_name and inherit_from))
        option_B = (not(path_to_bam and genome_name) and \
                (sample_name and inherit_from))
        if (option_A and option_B) or \
                (not(option_A) and not(option_B)):
            raise Exception('supply path_to_bam and genome_name or sample_name '\
                    'and inherit_from')
        
        if path_to_bam and genome_name:
            # this is not logged until super() so exceptions should be sent up to
            # CLI or whatever is calling module
            self.path_to_bam = path_to_bam
            self.reads = _pysam.Samfile(path_to_bam, 'rb')
            # baga only compatible with one reads group per bam currently
            self.reads_name = self.reads.header['RG'][0]['ID']
            sample_name = '{}__{}'.format(self.reads_name, genome_name)
            super(Checker, self).__init__(sample_name, module_name, 
                    inherit_from = inherit_from, **kwargs)
            
            # to retain crucial info without needing an open Samfile
            self.genome_lengths = {}
            for chrm_name,length in zip(self.reads.references, self.reads.lengths):
                self.genome_lengths[chrm_name] = length
            
            self.genome_name = genome_name
            
        elif inherit_from == 'self':
            super(Checker, self).__init__(sample_name, module_name, 
                    inherit_from = inherit_from, **kwargs)

    def regionProperProportion(self, start, end):
        '''of all fragments with at least two reads, what proportion does BWA consider "proper"?'''
        reads_iter = self.reads.fetch( self.reads.references[0], start, end)
        both_mapped_proper = set()
        both_mapped_notproper = set()
        for r in reads_iter:

            if not r.is_duplicate and not r.mate_is_unmapped:
                if r.is_proper_pair:
                    both_mapped_proper.add(r.query_name)
                else:
                    both_mapped_notproper.add(r.query_name)

        proportion = float(len(both_mapped_proper)) / \
                (len(both_mapped_proper)+len(both_mapped_notproper))
        #print('%s / %s = %.03f' % (len(both_mapped_proper), (len(both_mapped_proper)+len(both_mapped_notproper)), proportion))
        return(proportion)

    def getCoverageDepths(self, resolution = 10, 
                                min_mapping_quality = 30, 
                                force = False):
        '''
        collect per position along genome, aligned read coverage depth that pass a 
        quality standard and those that BWA considers additionally to be in a 
        "proper" pair.
        '''

        have_totals = []
        have_propers = []
        # check if needed data is present
        if not force:
            try:
                for seq_name,length in sorted(self.genome_lengths.items()):
                    if type(self.depths_totals[seq_name]) is _array and \
                            len(self.depths_totals[seq_name]) == \
                            int(length / self.depth_resolution) + 1:
                        have_totals += [True]
                        self.logger.info('Found total read depths for {} '\
                                'against {} from previous analysis'.format(
                                self.path_to_bam, seq_name))
            except (AttributeError, KeyError) as e:
                self.logger.debug('Could not find total read depths for {} '\
                        'against {} from a previous analysis: "{}"'.format(
                        self.path_to_bam, seq_name, e.message))
                pass
            
            # correct looking data found for all sequences?
            if len(have_totals):
                have_totals = all(have_totals)
            else:
                have_totals = False
            
            try:
                for seq_name,length in sorted(self.genome_lengths.items()):
                    if type(self.depths_propers[seq_name]) is _array and \
                            len(self.depths_propers[seq_name]) == \
                            int(length / self.depth_resolution) + 1:
                        have_propers += [True]
                        self.logger.info('Found proper-paired read depths for {} '\
                                'against {} from previous analysis'.format(
                                self.path_to_bam, seq_name))
            except (AttributeError, KeyError) as e:
                self.logger.debug('Could not find proper-paired read depths for '\
                        '{} against {} from a previous analysis: "{}"'.format(
                        self.path_to_bam, seq_name, e.message))
                pass
            
            # correct looking data found for all sequences?
            if len(have_propers):
                have_propers = all(have_propers)
            else:
                have_propers = False


        # check if need access to reads
        if (not all([have_totals, have_propers])) and not hasattr(self, 'reads'):
            try:
                self.reads = _pysam.Samfile(self.path_to_bam, 'rb')
            except AttributeError:
                raise Exception('Problems restoring data from disk, please '\
                        're-run with the force option')

        if not hasattr(self, 'depth_totals'):
            self.depths_totals = {}

        if not hasattr(self, 'depths_propers'):
            self.depths_propers = {}

        if force or not (have_totals and have_propers):
            if force and (have_totals or have_propers):
                self.logger.info('Recalculating read coverage depths '\
                        'because force is True')
            for h,header in enumerate(sorted(self.reads.header['SQ'])):
                num_ref_positions = header['LN']
                seq_name = header['SN']
                self.logger.info('Scanning {} in {} for read coverage depths'\
                        ''.format(seq_name, self.genome_name))
                depth_total = _array('i', (0,) * \
                        (int(num_ref_positions / resolution) + 1))
                depth_proper = _array('i', (0,) * \
                        (int(num_ref_positions / resolution) + 1))
                pl_iter = self.reads.pileup( seq_name )
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
                
                self.depths_totals[seq_name] = depth_total
                self.depths_propers[seq_name] = depth_proper
                # needed for plotting
                self.depth_resolution = resolution


    def getMeanInsertSize(self, upper_limit = 10000, 
                                min_mapping_quality = 30, 
                                force = False):
        '''
        omits pairs at either end of sequence break point of 
        circular chromosome (by ampplying an upper limit) that 
        would cause a few chromosome length insert sizes.
        '''

        if force or not hasattr(self, 'mean_insert_size'):
            if hasattr(self, 'mean_insert_size'):
                self.logger.debug('Recalculating mean insert size from BAM '\
                        'file because --force specified . . .')
            else:
                self.logger.debug('Calculating mean insert size from BAM '\
                        'file . . .')
            isizes = _array('i')
            for ref in sorted(self.reads.references):
                alignment = self.reads.fetch(ref)
                for read in alignment:
                    if read.is_proper_pair and \
                               read.is_read1 and \
                               not read.is_duplicate and \
                               not read.is_qcfail and \
                               read.mapping_quality >= min_mapping_quality:
                        
                        insert_length = abs(read.template_length)
                        if insert_length < upper_limit:
                            isizes.append(insert_length)
                
            try:
                mean_insert_size = sum(isizes) / float(len(isizes))
            except ZeroDivisionError:
                self.logger.warning('No paired, aligned reads found in {}'\
                        ''.format(self.reads.filename))
                self.logger.warning('Structure analysis for rearrangements is '\
                        'not possible')
                mean_insert_size = -1
                # raise exception?
            
            self.mean_insert_size = mean_insert_size
            
        elif hasattr(self, 'mean_insert_size'):
            if self.mean_insert_size == -1:
                self.logger.warning('Previously, no paired, aligned reads '\
                        'found in {}'.format(self.reads.filename))
                self.logger.warning('Structure analysis for rearrangements is '\
                        'not possible')
                # raise exception?
            else:
                self.logger.debug('Mean read pair insertsize found '\
                        'from previous analysis: use --force to recalculate')



    def getProperRatios(self, include_zeros = False, force = False):
        '''
        Calculates ratios non-proper pair to proper pair per position. Optionally 
        including regions with zero depths as zero ratio, defaults to excluding.
        Relatively more non-proper pair gives a higher ratio.
        '''
        # could pre-declare and use
        # n, total, proper in _izip(xrange(self.depth_total), self.depth_total, self.depth_proper)
        # ratios = _array('f', (-1,) * len(self.depth_total))

        # currently appending
        #a, b, c = 0, 0, 0

        if force or not hasattr(self, 'ratios'):
            if hasattr(self, 'ratios'):
                self.logger.debug('Recalculating proper-paired to total read '\
                        'depths ratios because --force specified.')
            else:
                self.logger.debug('Calculating proper-paired to total read '\
                        'depths ratios.')
            ratios = {}
            for seq_name in self.depths_totals:
                these_ratios = _array('f')
                for total, proper in _zip(self.depths_totals[seq_name], 
                        self.depths_propers[seq_name]):
                    if total >= proper > 0:
                        #a += 1
                        # both total aligned and proper pairs positive (present) and at least some nonpropers
                        # typically proper > nonproper, ratio 0.1-0.3
                        # as total non-proper approaches proper pairs, 
                        # within an insert length of disturbance,
                        # ratio approaches 1 (50:50 nonproper:proper)
                        # all proper, ratio == 0
                        these_ratios.append((total - proper) / float(proper))
                    elif total > proper == 0:
                        #b += 1
                        # no proper, set upper limit of ratio == total reads
                        # not important here because we are interested exceeding a 
                        # threshold of relatively many non-proper pair reads which this
                        # would.
                        these_ratios.append(total)
                    else:
                        #c += 1
                        # total == proper == 0
                        # zero aligned reads
                        # avoid zero division and distinguish from all proper pairs i.e., 0 / 100 == 0
                        these_ratios.append(-1)
                ratios[seq_name] = these_ratios
            self.ratios = ratios
        elif hasattr(self, 'ratios'):
            self.logger.debug('Proper-paired to total read depths ratios found '\
                    'from previous analysis: use --force to recalculate')



    def getSmoothedRatios(self, window = 500, resolution = 10, force = False):
        '''
        calculate:
            mean, 

        over window of specified width moving by specified step size. Input can have 
        a lower than 1:1 resolution e.g., 1 in 10 positions with resolution = 10
        '''
        # set zero depths to zero ratios for mean (~ omit NAs)
        # eventually omit from actual calculation (below)
        #use_ratios = _array('f',[r if r != -1 else 0 for r in self.ratios])

        # omit -1 no depths from ratios below (but include zeros)

        ### resolution == step size
        ### resolution here needs to match depth and/or ratios resolution previously?

        if force or not hasattr(self, 'smoothed_ratios'):
            if hasattr(self, 'smoothed_ratios'):
                self.logger.debug('Recalculating smoothed ratios of read '\
                        'depths because --force specified.')
            else:
                self.logger.debug('Calculating smoothed ratios of read depths.')
            smoothed_ratios = {}
            for seq_name in self.depths_totals:
                use_ratios = self.ratios
                these_means = _array('f')
                for i in range(0, len(use_ratios[seq_name]) * resolution - \
                        window, resolution):
                    # collect nt range allowing for resolution of ratios
                    start = int(i / resolution)
                    end = int((i + window) / resolution)
                    these = use_ratios[seq_name][start:end]
                    # exclude zero ratios from calculation of mean: 
                    # associated with areas of low quality read alignments
                    # also present at edges of large deletions so prevent
                    # them lowering window mean prematurely
                    # require minimum half window length 
                    these = [t for t in these if t >= 0]
                    if len(these) <= window / 2.0 / resolution:
                        # when length of non-zero window is less than half
                        # that specified, make up with zeros to decrease mean 
                        # moderately close to edge of reference chromosome
                        # region with reads mapped
                        these += [0] * int(window / 2.0 / resolution - len(these))
                    
                    these_mean = sum(these) / float(len(these))
                    these_means.append(these_mean)
                smoothed_ratios[seq_name] = these_means
            
            self.smoothed_ratios = smoothed_ratios
        elif hasattr(self, 'smoothed_ratios'):
            self.logger.debug('Smoothed ratios of read depths found '\
                    'from previous analysis: use --force to recalculate')




    def scan_threshold(self, values, threshold, offset, 
                             test = '>=', 
                             filter_name = 'threshold', 
                             resolution = 10, 
                             buffer_distance = 0):
        '''
        Reports ranges for regions that exceed a threshold.
        Values given as list, resolution and offset determines positions they 
        correspond to.
        Offset should be half the width of moving window, if used.
        Buffer_distance is excluded region at each end of chromosome sequence if 
        known to be circular.
        Threshold testing can be '>', '>=', '<' or '<='
        '''
        # sort out test function
        if test == '>=':
            
            def threshold_crossed(value):
                return(_operator.ge(value, threshold))
            
        elif test == '<=':
            
            def threshold_crossed(value):
                return(_operator.le(value, threshold))

        elif test == '<':
            
            def threshold_crossed(value):
                return(_operator.lt(value, threshold))

        elif test == '>':
            
            def threshold_crossed(value):
                return(_operator.gt(value, threshold))

        else:
            print('Warning: "test" can be ">", ">=", "<" or "<=", not "{}"'.format(test))
            
            def threshold_crossed(value):
                return(False)


        # if input is a (sparse) dict, make it list-like

        #for v in these_values[::10]: threshold_crossed(v)

        suspect_regions = {}
        for seq_name,these_values in values.items():
            # convert to array adding zeros if/when required
            if isinstance(these_values, dict):
                data_type = str(type(these_values.values()[0])).split("<type '")[1][0]
                use_values = _array(data_type)
                for pos1 in range(min(these_values), max(these_values), resolution):
                    try:
                        use_values.append(these_values[pos1])
                    except KeyError:
                        use_values.append(0)
            else:
                use_values = these_values
            
            these_suspect_regions = []
            positions = range(offset, len(these_values)*resolution + offset, 
                    resolution)
            # initialise collecting depending on state at beginning
            if threshold_crossed(use_values[0]):
                collecting = True
                these_suspect_regions += [positions[0]]
            else:
                collecting = False
            for pos,v in zip(positions, use_values):
                if threshold_crossed(v):
                    if not collecting:
                        # exceeds limit and not collecting so start
                        collecting = True
                        these_suspect_regions += [pos]
                else:
                    if collecting:
                        # below limit and collecting so stop
                        collecting = False
                        these_suspect_regions += [pos]
            
            if len(these_suspect_regions) % 2 != 0:
                # complete terminal range
                these_suspect_regions += [self.genome_lengths[seq_name]]
            
            these_suspect_regions = [these_suspect_regions[n:n+2] for n in \
                    range(0,len(these_suspect_regions),2)]
            
            # omit either end according to buffer_distance
            these_suspect_regions = [(s,e) for s,e in these_suspect_regions if \
                    s > buffer_distance and \
                    e < (len(use_values)*resolution + offset - buffer_distance)]
            
            suspect_regions[seq_name] = these_suspect_regions

        if hasattr(self, 'suspect_regions'):
            self.suspect_regions[filter_name] = suspect_regions
        else:
            self.suspect_regions = {filter_name : suspect_regions}


    def extend_regions(self,    filter_name, 
                                filter_to_extend, 
                                resolution = 10, 
                                threshold = 0.5):
        '''
        Given suspect regions for filtering, extend if adjacent windows have 
        majority positions with ratio -1 for no reads at default settings with
        threshold = 0.5.
        Lowering threshold so that fewer positions must have no reads makes filter
        more greedy.
        '''

        extensions = {}
        for seq_name,suspect_regions in self.suspect_regions[filter_to_extend].items():
            these_extensions = []
            for n in range(len(suspect_regions) - 1 ):
                right_edge = suspect_regions[n][1]
                next_left_edge = suspect_regions[n + 1][0]
                window_size = int(round(self.mean_insert_size))
                # extend from left to right
                join = True
                if right_edge < next_left_edge - window_size:
                    for p in range(right_edge, next_left_edge - window_size):
                        these_ratios = self.ratios[seq_name][p / resolution: (p + window_size) / resolution]
                        if len([r for r in these_ratios if r < 0]) < window_size / resolution * threshold:
                            join = False
                            break
                else:
                    these_ratios = self.ratios[seq_name][right_edge / resolution: next_left_edge / resolution]
                    effective_window_size = next_left_edge - right_edge
                    if len([r for r in these_ratios if r <= 0]) > effective_window_size / resolution * threshold:
                        these_extensions += [[right_edge,next_left_edge]]
                    
                    continue
                
                if join:
                    # got to next disrupted region: add joining non-aligned region
                    these_extensions += [[right_edge,next_left_edge]]
                else:
                    if p > right_edge:
                        # got past at least first window: store non-aligned region
                        extension_right_edge = p + int(round(window_size / 2.0))
                        these_extensions += [[right_edge, extension_right_edge]]
                        end = extension_right_edge
                        
                    else:
                        end = right_edge
                    
                    # extend from right to left (if necessary)
                    join = True
                    for p in range(end, next_left_edge - window_size)[::-1]:  #break
                        these_ratios = self.ratios[seq_name][p / resolution: (p + window_size) / resolution]
                        if len([r for r in these_ratios if r < 0]) < window_size / resolution * threshold:
                            join = False
                            break
                    if join:
                        # got to back next disrupted region: add joining non-aligned region
                        # (unlikely if join False above)
                        these_extensions += [[end, next_left_edge]]
                        
                    elif p < next_left_edge - window_size - 1:
                        # got past at least first window: store non-aligned region
                        these_extensions += [[p + int(round(window_size / 2.0)), next_left_edge]]
            
            extensions[seq_name] = these_extensions

        if hasattr(self, 'suspect_regions'):
            self.suspect_regions[filter_name] = extensions
        else:
            self.suspect_regions = {filter_name : extensions}

def plotRegions(BAMs, genome, 
        regions = False,
        force = False,
        # if part of a pipeline with existing logger:
        task_name = False,
        console_verbosity_lvl = False,
        log_folder = False):
    '''
    Plot all or some of the regions affected by rearrangements

    This function should be called on a per sample basis

    Parameters
    ----------
    regions : dict
        sequence names as keys and tuple integer ranges as values for regions to plot
    '''

    # this is a module level function calling other objects from module

    module_name = __name__
    # pass log_folder on to objects.__init__()
    # also console_verbosity_lvl but not sure when used <== missing .set_verbosity bug? <==
    # analysis_path needed?

    if task_name:
        # if this has been called via the CLI it is part of a task
        logger = _logging.getLogger(task_name)
        # need to add adaptor after getting (handlers etc are retained from configureLogger)
        # in this namespace logger is the adaptor but not beyond
        logger = _logging.LoggerAdapter(logger, {'task': task_name})
    else:
        # else conform to conventional Python logging and use module name
        logger = _logging.getLogger(module_name)


    checkers = {}
    for BAM in BAMs:
        reads = _pysam.Samfile(BAM, 'rb')
        # baga only compatible with one reads group per bam currently
        reads_name = reads.header['RG'][0]['ID']
        sample_name = '{}__{}'.format(reads_name, genome.sample_name)
        logger.info('Loading data for: "{}"'.format(sample_name))
        checker = Checker(sample_name = sample_name, 
                task_name = task_name, 
                console_verbosity_lvl = console_verbosity_lvl, 
                log_folder = log_folder, inherit_from = 'self')
        
        if not hasattr(checker, 'depths_propers'):
            raise Exception('Sample {} seems to be missing data for plotting: '\
                    'please check stucture first')
        else:
            checkers[reads_name] = checker

    # needed for plotting
    # checker.smoothed_resolution = smoothed_resolution
    # checker.threshold = ratio_threshold

    plot_folder = 'plots_structure'
    filter_name = 'high non-proper pairs'

    # assert checker_info['genome_name'] == genome.id, e
    outdir = _os.path.sep.join([plot_folder, genome.sample_name])
    if not _os.path.exists(outdir):
        _os.makedirs(outdir)

    if not regions:
        # all ranges for filtering for this sample
        # currently range selection is rearrangements filter, not including the extensions
        # do_ranges = checker_info['suspect_region']['rearrangements_extended']
        do_regions = {}
        for reads_name,checker in checkers.items():
            do_regions[reads_name] = {}
            for seq_name,these_ranges in checker.suspect_regions['rearrangements'].items():
                do_regions[reads_name][seq_name] = []
                for s,e in sorted(these_ranges):
                    # select consistant plotting region for comparison between samples
                    plot_chrom_start = int(round(s - 500 - 100, -3))
                    if s + 2500 > e:
                        plot_chrom_end = plot_chrom_start + 2500
                    else:
                        plot_chrom_end = int(round(e + 500 + 100, -3))
                    
                    do_regions[reads_name][seq_name] += [(plot_chrom_start,plot_chrom_end)]
    else:
        # specific region of specific sequence already provided
        # will plot for all BAMs provided
        for reads_name,checker in checkers.items():
            do_regions[reads_name] = regions

    for reads_name,checker in checkers.items():
        logger.info('Plotting filter regions for {} reads aligned to {}'\
                ''.format(reads_name, genome.sample_name))
        for seq_name,regions in do_regions[reads_name].items():
            for plot_chrom_start,plot_chrom_end in regions:
                plot_filename = '{:07d}_{:07d}_{}__{}__{}.svg'.format(
                        plot_chrom_start, plot_chrom_end, seq_name, 
                        checker.genome_name, checker.reads_name)
                plot_output_path = [outdir, plot_filename]
                plot_output_path = _os.path.sep.join(plot_output_path)
                logger.debug('Plotting to: {}'.format(plot_output_path))
                # one Plotter per plot . . would be better to supply plot_output_path
                # for each plotter.doPlot() ? would need moving code from __init__()
                # refactor if/wehen creating unified plotter for this and repeats
                plotter = Plotter(checker, genome, plot_output_path)
                plotter.doPlot(seq_name, plot_chrom_start, plot_chrom_end, 
                        panel = ((1,1),(1,1)), label = checker.reads_name)


    #return(checkers)

class Plotter:
    '''
    Plotter class of the Structure module contains methods to plot the regions likely 
    to have undergone structural rearrangements as found by an instance of the 
    Checker class.
    '''
    def __init__(self, checker, genome, plot_output_path,
        width_cm = 30, height_cm = 10, 
        viewbox_width_px = 1800, viewbox_height_px = 600,
        plot_width_prop = 0.8, plot_height_prop = 0.8, 
        white_canvas = True):
        '''
        Plot pairs of aligned homologous chromosome regions with percent 
        identity calculated over a moving window.
        
        genome: an instance of Genome from RepeatFinder for which repeats 
        have been inferred.
        
        plot_width_prop and plot_height_prop: Proportion of whole plot 
        area covered by actual plot, to allow space for labels.
        '''
        
        e = 'The provided genome ({}) does not seem to match the reference '\
                'sequence of the provided read alignment ({})'.format(
                genome.sample_name, checker.genome_name)
        assert genome.sample_name == checker.genome_name, e
        
        self.genome = genome
        
        import svgwrite as _svgwrite
        dwg = _svgwrite.Drawing(plot_output_path, width='%scm' % width_cm, height='%scm' % height_cm,
                                profile='full', debug=True)
        
        dwg.viewbox(width = viewbox_width_px, height = viewbox_height_px)
        
        if white_canvas:
            dwg.add(_svgwrite.shapes.Rect(insert=(0, 0), size=(viewbox_width_px, viewbox_height_px), fill = _svgwrite.rgb(100, 100, 100, '%')))
        
        self.viewbox_width_px = viewbox_width_px
        self.viewbox_height_px = viewbox_height_px
        self.plot_width_prop = plot_width_prop
        self.plot_height_prop = plot_height_prop
        self.checker = checker
        self.dwg = dwg
        self.rgb = _svgwrite.rgb

    def chrom2plot_x(self, pos_chrom, start, end):
        '''convert chromosome position to plotting position in canvas'''
        plotlen_chrom = end - start
        pos_plot = (pos_chrom - start) * ((self.viewbox_width_px * self.plot_width_prop) / plotlen_chrom)
        return(pos_plot)

    def plot_scale(self,    seq_name,
                            start, 
                            end, 
                            panel, 
                            num_ticks = 5,
                            tick_len_px = 20, 
                            colour=(0,0,0,'%'), 
                            stroke_width=3, 
                            font_size = '20pt', 
                            use_fontfamily = 'Nimbus Sans L',
                            plot_label = True):
        '''given the real position ordinate on chromosome, and range of plotting window, plot ticks with position'''

        # get tick positions on chromosomes to plot
        if end - start <= 2000:
            tick_rounding = 100
        elif end - start <= 20000:
            tick_rounding = 1000
        else:
            tick_rounding = 5000
            
        tick_dist = (end - start) / num_ticks
        if tick_rounding > tick_dist:
            tick_dist = tick_rounding
        else:
            rm = tick_dist % tick_rounding
            if (tick_rounding / 2.0) < rm:
                tick_dist += (tick_rounding - rm)
            else:
                tick_dist -= rm

        plotpositions_chrom = []
        for pos in range(0, len(self.genome.sequence[seq_name]), tick_dist):
            if start <= pos <= end:
                plotpositions_chrom += [pos]

        # convert to x-positions for plotting
        plotlen_chrom = end - start

        # plot into full viewbox or other specified panel
        ((this_x_panel, num_x_panels),(this_y_panel, num_y_panels)) = panel

        # calculate plotting area within viewbox
        # currently only single x panel implemented (full width of view box)
        plotstart_x = (self.viewbox_width_px - (self.viewbox_width_px * self.plot_width_prop)) / 2
        plotend_x = plotstart_x + (self.viewbox_width_px * self.plot_width_prop)
        if num_x_panels > 1:
            print('>1 x panel not implemented')
            return(False)

        # this area is just the data i.e., depth line plot
        # from which feature y positions calculated
        # start with pre-calculated
        plottop_y = self.plottop_y
        plotbottom_y = self.plotbottom_y
        if num_y_panels > 1:
            # update existing plot area in y direction
            # this uses num_y_panels for overall shape and this_y_panel for target
            y_per_panel = (plotbottom_y - plottop_y) / num_y_panels
            plottop_y = plotbottom_y - y_per_panel * this_y_panel
            plotbottom_y -= y_per_panel * (this_y_panel - 1)

        # get tick positions on x-axis to plot
        plotpositions_x = [self.chrom2plot_x(pos_chrom, start, end) for pos_chrom in plotpositions_chrom]
        # make tick font size slightly smaller than labels
        font_size_int = int(font_size.replace('pt',''))
        tick_font_size = '{}pt'.format(font_size_int * 0.85)

        for n,x in enumerate(plotpositions_x):
            MoveTo = "M %s %s" % (x + plotstart_x, plotbottom_y)
            Line = "L %s %s" % (x + plotstart_x, plotbottom_y + tick_len_px)
            # use 'z' to close path with a line
            self.dwg.add(
                self.dwg.path(
                    d="%s %s z" % (MoveTo, Line), 
                    stroke = self.rgb(*colour), 
                    stroke_linecap='round', 
                    fill = 'none', stroke_width = stroke_width
                    )
            )
            textanc = 'middle'
            if tick_rounding == 100:
                fmt = '%.01f kb'
            else:
                fmt = '%d kb'
            
            tiplabel = self.dwg.text(
                                fmt % (plotpositions_chrom[n]/1000.0), 
                                insert = (x + plotstart_x,plotbottom_y + tick_len_px * 2.1),
                                fill='black', 
                                font_family = use_fontfamily, 
                                text_anchor = textanc, 
                                font_size = tick_font_size)
            
            self.dwg.add(tiplabel)

        if plot_label:
            # label x-axis
            
            xaxislabel = self.dwg.text(
                                "Reference Chromosome Position", 
                                insert = ((plotstart_x + plotend_x)/2, plotbottom_y + tick_len_px * 2.1 * 2),
                                fill = 'black', 
                                font_family = use_fontfamily, 
                                text_anchor = "middle", 
                                font_size = font_size)
            
            self.dwg.add(xaxislabel)


    def plot_ORFs(self, seq_name, start, end, 
                        panel = ((1,1),(1,1)), 
                        stroke_width = 40, 
                        colour = (0,0,0,'%'), 
                        font_size = 15, 
                        use_fontfamily = 'Nimbus Sans L'):
        
        # upper < lower because SVG upside down <<<<<<<<<<<<<<<
        # plot into full viewbox or other specified panel
        ((this_x_panel, num_x_panels),(this_y_panel, num_y_panels)) = panel

        # calculate plotting area within viewbox
        # currently only single x panel implemented (full width of view box)
        plotstart_x = (self.viewbox_width_px - (self.viewbox_width_px * self.plot_width_prop)) / 2
        plotend_x = plotstart_x + (self.viewbox_width_px * self.plot_width_prop)
        if num_x_panels > 1:
            print('>1 x panel not implemented')
            return(False)

        # this area is just the data i.e., depth line plot
        # from which feature y positions calculated
        # start with pre-calculated
        plottop_y = self.plottop_y
        plotbottom_y = self.plotbottom_y
        if num_y_panels > 1:
            # update existing plot area in y direction
            # this uses num_y_panels for overall shape and this_y_panel for target
            y_per_panel = (plotbottom_y - plottop_y) / num_y_panels
            plottop_y = plotbottom_y - y_per_panel * this_y_panel
            plotbottom_y -= y_per_panel * (this_y_panel - 1)

        # collect names, strand and ranges of ORFs
        ORF_plot_info = []
        # ORFs; RNA; large
        for ID,(s,e,d,name) in self.genome.annotations[seq_name][0].items():
            status = False
            if start <= s and e < end:
                plot_x_s = self.chrom2plot_x(s, start, end) + plotstart_x
                plot_x_e = self.chrom2plot_x(e, start, end) + plotstart_x
                status = 'complete'
            elif s < start < e:
                plot_x_s = self.chrom2plot_x(start, start, end) + plotstart_x
                plot_x_e = self.chrom2plot_x(e, start, end) + plotstart_x
                status = 'left cut'
            elif s < end < e:
                plot_x_s = self.chrom2plot_x(s, start, end) + plotstart_x
                plot_x_e = self.chrom2plot_x(end, start, end) + plotstart_x
                status = 'right cut'
            
            if status:
                if len(name):
                    use_name = '{} ({})'.format(ID, name)
                else:
                    use_name = ID
                
                ORF_plot_info += [(plot_x_s, plot_x_e, d, use_name, status)]

        # some additional parameters here for tweaking layout of features
        feature_thickness = (plotbottom_y - plottop_y) * 0.07
        point_width = (plotend_x - plotstart_x) * 0.008
        # half a feature thickness above scale, one feature thickness for reverse strand, half a feature thickness above reverse strand
        forward_y_offset = feature_thickness * 0.5 + feature_thickness + feature_thickness * 0.5
        reverse_y_offset = feature_thickness * 0.5

        # plot feature lane guide lines
        commands = ["M %s %s" % (
                                            plotstart_x, 
                                            plotbottom_y - reverse_y_offset - feature_thickness * 0.5
                                            )]
        commands += ["L %s %s" % (
                                            plotend_x, 
                                            plotbottom_y - reverse_y_offset - feature_thickness * 0.5
                                            )]
        self.dwg.add(
            self.dwg.path(
                d=' '.join(commands), 
                stroke=self.rgb(70, 70, 70,'%'), stroke_linecap='round', 
                fill = 'none', stroke_width = 5
            )
        )

        commands = ["M %s %s" % (
                                            plotstart_x, 
                                            plotbottom_y - forward_y_offset - feature_thickness * 0.5
                                            )]
        commands += ["L %s %s" % (
                                            plotend_x, 
                                            plotbottom_y - forward_y_offset - feature_thickness * 0.5
                                            )]
        self.dwg.add(
            self.dwg.path(
                d=' '.join(commands), 
                stroke=self.rgb(70, 70, 70,'%'), stroke_linecap='round', 
                fill = 'none', stroke_width = 5
            )
        )

        # plot labels
        textanc = 'end'
        horizontal_offset = 20
        label = "Forward strand"
        feature_lane_label = self.dwg.text(
                                            label, 
                                            insert = (
                                                    plotstart_x - horizontal_offset, 
                                                    plotbottom_y - (forward_y_offset + feature_thickness * 0.5)
                                                    ), 
                                             fill = 'black', font_family = use_fontfamily, 
                                             text_anchor = textanc, font_size = '%dpt' % font_size, 
                                             baseline_shift='-50%'
                                             )
        self.dwg.add(feature_lane_label)

        label = "Reverse strand"
        feature_lane_label = self.dwg.text(
                                            label, 
                                            insert = (
                                                plotstart_x - horizontal_offset, 
                                                plotbottom_y - (reverse_y_offset + feature_thickness * 0.5)
                                                ), 
                                             fill='black', font_family=use_fontfamily, 
                                             text_anchor=textanc, font_size = '%dpt' % font_size,
                                             baseline_shift='-50%'
                                             )
        self.dwg.add(feature_lane_label)

        ORF_plot_info_forward_strand = sorted(
                                                [o for o in ORF_plot_info if o[2] == 1],
                                                reverse = True
                                                )
        ORF_plot_info_reverse_strand = sorted(
                                                [o for o in ORF_plot_info if o[2] == -1],
                                                reverse = False
                                                )

        for s,e,d,name,state in ORF_plot_info_forward_strand + ORF_plot_info_reverse_strand:
            if d == 1:
                start,end = s,e
                # forward strand
                #if start ==  plotstart_x:
                if state == 'left cut':
                    # 'left cut', start lower left plus a bit for cut angle
                    commands = ["M %s %s" % (
                                                start - point_width, 
                                                plotbottom_y - forward_y_offset
                                                )]
                else:
                    # start lower left at feature start
                    commands = ["M %s %s" % (
                                                start, 
                                                plotbottom_y - forward_y_offset
                                                )]
                
                #if end == plotend_x:
                if state == 'right cut':
                    # 'right cut', go to lower right then top right plus a bit for angle
                    commands += ["L %s %s" % (
                                                end, 
                                                plotbottom_y - forward_y_offset
                                                )]
                    commands += ["L %s %s" % (
                                                end + point_width, 
                                                plotbottom_y - (forward_y_offset + feature_thickness * 1)
                                                )]
                else:
                    # else point on right for forward strand
                    commands += ["L %s %s" % (
                                                end - point_width, 
                                                plotbottom_y - forward_y_offset
                                                )]
                    commands += ["L %s %s" % (
                                                end, 
                                                plotbottom_y - (forward_y_offset + feature_thickness * 0.5)
                                                )]
                    commands += ["L %s %s" % (
                                                end - point_width, 
                                                plotbottom_y - (forward_y_offset + feature_thickness * 1)
                                                )]
                
                commands += ["L %s %s z" % (
                                                start, 
                                                plotbottom_y - (forward_y_offset + feature_thickness * 1)
                                                )]
                use_y_offset = forward_y_offset
                
            else:
                # reverse strand
                start,end = e,s
                #if start == plotend_x:
                if state == 'right cut':
                    # 'right cut', go upper right plus a bit for cut angle
                    commands = ["M %s %s" % (
                                                start + point_width, 
                                                plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                                                )]
                else:
                    # start upper right plus a bit for cut angle
                    commands = ["M %s %s" % (
                                                start,
                                                plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                                                )]
                
                #if end == plotstart_x:
                if state == 'left cut':
                    # 'left cut', go to upper left
                    commands += ["L %s %s" % (
                                                end, 
                                                plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                                                )]
                    # plus a bit lower left
                    commands += ["L %s %s" % (
                                                end - point_width, 
                                                plotbottom_y - reverse_y_offset
                                                )]
                else:
                    # else point on left for reverse strand
                    commands += ["L %s %s" % (
                                                end + point_width, 
                                                plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                                                )]
                    commands += ["L %s %s" % (
                                                end, 
                                                plotbottom_y - (reverse_y_offset + feature_thickness * 0.5)
                                                )]
                    commands += ["L %s %s" % (
                                                end + point_width, 
                                                plotbottom_y - reverse_y_offset
                                                )]
                
                commands += ["L %s %s z" % (
                                                start, 
                                                plotbottom_y - reverse_y_offset
                                                )]
                
                use_y_offset = reverse_y_offset
                
            self.dwg.add(
                self.dwg.path(
                    d=' '.join(commands), stroke='white', stroke_linecap='round', 
                    fill = self.rgb(*colour), stroke_width = 1
                )
            )
            added_label = self.dwg.add(
                self.dwg.text(
                    name,
                    insert = (
                        (start+end)/2,
                        plotbottom_y - (use_y_offset + feature_thickness * 0.5)
                        ), 
                    fill='white', stroke='black', stroke_width=0.4, 
                    font_family=use_fontfamily, font_weight='bold',
                    text_anchor='middle', font_size = '%dpt' % font_size, 
                    baseline_shift='-50%'
                )
            )
            # this would vary with font and font size
            max_length_for_unrotated_label = self.viewbox_width_px * 1.0/18
            if max(start,end) - min(start,end) < max_length_for_unrotated_label:
                x = (start+end)/2
                y = plotbottom_y - (use_y_offset + feature_thickness * 0.5)
                added_label.rotate(-25, center = (x, y))

    def plot_LargeFeatures(self, seq_name, start, end, panel = ((1,1),(1,1)),
                          stroke_width = 40, colour = (20, 20, 20,'%'), 
                          font_size = 15, use_fontfamily = 'Nimbus Sans L'):
        
        # plot into full viewbox or other specified panel
        ((this_x_panel, num_x_panels),(this_y_panel, num_y_panels)) = panel

        # calculate plotting area within viewbox
        # currently only single x panel implemented (full width of view box)
        plotstart_x = (self.viewbox_width_px - (self.viewbox_width_px * self.plot_width_prop)) / 2
        plotend_x = plotstart_x + (self.viewbox_width_px * self.plot_width_prop)
        if num_x_panels > 1:
            print('>1 x panel not implemented')
            return(False)


        # this area is just the data i.e., depth line plot
        # from which feature y positions calculated
        # start with pre-calculated
        plottop_y = self.plottop_y
        plotbottom_y = self.plotbottom_y
        if num_y_panels > 1:
            # update existing plot area in y direction
            # this uses num_y_panels for overall shape and this_y_panel for target
            y_per_panel = (plotbottom_y - plottop_y) / num_y_panels
            plottop_y = plotbottom_y - y_per_panel * this_y_panel
            plotbottom_y -= y_per_panel * (this_y_panel - 1)


        Features_to_plot = []
        #for feat_chrom_start,feat_chrom_end,n in Features:
        for name, (feat_chrom_start, feat_chrom_end) in self.genome.annotations[seq_name][2].items():
            # don't plot unless part of feature is within plotting range
            if feat_chrom_start < end and feat_chrom_end > start:
                s = self.chrom2plot_x(max( start, feat_chrom_start), start, end ) + plotstart_x
                e = self.chrom2plot_x(min( end, feat_chrom_end), start, end ) + plotstart_x
                Features_to_plot += [(s, e, name)]


        # some more parameters for tweaking feature layout here
        # same set in ORF function
        feature_thickness = (plotbottom_y - plottop_y) * 0.07
        # used if feature cut at either end
        point_width = (plotend_x - plotstart_x) * 0.008
        # half a feature thickness above scale, 
        # two feature thickness for strands plus a half inbetween, 
        # half a feature thickness above forward strand
        y_offset = feature_thickness * 0.5 + feature_thickness + \
                         feature_thickness * 0.5 + feature_thickness + \
                         feature_thickness * 0.5

        # plot feature lane guide lines
        commands = ["M %s %s" % (
                                    plotstart_x, 
                                    plotbottom_y - y_offset - feature_thickness * 0.5
                                    )]
        commands += ["L %s %s" % (
                                    plotend_x, 
                                    plotbottom_y - y_offset - feature_thickness * 0.5
                                    )]
        self.dwg.add(
            self.dwg.path(
                d=' '.join(commands), stroke=self.rgb(70, 70, 70,'%'), 
                stroke_linecap='round', fill = 'none', stroke_width = 5
            )
        )

        # plot labels
        textanc = 'end'
        horizontal_offset = 20
        label = "Large Features"
        feature_lane_label = self.dwg.text(
                            label, 
                            insert = (
                                plotstart_x - horizontal_offset, 
                                plotbottom_y - (y_offset + feature_thickness * 0.5)
                            ), 
                            fill='black', font_family=use_fontfamily, 
                            text_anchor=textanc, 
                            font_size = '%dpt' % font_size, 
                            baseline_shift='-50%'
                        )
        self.dwg.add(feature_lane_label)

        # plot feature
        for s, e, name in Features_to_plot:
            start,end = s,e
            if start ==  plotstart_x:    # 'left cut'
                commands = ["M %s %s" % (
                                            start - point_width, 
                                            plotbottom_y - y_offset
                                            )]
            else:
                commands = ["M %s %s" % (
                                            start, 
                                            plotbottom_y - y_offset
                                            )]
            
            commands += ["L %s %s" % (
                                            end, 
                                            plotbottom_y - y_offset
                                            )]
            
            if end == plotend_x:   # 'right cut':
                commands += ["L %s %s" % (
                                            end + point_width, 
                                            plotbottom_y - (y_offset + feature_thickness * 1)
                                            )]
            else:
                commands += ["L %s %s" % (
                                            end, 
                                            plotbottom_y - (y_offset + feature_thickness * 1)
                                            )]
            
            commands += ["L %s %s z" % (
                                            start, 
                                            plotbottom_y - (y_offset + feature_thickness * 1)
                                            )]
            
            self.dwg.add(
                self.dwg.path(
                    d=' '.join(commands), 
                    stroke='white', stroke_linecap='round', 
                    fill = self.rgb(*colour), stroke_width = 1
                )
            )
            self.dwg.add(
                self.dwg.text(
                    name, 
                    insert = (
                        (start+end)/2, 
                        plotbottom_y - (y_offset + feature_thickness * 0.5)
                    ), 
                    fill='white', stroke='black', stroke_width=0.2, 
                    font_family=use_fontfamily, font_weight='bold',
                    text_anchor='middle', font_size = '%dpt' % font_size, 
                    baseline_shift='-50%'
                    ))

    def calc_plot_region(self, start, 
                          end, 
                          panel = ((1,1),(1,1)),
                          values_upper_prop = 0.75,
                          values_lower_prop = 0.4):
        
        '''
        given the real position ordinate on chromosome and panel for plotting
        calculate region for plotting.
        '''
        
        # plot into full viewbox or other specified panel
        ((this_x_panel, num_x_panels), (this_y_panel, num_y_panels)) = panel

        # calculate plotting area within viewbox
        # currently only single x panel implemented (full width of view box)
        # this area includes all: depth, scale, ORFs, features etc
        plotstart_x = (self.viewbox_width_px - (self.viewbox_width_px * self.plot_width_prop)) / 2
        plotend_x = plotstart_x + (self.viewbox_width_px * self.plot_width_prop)
        if num_x_panels > 1:
            print('>1 x panel not implemented')
            return(False)

        plottop_y = (self.viewbox_height_px - (self.viewbox_height_px * self.plot_height_prop)) / 2
        plotbottom_y = plottop_y + (self.viewbox_height_px * self.plot_height_prop)

        # to prevent x-axis label going off the bottom
        # probably a better way of achieving this
        shift_up_offset = self.viewbox_height_px * self.plot_height_prop * 0.1
        plottop_y -= shift_up_offset
        plotbottom_y -= shift_up_offset

        # this area is just the data i.e., depth line plot
        if num_y_panels > 1:
            # update existing plot area in y direction
            # this uses num_y_panels for overall shape and this_y_panel for target
            y_per_panel = (plotbottom_y - plottop_y) / num_y_panels
            plottop_y = plotbottom_y - y_per_panel * this_y_panel
            plotbottom_y -= y_per_panel * (this_y_panel - 1)

        # lane height . . for plotting multiple data sets over same region.
        upper_y = plotbottom_y - (plotbottom_y - plottop_y) * values_upper_prop
        lower_y = plotbottom_y - (plotbottom_y - plottop_y) * values_lower_prop
        # Single lane for comparing different regions in one plot.
        y_per_lane = (lower_y - upper_y) / 1

        self.plotstart_x = plotstart_x
        self.plotend_x = plotend_x
        self.upper_y = upper_y
        self.lower_y = lower_y
        self.y_per_lane = y_per_lane
        self.plottop_y = plottop_y
        self.plotbottom_y = plotbottom_y

    def plot_values(self, values, 
                          start, 
                          end, 
                          label, 
                          offset = 0,
                          resolution = 10,
                          panel = ((1,1),(1,1)),
                          colour = ('60', '60', '60', '%'), 
                          fill = True,
                          # 'Left' or 'Right'
                          plot_axis = False,
                          max_y_scale = False,
                          num_ticks = 5,
                          add_sample_label = True,
                          y_axis_label = 'Aligned Reads',
                          font_size = '20pt', 
                          use_fontfamily = 'Nimbus Sans L'):
        
        '''
        given the real position ordinate on chromosome with e.g., read depth data, 
        plot line. Also plot y-axis and scale
        Values must be either:
            a dictionary: chromosome position => value
        or:
            an ordered iterable like a list or array. If the list len is less than 
            the chromosome, the resolution must account for the difference e.g., 
            resolution = 10 if the list lenth is 1/10 of the chromosome length
        
        use offset as half moving window width if required
        '''
        
        ## could check for type tuple here for plotting ratios

        if isinstance(values, dict):
            #region_values = dict([(p,d) for p,d in use_depths.items() if start <= p <= end])
            region_values = {}
            # pileup doesn't return zero depths so add them here
            for pos1 in xrange(start - offset, end - offset):
                if pos1 % resolution == 0:
                    try:
                        region_values[pos1 + offset] = values[pos1]
                    except KeyError:
                        region_values[pos1 + offset] = 0
            
            plotpositions_chrom = sorted(region_values)
            
        elif isinstance(values, list) or isinstance(values, _array):
            region_values = {}
            plotpositions_chrom = []
            for pos1 in xrange(start - offset, end - offset):
                if pos1 % resolution == 0:
                    try:
                        region_values[pos1 + offset] = values[int(pos1/resolution)]
                        plotpositions_chrom += [pos1 + offset]
                    except IndexError:
                        # should be off the end of the chromosome
                        pass

        # get positions on x-axis to plot
        plotpositions_x = [self.chrom2plot_x(pos_chrom, start, end) \
                for pos_chrom in plotpositions_chrom]

        #### plot values ####

        # get corresponding depths over this chromosome region
        # normalise
        if not max_y_scale:
            # unless provided . . .
            max_plot_depth = max(region_values.values()) * 1
            # round up to the nearest 10 reads
            max_y_scale = round((max(region_values.values()) + 10) / 10.0) * 10
            
        plotdepths = []
        for pos0 in plotpositions_chrom:
            d = region_values[pos0]
            if d >= max_y_scale:
                use_d = 1.0
            else:
                use_d = d / float(max_y_scale)
            
            plotdepths += [use_d]

        # make paths for each isolate plot lane 
        lane_num = 0

        # establish lane plot area (always single for percent ID at duplicate regions)
        plot_spacing = 3
        lane_upper_y = self.lower_y - (lane_num + 1) * self.y_per_lane + plot_spacing
        lane_lower_y = self.lower_y - lane_num * self.y_per_lane - plot_spacing
        lane_plot_height = lane_lower_y - lane_upper_y
        #print(i, lane_lower_y, lane_lower_y)


        if fill:
            # start at lower left then move up to first position
            commands = ['M %s %s' % (
                                        plotpositions_x[0] + self.plotstart_x, 
                                        lane_lower_y
                                        )]
        else:
            commands = ['M %s %s' % (
                                        plotpositions_x[0] + self.plotstart_x, 
                                        lane_lower_y - lane_plot_height * plotdepths[0]
                                        )]

        for n,d in enumerate(plotdepths):
            # because everything is upside down in SVG minus goes up on page
            commands += ['L %s %s' % (
                                        plotpositions_x[n] + self.plotstart_x,
                                        lane_lower_y - lane_plot_height * d
                                        )]

        if fill:
            use_stroke_width = '0'
            use_fill_colour = self.rgb(*colour)
            # finish at lower right then close
            commands += ['L %s %s z' % (
                                        plotpositions_x[-1] + self.plotstart_x, 
                                        lane_lower_y
                                        )]
        else:
            use_stroke_width = '3'
            use_fill_colour = 'none'
            # plot ended at last point without closing

        plot_path = self.dwg.path(
                            d=' '.join(commands), stroke = self.rgb(*colour), 
                            stroke_linecap='round', stroke_width= use_stroke_width, 
                            fill = use_fill_colour, fill_rule='evenodd'
                            )

        self.dwg.add(plot_path)


        #### plot axis ####

        label_horizontal_offset = 0
        if plot_axis:
            if plot_axis == 'Right':
                use_x_offset = self.plotend_x
                direction = 1
                label_rot = 90
                textanc = 'start'
            else:
                # defaults to left side
                use_x_offset = self.plotstart_x
                direction = -1
                label_rot = 270
                textanc = 'end'
            
            tick_offset = 20 * direction
            # set 2% space between axis and plotting area
            use_x_offset += (self.plotend_x - self.plotstart_x) * 0.02 * direction
            
            commands = ['M %s %s' % (use_x_offset, self.upper_y), 'L %s %s' % (use_x_offset, self.lower_y)]
            plot_path = self.dwg.path(
                                d=' '.join(commands), stroke = 'black', 
                                stroke_linecap='round', stroke_width= '3',
                                fill = 'none', fill_rule='evenodd'
                                )
            
            self.dwg.add(plot_path)
            
            # ticks <== vary by divisibility of max_y_scale (implement auto feature?)
            # num_ticks = 5.0
            
            if max_y_scale > num_ticks**2:
                make_integers = True
            else:
                make_integers = False
            
            # make tick font size slightly smaller than labels
            font_size_int = int(font_size.replace('pt',''))
            tick_font_size = '{}pt'.format(font_size_int * 0.85)
            
            for n in range(int(num_ticks) + 1):
                commands = ['M %s %s' % (
                                            use_x_offset, 
                                            self.lower_y - (self.lower_y - self.upper_y) * (n/float(num_ticks))
                                        ), 
                                        'L %s %s' % (
                                            use_x_offset + tick_offset, 
                                            self.lower_y - (self.lower_y - self.upper_y) * (n/float(num_ticks))
                                        )]
                
                plot_path = self.dwg.path(
                                    d=' '.join(commands), 
                                    stroke = 'black', 
                                    stroke_linecap='round', 
                                    stroke_width= '3',
                                    fill = 'none', 
                                    fill_rule='evenodd'
                                    )
                
                self.dwg.add(plot_path)
                
                value = round(n*(max_y_scale/num_ticks), 2)
                if make_integers:
                    value = int(value)
                
                ticklabel = self.dwg.text('%s' % value, 
                                                  insert = (
                                                    use_x_offset + tick_offset * 1.5, 
                                                    self.lower_y - (self.lower_y - self.upper_y) * (n/float(num_ticks))
                                                  ),
                                                  fill='black', 
                                                  font_family=use_fontfamily, 
                                                  text_anchor=textanc, 
                                                  font_size = tick_font_size, 
                                                  baseline_shift='-50%'
                            )
                
                self.dwg.add(ticklabel)
            
            # label y-axis
            label_horizontal_offset = 0
            textanc = 'middle'
            x, y = use_x_offset - label_horizontal_offset, lane_lower_y - (self.lower_y - self.upper_y) * 0.5
            x_pos = x + 90 * direction
            yaxislabel = self.dwg.text(
                                    '', 
                                    insert = (x_pos, y), 
                                    fill='black', font_family=use_fontfamily, 
                                    text_anchor=textanc, font_size = font_size)  # , baseline_shift='-50%'
            
            if isinstance(y_axis_label, str):
                yaxislabel.add(self.dwg.tspan(y_axis_label))
            else:
                # must be tuple
                # InkScape seems not to support em (font height) unites . . . .
                # but knowning font size in pt seems just as good
                font_size_int = int(font_size.replace('pt',''))
                
                # offset away from axis if additional lines in label
                # should probably do a proper pt to px translation?
                extra_dist = len(y_axis_label) * font_size_int * 0
                
                x_pos += extra_dist * direction
                yaxislabel.add(self.dwg.tspan(y_axis_label[0], x = [x_pos], dy = ['-{}pt'.format(font_size_int * 1.2)]))
                yaxislabel.add(self.dwg.tspan(y_axis_label[1], x = [x_pos], dy = ['{}pt'.format(font_size_int * 1.2)]))
            
            added_yaxislabel = self.dwg.add(yaxislabel)
            added_yaxislabel.rotate(label_rot, center = (x_pos, y))

        if add_sample_label:
            # label duplicate (A, B, etc)
            textanc = 'middle'
            isolatelabel = self.dwg.text(label, 
                                        insert = (
                                            self.plotstart_x - label_horizontal_offset * 4, 
                                            lane_lower_y - (self.lower_y - self.upper_y) * 1.3
                                        ), 
                                        fill='black', font_family=use_fontfamily, 
                                        text_anchor=textanc, 
                                        font_size = font_size, baseline_shift='-50%')
            
            self.dwg.add(isolatelabel)

        return(max_y_scale)

    def plot_suspect_regions(self, 
                              seq_name,
                              start, end, 
                              suspect_regions, 
                              panel = ((1,1),(1,1)),
                              colour = ('80', '10', '10', '%')):

        '''given the real position ordinate on chromosome with suspicious regions, plot them.'''
        
        # establish lane plot area (always single for percent ID at duplicate regions)
        lane_num = 0
        plot_spacing = 3
        lane_upper_y = self.lower_y - (lane_num + 1) * self.y_per_lane + plot_spacing
        lane_lower_y = self.lower_y - lane_num * self.y_per_lane - plot_spacing
        lane_plot_height = lane_lower_y - lane_upper_y

        # draw ambiguous (translocated) regions
        for chrm_s,chrm_e in suspect_regions[seq_name]:
            if chrm_e > start and chrm_s < end:
                
                s = self.chrom2plot_x(max(chrm_s, start), start, end)
                e = self.chrom2plot_x(min(chrm_e, end), start, end)
                
                commands = ['M %s %s' % (
                                        s + self.plotstart_x, 
                                        lane_lower_y
                                        )]
                commands += ['L %s %s' % (
                                        s + self.plotstart_x, 
                                        lane_upper_y
                                        )]
                commands += ['L %s %s' % (
                                        e + self.plotstart_x, 
                                        lane_upper_y
                                        )]
                commands += ['L %s %s z' % (
                                        e + self.plotstart_x, 
                                        lane_lower_y
                                        )]
                plot_path = self.dwg.path(
                                    d=' '.join(commands), stroke = self.rgb(*colour), 
                                    stroke_linecap='round', stroke_width= '3', 
                                    fill = self.rgb(*colour), fill_rule='evenodd',
                                    fill_opacity = 0.2
                                    )
                
                self.dwg.add(plot_path)

    def plot_threshold(self, 
                              threshold,
                              max_y_scale,
                              panel = ((1,1),(1,1)),
                              colour = ('0', '100', '0', '%')):

        '''given the real position ordinate on chromosome with suspiscious regions, plot them.'''
        
        # establish lane plot area (always single for percent ID at duplicate regions)
        lane_num = 0
        plot_spacing = 3
        lane_upper_y = self.lower_y - (lane_num + 1) * self.y_per_lane + plot_spacing
        lane_lower_y = self.lower_y - lane_num * self.y_per_lane - plot_spacing
        lane_plot_height = lane_lower_y - lane_upper_y


        # a "y =" line
        # blue for now
        colour = ('80', '10', '10', '%')
        commands = ['M %s %s' % (
                                self.plotstart_x, 
                                self.lower_y - (self.lower_y - self.upper_y) * threshold * (1.0 / max_y_scale)
                                #lane_lower_y - (lane_lower_y - lane_upper_y) * threshold  <== this required if >1 lane?
                                )]
        commands += ['L %s %s' % (
                                self.plotend_x, 
                                self.lower_y - (self.lower_y - self.upper_y) * threshold * (1.0 / max_y_scale)
                                )]
        plot_path = self.dwg.path(
                            d=' '.join(commands), stroke = self.rgb(*colour), 
                            stroke_linecap='round', stroke_width= '3', stroke_opacity = 0.5, 
                            fill = self.rgb(*colour), fill_rule='evenodd'
                            )

        self.dwg.add(plot_path)

    def doPlot(self,    seq_name,
                        plot_chrom_start, 
                        plot_chrom_end, 
                        panel = ((1,1),(1,1)), 
                        label = False, 
                        axis_ratio_max = 1.2):
        
        # reuse max_y_scale for multiple plot on same axis
        self.calc_plot_region(              plot_chrom_start, 
                                            plot_chrom_end, 
                                            panel = panel)

        # plot scale
        self.plot_scale(seq_name, plot_chrom_start, plot_chrom_end, panel)
        # plot ORFs
        self.plot_ORFs(seq_name, plot_chrom_start, plot_chrom_end, panel)
        # plot large features
        self.plot_LargeFeatures(seq_name, plot_chrom_start, plot_chrom_end, panel)


        max_y_scale = self.plot_values(     self.checker.depths_totals[seq_name], 
                                            plot_chrom_start, 
                                            plot_chrom_end, 
                                            label, 
                                            panel = panel, 
                                            colour = ('55', '55', '55', '%'), 
                                            plot_axis = 'Left',
                                            y_axis_label = 'Aligned Reads',
                                            add_sample_label = True)

        self.plot_values(                   self.checker.depths_propers[seq_name], 
                                            plot_chrom_start, 
                                            plot_chrom_end, 
                                            label, 
                                            panel = panel, 
                                            colour = ('75', '75', '75', '%'), 
                                            plot_axis = False, 
                                            max_y_scale = max_y_scale,
                                            add_sample_label = False)

        ## plot raw ratios
        self.plot_values(  [r if r != -1 else 0 for r in self.checker.ratios[seq_name]], 
                              plot_chrom_start, 
                              plot_chrom_end, 
                              label, 
                              panel = panel,
                              colour = ('0', '100', '0', '%'), 
                              num_ticks = 4,
                              fill = False,
                              plot_axis = 'Right',
                              y_axis_label = ('Non-Proper Pair :','Proper Pair Reads'),
                              max_y_scale = axis_ratio_max,
                              add_sample_label = False)

        ## plot smoothed ratios
        #offset = int(round( self.checker.depth_resolution * \
        #        self.checker.smoothed_resolution / 2.0))

        offset = int(round(self.checker.smoothed_window / 2.0))

        ## this offset makes smoothed line match zones to exclude, but those zones must have been generated with this offset also
        ## so why aren't smoothered stored with this offset
        ## and finally, offset is too little (half?)

        self.plot_values(     self.checker.smoothed_ratios[seq_name], 
                              plot_chrom_start, 
                              plot_chrom_end, 
                              label, 
                              offset = offset,
                              panel = panel,
                              colour = ('0', '80', '70', '%'), 
                              num_ticks = 4,
                              fill = False,
                              plot_axis = False,
                              max_y_scale = axis_ratio_max,
                              add_sample_label = False)

        # plot threshold
        self.plot_threshold(self.checker.threshold, axis_ratio_max)

        # high proportion of non-proper pairs
        self.plot_suspect_regions(       seq_name,
                                         plot_chrom_start, 
                                         plot_chrom_end,
                                         self.checker.suspect_regions['rearrangements'],
                                         panel = panel,
                                         colour = ('80', '10', '0', '%'))

        # orange in inkscape
        #[float(int('e6aa00'[s:s+2], base=16))/255 for s in (0,2,4)]
        # extension near high proportion of non-proper pairs
        self.plot_suspect_regions(       seq_name,
                                         plot_chrom_start, 
                                         plot_chrom_end,
                                         self.checker.suspect_regions['rearrangements_extended'],
                                         panel = panel,
                                         colour = ('90', '60', '00', '%'))


        self.dwg.save()


def collectRegions(BAMs, genome, 
        regions = False,
        max_memory = False, 
        num_padding_positions = 5000, 
        min_align_region = 200, 
        force = False,
        # if part of a pipeline with existing logger:
        task_name = False,
        console_verbosity_lvl = False,
        log_folder = False):
    '''
    Plot all or some of the regions affected by rearrangements

    This function should be called on a per sample basis

    Parameters
    ----------
    regions : dict
        sequence names as keys and tuple integer ranges as values for regions to plot
    '''

    # this is a module level function calling other objects from module

    module_name = __name__
    # pass log_folder on to objects.__init__()
    # also console_verbosity_lvl but not sure when used <== missing .set_verbosity bug? <==
    # analysis_path needed?

    if task_name:
        # if this has been called via the CLI it is part of a task
        logger = _logging.getLogger(task_name)
        # need to add adaptor after getting (handlers etc are retained from configureLogger)
        # in this namespace logger is the adaptor but not beyond
        logger = _logging.LoggerAdapter(logger, {'task': task_name})
    else:
        # else conform to conventional Python logging and use module name
        logger = _logging.getLogger(module_name)


    out_folder = 'read_collections'
    # assert checker_info['genome_name'] == genome.id, e
    out_path = _os.path.sep.join([out_folder, genome.sample_name])
    if not _os.path.exists(out_path):
        _os.makedirs(out_path)

    # --collect always assembles
    from baga import AssembleReads
    if max_memory:
        use_mem_gigs = max_memory
    else:
        # round down available GBs
        use_mem_gigs = int(_get_available_memory())
        if use_mem_gigs == 0:
            logger.warning('Apparently <1GB memory available . . SPAdes will '\
                    'be allowed to use the minimum 1GB memory')
            use_mem_gigs = 1

    if regions:
        # target regions supplied, use those exactly
        use_num_padding_positions = 0
    else:
        # add some padding around the disrupted regions
        use_num_padding_positions = num_padding_positions

    #checkers = {}
    for BAM in BAMs:
        logger.info('Collecting short reads from {} aligned to {}'\
                ''.format(BAM, genome.sample_name))
        
        collector = Collector(BAM, genome_name = genome.sample_name, 
                    task_name = task_name, 
                    console_verbosity_lvl = console_verbosity_lvl, 
                    log_folder = log_folder)
        
        # could check sequence names of genome against those in BAM here?
        # collector.reads.references vs genome.sequences ?
        if set(collector.reads.references) != set(genome.sequence):
            raise Exception('sequence names in BAM and genome do not match')
        #set(collector.genome_lengths) == set()
        
        # always need unmapped and poorly mapped reads
        # extract from BAM and write fastq
        ##### update for all sequences:
        collector.getUnmapped(low_quality_mapping_threshold = 10)
        r1_out_path_um, r2_out_path_um, rS_out_path_um = \
                collector.writeUnmapped(out_path)
        
        ##### is this needed even for single assembly?
        # assemble poorly/unmapped alone first to subtract from other assemblies
        reads_path_unmapped = {} #?
        output_folder_um = \
                '_'.join(r1_out_path_um.split('_')[:-1]).split(_os.path.sep)[-1]
        reads_path_unmapped[output_folder_um] = (r1_out_path_um, r2_out_path_um, 
                rS_out_path_um)
        path_to_bad_unmapped_contigs = _os.path.sep.join([out_folder, 
                genome.sample_name, output_folder_um, 'contigs.fasta'])
        if _os.path.exists(path_to_bad_unmapped_contigs) and \
                _os.path.getsize(path_to_bad_unmapped_contigs) > 0 and \
                not force:
            logger.info('Found assembly at {}\nUse --force/-F to overwrite. '\
                    'Skipping . . .'.format(path_to_bad_unmapped_contigs))
        else:
            if not force:
                logger.info('Nothing found at {}. Doing assembly.'.format(
                        path_to_bad_unmapped_contigs))
            
            reads = AssembleReads.DeNovo(paths_to_reads = reads_path_unmapped)
            reads.SPAdes(output_folder = [out_folder, genome.sample_name], 
                    mem_num_gigs = use_mem_gigs, only_assembler = True, 
                    careful = False)
        
        # first decide regions to do for this BAM
        do_regions = {}
        if not regions:
            # no ranges provided so assemble each region detected when checking structure
            # assemble each region separately
            # currently range selection is rearrangements with the extensions
            # need to load checker
            single_assembly = False
            sample_name = '{}__{}'.format(collector.reads_name, 
                    genome.sample_name)
            logger.info('Loading data for: "{}"'.format(sample_name))
            checker = Checker(sample_name = sample_name, 
                    task_name = task_name, 
                    console_verbosity_lvl = console_verbosity_lvl, 
                    log_folder = log_folder, inherit_from = 'self')
            
            for seq_name,rearrangements in checker.suspect_regions['rearrangements'].items():
                rearrangements_extended = checker.suspect_regions['rearrangements_extended'][seq_name]
                # join main rearrangement zones with extended regions if found
                # to get contiguous blocks for investigation
                # extensions start on same pos so just remove ordinates listed twice
                all_regions = sorted(rearrangements + rearrangements_extended)
                use_regions = [a for b in all_regions for a in b]
                c = _Counter(use_regions)
                use_regions = [a for a in use_regions if c[a] == 1]
                use_regions = zip(use_regions[::2],use_regions[1::2])
                do_regions[seq_name] = use_regions
        else:
            # specific region(s) of specific sequence(s) already provided
            # will assemble together for all BAMs provided
            single_assembly = True
            do_regions = regions
            # single assembly of reads aligned to one or more ranges in
            # reference(s) suppied by user. (not separate assemblies of
            # multiple ranges i.e., those with putative rearrangements)
        
        # then do the collecting and assembling
        
        # assemble read from each region with poorly/unmapped or as requested
        reads_paths = {}
        # make a second dict of reads for assembly, all values for unmapped reads
        # that need to be included in each assembly
        reads_path_unmapped = {}
        assemblies_by_region = {}
        for seq_name,these_regions in sorted(do_regions.items()):
            for (s,e) in these_regions:
                collector.makeCollection(seq_name, s, e, 
                        use_num_padding_positions)
                r1_out_path, r2_out_path, rS_out_path = collector.writeCollection(out_path)
                if not r1_out_path:
                    # if no reads found, False returned
                    # do not add to reads_paths dict for assembly
                    logger.debug('no reads found by collector')
                    continue
                # put assembly in folder with same name as read files
                output_folder = \
                        '_'.join(r1_out_path.split('_')[:-1]).split(_os.path.sep)[-1]
                path_to_contigs = _os.path.sep.join([out_path, output_folder, 'contigs.fasta'])
                assemblies_by_region[s,e] = path_to_contigs
                if _os.path.exists(path_to_contigs) and \
                        _os.path.getsize(path_to_contigs) > 0 and \
                        not args.force:
                    logger.info('Found assembly at {}\nUse --force/-F to '\
                            'overwrite. Skipping . . .'.format(path_to_contigs))
                else:
                    reads_paths[output_folder] = (r1_out_path, r2_out_path, rS_out_path)
                    reads_path_unmapped[output_folder] = (r1_out_path_um, 
                            r2_out_path_um, rS_out_path_um)
            
            logger.debug('len(assemblies_by_region) == {}'.format(
                    len(assemblies_by_region)))
            reads = AssembleReads.DeNovo(paths_to_reads = reads_paths, 
                    paths_to_reads2 = reads_path_unmapped)
            reads.SPAdes(output_folder = out_path, mem_num_gigs = use_mem_gigs, 
                    single_assembly = single_assembly, only_assembler = True, 
                    careful = False)
            if False:
                ### not yet implemented for multi-chromsomes
                # a dict of paths to contigs per region
                aligner = Aligner(genome)
                unmappedfasta = _os.path.sep.join([out_path, 
                        output_folder_um, 'contigs.fasta'])
                if _os.path.exists(unmappedfasta) and \
                        _os.path.getsize(unmappedfasta) > 0:
                    if len(assemblies_by_region) > 0:
                        # provide dict of range tuples
                        aligner.alignRegions(assemblies_by_region, 
                                use_num_padding_positions, 
                                path_to_omit_sequences = unmappedfasta, 
                                single_assembly = single_assembly,
                                min_region_length = min_align_region)
                        aligner.reportAlignments()
                    else:
                        logger.warning('no assembled regions found. Either there are '\
                                'none which is fine, or SPAdes assemblies failed '\
                                'to finish. You could check SPAdes log files in '\
                                'folders in {}'.format(os.path.sep.join([
                                'read_collections', genome_name])))
                else:
                    logger.warning('no assembled unmapped and poorly mapped reads '\
                            'found at: {}'.format(unmappedfasta))
                    try:
                        r1_size = os.path.getsize(r1_out_path_um)
                        r2_size = os.path.getsize(r2_out_path_um)
                        logger.warning('but reads, {} ({:,} bytes) and {} ({:,} bytes), '\
                                'exist . . check SPAdes assembly log in {}'.format(
                                r1_out_path_um, r1_size, r2_out_path_um, r2_size, 
                                unmappedfasta.replace('contigs.fasta','')))
                    except IOError:
                        logger.warning('could not find unmapped and poorly aligned '\
                                'reads at:\n{}\n{}\nthis is unexpected but '\
                                'conceivable (if ALL reads really did map to '\
                                'reference!).'.format(r1_out_path_um, r2_out_path_um))
                    logger.info('proceeding with alignment of assembled putatively '\
                            'rearranged regions to reference nonetheless')
                    if len(assemblies_by_region) > 0:
                        aligner.alignRegions(assemblies_by_region, 
                                use_num_padding_positions, 
                                single_assembly = single_assembly)
                        aligner.reportAlignments()
                    else:
                        logger.warning('no assembled regions found. Either there '\
                                'are none which is fine, or SPAdes assemblies '\
                                'failed to finish. You could check SPAdes log files '\
                                'in folders in {}'.format(out_path))


class Collector(_MetaSample):
    '''
    Collector class of the Structure module contains methods to extract aligned 
    short reads around regions of probably rearrangements (as found by an instance 
    of the Checker class). The reads in the resulting fastq files can be used for 
    de novo reconstruction along with the unmapped reads to investigate sequences 
    in regions of interest.
    '''

    def __init__(self, path_to_bam, genome_name = False, inherit_from = False, 
            **kwargs):
        '''
        Extract paired reads mapped to a sequence described in a BAM file within 
        a single range and write to a fastq file.
        '''
        module_name = __name__
        
        if not _os.path.exists(path_to_bam):
            raise Exception('Could not find {}. Please ensure all BAM files '\
                    'exist'.format(path_to_bam))
        
        # this is not logged until super() so exceptions should be sent up to
        # CLI or whatever is calling module
        self.path_to_bam = path_to_bam
        self.reads = _pysam.Samfile(path_to_bam, 'rb')
        # baga only compatible with one reads group per bam currently
        self.reads_name = self.reads.header['RG'][0]['ID']
        if genome_name:
            self.genome_name = genome_name
            sample_name = '{}__{}'.format(self.reads_name, genome_name)
        else:
            sample_name = self.reads_name
        
        super(Collector, self).__init__(sample_name, module_name, 
                inherit_from = inherit_from, **kwargs)
        
        self.genome_lengths = {}
        for chrm_name,length in zip(self.reads.references, self.reads.lengths):
            self.genome_lengths[chrm_name] = length
        
        # needed for complementing read sequences in BAMs
        self.transtable = _string.maketrans('ACGT','TGCA')

    def getUnmapped(self, low_quality_mapping_threshold = 10):
        '''
        Using pySAM, fetch all unmapped and poorly mapped (aligned) paired end reads
        default for poorly mapped is 10% (0.1) chance of being wrong or worse:
        -10*math.log(0.1)/math.log(10) == 10

        If one read is unmapped or poorly mapped, both reads are collected.
        '''

        c = 0
        read_pairs = {}
        for seq_name in self.reads.references:
            these_pairs = {}
            self.logger.info('Searching for unmapped and poorly mapped '\
                    '(aligned) reads for {} against {}'.format(
                    self.reads.header['RG'][0]['ID'], seq_name))
            # a pos0 slice from a pos1 index should be s-1,e
            # so this is approximate as pos1 or pos0 not specified
            reads_iter = self.reads.fetch(seq_name)
            for r in reads_iter:
                c += 1
                if r.is_unmapped or r.mapping_quality <= \
                        low_quality_mapping_threshold:
                    if r.is_reverse:
                        # SEQ being reverse complemented
                        # not sure why unmapped reads would be reverse
                        # complemented from fastq sequence
                        read_info = (r.query_sequence[::-1].translate(self.transtable), 
                                r.qual[::-1])
                    else:
                        read_info = (r.query_sequence, r.qual)
                    try:
                        these_pairs[r.query_name][int(r.is_read2)+1] = \
                                (r.query_sequence, r.qual)
                    except KeyError:
                        these_pairs[r.query_name] = {}
                        these_pairs[r.query_name][int(r.is_read2)+1] = \
                                (r.query_sequence, r.qual)
            
            self.logger.info('Found {} pairs (of {:,} total) with at least one '\
                    'read unmapped against {}'.format(len(these_pairs), 
                    c, seq_name))
            read_pairs[seq_name] = these_pairs

        self.unmapped_read_pairs = read_pairs


    def writeUnmapped(self, path_to_fastq_folder):
        '''
        Write reads to fastq files returning a tuple of paths to: reads 1, reads 2 and singletons
        '''


        prefix = '{}__{}'.format(self.reads_name, self.genome_name)

        r1_fastq_filename = '{0}_unmapped_R1.fastq'.format(
                        prefix)

        r2_fastq_filename = '{0}_unmapped_R2.fastq'.format(
                        prefix)

        rS_fastq_filename = '{0}_unmapped_S.fastq'.format(
                        prefix)

        print('Writing to:\n{}/{}\n{}/{}\n{}/{}\n'.format(
                        path_to_fastq_folder,
                        r1_fastq_filename,
                        path_to_fastq_folder,
                        r2_fastq_filename,
                        path_to_fastq_folder,
                        rS_fastq_filename))

        r1_out_path = _os.path.sep.join([path_to_fastq_folder,r1_fastq_filename])
        r1_out = open(r1_out_path, 'w')

        r2_out_path = _os.path.sep.join([path_to_fastq_folder,r2_fastq_filename])
        r2_out = open(r2_out_path, 'w')

        rS_out_path = _os.path.sep.join([path_to_fastq_folder,rS_fastq_filename])
        rS_out = open(rS_out_path, 'w')

        num_written = 0
        for seq_name,read_pairs in self.unmapped_read_pairs.items():
            if len(self.unmapped_read_pairs) == 0:
                self.logger.warning('no unmapped or poorly mapped reads found '\
                        'between {} and {} in {} (this might not be a problem).'\
                        ''.format(self.reads_name, seq_name, self.genome_name))
            for read_id,reads in read_pairs.items():
                num_written += 1
                if len(reads) == 2:
                    r1seq,r1qual = reads[1]
                    r1_out.write('@{}/1\n{}\n+\n{}\n'.format(
                                read_id,
                                r1seq,
                                r1qual))
                    r2seq,r2qual = reads[2]
                    r2_out.write('@{}/1\n{}\n+\n{}\n'.format(
                                read_id,
                                r2seq,
                                r2qual))
                elif len(reads) == 1:
                    n,(rseq,rqual) = reads.items()[0]
                    rS_out.write('@{}/{}\n{}\n+\n{}\n'.format(
                                read_id,
                                n,
                                rseq,
                                rqual))

        r1_out.close()
        r2_out.close()
        rS_out.close()

        if not num_written:
            self.logger.debug('No reads written so removing empty fastq files '\
                    'at: {}, {} and {}'.format(r1_out_path, r2_out_path, 
                    rS_out_path))
            _os.unlink(r1_out_path)
            _os.unlink(r2_out_path)
            _os.unlink(rS_out_path)

        # return output destination for assembly
        return(r1_out_path, r2_out_path, rS_out_path)


    def makeCollection(self, seq_name,
                             seq_start, 
                             seq_end, 
                             num_padding_positions,
                             circular = True):
        '''
        Using pySAM, fetch all aligned reads in a specified region . . .
        '''

        ## this bit accounts for padding off one end requiring reads to be
        ## collected from the other end. Set circular = False for assembly
        ## contigs or linear chromosomes
        # a pos0 slice from a pos1 index should be s-1,e
        # so this is approximate as pos1 or pos0 not specified
        range_min = 0
        range_max = self.reads.lengths[self.reads.references.index(seq_name)]
        ranges = []
        range_start = seq_start - num_padding_positions
        range_end = seq_end + num_padding_positions

        if range_start < range_min and range_end > range_max:
            self.logger.warning('requested alignment range ({}-{}bp), which '\
                    'includes additional padding ({}bp) exceeds alignment '\
                    '(reference sequence) length: {}. Using all reads in this '\
                    'alignment.'.format(range_start, range_end, 
                    num_padding_positions, self.reads.lengths[0]))
            range_start = range_min
            range_end = range_max
        elif range_start < range_min:
            if circular:
                self.logger.debug('Assuming sequences are circlur and not '\
                        'contigs: collecting some reads from far end of a '\
                        'sequence ({})'.format(seq_name))
                ranges += [(range_max + range_start, range_max)]
            range_start = range_min
        elif range_end > range_max:
            if circular:
                self.logger.debug('Assuming sequences are circlur and not '\
                        'contigs: collecting some reads from far end of a '\
                        'sequence ({})'.format(seq_name))
                ranges += [(range_min, range_end - range_max)]
            range_end = range_max

        if len(ranges) == 1:
            self.logger.warning('Alignment range ({} to {} bp) spans start or '\
                    'end. Chromosome assumed to be circular so some reads from '\
                    'other end included in reads collection'.format(
                    seq_start - num_padding_positions, 
                    seq_end + num_padding_positions))

        ranges += [(range_start, range_end)]

        read_pairs = {}
        for start,end in ranges:
            self.logger.info('Collecting from {} to {} bp in {}'.format(start, 
                    end, seq_name))
            reads_iter = self.reads.fetch(str(seq_name), start, end)
            for r in reads_iter:
                if r.is_reverse:
                    # SEQ being reverse complemented
                    read_info = (r.query_sequence[::-1].translate(self.transtable), 
                            r.qual[::-1])
                else:
                    read_info = (r.query_sequence, r.qual)
                
                try:
                    read_pairs[r.query_name][int(r.is_read2)+1] = read_info
                except KeyError:
                    read_pairs[r.query_name] = {}
                    read_pairs[r.query_name][int(r.is_read2)+1] = read_info

        self.collected_read_pairs = read_pairs
        self.collected_seq_name = seq_name
        self.collected_start = seq_start
        self.collected_end = seq_end
        self.collected_num_padding_positions = num_padding_positions

        self.logger.info('Found {} pairs with at least one read mapped to this '\
                'region'.format(len(read_pairs)))

    def writeCollection(self, path_to_fastq_folder):
        '''
        Write reads to a fastq file
        '''

        if len(self.collected_read_pairs) == 0:
            self.logger.warning('no reads found at {}-{} bp between {} and {}. '\
                    'Not writing fastq files.'.format(self.collected_start,
                    self.collected_end, self.reads_name, 
                    self.collected_seq_name))
            return(None, None, None)
        else:
            prefix = '{}__{}__{}'.format(self.reads_name, self.genome_name, 
                    self.collected_seq_name)
            seq_name_i = self.reads.references.index(self.collected_seq_name)
            zeropadding = len(str(self.reads.lengths[seq_name_i]))
            
            r1_fastq_filename = '{0}_{1:0{3}d}-{2:0{3}d}+{4}_R1.fastq'.format(
                            prefix,
                            self.collected_start,
                            self.collected_end,
                            zeropadding,
                            self.collected_num_padding_positions)
            
            r2_fastq_filename = '{0}_{1:0{3}d}-{2:0{3}d}+{4}_R2.fastq'.format(
                            prefix,
                            self.collected_start,
                            self.collected_end,
                            zeropadding,
                            self.collected_num_padding_positions)
            
            rS_fastq_filename = '{0}_{1:0{3}d}-{2:0{3}d}+{4}_S.fastq'.format(
                            prefix,
                            self.collected_start,
                            self.collected_end,
                            zeropadding,
                            self.collected_num_padding_positions)
            
            self.logger.info('Writing to pair member 1 to: {}{}{}'.format(
                    path_to_fastq_folder, _os.path.sep, r1_fastq_filename))
            self.logger.info('Writing to pair member 2 to: {}{}{}'.format(
                    path_to_fastq_folder, _os.path.sep, r2_fastq_filename))
            self.logger.info('Writing to unpaired to: {}{}{}'.format(
                    path_to_fastq_folder, _os.path.sep, rS_fastq_filename))
            
            r1_out_path = _os.path.sep.join([path_to_fastq_folder,r1_fastq_filename])
            r1_out = open(r1_out_path, 'w')
            
            r2_out_path = _os.path.sep.join([path_to_fastq_folder,r2_fastq_filename])
            r2_out = open(r2_out_path, 'w')
            
            rS_out_path = _os.path.sep.join([path_to_fastq_folder,rS_fastq_filename])
            rS_out = open(rS_out_path, 'w')
            
            
            for read_id,reads in self.collected_read_pairs.items():
                if len(reads) == 2:
                    r1seq,r1qual = reads[1]
                    r1_out.write('@{}/1\n{}\n+\n{}\n'.format(
                                read_id,
                                r1seq,
                                r1qual))
                    r2seq,r2qual = reads[2]
                    r2_out.write('@{}/2\n{}\n+\n{}\n'.format(
                                read_id,
                                r2seq,
                                r2qual))
                elif len(reads) == 1:
                    n,(rseq,rqual) = reads.items()[0]
                    rS_out.write('@{}/{}\n{}\n+\n{}\n'.format(
                                read_id,
                                n,
                                rseq,
                                rqual))
            
            r1_out.close()
            r2_out.close()
            rS_out.close()
            
            # return output destination for assembly
            return(r1_out_path, r2_out_path, rS_out_path)


class Aligner:
    '''
    Aligner class of the Structure module contains methods to align contigs, 
    assembled de novo from regions likely to be affected by rearrangements (as 
    found by an instance of the Checker class), back to those regions in the reference chromosome.
    '''
    def __init__(self, genome):
        '''
        A region of interest aligner be instantiated with a CollectData.Genome 
        instance to align to.
        '''
        
        self.genome = genome
        self.ORFs_ordered = {}
        for seq_name,annotations in genome.annotations.items():
            self.ORFs_ordered[seq_name] = sorted(annotations[0], 
                    key = annotations[0].get)
        self.exe_aligner = _get_exe_path('seq-align')
    def seqalign(self, seqA, seqB, algorithm = 'needleman_wunsch', protein = 'no'):
        
        '''
        Call seq-align as a subprocess to do a gapped pairwise optimal alignment
        algorithm = 'needleman_wunsch'|'smith_waterman'
        '''

        # update aligner executable
        # 'smith_waterman' or 'needleman_wunsch'
        exe = _os.path.sep.join(self.exe_aligner.split(_os.path.sep)[:-1]+[algorithm])

        if protein in (True, 'yes'):
            seqA = seqA.translate()
            seqB = seqB.translate()

        # seqA = _Seq('cgcacgatttgtagtcccgtgtctgcatggacatagcca')
        # seqB = _Seq('aacacaaacaacacaccgggcctagtagagagttggcggcgcc')
        # seqA = _Seq('AMKINVFYAEEEESRCSRPFIISPSLVVGLREQRNLKLDSKKASLV')
        # seqB = _Seq('AMKIIKIKNVFYAESRCSRPFIISPSLVVGVQRREQRNLKLDSKKASLV')
        if algorithm == 'needleman_wunsch':
            out_handle = _StringIO()
            seqrecords = [_SeqRecord(seqA, id = 'A'), _SeqRecord(seqB, id = 'B')]
            _SeqIO.write(seqrecords, out_handle, 'fasta')
            cmd = [exe]
            cmd += ['--stdin']
            cmd += ['--freestartgap']
            cmd += ['--freeendgap']
            # print(' '.join(cmd))
            proc = _subprocess.Popen(
                cmd,
                stdout = _subprocess.PIPE,
                stderr = _subprocess.PIPE,
                stdin = _subprocess.PIPE)
            
            result, err = proc.communicate(out_handle.getvalue())
            try:
                A, B = result.split('\n')[:2]
                return(A, B)
            except ValueError as e:
                print(seqrecords)
                _SeqIO.write(seqrecords, 'temp.fasta', 'fasta')
                _sys.exit('There seems to be a problem with the seqalign output (return code {}):\n{} ({})'.format(proc.returncode,err,e))

        else:
            print('WARNING: smith-waterman wrapper is unimplemented: output much more complex and is simply printed to stdout for now')
            out_handle = _StringIO()
            #seqrecords = [_SeqRecord(seqA, id = 'A'), _SeqRecord(seqB, id = 'B')]
            #_SeqIO.write(seqrecords, out_handle, 'fasta')
            # smith_waterman seems to required raw sequences, not fastas
            seqrecords = '\n'.join([str(seqA), str(seqB)])
            out_handle.write(seqrecords)
            cmd = [exe]
            cmd += ['--stdin']
            # because expected to deteriorate?
            # cmd += ['--freeendgap']
            # print(' '.join(cmd))
            proc = _subprocess.Popen(
                cmd,
                stdout = _subprocess.PIPE,
                stderr = _subprocess.PIPE,
                stdin = _subprocess.PIPE)
            
            # proc.stdin.write(out_handle.getvalue())
            # proc.stdin.close()
            #result = proc.stdout.read()
            result, err = proc.communicate(out_handle.getvalue())
            #proc.wait()
            print(result)
            # try:
                # A, B = result.split('\n')[:2]
                # print(A,B)
                # return(A, B)
            # except ValueError as e:
                # _sys.exit('There seems to be a problem with the seqalign output:\n{} ({})'.format(err,e))
            



    def countEndGaps(self, seq, pre = True):
        if pre:
            use = seq
        else:
            use = reversed(seq)
        num_gaps = 0
        for c in use:
            if c == '-':
                num_gaps += 1
            else:
                break

        return(num_gaps)

    def alignRegions(self, assemblies_by_region, 
                           num_padding_positions = 5000, 
                           pID_window = 100, 
                           pID_step = 10, 
                           min_region_length = 200,
                           path_to_omit_sequences = False,
                           min_pID_aligned = 0.9,
                           single_assembly = False,
                           force = False):
        
        '''
        Align contigs to reference chromosome regions

        Provided with a dict of chromosome range tuples to contig file paths, 
        align each range to the contigs using seq-align as a subprocess to do 
        a Needleman-Wunch gapped pairwise globl alignment.

        path_to_omit_sequences can be a fasta file of contigs that should be ignored if
        encountered (e.g. contigs from unmapped reads to subtract from assemblies
        that include reads from regions of interest).

        min_pID_aligned the minimum percent identity required at one or more windows 
        along the alignment.
        '''

        if path_to_omit_sequences:
            unmapped_read_contigs = set([str(rec.seq) for rec in _SeqIO.parse(path_to_omit_sequences, 'fasta')])
        else:
            unmapped_read_contigs = set()

        assert len(assemblies_by_region) > 0, 'no assemblies by region were provided: '\
                'cannot align them'

        aligned = {}
        too_small = []
        for (s,e),contigfile in sorted(assemblies_by_region.items()):
            if e - s < min_region_length:
                too_small += [(contigfile,e - s)]
            else:
                ref_chrom_region = _Seq(self.genome.sequence[s-num_padding_positions:e+num_padding_positions].tostring())
                ref_region_id = 'ref_{:07d}_{:07d}'.format(s-num_padding_positions,e+num_padding_positions)
                aligned[ref_region_id] = {}
                # get appropriate filename
                if single_assembly:
                    # contigfile == .../<sample>__<genome>_<start>-<end>+<padding>/contigs.fasta
                    use_contigfile = '_'.join(contigfile.split('_')[:-1] + ['multi_region']) + _os.path.sep + 'contigs.fasta'
                    # use_contigfile == .../<sample>__<genome>_multi_region/contigs.fasta
                else:
                    use_contigfile = contigfile
                
                # check if already did alignment
                pattern = _os.path.sep.join(use_contigfile.split(_os.path.sep)[:-1]) + '_{}_vs_contig*.fna'.format(ref_region_id)
                previous_alignments = _glob(pattern)
                if len(previous_alignments) > 0 and not force:
                    print('Found previous alignments for region {}:\n{}'.format(ref_region_id,'\n'.join(previous_alignments)))
                    print('Use --force/-F to realign and overwrite')
                    for aln in previous_alignments:
                        print('contig'+aln[:-4].split('contig')[-1])
                        aligned[ref_region_id]['contig'+aln[:-4].split('contig')[-1]] = tuple([a.seq for a in _SeqIO.parse(aln, 'fasta')])
                else:
                    try:
                        num_contigs = 0
                        for n,rec in enumerate(_SeqIO.parse(use_contigfile,'fasta')):
                            num_contigs += 1
                    except IOError:
                        print('WARNING: cannot access {}'.format(use_contigfile))
                        print('This assembly may have failed . . .')
                        continue
                    
                    print('Found {} contigs in {} for analysis.'.format(num_contigs, use_contigfile))
                    for n,rec in enumerate(_SeqIO.parse(use_contigfile,'fasta')):
                        if str(rec.seq) not in unmapped_read_contigs:
                            print('Aligning novel contig: {} to chromosome region {}-{} bp'.format(rec.id, s, e))
                            use_seq = rec
                            #ref_alnd,contig_alnd = self.seqalign(ref_chrom_region, rec.seq, 'smith_waterman')
                            ref_alnd,contig_alnd = self.seqalign(ref_chrom_region, rec.seq, 'needleman_wunsch')
                            #print(ref_alnd,contig_alnd)
                            pIDs = dict(self.get_percent_ID(ref_alnd, contig_alnd, window = pID_window, step = pID_step))
                            #ref_alnd_rc,contig_alnd_rc = self.seqalign(ref_chrom_region, rec.seq.reverse_complement(), 'smith_waterman')
                            ref_alnd_rc,contig_alnd_rc = self.seqalign(ref_chrom_region, rec.seq.reverse_complement(), 'needleman_wunsch')
                            #print(ref_alnd,contig_alnd)
                            pIDs_rc = dict(self.get_percent_ID(ref_alnd_rc, contig_alnd_rc, window = pID_window, step = pID_step))
                            retained = False
                            if sum(pIDs.values()) > sum(pIDs_rc.values()):
                                if max(pIDs.values()) >= min_pID_aligned:
                                    # only print if contains a pID_window bp window with > 90% identity
                                    fout = _os.path.sep.join(use_contigfile.split(_os.path.sep)[:-1]) + '_{}_vs_contig{:02d}.fna'.format(ref_region_id,n+1)
                                    seqs = []
                                    seqs += [_SeqRecord(seq = _Seq(ref_alnd), 
                                                        id = ref_region_id)]
                                    seqs += [_SeqRecord(seq = _Seq(contig_alnd), 
                                                        id = rec.id)]
                                    _SeqIO.write(seqs, fout, 'fasta')
                                    print('Writing: {}'.format(fout))
                                    retained = True
                                    aligned[ref_region_id]['contig{:02d}'.format(n+1)] = (ref_alnd, contig_alnd)
                            else:
                                try:
                                    if max(pIDs_rc.values()) >= min_pID_aligned:
                                        # only print if contains a pID_window bp window with > 90% identity
                                        fout = _os.path.sep.join(
                                                use_contigfile.split(_os.path.sep)[:-1]) + \
                                                '_{}_vs_contig{:02d}_rc.fna'.format(ref_region_id,n+1)
                                        seqs = []
                                        seqs += [_SeqRecord(seq = _Seq(ref_alnd_rc), 
                                                            id = ref_region_id)]
                                        seqs += [_SeqRecord(seq = _Seq(contig_alnd_rc), 
                                                            id = rec.id)]
                                        _SeqIO.write(seqs, fout, 'fasta')
                                        print('Writing: {}'.format(fout))
                                        retained = True
                                        # only thr contig is reverse complemented here
                                        aligned[ref_region_id]['contig{:02d}_rc'.format(n+1)] = (ref_alnd_rc, contig_alnd_rc)
                                except ValueError:
                                    # no regions found
                                    pass
                            if not retained:
                                print('No alignment with a percent identity >= {:.0%} over a window of {} bp'.format(min_pID_aligned, pID_window))
                        else:
                            print('Omitting contig from unmapped/poorly mapped reads: {}'.format(rec.id))
            
                print(ref_region_id,len(aligned[ref_region_id]))
                if len(aligned[ref_region_id]) == 0:
                    del aligned[ref_region_id] 

        # retain info useful for summarising
        self.aligned = aligned
        self.aligned_from_single_assembly = single_assembly
        if len(too_small):
            print('These regions were shorter than "min_region_length" ({} bp) and not aligned:'.format(min_region_length))
            for r,l in too_small:
                print('{}: {} bp'.format(r,l))
        if len(too_small) == len(assemblies_by_region):
            print('All {} regions too short . . . leaving nothing to align.'.format(len(assemblies_by_region)))
            self.aligned_to_sample = False
            self.aligned_to_genome = False
            self.path_to_alignments = False
        else:
            # use_contigfile == .../<sample>__<genome>_multi_region/contigs.fasta
            s, x = use_contigfile.split(_os.path.sep)[-2].split('__')
            g = x.split('_')[0]
            self.aligned_to_sample = s
            self.aligned_to_genome = g
            self.path_to_alignments = _os.path.sep.join(use_contigfile.split(_os.path.sep)[:-2])

    def reportAlignments(self):
        '''
        Summarise variants in alignments of contigs to reference chromosome regions

        Only SNPs are supported: indels are reported for each inserted or deleted 
        positions treating gaps as an additional character

        Generate a pseudo multiple alignment including each aligned contig and the 
        referance region. Write out variants to a table as a comma separated values 
        text file.
        '''

        assert hasattr(self, 'aligned'), '"aligned" attribute not found. .alignRegions() not run?'

        if len(self.aligned) == 0:
            print('No alignments to report')

        print(self.aligned_from_single_assembly, self.aligned_to_sample, self.aligned_to_genome, self.path_to_alignments)

        variants_by_contig = {}
        alnd_range_by_contig = {}
        aligned_regions = {}
        for ref_region_id, alnd_contigs in self.aligned.items():
            these_aligned_regions = []
            print(ref_region_id, len(alnd_contigs))
            ref_region_start = int(ref_region_id.split('_')[1])
            ref_region_end = int(ref_region_id.split('_')[2])
            # update if near end of chromosome for calcs below
            # (perhaps this should happen earlier?)
            ref_region_end = min(len(self.genome.sequence),ref_region_end)
            ref_region_len = ref_region_end - ref_region_start
            variants_by_contig[ref_region_id] = {}
            alnd_range_by_contig[ref_region_id] = {}
            for contig_id, (ref_alnd_seq, contig_alnd_seq) in alnd_contigs.items():
                print(contig_id)
                ref_pre_gaps = self.countEndGaps(ref_alnd_seq, pre = True)
                ref_post_gaps = self.countEndGaps(ref_alnd_seq, pre = False)
                contig_pre_gaps = self.countEndGaps(contig_alnd_seq, pre = True)
                contig_post_gaps = self.countEndGaps(contig_alnd_seq, pre = False)
                pre_gaps = max(ref_pre_gaps,contig_pre_gaps)
                post_gaps = max(ref_post_gaps,contig_post_gaps)
                print("ref_pre_gaps: {}, ref_post_gaps: {}, contig_pre_gaps: {}, "\
                        "contig_post_gaps: {}, pre_gaps: {} ,post_gaps: {}".format(
                        ref_pre_gaps, ref_post_gaps, contig_pre_gaps, contig_post_gaps, 
                        pre_gaps, post_gaps))
                alnd_len = len(ref_alnd_seq)
                query_insertions = '-' in ref_alnd_seq[pre_gaps:alnd_len-post_gaps]
                A_alnd = str(ref_alnd_seq[ref_pre_gaps:alnd_len-ref_post_gaps])
                B_orig = self.genome.sequence[ref_region_start:ref_region_end].tostring()
                # confirm, minus gaps in aligned contig, we have the correct region
                assert A_alnd.replace('-','') == B_orig, '{}\n!=\n{}'\
                        ''.format(A_alnd.replace('-',''),B_orig)
                
                # identify variants
                these_variants = {}
                refpos0 = ref_region_start
                alnpos0 = 0
                for ref in str(ref_alnd_seq[pre_gaps:alnd_len-post_gaps]):
                    if ref != contig_alnd_seq[pre_gaps:][alnpos0]:
                        # save SNPs and indels, but for indels, each contiguous gap is saved separately
                        #print(alnpos0, refpos0, ref, contig_alnd_seq[pre_gaps:][alnpos0], self.genome.sequence[refpos0])
                        these_variants[refpos0] = (ref, contig_alnd_seq[pre_gaps:][alnpos0])
                    
                    alnpos0 += 1
                    if ref != '-':
                        refpos0 += 1
                
                variants_by_contig[ref_region_id][contig_id] = these_variants
                alnd_range_by_contig[ref_region_id][contig_id] = (
                        ref_region_start + contig_pre_gaps, 
                        ref_region_end - contig_post_gaps)
                
                ref_sequence_ends_trimmed = []
                contig_sequence_ends_trimmed = []
                del_warning = 0
                for ref,cont in zip(ref_alnd_seq[pre_gaps:len(ref_alnd_seq)-post_gaps], 
                        contig_alnd_seq[pre_gaps:len(ref_alnd_seq)-post_gaps]):
                    if '-' == ref:
                        del_warning += 1
                    else:
                        ref_sequence_ends_trimmed += [ref]
                        contig_sequence_ends_trimmed += [cont]
                if del_warning:
                    print('WARNING: deletions in reference/insertions in sample are not '\
                            'supported: {} gap(s) in reference chromsome region was omitted'.format(del_warning))
                effective_region_end = min(ref_region_end, ref_region_end - contig_post_gaps)
                effective_region_start = max(ref_region_start, ref_region_start + contig_pre_gaps)
                e = '{} vs {} - {} == {}'.format(len(ref_sequence_ends_trimmed), 
                        effective_region_end, effective_region_start, 
                        effective_region_end - effective_region_start)
                assert len(ref_sequence_ends_trimmed) == effective_region_end - effective_region_start, e
                these_aligned_regions += [(
                        _SeqRecord(_Seq(''.join(ref_sequence_ends_trimmed)), id = ref_region_id), 
                        _SeqRecord(_Seq(''.join(contig_sequence_ends_trimmed)), id = contig_id))]
            
            aligned_regions[ref_region_id] = these_aligned_regions


        if self.aligned_from_single_assembly:
            alnfilename_variants = '{}/{}__{}_multi_region_alnd_variants.csv'.format(
                    self.path_to_alignments, self.aligned_to_sample, 
                    self.aligned_to_genome)
            alnfilename_summary = '{}/{}__{}_multi_region_alnd_summary.csv'.format(
                    self.path_to_alignments, self.aligned_to_sample, 
                    self.aligned_to_genome)
            with open(alnfilename_variants, 'w') as foutvar, \
                    open(alnfilename_summary, 'w') as foutsumm:
                foutvar.write('"assembly start","assembly end","contig","position",'\
                        '"reference","variant"\n')
                foutsumm.write('"assembly start","assembly end","contig name",'\
                        '"contig start","contig end","aligned length (bp)",'\
                        '"total variants"\n')
                for ref_region_id, alnd_contigs in alnd_range_by_contig.items():
                    # ref_0000100_0000200
                    ref_start, ref_end = map(int,ref_region_id[4:].split('_'))
                    alnfilename_msa = '{}/{}__{}_multi_region_{}.fna'.format(
                            self.path_to_alignments, self.aligned_to_sample, 
                            self.aligned_to_genome, ref_region_id[4:])
                    # write reference once then each aligned contig sequence
                    with open(alnfilename_msa, 'w') as foutalnd:
                        _SeqIO.write(aligned_regions[ref_region_id][0][0], foutalnd, 'fasta')
                        for ref_aligned,contig_aligned in aligned_regions[ref_region_id]:
                            _SeqIO.write(contig_aligned, foutalnd, 'fasta')
                    for contig_id, (contig_start, contig_end) in alnd_contigs.items():
                        print('{} to {} is {} to {}'.format(ref_region_id, contig_id, contig_start+1, contig_end))
                        these_variants = variants_by_contig[ref_region_id][contig_id]
                        foutsumm.write('{}, {}, "{}", {}, {}, {}\n'.format(
                                ref_start+1, ref_end, contig_id, contig_start+1, 
                                contig_end, contig_end - contig_start, len(these_variants)))
                        for refpos0,(ref_char,alnd_char) in these_variants.items():
                            # assembly start, assembly end, contig, refpos1, ref_char, aligned_char
                            # base-1 counts, position inclusive ranges i.e., not Python slices
                            foutvar.write('{},{},"{}",{},"{}","{}"\n'.format(
                                    ref_start+1, ref_end, contig_id, refpos0+1, ref_char, alnd_char))
        else:
            for ref_region_id, alnd_contigs in alnd_range_by_contig.items():
                if len(aligned_regions[ref_region_id]) == 0:
                    print('No alignment to summarise for {}'.format(ref_region_id))
                    continue
                # ref_0000100_0000200
                ref_start, ref_end = map(int,ref_region_id[4:].split('_'))
                # write reference once then each aligned contig sequence
                alnfilename_msa = '{}/{}__{}_{}_multi_alnd.fna'.format(
                        self.path_to_alignments, self.aligned_to_sample, 
                        self.aligned_to_genome, ref_region_id[4:])
                with open(alnfilename_msa, 'w') as foutalnd:
                    _SeqIO.write(aligned_regions[ref_region_id][0][0], foutalnd, 'fasta')
                    for ref_aligned,contig_aligned in aligned_regions[ref_region_id]:
                        _SeqIO.write(contig_aligned, foutalnd, 'fasta')
                        #_SeqIO.write(aligned_regions[ref_region_id][0][1], foutalnd, 'fasta')
                alnfilename_variants = '{}/{}__{}_{}_alnd_variants.csv'.format(
                        self.path_to_alignments, self.aligned_to_sample, 
                        self.aligned_to_genome, ref_region_id[4:])
                alnfilename_summary = '{}/{}__{}_{}_alnd_summary.csv'.format(
                        self.path_to_alignments, self.aligned_to_sample, 
                        self.aligned_to_genome, ref_region_id[4:])
                print('Summarising alignment for {} in {}'.format(ref_region_id, alnfilename_summary))
                with open(alnfilename_variants, 'w') as foutvar, open(alnfilename_summary, 'w') as foutsumm:
                    foutvar.write('"assembly start","assembly end","contig","position",'\
                            '"reference","variant"\n')
                    foutsumm.write('"assembly start","assembly end","contig name",'\
                            '"contig start","contig end","aligned length (bp)",'\
                            '"total variants"\n')
                    for contig_id, (contig_start, contig_end) in alnd_contigs.items():
                        print('{} to {} is {} to {}'.format(ref_region_id, contig_id, contig_start+1, contig_end))
                        these_variants = variants_by_contig[ref_region_id][contig_id]
                        foutsumm.write('{}, {}, "{}", {}, {}, {}\n'.format(
                                ref_start+1, ref_end, contig_id, contig_start+1, 
                                contig_end, contig_end - contig_start, len(these_variants)))
                        for refpos0,(ref_char,alnd_char) in these_variants.items():
                            # assembly start, assembly end, contig, refpos1, ref_char, aligned_char
                            # base-1 counts, position inclusive ranges i.e., not Python slices
                            foutvar.write('{},{},"{}",{},"{}","{}"\n'.format(
                                    ref_start+1, ref_end, contig_id, refpos0+1, ref_char, alnd_char))

        info = {'variants_by_contig':variants_by_contig, 
                'alnd_range_by_contig':alnd_range_by_contig, 
                'aligned_regions':aligned_regions}


        return(info)

    def get_percent_ID(self, A, B, window = 100, step = 20):
        pID_per_window = []
        for i in range(0, len(A)-window, step):
            Achunk = A[i:i+window]
            Bchunk = B[i:i+window]
            pID = sum([a == b for a,b in zip(Achunk,Bchunk)])/float(window)
            pID_per_window += [(i+window/2,pID)]

        return(pID_per_window)

if __name__ == '__main__':
    main()
