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
AssembleReads module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module a wrapper to a de novo short read assembler.

The wrapped short read assembler is SPAdes 3.6.1.
http://spades.bioinf.spbau.ru/
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
class DeNovo:
    '''
    Methods to wrap third party genome read assemblers: currently only SPAdes 3.6.1.
    '''
    def __init__(self, reads):
        '''
        Initialise with a baga.PrepareReads.Reads object
        '''
        try:
            self.read_files = reads.trimmed_read_files
        except AttributeError:
            text = 'WARNING: baga was not used to quality-score trim these reads. Read trimming is recommended for most types of analysis. This can be achieved with the "trim()" method of the Reads class in the PrepareReads module.'
            print(text)
            try:
                self.read_files = reads.adaptorcut_read_files
            except AttributeError:
                text = 'WARNING: baga was not used to remove library preparation adaptor sequences from these reads. Adaptor removal is highly recommended so hopefully you already removed adaptor sequences! This can be achieved with the "cutAdaptors()" method of the Reads class in the PrepareReads module.'
                try:
                    self.read_files = reads.read_files
                    print(text)
                    print('continuing with these reads . . .')
                except AttributeError:
                    print('No reads info . . . baga.Assemble.DeNovo needs a baga.PrepareReads.Reads object with a "read_files" attribute.')    
        # currently baga CollectData includes path to reads in pairname keys to read file pair values
        # check and remove here
        for pairname, files in self.read_files.items():
            if _os.path.sep in pairname:
                self.read_files[pairname.split(_os.path.sep)[-1]] = files
                del self.read_files[pairname]


    def SPAdes(self, 
            exe = [], 
            output_folder = ['assemblies','SPAdes'],
            mem_num_gigs = 8, 
            max_cpus = -1):
        '''
    http://spades.bioinf.spbau.ru/release3.6.1/manual.html
    relevent inputs:
    -o <output_dir> Specify the output directory. Required option.
    --sc required for MDA (single-cell) data.
    --only-error-correction
    --only-assembler
    --careful reduce the number of mismatches and short indels. Run MismatchCorrector – a post processing tool. Recommended.
    --continue from the specified output folder starting from the last available check-point
    --restart-from <check_point>
        ec start from error correction
        as restart assembly module from the first iteration
        k<int> restart from the iteration with specified k values, e.g. k55
        mc restart mismatch correction
    --pe1-12 <file_name> interlaced forward and reverse paired-end reads.
    --pe1-1 <file_name> File with forward reads.
    --pe1-2 <file_name> File with reverse reads.
    --pe1-s <file_name> File with unpaired reads . . use --pe2-... for next library
    --threads <int>
    --memory <int> max memory in Gb
    -k <int,int,...>  Comma-separated list of odd ascending k-mers
    If --sc is set the default value are 21,33,55, for multicell data sets it is auto
    --cov-cutoff <float> positive float value, or 'auto', or 'off'. Default value is 'off'
        '''

        output_path = _os.path.sep.join(output_folder)
        if not _os.path.exists(output_path):
            _os.makedirs(output_path)

        # max threads is slightly different to cpus
        # . . can probably use more
        max_processes = _decide_max_processes( max_cpus )

        start_time = _time.time()

        if len(exe) == 0:
            from baga import Dependencies
            try:
                # for e.g. python2 when python3.3 is last python3 supported
                use_exe = Dependencies.dependencies['spades']['checker']['arguments']['extra_pre'] + \
                                [_os.path.sep.join(
                                    [Dependencies.destination_programs] + \
                                    Dependencies.dependencies['spades']['checker']['arguments']['path']
                                )]
            except KeyError:
                use_exe = [_os.path.sep.join(
                                    [Dependencies.destination_programs] + \
                                    Dependencies.dependencies['spades']['checker']['arguments']['path']
                                )]
        else:
            use_exe = _os.path.sep.join(exe)

        contigs = []
        for pairname, files in self.read_files.items():
            # check for singletons?
            print(files)
            cmd = use_exe
            cmd += ['--pe1-1', files[1]]
            cmd += ['--pe1-2', files[2]]
            cmd += ['-o', output_path]
            cmd += ['--threads', str(max_processes)]
            cmd += ['--memory', str(mem_num_gigs)]
            cmd += ['--careful']
            print(' '.join(cmd))
            _subprocess.call(cmd)

        self.paths_to_contigs = contigs


if __name__ == '__main__':
    main()
