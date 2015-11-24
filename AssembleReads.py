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
# Dr Michael A Brockhurst (The University of York, UK)
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
    def __init__(self, baga = False, paths_to_reads = False, paths_to_reads2 = False):
        '''
        Initialise with a baga.PrepareReads.Reads object or, 
        dict of paths to read files. For each element either dict 
        (keys: 1,2,'S') or tuple:
            Length two: read 1 and read 2 of pairs
            Length three: read 1, read 2 of pairs and singletons
        Keys are path to folder for output.
        Optionally, provide a list of secondary reads, same length as main reads list.
        '''
        assert sum([baga is False, paths_to_reads is False]), 'initiatiate with a baga.PrepareReads.Reads object or a list of paths to paired reads'
        if baga:
            try:
                self.read_files = baga.trimmed_read_files
            except AttributeError:
                text = 'WARNING: baga was not used to quality-score trim these reads. Read trimming is recommended for most types of analysis. This can be achieved with the "trim()" method of the Reads class in the PrepareReads module.'
                print(text)
                try:
                    self.read_files = baga.adaptorcut_read_files
                except AttributeError:
                    text = 'WARNING: baga was not used to remove library preparation adaptor sequences from these reads. Adaptor removal is highly recommended so hopefully you already removed adaptor sequences! This can be achieved with the "cutAdaptors()" method of the Reads class in the PrepareReads module.'
                    try:
                        self.read_files = baga.read_files
                        print(text)
                        print('continuing with these reads . . .')
                    except AttributeError:
                        print('No reads info . . . baga.Assemble.DeNovo needs a baga.PrepareReads.Reads object with a "read_files" attribute.')    
            # currently baga CollectData includes path to reads in pairname keys to read file pair values
            # check and remove here, use just filename
            for pairname, files in self.read_files.items():
                if _os.path.sep in pairname:
                    self.read_files[pairname.split(_os.path.sep)[-1]] = files
                    del self.read_files[pairname]
        else:
            self.read_files = paths_to_reads
            if paths_to_reads2:
                assert len(paths_to_reads) == len(paths_to_reads2), 'List of secondary read libraries must be same length as first'
                self.read_files2 = paths_to_reads2


            


    def SPAdes(self, 
            exe = [], 
            output_folder = ['assemblies','SPAdes'],
            mem_num_gigs = 8, 
            max_cpus = -1,
            single_assembly = False):
        '''
        de novo assembly of short reads using SPAdes

        By default, the provided short reads in dictionary: self.paths_to_reads
        will be assembled separately, unless single_assembly set to True in 
        which case each set of paired read fastq files will be used in a 
        single assembly.

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

        assert isinstance(output_folder, list), 'Provide output folder as list of folders forming path'

        base_output_path = _os.path.sep.join(output_folder)

        if not _os.path.exists(base_output_path):
            _os.makedirs(base_output_path)

        # max threads is slightly different to cpus
        # . . can probably use more
        max_processes = _decide_max_processes( max_cpus )

        # if an exe is not provided, use that stored in Dependencies
        if len(exe):
            use_exe = _os.path.sep.join(exe)
        else:
            from baga import Dependencies
            use_exe = _get_exe_path('spades')

        def run_SPAdes(cmd):
            proc = _subprocess.Popen(cmd, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)
            # allow for failed SPAdes runs (possibly caused by small fastq files) <== but also check they were actually built properly
            try:
                stdout_value, stderr_value = proc.communicate()
                checkthese = []
                getline = False
                for line in stdout_value.split('\n'):
                    if 'Warnings saved to' in line:
                        getline = False
                    if getline:
                        l = line.rstrip()
                        if len(l):
                            checkthese += [l]
                    if 'SPAdes pipeline finished WITH WARNINGS!' in line:
                        getline = True
                
                if len(checkthese):
                    print('SPAdes completed with warnings:\n{}\n'.format('\n'.join(checkthese)))
                else:
                    print('SPAdes completed without warnings')
                
                # with open('___SPAdes_{}_good_{}.log'.format(cnum, thetime), 'w') as fout:
                    # fout.write(stdout_value)
                path2contigs = _os.path.sep.join([this_output_path,'contigs.fasta'])
            except _subprocess.CalledProcessError as e:
                print('SPAdes probably did not complete: error returned ({})'.format(proc.returncode))
                print('Error: {}'.format(e))
                print('Writing some info relevent to SPAdes crash to ___SPAdes_{}_bad_{}.log'.format(cnum, thetime))
                with open('___SPAdes_{}_bad_{}.log'.format(cnum, thetime), 'w') as fout:
                    fout.write(dir(proc))
                    fout.write('\n' + str(e.returncode) + '\n')
                    fout.write(_os.path.sep.join([this_output_path,'contigs.fasta']))
                
                path2contigs = None
            
            return(path2contigs)

        if isinstance(use_exe, list):
            # allow for use of prepended executable with script to run
            cmd = list(use_exe)
        else:
            # or just executable
            cmd = [use_exe]

        contigs = {}
        if single_assembly:
            print('Combining reads aligned at multiple regions into single assembly')
            if isinstance(use_exe, list):
                # allow for use of prepended executable with script to run
                cmd = list(use_exe)
            else:
                # or just executable
                cmd = [use_exe]
            for cnum, (pairname, files) in enumerate(self.read_files.items()):
                # allow use of tuples or dicts by converting dicts to lists
                if isinstance(files, dict):
                    use_files = []
                    for k,v in sorted(files.items()):
                        use_files += [v]
                else:
                    use_files = files
                
                cmd += ['--pe{}-1'.format(cnum+1), use_files[0]]
                cmd += ['--pe{}-2'.format(cnum+1), use_files[1]]
                try:
                    # use unpaired reads if available
                    cmd += ['--pe{}-s'.format(cnum+1), use_files[2]]
                except IndexError:
                    pass
            try:
                # add a second library if provided
                if isinstance(self.read_files2[pairname], dict):
                    # if a dict supplied, make it a list
                    use_files2 = []
                    for k,v in sorted(self.read_files2[pairname].items()):
                        use_files2 += [v]
                else:
                    use_files2 = self.read_files2[pairname]
                cmd += ['--pe{}-1'.format(cnum+2), use_files2[0]]
                cmd += ['--pe{}-2'.format(cnum+2), use_files2[1]]
                try:
                    cmd += ['--pe{}-s'.format(cnum+2), use_files2[2]]
                except IndexError:
                    pass
            except AttributeError:
                pass
            
            ## this isn't very flexible:
            # retain <sample>__<genome> from pairname:
            # pairname == <sample>__<genome>_<start>-<end>+<padding>
            # and replace with multiregion
            folder = '{}__{}_{}'.format(pairname.split('__')[0],
                                        pairname.split('__')[1].split('_')[0],
                                        'multi_region')
            this_output_path = _os.path.sep.join(output_folder + [folder])
            if not _os.path.exists(this_output_path):
                _os.makedirs(this_output_path)
            
            cmd += ['-o', this_output_path]
            cmd += ['--threads', str(max_processes)]
            cmd += ['--memory', str(mem_num_gigs)]
            cmd += ['--careful']
            thetime = _time.asctime( _time.localtime(_time.time()) )
            print('about to launch SPAdes . . . at {}'.format(thetime))
            print(' '.join(cmd))
            contigs['multi_region'] = run_SPAdes(cmd)
        else:
            start_time = _time.time()
            # prepare commandline and launch each SPAdes assembly
            contigs = {}
            for cnum, (pairname, files) in enumerate(sorted(self.read_files.items())):
                if isinstance(use_exe, list):
                    # allow for use of prepended executable with script to run
                    cmd = list(use_exe)
                else:
                    # or just executable
                    cmd = [use_exe]
                # allow use of tuples or dicts by converting dicts to lists
                if isinstance(files, dict):
                    use_files = []
                    for k,v in sorted(files.items()):
                        use_files += [v]
                else:
                    use_files = files
                
                cmd += ['--pe1-1', use_files[0]]
                cmd += ['--pe1-2', use_files[1]]
                try:
                    # use unpaired reads if available
                    cmd += ['--pe1-s', use_files[2]]
                except IndexError:
                    pass
                try:
                    # add a second library if provided
                    if isinstance(self.read_files2[pairname], dict):
                        # if a dict supplied, make it a list
                        use_files2 = []
                        for k,v in sorted(self.read_files2[pairname].items()):
                            use_files2 += [v]
                    else:
                        use_files2 = self.read_files2[pairname]
                    cmd += ['--pe2-1', use_files2[0]]
                    cmd += ['--pe2-2', use_files2[1]]
                    try:
                        cmd += ['--pe2-s', use_files2[2]]
                    except IndexError:
                        pass
                except AttributeError:
                    pass
                
                this_output_path = _os.path.sep.join(output_folder + [pairname])
                if not _os.path.exists(this_output_path):
                    _os.makedirs(this_output_path)
                
                cmd += ['-o', this_output_path]
                cmd += ['--threads', str(max_processes)]
                cmd += ['--memory', str(mem_num_gigs)]
                cmd += ['--careful']
                thetime = _time.asctime( _time.localtime(_time.time()) )
                print('about to launch SPAdes . . . at {}'.format(thetime))
                print(' '.join(cmd))
                contigs[pairname] = run_SPAdes(cmd)
                if len(self.read_files) > 1:
                    # report durations, time left etc
                    _report_time(start_time, cnum, len(self.read_files))

        self.paths_to_contigs = contigs


if __name__ == '__main__':
    main()
