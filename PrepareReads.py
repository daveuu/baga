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
PrepareReads module from the Bacterial and Archaeal Genome (BAG) Analyzer.

This module contains functions to subsample large read sets, and wrappers 
around tools remove adaptor sequences and trim low quality positions.
'''

# stdlib
from time import sleep as _sleep
from baga import _subprocess
from baga import _os
from baga import _multiprocessing
from baga import _cPickle
from baga import _gzip
from baga import _time
from cStringIO import StringIO as _StringIO
from random import sample as _sample

# external Python modules
from Bio import SeqIO as _SeqIO

# package functions
from baga import decide_max_processes as _decide_max_processes
from baga import get_exe_path as _get_exe_path
from baga import report_time as _report_time
def main():
    pass
def cpu_count():
    try:
        return _multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass
    
    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')
        
        if res > 0:
            return res
    except IOError:
        pass
    
    # Windows
    try:
        res = int(_os.environ['NUMBER_OF_PROCESSORS'])
        
        if res > 0:
            return res
    except (KeyError, ValueError):
        pass
class Reads:
    '''
    Prepare reads for alignment to genome sequence by removing adaptor sequences 
    and trimming by position specific scores. Align reads to a reference genome 
    sequence.
    '''
    def __init__(self, reads):
        '''
        Initialise with a baga.CollectData.Reads object
        '''
        try:
            self.read_files = reads.read_files
        except AttributeError:
            print('''
    baga.CallVariants.Reads needs a baga.CollectData.Reads object 
    with a "read_files" attribute. This can be obtained with the 
    "getFromENA()" or "getFromPath()" methods.''')
    def saveLocal(self, name):
        '''
        Save processed read info to a local compressed pickle file.
        'name' can exclude extension: .baga will be added
        '''
        fileout = 'baga.PrepareReads.Reads-%s.baga' % name
        print('Saving to %s' % fileout)
        _cPickle.dump(self, _gzip.open(fileout, 'wb'))
    def subsample(self, genome_size = 6601757, read_cov_depth = 80, pc_loss = 0.2, force = False):
        '''
        Given the size in basepairs of a genome sequence, downsample fastq files to a 
        desired average read coverage depth predicted after read alignment. Read lengths
        are taken from the file. By default, 20% are assumed to be lost at downstream 
        quality control stages (e.g. quality score based trimming). The percent loss is 
        used in coverage depth estimation.
        '''

        subsampled_read_files = {}
        start_time = _time.time()
        for cnum,(run_acc,files) in enumerate(self.read_files.items()):
            #files = {1: 'reads/ERR_______1.fastq.gz', 2: 'reads/ERR_______2.fastq.gz'}
            #for run_acc,files in read_files.items(): #break
            original_name_1 = files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            original_name_2 = files[2].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            #extensions = _os.path.extsep.join(files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[1:])
            path = files[1].split(_os.path.sep)[:-1]
            processed_name_1 = original_name_1+'_subsmp.fastq.gz'
            processed_name_2 = original_name_2+'_subsmp.fastq.gz'
            processed_path_1 = _os.path.sep.join(path + [processed_name_1])
            processed_path_2 = _os.path.sep.join(path + [processed_name_2])
            
            if not all([_os.path.exists(processed_path_1), 
                        _os.path.exists(processed_path_2)]) \
                    or force:
                
                if files[1][-2:] == 'gz':
                    fh1 = _gzip.open(files[1])
                else:
                    fh1 = open(files[1])
                
                aread = _SeqIO.parse(fh1, 'fastq').next()
                read_len = len(aread.seq)
                
                print('Counting reads in %s' % original_name_1)
                fh1.seek(0)
                lines = 0
                # report per half million reads
                interval = 2000000
                nextreport = interval
                for line in fh1:
                    lines += 1
                    if lines == nextreport:
                        print('{:,} reads'.format(lines/4))
                        nextreport += interval
                
                totalreads = lines / 4.0
                print('Found %s reads' % totalreads)
                full_depth_coverage = read_len * 2 * totalreads * (1 - pc_loss) / genome_size
                print('These paired read files would provide approximately {:.1f}x coverage depth'.format(full_depth_coverage))
                numreads2keep = int( round(genome_size * read_cov_depth / (read_len * 2) /  (1 - pc_loss), 0) )
                
                if numreads2keep >= totalreads:
                    print('This pair of read files is estimated to provide only {:.1f}x coverage, but {}x requested.'.format(full_depth_coverage, read_cov_depth))
                    print('No sampling performed. Original files will be used')
                    # pass original files over with subsampled
                    subsampled_read_files[run_acc] = {}
                    subsampled_read_files[run_acc][1] = files[1]
                    subsampled_read_files[run_acc][2] = files[2]
                    fh1.close()
                    if len(self.read_files) > 1:
                        # report durations, time left etc
                        _report_time(start_time, cnum, len(self.read_files))
                    
                    continue
                else:
                    print('For approximately {}x read coverage, will retain {} of {} {}bp read pairs'.format(
                                    read_cov_depth, numreads2keep, totalreads, read_len))
                    
                    fh1.seek(0)
                    if files[2][-2:] == 'gz':
                        fh2 = _gzip.open(files[2])
                    else:
                        fh2 = open(files[2])
                    
                    fout1 = _gzip.open(processed_path_1, 'wb')
                    fout2 = _gzip.open(processed_path_2, 'wb')
                    
                    batch_size = 200000
                    keep_per_pop = int(numreads2keep / float(totalreads) * batch_size) + 1
                    nextwrite = batch_size
                    written = 0
                    n1 = 0
                    n2 = 0
                    these_lines1 = []
                    these_lines2 = []
                    reportfreq = 10
                    thisreport = 0
                    print('Subsampling . . .')
                    for line in fh1:
                        these_lines1 += [line]
                        if len(these_lines1) % 4 == 0:
                            n1 += 1
                            
                        if n1 == nextwrite:
                            keep_indices = sorted(_sample(xrange(batch_size), keep_per_pop))
                            keep_these = []
                            for i in keep_indices:
                                i1 = i * 4
                                i2 = i * 4 + 4
                                keep_these += these_lines1[i1:i2]
                            
                            # try parsing a read for QC
                            assert _SeqIO.read(_StringIO(''.join(keep_these[:4])), 'fastq')
                            fout1.write(''.join(keep_these))
                            these_lines1 = []
                            written += keep_per_pop
                            thisreport += 1
                            if thisreport == reportfreq or written == keep_per_pop:
                                # report first time and at intevals
                                print('Written {:,} reads ({:.1%}) to {}'.format(written,
                                                                                 written/float(numreads2keep),
                                                                                 processed_path_1))
                            
                            for line2 in fh2:
                                these_lines2 += [line2]
                                if len(these_lines2) % 4 == 0:
                                    n2 += 1
                                
                                if n2 == nextwrite:
                                    keep_these = []
                                    for i in keep_indices:
                                        i1 = i * 4
                                        i2 = i * 4 + 4
                                        keep_these += these_lines2[i1:i2]
                                    
                                    assert _SeqIO.read(_StringIO(''.join(keep_these[:4])), 'fastq')
                                    fout2.write(''.join(keep_these))
                                    these_lines2 = []
                                    if thisreport == reportfreq or written == keep_per_pop:
                                        thisreport = 0
                                        print('Written {:,} reads ({:.1%}) to {}'.format(written,
                                                                                         written/float(numreads2keep),
                                                                                         processed_path_2))
                                    nextwrite += batch_size
                                    break
                    
                    # write remainder
                    remainder = nextwrite - n1
                    keep_in_remainder = int(keep_per_pop * (remainder / float(batch_size))) + 1
                    keep_indices = sorted(_sample(xrange(remainder), keep_in_remainder))
                    keep_these = []
                    for i in keep_indices:
                        i1 = i * 4
                        i2 = i * 4 + 4
                        keep_these += these_lines1[i1:i2]
                    
                    # try parsing a read for QC
                    assert _SeqIO.read(_StringIO(''.join(keep_these[:4])), 'fastq')
                    fout1.write(''.join(keep_these))
                    written += keep_in_remainder
                    print('Written {:,} reads ({:.1%}) to {}'.format(written,
                                                                             written/float(numreads2keep),
                                                                             processed_path_1))
                    
                    # get remainder
                    for line2 in fh2:
                        these_lines2 += [line2]
                    
                    # write remainder
                    keep_these = []
                    for i in keep_indices:
                        i1 = i * 4
                        i2 = i * 4 + 4
                        keep_these += these_lines2[i1:i2]
                    
                    assert _SeqIO.read(_StringIO(''.join(keep_these[:4])), 'fastq') ###### check why keep_these was empty
                    fout2.write(''.join(keep_these))
                    print('Written {:,} reads ({:.1%}) to {}'.format(written,
                                                                             written/float(numreads2keep),
                                                                             processed_path_2))
                    
                    # not sure if this is quicker/slower (more calls to .join())
                    # this_read = []
                    # for line in fh1:
                        # this_read += [line]
                        # if len(this_read) == 4:
                            # these_reads1 += [''.join(this_read)]
                            # #these_reads1 += this_read
                            # this_read = []
                            # n1 += 1
                            
                        # if n1 == nextwrite:
                            # keep_indices = sorted(_sample(xrange(batch_size), keep_per_pop))
                            # # try parsing a read for QC
                            # assert _SeqIO.read(_StringIO(these_reads1[0]), 'fastq')
                            # fout1.write(''.join([these_reads1[i] for i in keep_indices]))
                            # these_reads1 = []
                            # written += keep_per_pop
                            # print('Written {:,} reads ({:.2%}) to {}'.format(written,
                                                                             # written/float(numreads2keep),
                                                                             # processed_path_1))
                            # for line2 in fh2:
                                # this_read += [line2]
                                # if len(this_read) == 4:
                                    # these_reads2 += [''.join(this_read)]
                                    # this_read = []
                                    # n2 += 1
                                
                                # if n2 == nextwrite:
                                    # assert _SeqIO.read(_StringIO(these_reads2[0]), 'fastq')
                                    # fout2.write(''.join([these_reads2[i] for i in keep_indices]))
                                    # these_reads2 = []
                                    # print('Written {:,} reads ({:.2%}) to {}'.format(written,
                                                                                     # written/float(numreads2keep),
                                                                                     # processed_path_2))
                                    # nextwrite += batch_size
                                    # break
                    
                    fout1.close()
                    fout2.close()
                    fh1.close()
                    fh2.close()
                
            else:
                print('Found:')
                print(processed_path_1)
                print(processed_path_2)
                print('use "force = True" to overwrite')
            
            if len(self.read_files) > 1:
                # report durations, time left etc
                _report_time(start_time, cnum, len(self.read_files))
            
            subsampled_read_files[run_acc] = {}
            subsampled_read_files[run_acc][1] = processed_path_1
            subsampled_read_files[run_acc][2] = processed_path_2

        # replace here as this step is optional
        self.fullsized_read_files = list(self.read_files)
        self.read_files = subsampled_read_files

    def cutAdaptors_old(self, path_to_exe = False, force = False, max_cpus = -1):

        if not path_to_exe:
            path_to_exe = _get_exe_path('cutadapt')

        adaptorcut_read_files = {}
        adaptor_seqs = [
            'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
            'AGATCGGAAGAGCACACGTCT',
            'AGATCGGAAGAGC',
            'GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG',
            'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        ]

        start_time = _time.time()

        for cnum,(run_acc,files) in enumerate(self.read_files.items()):
            original_name_1 = files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            original_name_2 = files[2].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            path = files[1].split(_os.path.sep)[:-1]
            # add a suffix for this analysis step
            extensions = _os.path.extsep.join([''] + files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[1:])
            processed_name_1 = original_name_1 + '_adpt' + extensions
            processed_name_2 = original_name_2 + '_adpt' + extensions
            processed_path_1 = _os.path.sep.join(path + [processed_name_1])
            processed_path_2 = _os.path.sep.join(path + [processed_name_2])
            # single end
            cmd = [path_to_exe] + \
                  [a for b in [('-a', a) for a in adaptor_seqs] for a in b] + \
                  ['-o', processed_path_1, files[1]]
            # paired end
            cmd = [path_to_exe] + \
                  [a for b in [('-a', a) for a in adaptor_seqs] for a in b] + \
                  [a for b in [('-A', a) for a in adaptor_seqs] for a in b] + \
                  ['-o', processed_path_1, '-p', processed_path_2] + \
                  [files[1], files[2]]
            
            if not all([_os.path.exists(processed_path_1), 
                        _os.path.exists(processed_path_2)]) \
                    or force:
                
                print('Called: "%s"' % ' '.join(cmd))
                with open(_os.path.sep.join(path + [run_acc+'_cutadapt.log']),"w") as out:
                    _subprocess.call(cmd, stdout = out)
            else:
                print('Found:')
                print(processed_path_1)
                print(processed_path_2)
                print('use "force = True" to overwrite')
            
            if len(self.read_files) > 1:
                # report durations, time left etc
                _report_time(start_time, cnum, len(self.read_files))
            
            adaptorcut_read_files[run_acc] = {}
            adaptorcut_read_files[run_acc][1] = processed_path_1
            adaptorcut_read_files[run_acc][2] = processed_path_2

        self.adaptorcut_read_files = adaptorcut_read_files
    def cutAdaptors(self, path_to_exe = False, force = False, max_cpus = -1):

        if not path_to_exe:
            path_to_exe = _get_exe_path('cutadapt')

        adaptorcut_read_files = {}
        adaptor_seqs = [
            'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
            'AGATCGGAAGAGCACACGTCT',
            'AGATCGGAAGAGC',
            'GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG',
            'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        ]

        cmds = []
        processed_paths_to_do = []
        for cnum,(run_acc,files) in enumerate(self.read_files.items()):
            original_name_1 = files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            original_name_2 = files[2].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            path = files[1].split(_os.path.sep)[:-1]
            # add a suffix for this analysis step
            extensions = _os.path.extsep.join([''] + files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[1:])
            processed_name_1 = original_name_1 + '_adpt' + extensions
            processed_name_2 = original_name_2 + '_adpt' + extensions
            processed_path_1 = _os.path.sep.join(path + [processed_name_1])
            processed_path_2 = _os.path.sep.join(path + [processed_name_2])
            # single end
            cmd = [path_to_exe] + \
                  [a for b in [('-a', a) for a in adaptor_seqs] for a in b] + \
                  ['-o', processed_path_1, files[1]]
            # paired end
            cmd = [path_to_exe] + \
                  [a for b in [('-a', a) for a in adaptor_seqs] for a in b] + \
                  [a for b in [('-A', a) for a in adaptor_seqs] for a in b] + \
                  ['-o', processed_path_1, '-p', processed_path_2] + \
                  [files[1], files[2]]
            
            if not all([_os.path.exists(processed_path_1), 
                        _os.path.exists(processed_path_2)]) \
                    or force:
                
                # collect expected outputs
                processed_paths_to_do += [(processed_path_1,processed_path_2)]
                # collect all the commands to be issued
                cmds += [(run_acc,cmd)]
                
            else:
                print('Found:')
                print(processed_path_1)
                print(processed_path_2)
                print('use "force = True" to overwrite')
                
                adaptorcut_read_files[run_acc] = {}
                adaptorcut_read_files[run_acc][1] = processed_path_1
                adaptorcut_read_files[run_acc][2] = processed_path_2



        if len(cmds):
            max_processes = _decide_max_processes(max_cpus)
            
            processes = {}
            
            ### how to combine this which hangs on _os.wait()
            for run_acc,cmd in cmds:
                
                print('Called: "%s"' % ' '.join(cmd))
                # process is key, open file being piped to is value
                this_stdout_file = open(_os.path.sep.join(path + [run_acc+'_cutadapt.log']),"w")
                thisprocess = _subprocess.Popen(cmd, shell=False, stdout = this_stdout_file)
                processes[thisprocess] = this_stdout_file
                
                if len(processes) >= max_processes:
                    _os.wait()
                    finished = dict([(p,f) for p,f in processes.items() if p.poll() is not None])
                    # close files for finished processes
                    for process,stdout_file in finished.items():
                        stdout_file.close()
                        # update active processes
                        del processes[process]
            
            # Check if all the child processes were closed
            for p in processes:
                if p.poll() is None:
                    p.wait()


        fails = []
        for (run_acc,cmd),(processed_path_1,processed_path_2) in zip(cmds,processed_paths_to_do):
            if _os.path.exists(processed_path_1) and _os.path.exists(processed_path_2):
                print('Found:')
                print(processed_path_1)
                print(processed_path_2)
                adaptorcut_read_files[run_acc] = {}
                adaptorcut_read_files[run_acc][1] = processed_path_1
                adaptorcut_read_files[run_acc][2] = processed_path_2
            else:
                print('Processing of the following pair seems to have failed')
                print(processed_path_1)
                print(processed_path_2)
                fails += [(processed_path_1,processed_path_2)]

        assert len(fails) == 0, 'There was a problem finding all of the output from cutadapt. Try repeating this or an eralier step with the --force option to overwite previous, possibly incomplete, files'

        self.adaptorcut_read_files = adaptorcut_read_files

    def trim(self, path_to_exe = False, 
                   force = False, 
                   max_cpus = -1):

        if not path_to_exe:
            exe_sickle = _get_exe_path('sickle')
        else:
            exe_sickle = _os.path.sep.join(path_to_exe)

        e1 = 'Could not find "adaptorcut_read_files" attribute. \
        Before quality score trimming, reads must be cleaned of \
        library preparation sequences. Please run cutAdaptors() \
        method on this Reads instance.'

        assert hasattr(self, 'adaptorcut_read_files'), e1

        e2 = 'Could not find %s. Either run cutAdaptors() again \
        or ensure file exists'

        for run_acc, files in self.adaptorcut_read_files.items():
            assert _os.path.exists(files[1]), e2 % files[1]
            assert _os.path.exists(files[1]), e2 % files[1]

        trimmed_read_files = {}

        cmds = []
        processed_paths_to_do = []
        for run_acc,files in self.adaptorcut_read_files.items():
            original_name_1 = files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            original_name_2 = files[2].split(_os.path.sep)[-1].split(_os.path.extsep)[0]
            path = files[1].split(_os.path.sep)[:-1]
            # add a suffix for this analysis step
            extensions = _os.path.extsep.join([''] + files[1].split(_os.path.sep)[-1].split(_os.path.extsep)[1:])
            processed_name_1 = original_name_1 + '_qual' + extensions
            processed_name_2 = original_name_2 + '_qual' + extensions
            processed_name_s = original_name_1 + '_singletons_qual' + extensions
            processed_path_1 = _os.path.sep.join(path + [processed_name_1])
            processed_path_2 = _os.path.sep.join(path + [processed_name_2])
            processed_path_s = _os.path.sep.join(path + [processed_name_s])
            # Illumina quality using CASAVA >= 1.8 is Sanger encoded
            QSscore_scale = 'sanger'
            cmd = [exe_sickle, 'pe',
            '-f', files[1] ,'-r', files[2], 
            '-t', QSscore_scale,
            '-o', processed_path_1,
            '-p', processed_path_2, 
            '-s', processed_path_s, 
            # quality 25, length 50 (of 150)
            '-q','25','-l','50']
            if not all([_os.path.exists(processed_path_1), 
                        _os.path.exists(processed_path_2),
                        _os.path.exists(processed_path_s)]) \
                    or force:
                
                # collect expected outputs
                processed_paths_to_do += [(processed_path_1,processed_path_2,processed_path_s)]
                # collect all the commands to be issued
                cmds += [(run_acc,cmd)]
                
            else:
                print('Found:')
                print(processed_path_1)
                print(processed_path_2)
                print(processed_path_s)
                print('use "force = True" to overwrite')
                
                trimmed_read_files[run_acc] = {}
                trimmed_read_files[run_acc][1] = processed_path_1
                trimmed_read_files[run_acc][2] = processed_path_2



        if len(cmds):
            max_processes = _decide_max_processes(max_cpus)
            
            processes = {}
            
            ### how to combine this which hangs on _os.wait()
            for run_acc,cmd in cmds:
                
                print('Called: "%s"' % ' '.join(cmd))
                # process is key, open file being piped to is value
                this_stdout_file = open(_os.path.sep.join(path + [run_acc+'_sickle.log']),"w")
                thisprocess = _subprocess.Popen(cmd, shell = False, stdout = this_stdout_file)
                processes[thisprocess] = this_stdout_file
                
                if len(processes) >= max_processes:
                    _os.wait()
                    finished = dict([(p,f) for p,f in processes.items() if p.poll() is not None])
                    # close files for finished processes
                    for process,stdout_file in finished.items():
                        stdout_file.close()
                        # update active processes
                        del processes[process]
            
            # Check if all the child processes were closed
            for p in processes:
                if p.poll() is None:
                    p.wait()


        fails = []
        for (run_acc,cmd),(processed_path_1,processed_path_2,processed_path_s) in zip(cmds,processed_paths_to_do):
            if _os.path.exists(processed_path_1) and _os.path.exists(processed_path_2):
                print('Found:')
                print(processed_path_1)
                print(processed_path_2)
                trimmed_read_files[run_acc] = {}
                trimmed_read_files[run_acc][1] = processed_path_1
                trimmed_read_files[run_acc][2] = processed_path_2
            else:
                print('Processing of the following pair seems to have failed')
                print(processed_path_1)
                print(processed_path_2)
                fails += [(processed_path_1,processed_path_2)]

        assert len(fails) == 0, 'There was a problem finding all of the output from sickle. Try repeating this or an eralier step with the --force option to overwite previous, possibly incomplete, files'


        self.trimmed_read_files = trimmed_read_files

if __name__ == '__main__':
    main()
