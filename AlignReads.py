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
AlignReads module from the Bacterial and Archaeal Genome (BAG) Analyzer.

This module contains wrappers around tools to align reads to a reference 
genome, including realignment around possible insertions and deletions.
'''

# stdlib
from time import sleep as _sleep
from baga import _subprocess
from baga import _os
from baga import _multiprocessing
from baga import _cPickle
from baga import _gzip
from baga import _StringIO
from baga import _tarfile
from baga import _array
from baga import _json

# external Python modules
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord
from Bio import SeqIO as _SeqIO
from random import sample as _sample
# import pysam
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord

# package functions
from baga import decide_max_processes as _decide_max_processes
from baga import get_exe_path as _get_exe_path
from baga import get_jar_path as _get_jar_path
def main():
    pass
class SAMs:
    '''
    A collection of short read datasets that have been aligned to the same genome 
    sequence in Sequence Alignment/Map (SAM) format for variant calling using GATK.
    '''

    def __init__(self, reads = False, genome = False, baga = False):
        '''
        Initialise with:
        a baga.PrepareReads.Reads object and,
        a baga.CollectData.Genome object.
        
        OR
        
        a baga.AlignReads.SAMs (like this one) object that 
        was previously saved.
        '''
        
        if (reads and genome) and not baga:
            try:
                self.read_files = reads.trimmed_read_files
            except AttributeError:
                text = 'WARNING: baga was not used to quality-score trim these reads. Read trimming is recommended for most types of analysis. This can be achieved with the "trim()" method of the Reads class in the PrepareReads module.'
                print(text)
                try:
                    self.read_files = reads.adaptorcut_read_files
                except AttributeError:
                    text = 'WARNING: baga was not used to remove library preparation adaptor sequences from these reads. Adaptor removal is highly recommended so hopefully you already removed adaptor sequences! This can be achieved with the "cutAdaptors()" method of the Reads class in the PrepareReads module.'
                    self.read_files = reads.read_files
                    print(text)
                    print('continuing with these reads . . .')
            
            self.genome_sequence = genome.sequence
            self.genome_id = genome.id
        
        elif baga and not (reads and genome):
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
            
        else:
            raise NameError('instantiate baga.AlignReads.SAMs with either genome and reads or previous saved alignments (baga.AlignReads.SAMs-*.baga)')


    def saveLocal(self, name):
        '''
        Save processed SAM file info to a local compressed pickle file.
        'name' can exclude extension: .baga will be added
        '''
        fileout = 'baga.AlignReads.SAMs-{}.baga'.format(name)
        
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

    def align(self, insert_size = False, 
                    path_to_exe = False, 
                    local_alns_path = ['alignments'], 
                    force = False, 
                    max_cpus = -1):
        if not path_to_exe:
            path_to_exe = _get_exe_path('bwa')

        # write genome sequence to a fasta file
        try:
            _os.makedirs('genome_sequences')
        except OSError:
            pass

        genome_fna = 'genome_sequences/%s.fna' % self.genome_id

        _SeqIO.write(_SeqRecord(_Seq(self.genome_sequence.tostring()), id = self.genome_id), 
                    genome_fna, 
                    'fasta')

        # make folder for alignments (BAMs)
        local_alns_path = _os.path.sep.join(local_alns_path)
        if not _os.path.exists(local_alns_path):
            _os.makedirs(local_alns_path)

        # make a subdir for this genome
        local_alns_path_genome = _os.path.sep.join([
                                local_alns_path, 
                                self.genome_id])
        if not _os.path.exists(local_alns_path_genome):
            _os.makedirs(local_alns_path_genome)


        max_processes = _decide_max_processes( max_cpus )


        e1 = 'Could not find "read_files" attribute. Before aligning to genome, reads must be quality score trimmed. Please run trim() method on this Reads instance.'

        assert hasattr(self, 'read_files'), e1

        e2 = 'Could not find %s. Either run trim() again or ensure file exists'

        for run_acc, files in self.read_files.items():
            assert _os.path.exists(files[1]), e2 % files[1]
            assert _os.path.exists(files[1]), e2 % files[1]

        have_index_files = [_os.path.exists(genome_fna + '.' + a) for a in ('ann','pac','amb','bwt','sa')]

        if not all(have_index_files):
            print('Writing BWA index files for %s' % genome_fna)
            _subprocess.call([path_to_exe, 'index', genome_fna])

        aligned_read_files = {}
        for run_acc,files in self.read_files.items():
            RGinfo = r"@RG\tID:%s\tSM:%s\tPL:ILLUMINA" % (run_acc,run_acc)
            if insert_size:
                cmd = [path_to_exe, 'mem', '-t', str(max_processes), '-M', '-a', '-I', insert_size, '-R', RGinfo, genome_fna, files[1], files[2]]
            else:
                # BWA can estimate on-the-fly
                cmd = [path_to_exe, 'mem', '-t', str(max_processes), '-M', '-a', '-R', RGinfo, genome_fna, files[1], files[2]]
            
            out_sam = _os.path.sep.join([local_alns_path_genome, '%s__%s.sam' % (run_acc, self.genome_id)])
            
            if not _os.path.exists(out_sam) or force:
                print('Called: "%s"' % ' '.join(cmd))
                with open(out_sam, "wb") as out:
                    _subprocess.call(cmd, stdout = out)
                
            else:
                print('Found:')
                print(out_sam)
                print('use "force = True" to overwrite')
            
            print(' '.join(cmd))
            
            aligned_read_files[run_acc] = out_sam

        self.aligned_read_files = aligned_read_files

    def toBAMs(self, path_to_exe = False, force = False, max_cpus = -1):
        if not path_to_exe:
            path_to_exe = _get_exe_path('samtools')

        processes = set()
        max_processes = _decide_max_processes( max_cpus )

        paths_to_BAMs = []

        for run_acc,SAM in self.aligned_read_files.items():
            BAM_out = SAM[:-4] + '.bam'
            if not _os.path.exists(BAM_out) or force:
                cmd = '{0} view -buh {1} | {0} sort - {2}; {0} index {2}.bam'.format(path_to_exe, SAM, SAM[:-4])
                print('Called: %s' % cmd)
                processes.add( _subprocess.Popen(cmd, shell=True) )
                if len(processes) >= max_processes:
                    (pid, exit_status) = _os.wait()
                    processes.difference_update(
                        [p for p in processes if p.poll() is not None])
            else:
                print('Found:')
                print(BAM_out)
                print('use "force = True" to overwrite')
            
            paths_to_BAMs += [BAM_out]

        # Check if all the child processes were closed
        for p in processes:
            if p.poll() is None:
                p.wait()

        self.paths_to_BAMs = paths_to_BAMs

    def removeDuplicates(self, path_to_jar = False, force = False, mem_num_gigs = 2, max_cpus = -1):
        if not path_to_jar:
            path_to_jar = _get_jar_path('picard')

        processes = set()
        max_processes = _decide_max_processes( max_cpus )

        print('Print: will use %s cpus for picard' % max_processes)

        paths_to_BAMs_dd = []
        for BAM in self.paths_to_BAMs:
            BAM_out = BAM[:-4] + '_dd.bam'
            if not _os.path.exists(BAM_out) or force:
                # an alternative parallel code with polling not waiting
                # while len(processes) >= max_processes:
                    # _sleep(1)
                    # print('processes: %s' % (len(processes)))
                    # processes.difference_update(
                        # [p for p in processes if p.poll() is not None]
                    # )
                picard_command = ['MarkDuplicates', 'I=', BAM, 'O=', BAM_out, 'M=', BAM[:-4] + '_dd.log'] #, 'VALIDATION_STRINGENCY=','LENIENT']
                cmd = ['java', '-Xmx%sg' % mem_num_gigs, '-jar', path_to_jar] + picard_command
                processes.add( _subprocess.Popen(cmd, shell=False) )
                print('Called: %s' % (' '.join(map(str, cmd))))
                #_subprocess.call(cmd, shell=False)
                print('processes: %s, max_processes: %s' % (len(processes),max_processes))
                if len(processes) >= max_processes:
                    (pid, exit_status) = _os.wait()
                    processes.difference_update(
                        [p for p in processes if p.poll() is not None])
                
            else:
                print('Found:')
                print(BAM_out)
                print('use "force = True" to overwrite')
            
            paths_to_BAMs_dd += [BAM_out]

        # Check if all the child processes were closed
        for p in processes:
            if p.poll() is None:
                p.wait()

        self.paths_to_BAMs_dd = paths_to_BAMs_dd


    def sortIndexBAMs(self, path_to_exe = False, force = False, max_cpus = -1):
        if not path_to_exe:
            path_to_exe = _get_exe_path('samtools')

        processes = set()
        max_processes = _decide_max_processes( max_cpus )

        paths_to_BAMs_dd_si = []
        for SAM in self.paths_to_BAMs_dd:
            BAM_out = SAM[:-4] + '_si.bam'
            if not _os.path.exists(BAM_out) or force:
                cmd = '{0} sort {1} {2}_si; {0} index {2}_si.bam'.format(path_to_exe, SAM, SAM[:-4])
                print('Called: %s' % cmd)
                processes.add( _subprocess.Popen(cmd, shell=True) )
                if len(processes) >= max_processes:
                    (pid, exit_status) = _os.wait()
                    processes.difference_update(
                        [p for p in processes if p.poll() is not None])
            else:
                print('Found:')
                print(BAM_out)
                print('use "force = True" to overwrite')
            
            paths_to_BAMs_dd_si += [BAM_out]

        # Check if all the child processes were closed
        for p in processes:
            if p.poll() is None:
                p.wait()

        self.paths_to_BAMs_dd_si = paths_to_BAMs_dd_si

    def IndelRealignGATK(self, 
            jar = ['external_programs', 'GenomeAnalysisTK', 'GenomeAnalysisTK.jar'], 
            picard_jar = False, 
            samtools_exe = False,
            force = False,
            mem_num_gigs = 2,
            max_cpus = -1):
        
        # GATK is manually downloaded by user and placed in folder of their choice
        jar = _os.path.sep.join(jar)

        if not picard_jar:
            picard_jar = _get_jar_path('picard')

        if not samtools_exe:
            samtools_exe = _get_exe_path('samtools')

        genome_fna = 'genome_sequences/{}.fna'.format(self.genome_id)

        e1 = 'Could not find "paths_to_BAMs_dd_si" attribute. \
        Before starting GATK analysis, read alignments must be \
        duplicates removed. Please run:\n\
        toBAMS()\n\
        removeDuplicates()\n\
        sortIndexBAMs()\n\
        method on this SAMs instance.'

        assert hasattr(self, 'paths_to_BAMs_dd_si'), e1

        e2 = 'Could not find %s. Please ensure file exists'

        for BAM in self.paths_to_BAMs_dd_si:
            assert _os.path.exists(BAM), e2 % BAM

        if not _os.path.exists(genome_fna[:-4] + '.dict'):
            print('Creating sequence dictionary for %s' % genome_fna)
            _subprocess.call(['java', '-jar', picard_jar, 'CreateSequenceDictionary', 'R=', genome_fna, 'O=', genome_fna[:-4] + '.dict'])

        #have_index_files = [_os.path.exists(genome_fna + '.' + a) for a in ('ann','pac','amb','bwt','sa','fai')]
        have_index_files = [_os.path.exists(genome_fna + '.' + a) for a in ('fai',)]

        if not all(have_index_files):
            print('Writing index files for %s' % genome_fna)
            _subprocess.call([samtools_exe, 'faidx', genome_fna])

        processes = set()
        max_processes = _decide_max_processes( max_cpus )

        for BAM in self.paths_to_BAMs_dd_si:
            intervals = BAM[:-4] + '.intervals'
            if not _os.path.exists(intervals) or force:
                cmd = ['java', '-Xmx%sg' % mem_num_gigs, '-jar', jar, 
                        '-T', 'RealignerTargetCreator', 
                        '-R', genome_fna, 
                        '-I', BAM, 
                        '-o', intervals] #, '--validation_strictness', 'LENIENT']
                print(' '.join(map(str, cmd)))
                processes.add( _subprocess.Popen(cmd, shell=False) )
                if len(processes) >= max_processes:
                    (pid, exit_status) = _os.wait()
                    processes.difference_update(
                        [p for p in processes if p.poll() is not None])
            else:
                print('Found:')
                print(intervals)
                print('use "force = True" to overwrite')

        # Check if all the child processes were closed
        for p in processes:
            if p.poll() is None:
                p.wait()


        paths_to_BAMs_dd_si_ra = []
        for BAM in self.paths_to_BAMs_dd_si:
            intervals = BAM[:-4] + '.intervals'
            bam_out = BAM[:-4] + '_realn.bam'
            if not _os.path.exists(bam_out) or force:
                cmd = ['java', '-Xmx4g', '-jar', jar, 
                        '-T', 'IndelRealigner', 
                        '-R', genome_fna, 
                        '-I', BAM, 
                        '-targetIntervals', intervals, 
                        '-o', bam_out,
                        '--filter_bases_not_stored']
                print(' '.join(map(str, cmd)))
                processes.add( _subprocess.Popen(cmd, shell=False) )
                if len(processes) >= max_processes:
                    _os.wait()
                    processes.difference_update(
                        [p for p in processes if p.poll() is not None])
            else:
                print('Found:')
                print(bam_out)
                print('use "force = True" to overwrite')
            
            paths_to_BAMs_dd_si_ra += [bam_out]

        for p in processes:
            if p.poll() is None:
                p.wait()

        # the last list of BAMs in ready_BAMs is input for CallgVCFsGATK
        # both IndelRealignGATK and recalibBaseScoresGATK put here
        self.ready_BAMs = [paths_to_BAMs_dd_si_ra]

if __name__ == '__main__':
    main()
