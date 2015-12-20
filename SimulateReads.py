#! /usr/bin/env python2
# -*- coding: utf-8 -*-
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
SimulateReads module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains wrappers around tools to generate Illumina-type short reads.
These are useful for benchmarking and validation of variant calling pipelines.
'''

# stdlib
from baga import _subprocess
from baga import _os
from baga import _cPickle
from baga import _gzip
from baga import _tarfile
from baga import _array
from baga import _json

# external Python modules
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord
from Bio import SeqIO as _SeqIO
import random as _random

# package functions
from baga import decide_max_processes as _decide_max_processes
from baga import get_exe_path as _get_exe_path
def main():
    pass
class Simulator:
    '''
    A simulator of short read datasets.
    '''

    def __init__(self, genome = False, 
                       baga = False, 
                       num_individuals = 1,
                       large_deletions = {},
                       large_del_padding = 1000,
                       random_seed = False):
        '''
        Initialise with:
        a baga.CollectData.Genome object.
        
        OR
        
        a path to baga.SimulateReads.Reads (like this one) object that 
        was previously saved.
        
        Large deletions can be included to simulate e.g. missing genomic islands or prophage.
        Currently, a set of genomes are generated with and a set without the large deletions if 
        specified.
        
        large_deletions should be a dict with arbitrary names of deletions as keys and
        tuples of (start,end) for python slices delineating each deletion. If supplied a
        set of genomes with and without deletions are generated, each set consisting of 
        num_individuals members.
        
        No variants will be generated within large_del_padding: a 'safe' distance around large 
        deletions outside of which variant calling should be reliable. Small deletions might run
        over into this zone.
        '''
        
        if random_seed:
            if not isinstance(random_seed, int):
                raise ValueError('random_seed must be integer')
            _random.seed(random_seed) # 684651
        
        if genome and not baga:
            self.genome = genome
        elif baga and not genome:
            # for reloading a previous instantiation
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
            raise NameError('instantiate baga.SimulateReads.Reads with a loaded baga.CollectData.Genome-*.baga object or previously saved baga.SimulateReads.Reads object')
        
        omit = set()
        for name,(start,end) in large_deletions.items():
            omit.update(range(start-large_del_padding, end+large_del_padding))
        
        samplefrom = sorted(set(range(len(self.genome.sequence))) - omit)
        self.samplefrom = samplefrom
        
        self.large_deletions = large_deletions
        self.num_individuals = num_individuals
        
        # to be optionally populated with methods
        self.SNPs_per_genome = []
        self.indel_dict_by_pos_pergenome = []
        
        

    def saveLocal(self, name):
        '''
        Save baga object
        '''
        fileout = 'baga.SimulateReads.Reads-{}.baga'.format(name)
        
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

    def generateSNPs(self, num_SNPs = 0, 
                           distribution = 'pseudo_tree'):
        '''
        Generate specified number of variants, randomly distributed among genome sites.

        Variants are shared among individual genomes. Variants can be distributed among 
        multiple genomes using a simple "pseudo tree" algorithm whereby a proportion are 
        shared i.e., ancestral, a proportion are unique and a proportion are common among 
        subgroups. This approximates an almost star-like phylogeny with some homoplasy 
        and several shared variants as if isolates from the same population were called 
        against an out-group reference chromosome.

        More realistic models will be implemented in the future.
        '''

        if len(self.large_deletions):
            # make one set with, one without large deletions
            # share total variants among double the specified population size
            use_num_individuals = self.num_individuals*2
        else:
            use_num_individuals = self.num_individuals


        # random SNPs (does same as GemHaps.py)
        SNPs = []
        for pos0 in _random.sample(self.samplefrom, num_SNPs):
            ref = self.genome.sequence[pos0]
            var = _random.choice(list(set(['T','A','G','C']) - set([ref])))
            SNPs += [(pos0, var)]

        SNPs_per_genome = []
        if distribution == 'pseudo_tree':
            # they all share a random 30%
            num_ancestral = int(round(0.3*num_SNPs))
            ancestral = set(_random.sample(SNPs, num_ancestral))
            # plus an additional 30% per genome, some shared
            num_pergenome = int(round(0.3*num_SNPs))
            for i in range(use_num_individuals):
                SNPs_per_genome += [sorted(set(_random.sample(SNPs, num_pergenome)) | ancestral)]
        elif distribution == 'random':
            # round down
            num_per_genome = int(num_SNPs / use_num_individuals)
            # some overlap . . cold allow variable amount of overlap . . .
            for i in range(use_num_individuals):
                SNPs_per_genome += [sorted(set(_random.sample(SNPs, num_pergenome)))]
        elif distribution == 'discrete':
            # round down
            num_per_genome = int(num_SNPs / use_num_individuals)
            # no overlap
            for i in range(use_num_individuals):
                SNPs_per_genome += [sorted(SNPs[i:(i+use_num_individuals)])]
        else:
            raise ValueError('Need a distribution!')

        self.SNPs_per_genome = SNPs_per_genome

        filename = 'SNPs_per_genome.csv'
        with open(filename, 'w') as fout:
            print('Writing to {}'.format(filename))
            fout.write('"Replicate","Position","Reference","Variant"\n')
            for gn,positions in enumerate(SNPs_per_genome):
                for pos0,var in positions:
                    fout.write('{},{},"{}","{}"\n'.format(gn+1,pos0,self.genome.sequence[pos0],var))


    def generateInDels(self, num_deletions = 0, 
                             num_insertions = 0,
                             distribution = 'pseudo_tree',
                             sizes = [2,2,2,2,3,3,3,4,4,5,6,7,8,9]):
        '''
        Generate specified number of small indels, randomly distributed among genome sites.

        sizes are the lengths to sample from
        '''

        if len(self.large_deletions):
            # make one set with, one without large deletions
            # share total variants among double the specified population size
            use_num_individuals = self.num_individuals*2
        else:
            use_num_individuals = self.num_individuals

        # random InDels
        include = list(self.samplefrom)
        # randomise first
        _random.shuffle(include)

        indel_dict_by_pos_pergenome = [dict() for i in range(use_num_individuals)]

        ## currently hardwired for 1/3 ancestral, 2/3 per genome exclusive
        if distribution == 'pseudo_tree':
            c = 0  # increment along randomised list of positions
            ## deletions
            num_ancestral = int(round(1.0/3 * num_deletions))
            num_exclusive = int(round(2.0/3 * num_deletions))
            anc_dels = {}
            for i in range(num_ancestral):
                anc_dels[include[c]] = _random.choice(sizes)
                c += 1
            for g in range(use_num_individuals):
                # add ancestors to each
                for pos0,length in anc_dels.items():
                    indel_dict_by_pos_pergenome[g][pos0] = length
                for i in range(num_exclusive):
                    indel_dict_by_pos_pergenome[g][include[c]] = _random.choice(sizes)
                    c += 1
            ## insertions
            num_ancestral = int(round(1.0/3 * num_insertions))
            num_exclusive = int(round(2.0/3 * num_insertions))
            anc_ins = {}
            for i in range(num_ancestral):
                anc_ins[include[c]] = ''.join(_random.sample(['A','T','C','G']*100, _random.choice(sizes)))
                c += 1
            for g in range(use_num_individuals):
                # add ancestors to each
                for pos0,seq in anc_ins.items():
                    indel_dict_by_pos_pergenome[g][pos0] = seq
                for i in range(num_exclusive):
                    indel_dict_by_pos_pergenome[g][include[c]] = ''.join(_random.sample(['A','T','C','G']*100, _random.choice(sizes)))
                    c += 1

        self.indel_dict_by_pos_pergenome = indel_dict_by_pos_pergenome

        filename = 'indel_dict_by_pos_pergenome.csv'
        with open(filename, 'w') as fout:
            print('Writing to {}'.format(filename))
            fout.write('"Replicate","Position","Reference","Variant"\n')
            for gn,positions in enumerate(indel_dict_by_pos_pergenome):
                for pos0,var in sorted(positions.items()):
                    fout.write('{},{},"{}","{}"\n'.format(gn+1,pos0,self.genome.sequence[pos0],var))

    def generateSequences(self):
        '''
        Create full length sequences with generated variants applied to the reference sequence.

        Generated variants are saved to a csv.

        If large deletions are present, variant positions appropriately are corrected when applied.
        '''

        save_these = {}
        save_these['SNPs'] = self.SNPs_per_genome
        save_these['InDels'] = self.indel_dict_by_pos_pergenome

        if len(self.large_deletions):
            # generate a version of reference genome with large deletions
            ranges = sorted(self.large_deletions.values())
            genome_large_deletions = self.genome.sequence[:ranges[0][0]]
            for n,(s,e) in enumerate(ranges[:-1]):
                genome_large_deletions.extend(self.genome.sequence[e:ranges[n+1][0]])
            
            genome_large_deletions.extend(self.genome.sequence[ranges[n+1][1]:])
            
            # adjust generated variant positions for geneome with deletions
            def adjust(pos0):
                offset = 0
                for s,e in ranges:
                    if pos0 > e:
                        offset += (e-s)
                return(pos0 - offset)
            
            # adjust the second half of the generated variants
            SNPs_per_genome_adjusted = self.SNPs_per_genome[:self.num_individuals]
            for SNPs in self.SNPs_per_genome[self.num_individuals:]:
                adjusted = []
                for pos0,variant in SNPs:
                    adjusted += [(adjust(pos0),variant)]
                SNPs_per_genome_adjusted += [adjusted]
            
            indel_dict_by_pos_pergenome_adjusted = self.indel_dict_by_pos_pergenome[:self.num_individuals]
            for indels in self.indel_dict_by_pos_pergenome[self.num_individuals:]:
                adjusted = {}
                for pos0,indel in indels.items():
                    adjusted[adjust(pos0)] = indel
                
                indel_dict_by_pos_pergenome_adjusted += [adjusted]
            
            save_these['SNPs_adjusted'] = SNPs_per_genome_adjusted
            save_these['InDels_adjusted'] = indel_dict_by_pos_pergenome_adjusted
            SNPs_per_genome = SNPs_per_genome_adjusted
            indel_dict_by_pos_pergenome = indel_dict_by_pos_pergenome_adjusted
        else:
            SNPs_per_genome = self.SNPs_per_genome
            indel_dict_by_pos_pergenome = self.indel_dict_by_pos_pergenome
            

        # adjusted are needed to apply the variants
        # unadjusted are needed to check the calling
        _cPickle.dump(save_these, open('baga.GemSIM_known_variants.p','w'))
        # save_these = cPickle.load(open('baga.GemSIM_known_variants.p','r'))


        ### generate genotypes (apply variants) ###
        genotypes = []
        for gn,SNPs in enumerate(SNPs_per_genome):
            if len(self.large_deletions) and gn >= self.num_individuals:
                # use genome with large deletions for second batch
                orig_genome = genome_large_deletions
                genome = _array('c',genome_large_deletions)
            else:
                # else use original
                orig_genome = self.genome.sequence
                genome = _array('c',self.genome.sequence)
            
            # first SNPs
            for pos0,SNP in SNPs:
                assert genome[pos0] != SNP
                genome[pos0] = SNP
            
            # check it worked
            changed = [pos0 for pos0,(new,old) in enumerate(zip(genome,list(orig_genome))) if new != old]
            assert changed == [pos0 for pos0,var in SNPs], 'SNPs were not correctly applied . . .'
            # then indels
            newgenome = _array('c')
            last_pos0 = 0
            for pos0,indel in sorted(indel_dict_by_pos_pergenome[gn].items()):
              if isinstance(indel,str):
                # insertion
                newgenome.extend(genome[last_pos0:pos0])
                newgenome.extend(indel)
                last_pos0 = pos0
              else:
                # deletion
                newgenome.extend(genome[last_pos0:pos0])
                last_pos0 = pos0 + indel
            
            newgenome.extend(genome[last_pos0:])
            genome_seqrecord = _SeqRecord(_Seq(newgenome.tostring()), 
                    id = self.genome.id+'_sim{:02d}'.format(gn+1), name = '', description = '')
            genotypes += [genome_seqrecord]
            print(len(self.genome.sequence),len(genotypes[-1]),genotypes[-1].id)

        self.genotypes = genotypes

    def writeSequences(self, schema = 'fasta'):
        '''
        Write out fasta files with generated variants applied to the reference sequence.
        '''

        written_genomes = []
        for i,genotype in enumerate(self.genotypes):
            genome_out = '{}_for_GemSim_{:02d}.fasta'.format(self.genome.id, i+1)
            print('Writing {}'.format(genome_out))
            _SeqIO.write(genotype,genome_out,schema)
            written_genomes += [genome_out]

        self.written_genomes = written_genomes

    def generateReads(self, path_to_exe = False, 
                            paths_to_genomes = False,
                            readcov = 60,
                            readlen = 100,
                            fraglen = 350,
                            sterrfraglen = 20,
                            model = 4,
                            max_cpus = -1):
        '''
        Call GemSIM to generate reads

        Need to have written genome sequences to generate from, possibly with 
        generated SNPs, small indels and large deletions.
        '''

        #max_cpus etc

        if paths_to_genomes:
            use_genomes = sorted(paths_to_genomes)
        elif hasattr(self, 'written_genomes'):
            use_genomes = sorted(self.written_genomes)
        else:
            raise ValueError('provide either paths_to_genomes or generate some then .writeSequences()')

        if not path_to_exe:
            path_to_exe = _get_exe_path('gemsim')

        comment2 = '''
        to generate reads put GemSIM v1.6 into subfolder GemSIM_v1.6 and issue these commands:
        GemSIM_v1.6/GemReads.py -r LESB58_for_GemSim_01.fasta -n 1980527 -l d -u 350 -s 20 -m GemSIM_v1.6/models/ill100v4_p.gzip -c -q 33 -p -o GemSimLESB58_01
        '''

        num_pairs = len(self.genome.sequence) * readcov / (readlen*2)

        if model == 4:
            path_to_model = _os.path.sep.join(path_to_exe.split(_os.path.sep)[:-1] + ['models','ill100v4_p.gzip'])
        elif model == 5:
            path_to_model = _os.path.sep.join(path_to_exe.split(_os.path.sep)[:-1] + ['models','ill100v5_p.gzip'])

        print('Using error model: {}'.format(path_to_model))
        print('Generating {:,} {}bp read pairs for {}x coverage depth of a {}bp genome ({})'.format(
                num_pairs, readlen, readcov, len(self.genome.sequence), self.genome.id))

        processes = set()
        max_processes = _decide_max_processes( max_cpus )

        import time
        start = time.time()
        out_raw = []
        for i,genome_in in enumerate(use_genomes):
            # could use per genome length . . less consistent than using reference
            # genome_len = len(_SeqIO.read(genome_in,'fasta').seq)
            # num_pairs = genome_len * readcov / (readlen*2)
            outprefix = 'GemSim_{}_{:02d}'.format(self.genome.id, i+1)
            cmd = [path_to_exe, 
                        '-r', genome_in,
                        '-n', num_pairs, 
                        '-l', 'd', '-u', fraglen, '-s', sterrfraglen, 
                        '-m', path_to_model, 
                        '-c', 
                        '-q', 33, '-p',
                        '-o', outprefix]
            out_raw += [outprefix+'_fir.fastq', outprefix+'_sec.fastq']
            # this would be better to rename and compress all in one
            # maybe as a shell script? Then resuming (--force) would be easier.
            if _os.path.exists(outprefix+'_fir.fastq') and \
                    _os.path.exists(outprefix+'_sec.fastq'):
                print('Found output for {}_fir.fastq (and sec), not regenerating, '\
                'delete these to start from scratch'.format(outprefix))
            else:
                cmd = map(str,cmd)
                print(' '.join(cmd))
                processes.add( _subprocess.Popen(cmd, shell=False) )
            if len(processes) >= max_processes:
                (pid, exit_status) = _os.wait()
                processes.difference_update(
                    [p for p in processes if p.poll() is not None])
            

        # Check if all the child processes were closed
        for p in processes:
            if p.poll() is None:
                p.wait()

        missing = []
        for o in out_raw:
            if not _os.path.exists(o):
                missing += [o]

        assert len(missing) == 0, 'Could not find:\n{}'.format('\n'.join(missing))
        print('all finished after {} minutes'.format(int(round((time.time() - start)/60.0))))

        outdir = _os.path.sep.join(['simulated_reads',self.genome.id])
        try:
            _os.makedirs(outdir)
        except OSError:
            pass

        for o in range(out_raw):
            new = _os.path.sep.join([outdir, o.replace('fir','R1').replace('sec','R2')])
            print('{} ==> {}'.format(o, new))
            _os.rename(o, new)
            cmd = ['gzip', new]
            print(' '.join(cmd))
            _subprocess.call(cmd)

    def do(self, num_SNPs = 0, num_deletions = 0, num_insertions = 0):


        if num_SNPs:
            self.generateSNPs(num_SNPs)

        if num_deletions or num_insertions:
            self.generateInDels(num_deletions, num_insertions)

        self.generateSequences()
        self.writeSequences()

if __name__ == '__main__':
    main()
