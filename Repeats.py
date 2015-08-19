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
Repeats module from the Bacterial and Archaeal Genome (BAG) Analyzer.

This module contains functions to detect repetitive regions in a genome 
sequence and plot them. The intended use is to mark the repetitive regions for 
exclusion from a short read mapping experiment. Repetitive regions are more 
likely to contain ambiguous variants that may be caused by divergence between 
duplicated regions in the reference genome or sample and not by variation at 
orthologous regions.
'''

# stdlib
from baga import _subprocess
from baga import _os
from baga import _cPickle
from baga import _gzip
from baga import _tarfile
from baga import _json
from baga import _StringIO

from collections import defaultdict as _defaultdict
from cStringIO import StringIO as _StringIO

# external Python modules
import pysam as _pysam
from Bio import SeqIO as _SeqIO
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord

# svgwrite imported within Plotter
def main():
    pass
class Finder:
    '''
    The Finder class of the Repeats module contains the methods to infer 
    sequence repeats up to many ORFs long in a genome sequence.
    '''
    def __init__(self, genome):
        '''
        A repeat finder object must be instantiated with a genome prepared by the CollectData option
        '''
        
        self.genome = genome
        self.ORFs_ordered = sorted(genome.ORF_ranges, key = genome.ORF_ranges.get)

    def parseBamForDoubleHits(self, filename):
        '''
        Parse a BWA generated BAM file for ORF queries that align to elsewhere
        than themselves
        '''
        alns = _pysam.Samfile(filename)
        ORF_hit_ranges = _defaultdict(dict)
        non_self_non_secondary = _defaultdict(list)
        for aln in alns:
            # get ORF that was hit
            for ORF,(s,e,strand,gene_name) in self.genome.ORF_ranges.items():
                if s < aln.reference_end and aln.reference_start < e:
                    # exclude self
                    if ORF != aln.query_name:
                        s_q, e_q, strand, gene_name = self.genome.ORF_ranges[aln.query_name]
                        # check how much of ORF was aligned,
                        # only accept if >50% covered
                        aligned_positions_in_hit_ORF = sorted(
                                        set(range(s,e)) & set(aln.positions)
                                                              )
                        if len(aligned_positions_in_hit_ORF) > (e-s) * 0.5 and \
                           len(aligned_positions_in_hit_ORF) > (e_q-s_q) * 0.5:
                            # print(  aln.query_name, 
                                    # e_q-s_q, ORF, e-s, 
                                    # len(aligned_positions_in_hit_ORF), 
                                    # aln.is_secondary
                                    # )
                            ORF_hit_ranges[aln.query_name][ORF] = \
                                aln.reference_start, aln.reference_end
                            if not aln.is_secondary:
                                non_self_non_secondary[aln.query_name] += [ORF]

        self.ORF_hit_ranges = dict(ORF_hit_ranges)

        self.non_self_non_secondary = dict(non_self_non_secondary)

    def followHitsAndAdd(self):
        '''
        make all hits symmetric:
            some hits are one way which is good enough evidence 
            for homology for now
        
        and generate ORF_hits dict from ORF_hit_ranges or update it if it 
        already exists
        '''
        if hasattr(self, 'ORF_hits'):
            ORFdict_in = dict(self.ORF_hits.items())
        else:
            ORFdict_in = dict(self.ORF_hit_ranges.items())

        ORFdict_out = {}
        # first fix one-way asymmetries
        hit_no_hits = set([a for b in ORFdict_in.values() for a in b]) - set(ORFdict_in)
        for query,hits in ORFdict_in.items():
            for hit in hit_no_hits.intersection(hits):
                try:
                    ORFdict_out[hit].add(query)
                except KeyError:
                    ORFdict_out[hit] = set([query])

        # then fix multi-way asymmetries
        for ORF in ORFdict_in:
            to_add = set([a for a,b in ORFdict_in.items() if len(set(b) & set(ORFdict_in[ORF])) > 0])
            to_add.update([a for a,b in ORFdict_in.items() if ORF in b])
            try:
                ORFdict_out[ORF].update(to_add)
            except KeyError:
                ORFdict_out[ORF] = to_add

        self.ORF_hits = ORFdict_out

    def orderORFHits(self):

        all_ORF_hits = set([a for b in self.ORF_hits.values() for a in b])

        self.ORFs_with_hits_ordered = sorted(all_ORF_hits, key = self.genome.ORF_ranges.get)
    def getAdjacentORFsWithHit(self, ORF, maxdist = 3, direction = 1, getall = False):
        '''return next ORF with a paralog and num ORFs to it within maxdist'''
        thisORFn_hits = self.ORFs_with_hits_ordered.index(ORF)
        thisORFn_all = self.ORFs_ordered.index(ORF)
        if direction == -1:
            s = thisORFn_all + maxdist * direction
            e = thisORFn_all
        else:
            s = thisORFn_all
            e = thisORFn_all + maxdist * direction

        ORFS_to_check = self.ORFs_ordered[s:e][::direction]
        next_hit_ORFs = {}
        for n in range(1,maxdist):
            if thisORFn_hits + n * direction == len(self.ORFs_with_hits_ordered):
                continue
            next_hit_ORF = self.ORFs_with_hits_ordered[thisORFn_hits + n * direction]
            try:
                next_hit_ORF_index = ORFS_to_check.index(next_hit_ORF)
                next_hit_ORFs[next_hit_ORF] = next_hit_ORF_index
            except ValueError:
                pass

        if len(next_hit_ORFs) == 0:
            if getall:
                return([])
            else:
                return(False,False)
            
        elif getall:
            return(next_hit_ORFs)
        else:
            return(sorted(next_hit_ORFs.items(), key = lambda x: x[1])[0])
    def getHomologousContiguousBlocks(self):
        blocks = []
        merging = False
        block_to_homoblocks = {}  ## need to add singletons
        for n,thisORF in enumerate(self.ORFs_with_hits_ordered):
            #n += 1
            #thisORF = ORFs_with_hits_ordered[n]
            #if thisORF in ('PLES_29041', 'PLES_29011', 'PLES_28981'): break
            #if thisORF == 'PLES_58641': break
            print('\n>>Checking %s (%s) for adjacent ORFs with homologs within the genome' % (thisORF,n))
            if n + 1 == len(self.ORFs_with_hits_ordered):
                thisORFnearestWithHit, dist_to_next = False, False
            else:
                thisORFnearestWithHit, dist_to_next = self.getAdjacentORFsWithHit(thisORF)
            
            if dist_to_next is False:
                # print('None found (within maxdist)')
                # so end block
                if not merging:
                    # add this ORF as a single-ORF block as complete
                    blocks += [[thisORF]]
                    # no homologous contiguous blocks found? contents of ORF_hit_ranges[thisORF] must be singletons
                    block_to_homoblocks[tuple([thisORF])] = sorted([[a] for a in self.ORF_hits[thisORF]])
                    # leadORF_to_homoblocks[(thisORF,)] = sorted(set(homoblocks))
                else:
                    # add this multi-ORF block as complete
                    blocks += [this_block]
                    # also save homologous contiguous blocks
                    # remove contiguous blocks which lost contiguity before others (must be at least one other same length as this_block though)
                    block_to_homoblocks[tuple(this_block)] = [b for b in sorted(map(sorted, homoblocks.values())) if len(b) == len(this_block)]
                    if len(block_to_homoblocks[tuple(this_block)]) == 1:
                        print('ERROR: no contiguous homologous blocks retained (except self) after removing short ones')
                    
                    # halt merging
                    merging = False
                
            # check for tandem repeats
            elif not merging and thisORFnearestWithHit in self.ORF_hits[thisORF]:
                # single ORF tandem repeat: NOT same block (two adjacent blocks probably of length 1)
                # print('tandem warning', thisORF, thisORFnearestWithHit, self.ORF_hits[thisORF])
                # due to tandem repeat, break blocks appropriately
                # add this ORF as a single-ORF block as complete
                blocks += [[thisORF]]
                # also save homologous contiguous blocks
                block_to_homoblocks[tuple(this_block)] = sorted(map(sorted, homoblocks.values()))
                
            elif merging and thisORFnearestWithHit in [a for b in map(self.ORF_hits.get, this_block) for a in b]:
                ## multi-ORF tandem repeat
                # this next adjacent is a homolog to an ORF in current contiguous block, so block should be terminated
                # add this multi-ORF block as complete
                blocks += [this_block]
                # remove contiguous blocks which lost contiguity before others (must be at least one other same length as this_block though)
                block_to_homoblocks[tuple(this_block)] = [b for b in sorted(map(sorted, homoblocks.values())) if len(b) == len(this_block)]
                # halt merging
                if len(block_to_homoblocks[tuple(this_block)]) == 1:
                    print('ERROR: no contiguous homologous blocks retained (except self) after removing short ones')
                
                merging = False
                
            else:
                # try:
                    # print('%s (%s, %s) is adjacent to %s which has a paralog to check' % (thisORF, ORFgene[locus_tag2ORFfeatID[thisORF]], ORFproduct[locus_tag2ORFfeatID[thisORF]], thisORFnearestWithHit))
                # except KeyError:
                    # if thisORF in locus_tag2ORFfeatID:
                        # print('%s (%s) is adjacent to %s which has a paralog to check' % (thisORF, ORFproduct[locus_tag2ORFfeatID[thisORF]], thisORFnearestWithHit))
                    # else:
                        # print('%s is adjacent to %s which has a paralog to check' % (thisORF, thisORFnearestWithHit))
                
                # did find an adjacent ORF within maxdist with a homolog elsewhere (and not self i.e., not tandem repeat): are homologs part of a homologous block?
                print('%s hits (para)logs %s' % (thisORF, ', '.join(self.ORF_hits[thisORF])))
                
                if not merging:
                    # assuming these share homology to adjacent blocks elsewhere . . .
                    # start a new block because we weren't merging
                    this_block = [thisORF, thisORFnearestWithHit]
                    homoblocks = {}
                    merging = True
                else:
                    # continue existing block because we were merging
                    this_block += [thisORFnearestWithHit]
                
                thisORFstrand = self.genome.ORF_ranges[thisORF][2]
                
                # do both have paralogs within maxdist of each other? (and is strand/orientation conserved?)
                for thisORFhit in self.ORF_hits[thisORF]:
                    # this time collect all ORFs within maxdist with homologs (not just the nearest)
                    thisORFstrandhit = self.genome.ORF_ranges[thisORFhit][2]
                    if n + 1 == len(self.ORFs_with_hits_ordered):
                        thisORFhit_nearestWithHits = {}
                    elif thisORFstrand == thisORFstrandhit:
                        thisORFhit_nearestWithHits = self.getAdjacentORFsWithHit(thisORFhit, direction = 1, getall = True)
                    else:
                        thisORFhit_nearestWithHits = self.getAdjacentORFsWithHit(thisORFhit, direction = -1, getall = True)
                    
                    if len(thisORFhit_nearestWithHits) > 0:
                        # print('checking this ORF homolog* (%s) for adjacents** also with homologs' % thisORFhit)
                        got = False
                        # got adjacents of one of this ORF's paralogs, now check if this ORF's adjacent hit any of these ORF's paralog adjacents
                        # i.e., same continguous block
                        for thisORFnearestWithHit_hit in self.ORF_hits[thisORFnearestWithHit]:
                            if thisORFnearestWithHit_hit in thisORFhit_nearestWithHits:
                                dist_to_thisORF_nearestWithHit = thisORFhit_nearestWithHits[thisORFnearestWithHit_hit]
                                # A = 'Focus:    %s  has %s   at %s dist\n' % (thisORF, thisORFnearestWithHit, dist_to_next)
                                # B = 'Paralogs: %s* has %s** at %s dist\n' % (thisORFhit, thisORFnearestWithHit_hit, dist_to_thisORF_nearestWithHit)
                                # print(A+B)
                                got = True
                                if thisORFhit in homoblocks:
                                    homoblocks[thisORFhit] += [thisORFnearestWithHit_hit]
                                    homoblocks[thisORFnearestWithHit_hit] = homoblocks[thisORFhit]
                                    del homoblocks[thisORFhit]
                                else:
                                    homoblocks[thisORFnearestWithHit_hit] = [thisORFhit, thisORFnearestWithHit_hit]
                        
                        if not got:
                            # print('Nothing for %s' % thisORFhit)
                            pass
                        
                    else:
                        # print('This ORFs homolog (%s) has no adjacents with hits (within maxdist)' % thisORFhit)
                        pass

        # add tandem singletons in both directions (only collected in one direction above)
        tandems = []
        for a,b in block_to_homoblocks.items():
            if len(b) > 1:
                if self.ORFs_ordered.index(b[0][0]) + 1 == self.ORFs_ordered.index(b[1][0]):
                    print(b)
                    tandems += [b]

        for t in tandems:
            for o in t:
                if tuple(o) in block_to_homoblocks:
                    print(block_to_homoblocks[tuple(o)])
                else:
                    print(tuple(o))
                    block_to_homoblocks[tuple(o)] = t

        homologous_groups = set()
        for this,those in sorted(block_to_homoblocks.items()):
            homologous_groups.add(tuple(sorted(map(tuple,those))))
            if len(those) > 0:
                for these in those:
                    if this != tuple(these):
                        print('%s => %s\n' % ('+'.join(this), '+'.join(these)))

        self.homologous_groups = sorted(homologous_groups)

        self.block_to_homoblocks = block_to_homoblocks
    def NeedlemanWunch_align(self, sA_info, sB_info, \
                             model = 'affine:global', exhaustive = 'no', protein = 'no', \
                             exe_exonerate = ['external_programs','exonerate-2.2.0-x86_64', 'bin','exonerate']):
        
        '''
        Call exonerate as a subprocess to do a Needleman-Wunch-like gapped pairwaise globl alignment
        http://www.ebi.ac.uk/~guy/exonerate/exonerate.man.html
        '''

        # optionally over-ride default external executable for this method, 
        # if set upstream
        if hasattr(self, '_exe_exonerate'):
            exe_exonerate = _os.path.sep.join(self._exe_exonerate)
        else:
            exe_exonerate = _os.path.sep.join(exe_exonerate)

        (sA, eA, strandA) = sA_info
        (sB, eB, strandB) = sB_info

        if strandA == 1:
            seqA = self.genome_genbank_record.seq[sA:eA]
        else:
            seqA = self.genome_genbank_record.seq[sA:eA].reverse_complement()

        if strandB == 1:
            seqB = self.genome_genbank_record.seq[sB:eB]
        else:
            seqB = self.genome_genbank_record.seq[sB:eB].reverse_complement()

        if protein in (True, 'yes'):
            seqA = seqA.translate()
            seqB = seqB.translate()

        _SeqIO.write(_SeqRecord(seqA, id = 'A'), 'A.fa', 'fasta')
        _SeqIO.write(_SeqRecord(seqB, id = 'B'), 'B.fa', 'fasta')

        cmd = [exe_exonerate, 
               '--model', model, 
               '--verbose', '0', 
               '--exhaustive', exhaustive, 
               #'--refine', 'full',
               # report the single best alignment, however bad it is
               '--showvulgar', 'yes', 
               '--showalignment', 'no', 
               '-n', '1', 
               '-s', '0',
               'A.fa', 'B.fa']

        try:
            aln = _subprocess.check_output(cmd).split('\n')[0]
            if _os.path.exists('A.fa'):
                _os.unlink('A.fa')
            if _os.path.exists('B.fa'):
                _os.unlink('B.fa')

        except _subprocess.CalledProcessError as e:
            print("Exonerate pairwise alignment of duplicated regions failed:\n"), e.output
            aln = ''

        return(aln)
    def getAlnFromVULGAR(self, sA_info, sB_info, vulgar, protein = 'yes'):
        '''convert VULGAR output of exonerate to aligned sequences'''

        (sA, eA, strandA) = sA_info
        (sB, eB, strandB) = sB_info

        bits = vulgar.rstrip().split(' ')[1:]
        Aid, sAaln, eAaln, Astrnd_aln, Bid, sBaln, eBaln, Bstrnd_aln, score = bits[:9]
        sAaln, eAaln, sBaln, eBaln = map(int, [sAaln, eAaln, sBaln, eBaln])
        alninfo = bits[9:]

        if '-' in (Astrnd_aln, Bstrnd_aln):
            # not alignable if exonerate had to resort to reverse strand
            return('X','X')
            
        if strandA == 1:
            seqA = self.genome_genbank_record.seq[sA:eA]
        else:
            seqA = self.genome_genbank_record.seq[sA:eA].reverse_complement()

        if strandB == 1:
            seqB = self.genome_genbank_record.seq[sB:eB]
        else:
            seqB = self.genome_genbank_record.seq[sB:eB].reverse_complement()

        if protein in (True, 'yes'):
            seqA = seqA.translate()
            seqB = seqB.translate()

        Aalnd = []
        Balnd = []
        Ai = 0
        Bi = 0
        for e in range(0,len(alninfo),3):  #break
            #print(vulgar[e])
            forA = int(alninfo[e+2])
            forB = int(alninfo[e+1])
            if alninfo[e] == 'M' or alninfo[e] == 'C':
                Aalnd += list(seqA[sAaln:eAaln][Ai:Ai+forA])
                Balnd += list(seqB[sBaln:eBaln][Bi:Bi+forB])
                Ai += forA
                Bi += forB
            elif alninfo[e] == 'G' or alninfo[e] == 'F':
                if forB == 0:
                    Balnd += list(seqB[sBaln:eBaln][Bi:Bi+forA])
                    Aalnd += ['-'] * forA
                    Bi += forA
                else:
                    Aalnd += list(seqA[sAaln:eAaln][Ai:Ai+forB])
                    Balnd += ['-'] * forB
                    Ai += forB

        return(''.join(Aalnd), ''.join(Balnd))
    def NeedlemanWunch_seqalign(self, sA_info, sB_info, protein = 'no'):
        
        '''
        Call seq-align as a subprocess to do a Needleman-Wunch gapped pairwise globl alignment
        '''
        (sA, eA, strandA) = sA_info
        (sB, eB, strandB) = sB_info

        if strandA == 1:
            #seqA = self.genome_genbank_record.seq[sA:eA]
            seqA = _Seq(self.genome.sequence[sA:eA].tostring())
        else:
            #seqA = self.genome_genbank_record.seq[sA:eA].reverse_complement()
            seqA = _Seq(self.genome.sequence[sA:eA].tostring()).reverse_complement()

        if strandB == 1:
            #seqB = self.genome_genbank_record.seq[sB:eB]
            seqB = _Seq(self.genome.sequence[sB:eB].tostring())
        else:
            #seqB = self.genome_genbank_record.seq[sB:eB].reverse_complement()
            seqB = _Seq(self.genome.sequence[sB:eB].tostring()).reverse_complement()

        if protein in (True, 'yes'):
            seqA = seqA.translate()
            seqB = seqB.translate()

        # seqA = _Seq('cgcacgatttgtagtcccgtgtctgcatggacatagcca')
        # seqB = _Seq('aacacaaacaacacaccgggcctagtagagagttggcggcgcc')
        # seqA = _Seq('AMKINVFYAEEEESRCSRPFIISPSLVVGLREQRNLKLDSKKASLV')
        # seqB = _Seq('AMKIIKIKNVFYAESRCSRPFIISPSLVVGVQRREQRNLKLDSKKASLV')
        seqrecords = [_SeqRecord(seqA, id = 'A'), _SeqRecord(seqB, id = 'B')]
        out_handle = _StringIO()
        _SeqIO.write(seqrecords, out_handle, 'fasta')
        cmd = [self.exe_aligner]
        cmd += ['--stdin']
        # because expected to deteriorate?
        # cmd += ['--freeendgap']
        proc = _subprocess.Popen(
            cmd, stdout = _subprocess.PIPE,
            stdin = _subprocess.PIPE)

        proc.stdin.write(out_handle.getvalue())
        proc.stdin.close()
        result = proc.stdout.read()
        proc.wait()
        A, B = result.split('\n')[:2]
        return(A, B)

    def make_non_alignment(self, sA_info, sB_info):

        (sA, eA, strandA) = sA_info
        (sB, eB, strandB) = sB_info

        if strandA == 1:
            A = str(self.genome_genbank_record.seq[sA:eA])
        else:
            A = str(self.genome_genbank_record.seq[sA:eA].reverse_complement())

        if strandB == 1:
            B = str(self.genome_genbank_record.seq[sB:eB])
        else:
            B = str(self.genome_genbank_record.seq[sB:eB].reverse_complement())

        A += '-' * (eB - sB)
        B = '-' * (eA - sA) + B

        return(A,B)
    def getORFsInTandemRepeats(self):
        ORFs_in_tandem_repeats = set()
        for groups in self.homologous_groups:
            # if any of these ORFs are involved in any tandem repeats
            tandem = False
            for n,group1 in enumerate(groups[:-1]):
                for group2 in groups[(n+1):]:
                    # allow up to 2 ORFs between tandem repeated contiguous groups of ORFs
                    if self.ORFs_ordered.index(group1[-1]) + 1 == self.ORFs_ordered.index(group2[0]) \
                    or self.ORFs_ordered.index(group2[-1]) + 1 == self.ORFs_ordered.index(group1[0]) \
                    or self.ORFs_ordered.index(group1[-1]) + 2 == self.ORFs_ordered.index(group2[0]) \
                    or self.ORFs_ordered.index(group2[-1]) + 2 == self.ORFs_ordered.index(group1[0]) \
                    or self.ORFs_ordered.index(group1[-1]) + 3 == self.ORFs_ordered.index(group2[0]) \
                    or self.ORFs_ordered.index(group2[-1]) + 3 == self.ORFs_ordered.index(group1[0]):
                        ORFs_in_tandem_repeats.update(group1)
                        ORFs_in_tandem_repeats.update(group2)

        self.ORFs_in_tandem_repeats = ORFs_in_tandem_repeats


    def align_blocks(self, extend_len = 200, max_extensions = 10, min_pID = 0.95, num_terminal_window_steps = 5):
        '''
        Align homologous blocks using the Needleman Wunch pairwise global alignment
        
        extension options:
            extend_len: initial extension length
            max_extensions: maximum extension iterations
            min_pID: min percent ID over . . .
            num_terminal_window_steps: . . . this numer of moving window steps, 
        to continue extension
        '''

        unexpecteds = []

        homologous_groups_alnd = []
        for g_n,groups in enumerate(self.homologous_groups):
            alignment_combos = {}
            for n,group1_forward in enumerate(groups[:-1]):
                strandA = self.genome.ORF_ranges[group1_forward[0]][2]
                group1 = group1_forward[::strandA]
                for group2_forward in groups[(n+1):]:   # break
                    strandB = self.genome.ORF_ranges[group2_forward[0]][2]
                    group2 = group2_forward[::strandB]
                    #print(group1,group2)
                    alndA = []
                    alndB = []
                    delimitersA = []
                    delimitersB = []
                    # same strand
                    # pre alignment: before first ORF
                    if strandA == 1:
                        sA = self.genome.ORF_ranges[group1[0]][0] - extend_len
                        eA = self.genome.ORF_ranges[group1[0]][0]
                    else:
                        sA = self.genome.ORF_ranges[group1[0]][1]
                        eA = self.genome.ORF_ranges[group1[0]][1] + extend_len
                    
                    if strandB == 1:
                        sB = self.genome.ORF_ranges[group2[0]][0] - extend_len
                        eB = self.genome.ORF_ranges[group2[0]][0]
                    else:
                        sB = self.genome.ORF_ranges[group2[0]][1]
                        eB = self.genome.ORF_ranges[group2[0]][1] + extend_len
                    
                    delimitersA += [(sA,eA)]
                    delimitersB += [(sB,eB)]
                    #aln_raw = self.NeedlemanWunch_align((sA,eA,strandA), (sB,eB,strandB), model = 'affine:global', exhaustive = 'yes', protein = 'no')
                    A,B = self.NeedlemanWunch_seqalign((sA,eA,strandA), (sB,eB,strandB), protein = 'no')
                    
                    # if 'vulgar' in aln_raw:
                        # A,B = self.getAlnFromVULGAR((sA,eA,strandA), (sB,eB,strandB), aln_raw, protein = 'no')
                        # if (A, B) == ('X','X'):
                            # # no alignment here, make alignment with A and B matching all gaps for percent identity i.e. zero
                            # A, B = self.make_non_alignment((sA, eA, strandA), (sB, eB, strandB))
                        # check for near-100% identity at beginning and extend further as necessary
                        #elif len(self.ORFs_in_tandem_repeats.intersection(group1 + group2)) == 0:
                    if len(self.ORFs_in_tandem_repeats.intersection(group1 + group2)) == 0:
                            # only if ORFs not in known tandem repeats
                        pIDs = self.get_percent_ID(A, B, window = 100, step = 20)
                        # if beginning of extended bit is high ID, extend and test again until pID drops off
                        for extend_multiplier in range(2,max_extensions):  #break
                            if sum([pID for pos,pID in pIDs][:num_terminal_window_steps])/num_terminal_window_steps > min_pID:
                                # print('extending beginning of contiguity: %s bp' % (extend_multiplier * extend_len))
                                # on reverse strand its like extending at the end WRT getting the target sequence out of the chromosome sequence
                                if strandA == 1:
                                    sA = self.genome.ORF_ranges[group1[0]][0] - extend_len * extend_multiplier
                                else:
                                    eA = self.genome.ORF_ranges[group1[0]][1] + extend_len * extend_multiplier
                                
                                if strandB == 1:
                                    sB = self.genome.ORF_ranges[group2[0]][0] - extend_len * extend_multiplier
                                else:
                                    eB = self.genome.ORF_ranges[group2[0]][1] + extend_len * extend_multiplier
                                
                                A ,B = self.NeedlemanWunch_seqalign((sA,eA,strandA), (sB,eB,strandB), protein = 'no')
                                    #aln_raw = self.NeedlemanWunch_align((sA,eA,strandA), (sB,eB,strandB), model = 'affine:global', exhaustive = 'yes', protein = 'no')
                                    # if 'vulgar' in aln_raw:
                                        # A, B = self.getAlnFromVULGAR((sA, eA, strandA), (sB, eB, strandB), aln_raw, protein = 'no')
                                        # if (A, B) == ('X','X'):
                                            # # no alignment here, make alignment with A and B matching all gaps for percent identity i.e. zero
                                            # A, B = self.make_non_alignment((sA, eA, strandA), (sB, eB, strandB))
                                pIDs = self.get_percent_ID(A, B, window = 100, step = 20)
                                    # else:
                                        # # unexpected: previously alignable but not after extension . . .
                                        # # print('no alignment!?')
                                        # unexpecteds += [g_n]
                                        # break
                            else:
                                # divergence at end of this contiguous block is wider than set by min_pID
                                # so end of block reached: stop extension
                                # update extended iterators, include the final extension as pID fall-off is somewhere within
                                delimitersA[-1] = (sA,eA)
                                delimitersB[-1] = (sB,eB)
                                # print('end of extension at start of contiguities')
                                break
                    #else:
                        # no alignment here, make alignment with A and B matching all gaps for percent identity i.e. zero
                     #   A, B = self.make_non_alignment((sA, eA, strandA), (sB, eB, strandB))
                    
                    ## end of pre-ORF alignment and extension
                    
                    ## now align each ORF and inter-ORF gap
                    
                    alndA += [A]
                    alndB += [B]
                    for o in range(len(group1)):
                        # first align ORFs as coding
                        sA = self.genome.ORF_ranges[group1[o]][0]
                        eA = self.genome.ORF_ranges[group1[o]][1]
                        sB = self.genome.ORF_ranges[group2[o]][0]
                        eB = self.genome.ORF_ranges[group2[o]][1]
                        
                        delimitersA += [(sA,eA)]
                        delimitersB += [(sB,eB)]
                        
                        #aln_raw = self.NeedlemanWunch_align((sA,eA,strandA), (sB,eB,strandB), model = 'affine:global', exhaustive = 'yes', protein = 'yes')
                        #A,B = self.getAlnFromVULGAR((sA, eA, strandA), (sB, eB, strandB), aln_raw, protein = 'yes')
                        A, B = self.NeedlemanWunch_seqalign((sA,eA,strandA), (sB,eB,strandB), protein = 'yes')
                        
                        # collect unaligned nucleotides
                        #Anuc = self.genome_genbank_record[sA:eA]
                        #Anuc.id = 'A'
                        Anuc = _SeqRecord(_Seq(self.genome.sequence[sA:eA].tostring()), id = 'A')
                        if strandA == -1:
                            Anuc.seq = Anuc.seq.reverse_complement()
                        
                        #Bnuc = self.genome_genbank_record[sB:eB]
                        #Bnuc.id = 'B'
                        Bnuc = _SeqRecord(_Seq(self.genome.sequence[sB:eB].tostring()), id = 'B')
                        if strandB == -1:
                            Bnuc.seq = Bnuc.seq.reverse_complement()
                        
                        # put aligned amino acids into a SeqRecord
                        A_aa = _SeqRecord(_Seq(A), id = 'A')
                        B_aa = _SeqRecord(_Seq(B), id = 'B')
                        
                        # align nucleotides as codons to aligned amino acids
                        alnd = self.Nuc2AA([Anuc, Bnuc], [A_aa, B_aa])
                        
                        alndA += [str(alnd['A'].seq)]
                        alndB += [str(alnd['B'].seq)]
                        
                        # then align post-ORF as non-coding to reach next ORF or off the end
                        
                        if o + 1 != len(group1):
                            # if not the last one, just to next ORF
                            if strandA == 1:
                                sA = eA
                                eA = self.genome.ORF_ranges[group1[o + 1]][0]
                            else:
                                eA = sA
                                sA = self.genome.ORF_ranges[group1[o + 1]][1]
                            
                            if strandB == 1:
                                sB = eB
                                eB = self.genome.ORF_ranges[group2[o + 1]][0]
                            else:
                                eB = sB
                                sB = self.genome.ORF_ranges[group2[o + 1]][1]
                            
                            if sA > eA or sB > eB:
                                # ORF overlaps: no inter-ORF to align.
                                # and not the last ORF pair so no extension to do.
                                # Need to trim end of last aligned ORFs so
                                # alignment can map back to chromosome.
                                # Find amount to trim allowing for gaps
                                if sA > eA:
                                    # print('ORF overlap: %s (%s..%s)' % (sA - eA, sA, eA))
                                    ORF_overlap = sA - eA
                                    c = 0
                                    for n,i in enumerate(alndA[-1][::-1]):
                                        #print(n,i)
                                        if i != '-':
                                            c += 1
                                            #print(c)
                                        if c == ORF_overlap:
                                            break
                                    
                                    #print(alndA[-1][-(n+1):])
                                    # trim end of last alignment to remove overlap
                                    alndA[-1] = alndA[-1][:-(n+1)]
                                
                                if sB > eB:
                                    # print('ORF overlap: %s (%s..%s)' % (sB - eB, sB, eB))
                                    ORF_overlap = sB - eB
                                    c = 0
                                    for n,i in enumerate(alndB[-1][::-1]):
                                        #print(n,i)
                                        if i != '-':
                                            c += 1
                                            #print(c)
                                        if c == ORF_overlap:
                                            break
                                    
                                    #print(alndB[-1][-(n+1):])
                                    # trim end of last alignment to remove overlap
                                    alndB[-1] = alndB[-1][:-(n+1)]
                                
                                # Then continue to next ORF pair
                                continue
                                
                            delimitersA += [(sA, eA)]
                            delimitersB += [(sB, eB)]
                            
                            #aln_raw = self.NeedlemanWunch_align((sA,eA,strandA), (sB,eB,strandB), model = 'affine:global', exhaustive = 'yes', protein = 'no')
                            A, B = self.NeedlemanWunch_seqalign((sA,eA,strandA), (sB,eB,strandB), protein = 'no')
                            
                            # if 'vulgar' in aln_raw:
                                # A, B = self.getAlnFromVULGAR((sA, eA, strandA), (sB, eB, strandB), aln_raw, protein = 'no')
                                # if (A, B) == ('X', 'X'):
                                    # # no alignment here, make alignment with A and B matching all gaps for percent identity i.e. zero
                                    # A, B = self.make_non_alignment((sA, eA, strandA), (sB, eB, strandB))
                            # else:
                                # # no alignment here, make alignment with A and B matching all gaps for percent identity i.e. zero
                                # A, B = self.make_non_alignment((sA, eA, strandA), (sB, eB, strandB))
                            
                            alndA += [A]
                            alndB += [B]
                            
                        elif len(self.ORFs_in_tandem_repeats.intersection(group1 + group2)) == 0:
                            # check for near-100% identity at end and extend further as necessary
                            # but not if ORFs in tandem repeats because:
                            # alignments will be to other parts of repeat, not extending contiguities
                            # last one, extend to get to end of duplication
                            # get current (non-extended percent identity)
                            pIDs = self.get_percent_ID(alndA[-1], alndA[-1], window = 100, step = 20)
                            for extend_multiplier in range(1, max_extensions + 1):
                                # if end of extended bit is high ID, extend and test again until pID drops off
                                if sum([pID for pos,pID in pIDs][-num_terminal_window_steps:])/num_terminal_window_steps < min_pID:
                                    # currently the end of the alignment is divergent enough be at end of duplication
                                    # save and continue
                                    if extend_multiplier > 1:
                                        # only if some extending was done
                                        delimitersA += [(sA,eA)]
                                        delimitersB += [(sB,eB)]
                                        alndA += [A]
                                        alndB += [B]
                                        # print('end of extension at ends of contiguities')
                                    else:
                                        # print('no extension at ends of contiguities required')
                                        pass
                                    break
                                else:
                                    # print('extending end of contiguity: %s bp' % (extend_multiplier * extend_len))
                                    # on reverse strand its like extending at the beginning WRT getting the target sequence out of the chromosome sequence
                                    if strandA == 1:
                                        sA = delimitersA[-1][-1]
                                        eA = self.genome.ORF_ranges[group1[o]][1] + extend_len * extend_multiplier
                                    else:
                                        eA = delimitersA[-1][0]
                                        sA = self.genome.ORF_ranges[group1[o]][0] - extend_len * extend_multiplier
                                    
                                    if strandB == 1:
                                        sB = delimitersB[-1][-1]
                                        eB = self.genome.ORF_ranges[group2[o]][1] + extend_len * extend_multiplier
                                    else:
                                        eB = delimitersB[-1][0]
                                        sB = self.genome.ORF_ranges[group2[o]][0] - extend_len * extend_multiplier
                                    
                                    A,B = self.NeedlemanWunch_seqalign((sA,eA,strandA), (sB,eB,strandB), protein = 'no')
                                    #aln_raw = self.NeedlemanWunch_align((sA,eA,strandA), (sB,eB,strandB), model = 'affine:global', exhaustive = 'yes', protein = 'no')
                                    pIDs = self.get_percent_ID(A, B, window = 100, step = 20)
                                    # if 'vulgar' in aln_raw:
                                        # A, B = self.getAlnFromVULGAR((sA, eA, strandA), (sB, eB, strandB), aln_raw, protein = 'no')
                                        # if (A, B) == ('X','X'):
                                            # # no alignment here, make alignment with A and B matching all gaps for percent identity i.e. zero
                                            # A, B = self.make_non_alignment((sA, eA, strandA), (sB, eB, strandB))
                                        # # get pIDs for last extension
                                        # pIDs = self.get_percent_ID(A, B, window = 100, step = 20)
                                    # else:
                                        # if extend_multiplier == 1:
                                            # # first extension did not align: no homology immediatley after ORFs
                                            # A, B = self.make_non_alignment((sA, eA, strandA), (sB, eB, strandB))
                                            # pIDs = self.get_percent_ID(A, B, window = 100, step = 20)
                                        # else:
                                            # # unexpected: previously alignable but not after further extension . . .
                                            # # print('no alignment!?')
                                            # unexpecteds += [g_n]
                                            # break
                                
                                if extend_multiplier == max_extensions:
                                    # should have broken by now
                                    # print('WARNING: extended %s x %s bp but regions not yet divergent enough' % (extend_multiplier, extend_len))
                                    delimitersA += [(sA,eA)]
                                    delimitersB += [(sB,eB)]
                                    alndA += [A]
                                    alndB += [B]
                    
                    # keep sequences in 5-3 orientation as aligned
                    A = ''.join(alndA)
                    B = ''.join(alndB)
                    
                    # but store ordinates reversed for -ve strand
                    if strandA == 1:
                        delimitersA = delimitersA[0][0], delimitersA[-1][-1]
                    else:
                        # from start last aligned ORF pair member (with extension) to end of first aligned pair member
                        delimitersA = delimitersA[-1][0], delimitersA[0][-1]
                    
                    if strandB == 1:
                        delimitersB = delimitersB[0][0], delimitersB[-1][-1]
                    else:
                        # from start last aligned ORF pair member (with extension) to end of first aligned pair member
                        delimitersB = delimitersB[-1][0], delimitersB[0][-1]
                    
                    # delimiter relative position determine strand
                    alignment_combos[group1,group2] = ((delimitersA, A),(delimitersB, B))
            
            homologous_groups_alnd += [alignment_combos]

        self.homologous_groups_alnd = homologous_groups_alnd
    def get_percent_ID(self, A, B, window = 100, step = 20):
        pID_per_window = []
        for i in range(0, len(A)-window, step):
            Achunk = A[i:i+window]
            Bchunk = B[i:i+window]
            pID = sum([a == b for a,b in zip(Achunk,Bchunk)])/100.0
            pID_per_window += [(i+window/2,pID)]

        return(pID_per_window)
    def Nuc2AA(self, unalignedNuc, alignedAA, remove_stops = False):
        '''Given unaligned nucleotides and aligned amino acids, align the nucleotides'''
        # ensure sequences are in dictionary form
        if type(unalignedNuc) is list:
            unalignedNuc = dict([(a.id,a) for a in unalignedNuc])
        elif type(unalignedNuc) is not dict:
            print('Need input as list or dict of SeqRecords')
        if type(alignedAA) is not dict:
            alignedAA = dict([(a.id,a) for a in alignedAA])
        elif type(alignedAA) is not dict:
            print('Need input as list or dict of SeqRecords')

        omit = []
        alignedNuc = {}
        for ID in sorted(unalignedNuc):
            if ID not in alignedAA:
                omit += [ID]
                continue
            rec = alignedAA[ID]
            userec = unalignedNuc[ID]
            if len(str(rec.seq).replace('-',''))*3 != len(userec.seq):
                print('sequences differ in length')
                if len(str(rec.seq).replace('-',''))*3 < len(userec.seq):
                    newrec = userec[0:0]
                    for i in range(0,len(userec.seq),3):
                        if str(userec.seq[i:(i+3)].translate()) != '*':
                            newrec += userec.seq[i:(i+3)]
                    if len(str(rec.seq).replace('-',''))*3 == len(newrec.seq):
                        # print('Warning: omitted some stop codons to make sequences equal in length, check this alignment')
                        userec = newrec
                    else:
                        print('Fail: on seq lengths. Returning None')
                        return(None)
            if len(userec.seq) % 3 != 0:
                userec.seq = userec.seq[:-(len(userec.seq) % 3)]
                # print('extra nucleotides removed')
            if remove_stops:
                if str(userec.seq[-3:].translate()) == '*':
                    userec.seq = userec.seq[:-3]
                    # print('stop codon removed')
            indels = [n for n,a in enumerate(str(rec.seq)) if a == '-']
            st = str(userec.seq)
            ns = []
            i = 0
            n = 0
            g = 0
            while indels:
                if n == indels[0]:
                    ns += '---'
                    indels.pop(0)
                else:
                    ns += st[(g*3):(g*3)+3]
                    g += 1
                n += 1
            ns = ''.join(ns)
            ns += st[(g*3):]
            alndrec = _SeqRecord(_Seq(ns, unalignedNuc[ID].seq.alphabet), id=ID)
            alignedNuc[ID] = alndrec

        if len(omit) > 0:
            print('MISSING: %s' % ' + '.join(omit))

        return(alignedNuc)
    def aln_pos0_2_chrm_pos0_pIDs(self, aligned_seq, pIDs, window = 100):
        '''map pIDs to chromosome given aligned sequence, pIDs, and start and end points'''

        strand = aligned_seq['strand']

        if aligned_seq['start'] > aligned_seq['end']:
            print(
            'NOT IMPLEMENTED:\n\
            plotting reverse strand in this way using aln_pos0_2_chrm_pos0()'
            )
            # plot reverse complement strand of genome for easier comparison?
            return(False)

        pIDs_by_chrmpos0 = {}

        if strand == 1:
            # this needs to start at the end if strand == -1
            alnd_chrom_pos0 = aligned_seq['start']
        else:
            alnd_chrom_pos0 = aligned_seq['end'] - 1

        mismatch = False

        for aln_pos0, char in enumerate(aligned_seq['seq_str']):
            if char != '-':
                if strand == 1:
                    #chromchar = self.genome_genbank_record.seq[alnd_chrom_pos0]
                    chromchar = self.genome.sequence[alnd_chrom_pos0]
                else:
                    #chromchar = self.genome_genbank_record.seq[alnd_chrom_pos0:alnd_chrom_pos0+1].reverse_complement()[0]
                    chromchar = _Seq(self.genome.sequence[alnd_chrom_pos0:alnd_chrom_pos0+1].tostring()).reverse_complement()[0]
                
                #print(char, chromchar)
                if char != chromchar:
                    mismatch = True
                    print('x')
                    break
                try:
                    # account for width of window
                    # within which percent ID calculated
                    pIDs_by_chrmpos0[alnd_chrom_pos0 + (window / 2 * strand)] = pIDs[aln_pos0]
                except KeyError:
                    pass
            
                alnd_chrom_pos0 += strand   ############ this needs to decrement if strand == -1

        if mismatch:
            print('WARNING: mismatch detected when assigning percent identity between duplications to chromosome positions . . . :-(')
            print("Reference chromosome used for plotting may not match that used for aligning duplicate region pairs? (else there's a bug)")
            assert mismatch == False

        return(pIDs_by_chrmpos0)

    def map_alignments_to_chromosome(self):

        homologous_groups_mapped = []

        for g_n,group in enumerate(self.homologous_groups_alnd):
            # one pairwise: two homologs
            # three pairwise: three homologs
            these_pairwise = []
            #print(g_n+1, len(self.homologous_groups_alnd))
            for ORFs, ORFs_info in group.items():
                
                ORFsA, ORFsB = ORFs
                (sA, eA), A = ORFs_info[0]
                (sB, eB), B = ORFs_info[1]
                # print(ORFsA, self.genome.ORF_ranges[ORFsA[0]][:2], sA, eA)
                # print(ORFsB, self.genome.ORF_ranges[ORFsB[0]][:2], sB, eB)
                
                aligned_A = {}
                aligned_B = {}
                
                # start and end are base-0 chromosome sequence Python slices
                (aligned_A['start'], aligned_A['end']), aligned_A['seq_str'] = ORFs_info[0]
                (aligned_B['start'], aligned_B['end']), aligned_B['seq_str'] = ORFs_info[1]
                
                # strands of each region
                aligned_A['strand'] = self.genome.ORF_ranges[ORFsA[0]][2]
                aligned_B['strand'] = self.genome.ORF_ranges[ORFsB[0]][2]
                
                # percent identity over aligned region
                pIDs = dict(self.get_percent_ID(aligned_A['seq_str'], aligned_B['seq_str'], window = 100, step = 20))
                
                # map the pIDs to from alinmnet to chromosome
                aligned_A['aln_pos0_2_chrm_pos0_pIDs'] = self.aln_pos0_2_chrm_pos0_pIDs(aligned_A, pIDs, window = 100)
                aligned_B['aln_pos0_2_chrm_pos0_pIDs'] = self.aln_pos0_2_chrm_pos0_pIDs(aligned_B, pIDs, window = 100)
                
                # store
                aligned = {}
                aligned['A'] = aligned_A
                aligned['B'] = aligned_B
                aligned['pIDs'] = pIDs
                
                these_pairwise += [aligned]
            
            homologous_groups_mapped += [these_pairwise]

        self.homologous_groups_mapped = homologous_groups_mapped
    def makeRanges(self, disjoint_consecs):
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
    def identify_ambiguous_regions(self, 
                    minimum_percent_identity = 0.98, 
                    window = 100, step = 20):
        '''
        Identify genome regions with high identity to other parts of the genome
        i.e., repeats, probably due to recent homology.
        window and step must be set as for get_percent_ID() applied previously
        '''
        for homologous_group in self.homologous_groups_mapped:
            for pair in homologous_group:
                for label in ('A','B'):
                    ambiguous_positions_by_step = []
                    strand = pair[label]['strand']
                    for pos,pID in pair[label]['aln_pos0_2_chrm_pos0_pIDs'].items():
                        if pID >= minimum_percent_identity:
                            # divide by step for calculating contiguous regions
                            ambiguous_positions_by_step += [pos/step]
                    
                    ambiguous_ranges = []
                    
                    if len(ambiguous_positions_by_step) > 0:
                        ambiguous_positions_by_step.sort()
                        ambiguous_ranges_by_step = self.makeRanges(ambiguous_positions_by_step)
                        for s,e in ambiguous_ranges_by_step:
                            ambiguous_ranges += [(s * step, e * step)]
                    
                    pair[label]['ambiguous_ranges'] = ambiguous_ranges
    def merge_ambiguous_regions(self, minimum_repeat_length = 400):
        '''
        combine all of the pairwise alignment-specific ambiguous regions genome-wide
        to remove redundancy and make a single list to check variants against. Also
        remove short, probably resolvable repeats determined by minimum_repeat_length.
        '''

        all_positions = set()

        for homologous_group in self.homologous_groups_mapped:
            for pair in homologous_group:
                for label in ('A','B'):
                    for s,e in pair[label]['ambiguous_ranges']:
                        all_positions.update(range(s,e))

        all_ambiguous_ranges = self.makeRanges(sorted(all_positions))

        initial_num_repeats = len(all_ambiguous_ranges)
        all_ambiguous_ranges = [(s,e) for s,e in all_ambiguous_ranges if e - s > minimum_repeat_length]
        s = sum([(e - s) for s,e in all_ambiguous_ranges])
        print('Dropped {} repeats less than {} basepairs; {} remain spanning {:,} basepairs'.format(
                                                                        initial_num_repeats, 
                                                                        minimum_repeat_length, 
                                                                        len(all_ambiguous_ranges),
                                                                        s))

        self.ambiguous_ranges = all_ambiguous_ranges


    def saveLocal(self, name):
        '''
        Save a FinderInfo object (the findings of a Finder instance for a 
        particular genome including information needed for plotting) local 
        baga file. This will be required for plotting any identified repeats.
        'name' will be incorporated into a longer filename with extension: .baga
        '''

        if name:
            fileout = 'baga.Repeats.FinderInfo-{}.baga'.format(name)
        else:
            fileout = 'baga.Repeats.FinderInfo-{}.baga'.format(self.genome.id)

        def add(obj, name):
            io = _StringIO()
            _json.dump(obj, io)
            io.seek(0, _os.SEEK_END)
            length = io.tell()
            io.seek(0)
            thisone = _tarfile.TarInfo(name = name)
            thisone.size = length
            tar.addfile(tarinfo = thisone, fileobj = io)

        with _tarfile.open(fileout, "w:gz") as tar:
            print('Writing to {} . . . '.format(fileout))
            # genome characteristics
            add(self.genome.id, 'genome_id')
            add(len(self.genome.sequence), 'genome_length')
            add(self.genome.ORF_ranges, 'ORF_ranges')
            add(self.genome.large_mobile_element_ranges, 'large_mobile_element_ranges')
            # repeats
            add(self.homologous_groups_mapped, 'homologous_groups_mapped')
            add(self.ambiguous_ranges, 'ambiguous_ranges')
    def findRepeats(self, minimum_percent_identity = 0.95, 
                          minimum_repeat_length = 400,
                          exe_bwa = False, 
                          exe_samtools = False,
                          exe_exonerate = False,
                          local_repeats_path = ['repeats'], 
                          local_genomes_path = ['genome_sequences'], 
                          force = False):
        '''
        Find repeats!
        Initial assignment and alignment of homologous blocks requires >=95% nucleotide identity
        Final selection of ambiguous repeats for filtering selectes regions >=98% identity
        '''

        exe_bwa = False
        exe_samtools = False
        exe_aligner = False

        def get_exe(name):
            from baga.Dependencies import dependencies as _dependencies
            exe = []
            exe += [_dependencies[name]['destination']]
            exe += _dependencies[name]['checker']['arguments']['path']
            exe = _os.path.sep.join(exe)
            return(exe)


        if not exe_bwa:
            self.exe_bwa = get_exe('bwa')
        elif exe_bwa == 'system':
            self.exe_bwa = _dependencies['bwa']['checker']['arguments']['path'][-1]
        else:
            self.exe_bwa = exe_bwa


        if not exe_samtools:
            self.exe_samtools = get_exe('samtools')
        elif exe_samtools == 'system':
            self.exe_samtools = _dependencies['samtools']['checker']['arguments']['path'][-1]
        else:
            self.exe_samtools = 'samtools'
            

        if not exe_aligner:
            self.exe_aligner = get_exe('seq-align')
        elif exe_bwa == 'system':
            self.exe_aligner = _dependencies['seq-align']['checker']['arguments']['path'][-1]
        else:
            self.exe_aligner = 'seq-align'


        local_repeats_path = _os.path.sep.join(local_repeats_path)
        local_genomes_path = _os.path.sep.join(local_genomes_path)

        genome_fna = '{}/{}.fna'.format(local_genomes_path, self.genome.id)
        genome_fna_ORFs = '{}/{}_ORFs.fna'.format(local_genomes_path, self.genome.id)

        try:
            _os.makedirs(local_genomes_path)
        except OSError:
            pass

        try:
            _os.makedirs(local_repeats_path)
        except OSError:
            pass

        if not _os.path.isfile(genome_fna) or _os.path.getsize(genome_fna) == 0:
            print('Writing genome to FASTA')
            _SeqIO.write(_SeqRecord(_Seq(self.genome.sequence.tostring()), id = self.genome.id),
                        genome_fna,
                        'fasta')

        #### align ORFs to chromosome ####
        cmd = [self.exe_bwa, 'index', genome_fna]
        print('Called: {}'.format(' '.join(cmd)))
        try:
            _subprocess.call(cmd)
        except OSError:
            print('Problem running BWA at {}. Please use Dependencies module to install locally or check system path'.format(cmd[0]))

        # collect ORF sequences
        recs = []
        for ID,(s,e,strand,gene_name) in self.genome.ORF_ranges.items():
            if strand == 1:
                recs += [_SeqRecord(_Seq(self.genome.sequence[s:e].tostring()), id = ID)]
            else:
                recs += [_SeqRecord(_Seq(self.genome.sequence[s:e].tostring()).reverse_complement(), id = ID)]

        _SeqIO.write(recs, genome_fna_ORFs, 'fasta')

        ## relaxed alignments using BWA i.e. find modestly divergent ORFs
        cmd = [self.exe_bwa, 'mem',
            # Band width. Essentially, gaps longer than INT will not be found ... also affected by the scoring matrix and the hit length [100]
            '-w', '100',
            # Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is ... avoids unnecessary extension ... [100]
            '-d', '100',
            # Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [1.5]
            '-r', '1.5',
            # Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000]
            '-c', '10000',
            # Matching score. [1]
            '-A', '1',
            # Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4]
            '-B', '2',
            # Gap open penalty. [6]
            '-O', '3',
            # Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1]
            '-E', '1',
            # Clipping penalty. SW extension: best score reaching the end of query is larger than best SW score minus the clipping penalty, clipping not be applied. [5]
            '-L', '5',
            # Penalty for an unpaired read pair .....applied without -P? set to zero here? [9]
            #-U 9
            # Assume the first input query file is interleaved paired-end FASTA/Q.
            #-p
            # Complete read group header line.
            #-R STR
            # Don't output alignment with score lower than INT. This option only affects output. [30]
            '-T', '30',
            # Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
            '-a',
            # Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output.
            # -C
            # Use hard clipping 'H' in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
            # -H
            # Mark shorter split hits as secondary (for Picard compatibility).
            # -M
            # Control the verbose level of the output. [3]
            # -v INT
            genome_fna, 
            genome_fna_ORFs]


        sam_name = 'ORFs_2_{}.sam'.format(self.genome.id)
        print('Called: {}\nPiping to {}'.format(' '.join(cmd), sam_name))
        with open(_os.path.sep.join([local_repeats_path, sam_name]), "wb") as out:
            _subprocess.call(cmd, stdout = out)

        bam_name = sam_name[:-4] + '_unsorted.bam'
        cmd = [self.exe_samtools, 'view', '-bS', _os.path.sep.join([local_repeats_path, sam_name])]
        print('Called: {}\nPiping to {}'.format(' '.join(cmd), bam_name))
        with open(_os.path.sep.join([local_repeats_path, bam_name]), "wb") as out:
            _subprocess.call(cmd, stdout = out)

        cmd = [self.exe_samtools, 'sort', _os.path.sep.join([local_repeats_path, bam_name]), _os.path.sep.join([local_repeats_path, sam_name[:-4]])] ##### where does output end up (which path?)
        print('Called: {}'.format(' '.join(cmd)))
        _subprocess.call(cmd)

        bam_name_sorted = sam_name.replace('.sam', '.bam')
        cmd = [self.exe_samtools, 'index', _os.path.sep.join([local_repeats_path, bam_name_sorted])]
        print('Called: {}'.format(' '.join(cmd)))
        _subprocess.call(cmd)

        f = _os.path.sep.join([local_repeats_path, sam_name])
        print('Removing: {}'.format(f))
        _os.unlink(f)

        f = _os.path.sep.join([local_repeats_path, bam_name])
        print('Removing: {}'.format(f))
        _os.unlink(f)

        #### parse and analyse ORF-chromosome alignment ####
        print('Parsing and analysing the ORF-chromosome alignments . . .')
        self.parseBamForDoubleHits(_os.path.sep.join([local_repeats_path, bam_name_sorted]))
        # this would vary for genomes other than LESB58
        self.followHitsAndAdd()
        max_iterations = 10
        i = 0
        while len(self.ORF_hits) != len(set([a for b in self.ORF_hits.values() for a in b])):
            self.followHitsAndAdd()
            i += 1
            print(i)
            if i == max_iterations:
                # for general use this should be a warning that could be logged but ignored?
                # maybe use 'force' or a new 'ignore'?
                print('Followed too many one-way ORF alignments . . . ')
                raise

        #### get contiguous blocks and their homologs (other blocks) ####
        self.orderORFHits()
        self.getHomologousContiguousBlocks()

        #### align contiguous homologous blocks ####
        self.getORFsInTandemRepeats()
        # min_pID is used to find an initial set of repeatitive and possibly 
        # divergent regions and should be lower than the eventual minimum
        # threshold
        self.align_blocks(min_pID = 0.95)
        self.map_alignments_to_chromosome()
        self.homologous_groups_mapped

        #### identify 98% (default) identical regions
        self.identify_ambiguous_regions(minimum_percent_identity = minimum_percent_identity)
        # repeats of unit length less than the insert size of the paired end 
        # sequenced fragments should be resolvable so omit them
        # (although tandem repeats of total length > insert size
        # may still be a problem . . .)
        self.merge_ambiguous_regions(minimum_repeat_length = minimum_repeat_length)
class Plotter:
    '''
    Plotter class of the Repeats module contains methods to plot
    the repetitive regions found by an instance of the Finder class.
    '''
    def __init__(self, finder_info, plot_output_path, 
                    width_cm = 30, height_cm = 20, 
                    viewbox_width_px = 1800, viewbox_height_px = 1200,
                    plot_width_prop = 0.8, plot_height_prop = 0.8, 
                    white_canvas = True):
        '''
        Plot pairs of aligned homologous chromosome regions with percent 
        identity calculated over a moving window.

        genome: an instance of Genome from CollectData module for which 
        repeats have been inferred using the Finder class of the Repeats
        module.

        plot_width_prop and plot_height_prop: Proportion of whole plot 
        area covered by actual plot, to allow space for labels.
        '''
        
        self.finder_info = finder_info

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
        self.dwg = dwg
        self.rgb = _svgwrite.rgb

    def doPlot(self, pair):
        # plot lower item 'B' (SVG works bottom to top on y-axis)
        # panel = ((this x-axis, num x-axis axes), (this y-axis, num y-axis axes))
        panel = ((1,1),(1,2))

        # plot scale
        self.plot_scale('B', pair, panel, plot_label = True)

        # plot ORFs
        self.plot_ORFs('B', pair, panel)

        # plot large features
        self.plot_LargeFeatures('B', pair, panel)

        # plot the identity to the other duplicates
        self.plot_pIDs('B', pair, 'B', panel)

        # indicate high-identity ambiguous regions
        self.plot_ambiguous_ranges('B', pair, panel)


        # plot upper item 'A' (SVG works bottom to top on y-axis)
        panel = ((1,1),(2,2))

        # plot scale
        self.plot_scale('A', pair, panel, plot_label = False)

        # plot ORFs
        self.plot_ORFs('A', pair, panel)

        # plot large features
        self.plot_LargeFeatures('A', pair, panel)

        # plot the identity to the other duplicates
        self.plot_pIDs('A', pair, 'A', panel)

        # indicate high-identity ambiguous regions
        self.plot_ambiguous_ranges('A', pair, panel)


        self.dwg.save()

    def chrom2plot_x(self, pos_chrom, pair, seq_label):
        '''convert chromosome position to plotting position in canvas'''
        plotlen_chrom = pair[seq_label]['end'] - pair[seq_label]['start']
        pos_plot = (pos_chrom - pair[seq_label]['start']) * ((self.viewbox_width_px * self.plot_width_prop) / plotlen_chrom)
        return(pos_plot)

    def plot_scale(self,    seq_label, 
                            pair, 
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
        if pair[seq_label]['end'] - pair[seq_label]['start'] <= 2000:
            tick_rounding = 100
        elif pair[seq_label]['end'] - pair[seq_label]['start'] <= 20000:
            tick_rounding = 1000
        else:
            tick_rounding = 5000
            
        tick_dist = (pair[seq_label]['end'] - pair[seq_label]['start']) / num_ticks
        if tick_rounding > tick_dist:
            tick_dist = tick_rounding
        else:
            rm = tick_dist % tick_rounding
            if (tick_rounding / 2.0) < rm:
                tick_dist += (tick_rounding - rm)
            else:
                tick_dist -= rm

        plotpositions_chrom = []
        for pos in range(0, self.finder_info['genome_length'], tick_dist):
            if pair[seq_label]['start'] <= pos <= pair[seq_label]['end']:
                plotpositions_chrom += [pos]

        # convert to x-positions for plotting
        plotlen_chrom = pair[seq_label]['end'] - pair[seq_label]['start']

        # plot into full viewbox or other specified panel
        ((this_x_panel, num_x_panels),(this_y_panel, num_y_panels)) = panel

        # calculate plotting area within viewbox
        # currently only single x panel implemented (full width of view box)
        plotstart_x = (self.viewbox_width_px - (self.viewbox_width_px * self.plot_width_prop)) / 2
        plotend_x = plotstart_x + (self.viewbox_width_px * self.plot_width_prop)
        if num_x_panels > 1:
            print('>1 x panel not implemented')
            return(False)

        plottop_y = (self.viewbox_height_px - (self.viewbox_height_px * self.plot_height_prop)) / 2
        plotbottom_y = plottop_y + (self.viewbox_height_px * self.plot_height_prop)

        # this area is just the data i.e., depth line plot
        if num_y_panels > 1:
            # update existing plot area in y direction
            # this uses num_y_panels for overall shape and this_y_panel for target
            y_per_panel = (plotbottom_y - plottop_y) / num_y_panels
            plottop_y = plotbottom_y - y_per_panel * this_y_panel
            plotbottom_y -= y_per_panel * (this_y_panel - 1)

        # get tick positions on x-axis to plot
        plotpositions_x = [self.chrom2plot_x(pos_chrom, pair, seq_label) for pos_chrom in plotpositions_chrom]
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
                                font_size = tick_font_size) #, baseline_shift='-50%')
            
            self.dwg.add(tiplabel)

        if plot_label:
            # label x-axis
            
            xaxislabel = self.dwg.text(
                                "Reference Chromosome Position", 
                                insert = ((plotstart_x + plotend_x)/2, plotbottom_y + tick_len_px * 2.1 * 2),
                                fill='black', font_family=use_fontfamily, 
                                text_anchor = "middle", font_size = font_size) #, baseline_shift='-50%')
            
            self.dwg.add(xaxislabel)


    def plot_ORFs(self, seq_label, 
                        pair, 
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

        plottop_y = (self.viewbox_height_px - (self.viewbox_height_px * self.plot_height_prop)) / 2
        plotbottom_y = plottop_y + (self.viewbox_height_px * self.plot_height_prop)

        # this area is just the data i.e., depth line plot
        if num_y_panels > 1:
            # update existing plot area in y direction
            # this uses num_y_panels for overall shape and this_y_panel for target
            y_per_panel = (plotbottom_y - plottop_y) / num_y_panels
            plottop_y = plotbottom_y - y_per_panel * this_y_panel
            plotbottom_y -= y_per_panel * (this_y_panel - 1)

        ORF_plot_info = []
        for ID,(s,e,d,name) in self.finder_info['ORF_ranges'].items():
            status = False
            if pair[seq_label]['start'] <= s and e < pair[seq_label]['end']:
                plot_x_s = self.chrom2plot_x(s, pair, seq_label) + plotstart_x
                plot_x_e = self.chrom2plot_x(e, pair, seq_label) + plotstart_x
                status = 'complete'
            elif s < pair[seq_label]['start'] < e:
                plot_x_s = self.chrom2plot_x(pair[seq_label]['start'], pair, seq_label) + plotstart_x
                plot_x_e = self.chrom2plot_x(e, pair, seq_label) + plotstart_x
                status = 'left cut'
            elif s < pair[seq_label]['end'] < e:
                plot_x_s = self.chrom2plot_x(s, pair, seq_label) + plotstart_x
                plot_x_e = self.chrom2plot_x(pair[seq_label]['end'], pair, seq_label) + plotstart_x
                status = 'right cut'
            
            if status:
                if len(name):
                    use_name = '{} ({})'.format(ID, name)
                else:
                    use_name = ID
                
                ORF_plot_info += [(plot_x_s, plot_x_e, d, use_name, status)]

        # some additional parameters here for tweaking layout of features
        ##### this determines feature lanes relative to values plot
        ##### see also values_upper_prop and values_upper_lower
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

        # plot ORFs by strand in 3' 5' order so points are above with white border
        #ORF_plot_info = ORF_partial_left + ORFs_to_plot + ORF_partial_right
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
                # 'forward' strand above reverse
                start,end = s,e
                # establish ORF strand lane positions
                #height_above = (plotbottom_y - plottop_y) * 0 ####### should be on top of scale
                if state == 'left cut':
                    commands = ["M %s %s" % (
                                start - point_width, 
                                plotbottom_y - forward_y_offset
                                )]
                else:
                    commands = ["M %s %s" % (
                                start, 
                                plotbottom_y - forward_y_offset
                                )]
                if state == 'right cut':
                    # no point but angled to indicate off-the-edge
                    commands += ["L %s %s" % (
                                end, 
                                plotbottom_y - forward_y_offset
                                )]
                    commands += ["L %s %s" % (
                                end + point_width, 
                                plotbottom_y - (forward_y_offset + feature_thickness * 1)
                                )]
                else:
                    # point
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
                            plotbottom_y - (forward_y_offset + feature_thickness * 0.5)
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
                    y = plotbottom_y - (forward_y_offset + feature_thickness * 0.5)
                    added_label.rotate(-25, center = (x, y))
                
            else:
                start,end = e,s
                #height_above = (plotbottom_y - plottop_y) * 0
                commands = ["M %s %s" % (
                            start, 
                            plotbottom_y - reverse_y_offset
                            )]
                if state == 'left cut':
                    # no point but angled to indicate off-the-edge
                    commands += ["L %s %s" % (
                            end, 
                            plotbottom_y - reverse_y_offset
                            )]
                    commands += ["L %s %s" % (
                            end - point_width, 
                            plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                            )]
                else:
                    # point
                    commands += ["L %s %s" % (
                            end + point_width, 
                            plotbottom_y - reverse_y_offset
                            )]
                    commands += ["L %s %s" % (
                            end, 
                            plotbottom_y - (reverse_y_offset + feature_thickness * 0.5)
                            )]
                    commands += ["L %s %s" % (
                            end + point_width, 
                            plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                            )]
                if state == 'right cut':
                    commands += ["L %s %s z" % (
                            start + point_width, 
                            plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                            )]
                else:
                    commands += ["L %s %s z" % (
                            start, 
                            plotbottom_y - (reverse_y_offset + feature_thickness * 1)
                            )]
                
                self.dwg.add(
                    self.dwg.path(
                        d=' '.join(commands), stroke='white', 
                        stroke_linecap='round', fill = self.rgb(*colour), 
                        stroke_width = 1
                    )
                )
                added_label = self.dwg.add(
                    self.dwg.text(
                        name,
                        insert = (
                            (start+end)/2, 
                            plotbottom_y - (reverse_y_offset + feature_thickness * 0.5)
                        ),
                        fill='white', stroke='black', stroke_width=0.4, 
                        font_family=use_fontfamily, font_weight='bold',
                        text_anchor='middle', font_size = '%dpt' % font_size, baseline_shift='-50%'
                    )
                )
                # this would vary with font and font size
                max_length_for_unrotated_label = self.viewbox_width_px * 1.0/18
                if max(start,end) - min(start,end) < max_length_for_unrotated_label:
                    x = (start+end)/2
                    y = plotbottom_y - (reverse_y_offset + feature_thickness * 0.5)
                    added_label.rotate(-25, center = (x, y))


    def plot_LargeFeatures(self, seq_label, pair, panel = ((1,1),(1,1)),
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

        plottop_y = (self.viewbox_height_px - (self.viewbox_height_px * self.plot_height_prop)) / 2
        plotbottom_y = plottop_y + (self.viewbox_height_px * self.plot_height_prop)

        # this area is just the data i.e., depth line plot
        if num_y_panels > 1:
            # update existing plot area in y direction
            # this uses num_y_panels for overall shape and this_y_panel for target
            y_per_panel = (plotbottom_y - plottop_y) / num_y_panels
            plottop_y = plotbottom_y - y_per_panel * this_y_panel
            plotbottom_y -= y_per_panel * (this_y_panel - 1)

        # Features = []
        # for f in self.finder.genome_genbank_record.features:
            # if f.type == 'misc_feature':
                # #print(f.qualifiers)
                # if 'uplication' in f.qualifiers['note'][0]:
                    # # not plotting annotated duplications: but mention in legend
                    # continue
                # Features += [(f.location.nofuzzy_start + 1, f.location.nofuzzy_end + 1, f.qualifiers['note'][0])]

        # select features to plot and their x plotting ordinates
        Features_to_plot = []
        #for feat_chrom_start,feat_chrom_end,n in Features:
        for name, (feat_chrom_start, feat_chrom_end) in self.finder_info['large_mobile_element_ranges'].items():
            # don't plot unless part of feature is within plotting range
            if feat_chrom_start < pair[seq_label]['end'] and feat_chrom_end > pair[seq_label]['start']:
                s = self.chrom2plot_x(max( pair[seq_label]['start'], feat_chrom_start), pair, seq_label ) + plotstart_x
                e = self.chrom2plot_x(min( pair[seq_label]['end'], feat_chrom_end), pair, seq_label ) + plotstart_x
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

    def plot_pIDs(self, seq_label, 
                        pair, 
                        label, 
                        panel = ((1,1),(1,1)), 
                        font_size = '20pt', 
                        use_fontfamily = 'Nimbus Sans L',
                        values_upper_prop = 0.75,
                        values_lower_prop = 0.4):
        '''given the real position ordinate on chromosome with percent identity to a duplication, 
        plot line. Also plot y-axis and scale. lower values for values_lower_prop, pID ends 
        further down'''

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

        # get positions on chromosome to plot
        plotpositions_chrom = sorted(pair[seq_label]['aln_pos0_2_chrm_pos0_pIDs'])

        # get positions on x-axis to plot
        plotpositions_x = [self.chrom2plot_x(pos_chrom, pair, seq_label) for pos_chrom in plotpositions_chrom]

        # get corresponding depths over this chromosome region
        plotdepths = []
        # normalise so median values are close to middle of each plot
        # for percent identity which will be 100% over most of the plot, no normalisation
        max_plot_depth = max(pair[seq_label]['aln_pos0_2_chrm_pos0_pIDs'].values()) * 1

        for pos0 in plotpositions_chrom:
            d = pair[seq_label]['aln_pos0_2_chrm_pos0_pIDs'][pos0]
            if d >= max_plot_depth:
                use_d = 1
            else:
                use_d = d / max_plot_depth
            
            plotdepths += [use_d]

        # make paths for each isolate plot lane 
        lane_num = 0
        colour = ('60', '60', '60', '%')

        # establish lane plot area (always single for percent ID at duplicate regions)
        plot_spacing = 3
        lane_upper_y = lower_y - (lane_num + 1) * y_per_lane + plot_spacing
        lane_lower_y = lower_y - lane_num * y_per_lane - plot_spacing
        lane_plot_height = lane_lower_y - lane_upper_y
        #print(i, lane_lower_y, lane_lower_y)
        commands = ['M %s %s' % (
                                    plotpositions_x[0] + plotstart_x, 
                                    lane_lower_y
                                    )]

        for n,d in enumerate(plotdepths):
            # because everything is upside down in SVG minus goes up on page
            commands += ['L %s %s' % (
                                        plotpositions_x[n] + plotstart_x,
                                        lane_lower_y - lane_plot_height * d
                                        )]

        commands += ['L %s %s z' % (
                                    plotpositions_x[-1] + plotstart_x, 
                                    lane_lower_y
                                    )]

        plot_path = self.dwg.path(
                            d=' '.join(commands), stroke = self.rgb(*colour), 
                            stroke_linecap='round', stroke_width= '0', 
                            fill = self.rgb(*colour), fill_rule='evenodd'
                            )

        self.dwg.add(plot_path)

        # plot y-axis
        axis_horizontal_offset = 5
        commands = ['M %s %s' % (plotstart_x - axis_horizontal_offset, upper_y), 
                    'L %s %s' % (plotstart_x - axis_horizontal_offset, lower_y)]
        plot_path = self.dwg.path(
                            d=' '.join(commands), stroke = 'black', 
                            stroke_linecap='round', stroke_width= '3',
                            fill = 'none', fill_rule='evenodd'
                            )

        self.dwg.add(plot_path)

        # ticks
        tick_length = 20
        # make tick font size slightly smaller than labels
        font_size_int = int(font_size.replace('pt',''))
        tick_font_size = '{}pt'.format(font_size_int * 0.85)

        for n in range(6):
            commands = ['M %s %s' % (
                                        plotstart_x - axis_horizontal_offset, 
                                        lower_y - (lower_y - upper_y) * (n/5.0)
                                    ), 
                                    'L %s %s' % (
                                        plotstart_x - axis_horizontal_offset - tick_length, 
                                        lower_y - (lower_y - upper_y) * (n/5.0)
                                    )]
            
            plot_path = self.dwg.path(
                                d=' '.join(commands), 
                                stroke = 'black', stroke_linecap='round', stroke_width= '3',
                                fill = 'none', fill_rule='evenodd'
                                )
            
            self.dwg.add(plot_path)
            
            textanc = 'end'
            ticklabel = self.dwg.text('%s %%' % (int((n/5.0)*100)), 
                                                    insert = (
                                                            plotstart_x - 30, 
                                                            lower_y - (lower_y - upper_y) * (n/5.0)
                                                    ),
                                                    fill = 'black', 
                                                    font_family = use_fontfamily, 
                                                    text_anchor = textanc, 
                                                    font_size = tick_font_size, 
                                                    baseline_shift = '-50%'
                                                )
            
            self.dwg.add(ticklabel)

        # label y-axis
        label_horizontal_offset = 35
        textanc = 'middle'
        lab = 'Nucleotide Identity'
        x, y = plotstart_x - label_horizontal_offset, lane_lower_y - (lower_y - upper_y) * 0.5
        yaxislabel = self.dwg.text(
                                lab, 
                                insert = (x - 100, y), 
                                fill='black', 
                                font_family = use_fontfamily, 
                                text_anchor = textanc, 
                                font_size = font_size, 
                                transform="rotate(270 %s,%s)" % (x - 100,y))  # , baseline_shift='-50%'

        self.dwg.add(yaxislabel)

        # label duplicate (A, B, etc)
        textanc = 'middle'
        isolatelabel = self.dwg.text(label, 
                                    insert = (
                                        plotstart_x - label_horizontal_offset * 4, 
                                        lane_lower_y - (lower_y - upper_y) * 1.3
                                    ), 
                                    fill='black', font_family=use_fontfamily, 
                                    text_anchor=textanc, 
                                    font_size = font_size, baseline_shift='-50%')

        self.dwg.add(isolatelabel)

    def plot_ambiguous_ranges(self, seq_label, 
                                    pair, 
                                    panel = ((1,1),(1,1)), 
                                    font_size = '20pt', 
                                    use_fontfamily = 'Nimbus Sans L',
                                    values_upper_prop = 0.75,
                                    values_lower_prop = 0.4):
        '''given the real position ordinate on chromosome with percent identity to a duplication, plot line. Also plot y-axis and scale'''

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

        # establish lane plot area (always single for percent ID at duplicate regions)
        lane_num = 0
        plot_spacing = 3
        lane_upper_y = lower_y - (lane_num + 1) * y_per_lane + plot_spacing
        lane_lower_y = lower_y - lane_num * y_per_lane - plot_spacing
        lane_plot_height = lane_lower_y - lane_upper_y

        # draw ambiguous (high identity) regions
        # first indicate regions for exclusion due to high identity anywhere
        # pale green
        colour = ('20', '90', '20', '%')
        for chrm_s,chrm_e in self.finder_info['ambiguous_ranges']:
            if chrm_s < pair[seq_label]['end'] and chrm_e > pair[seq_label]['start']:
                
                s = self.chrom2plot_x(chrm_s, pair, seq_label)
                e = self.chrom2plot_x(chrm_e, pair, seq_label)
                s_plot = self.chrom2plot_x(pair[seq_label]['start'], pair, seq_label)
                e_plot = self.chrom2plot_x(pair[seq_label]['end'], pair, seq_label)
                
                commands = ['M %s %s' % (
                                        max(s_plot, s) + plotstart_x, 
                                        lane_lower_y
                                        )]
                commands += ['L %s %s' % (
                                        max(s_plot, s) + plotstart_x, 
                                        lane_upper_y
                                        )]
                commands += ['L %s %s' % (
                                        min(e_plot, e) + plotstart_x, 
                                        lane_upper_y
                                        )]
                commands += ['L %s %s z' % (
                                        min(e_plot, e) + plotstart_x, 
                                        lane_lower_y
                                        )]
                plot_path = self.dwg.path(
                                    d=' '.join(commands), stroke = self.rgb(*colour), 
                                    stroke_linecap='round', stroke_width= '3', 
                                    fill = self.rgb(*colour), fill_rule='evenodd',
                                    fill_opacity = 0.25
                                    )
                
                self.dwg.add(plot_path)

        # first indicate regions with high identity in the current plot
        colour = ('80', '10', '10', '%')
        for chrm_s,chrm_e in pair[seq_label]['ambiguous_ranges']:
                
                s = self.chrom2plot_x(chrm_s, pair, seq_label)
                e = self.chrom2plot_x(chrm_e, pair, seq_label)
                
                commands = ['M %s %s' % (
                                        s + plotstart_x, 
                                        lane_lower_y
                                        )]
                commands += ['L %s %s' % (
                                        s + plotstart_x, 
                                        lane_upper_y
                                        )]
                commands += ['L %s %s' % (
                                        e + plotstart_x, 
                                        lane_upper_y
                                        )]
                commands += ['L %s %s z' % (
                                        e + plotstart_x, 
                                        lane_lower_y
                                        )]
                plot_path = self.dwg.path(
                                    d=' '.join(commands), stroke = self.rgb(*colour), 
                                    stroke_linecap='round', stroke_width= '3', 
                                    fill = self.rgb(*colour), fill_rule='evenodd',
                                    fill_opacity = 0.25
                                    )
                
                self.dwg.add(plot_path)


def loadFinderInfo(filein):
    finder_info = {}
    with _tarfile.open(filein, "r:gz") as tar:
        for member in tar:
            contents = _StringIO(tar.extractfile(member).read())
            try:
                # either json serialised conventional objects
                contents = _json.loads(contents.getvalue())
                if member.name == 'homologous_groups_mapped':
                    # need to restore some dict keys to integers
                    # or json them as .items()
                    for groups in contents:
                        for pair in groups:
                            pair['A']['aln_pos0_2_chrm_pos0_pIDs'] = dict([(int(k),v) for k,v in pair['A']['aln_pos0_2_chrm_pos0_pIDs'].items()])
                            pair['B']['aln_pos0_2_chrm_pos0_pIDs'] = dict([(int(k),v) for k,v in pair['B']['aln_pos0_2_chrm_pos0_pIDs'].items()])
                            pair['pIDs'] = dict([(int(k),v) for k,v in pair['pIDs'].items()])
                    
            except ValueError:
                # or longer python array.array objects
                contents = _array('c', contents.getvalue())
            
            finder_info[member.name] = contents

    return(finder_info)
def plotRepeats(genome_name, outdir = ['plots_repeats'], force = False):
    '''Plot repeats!'''

    # first load Repeats.Finder analysis results and related information

    from baga import Repeats

    filein = 'baga.Repeats.FinderInfo-{}.baga'.format(genome_name)
    finder_info = Repeats.loadFinderInfo(filein)

    outdir = _os.path.sep.join(outdir + [finder_info['genome_id']])
    if not _os.path.exists(outdir):
      _os.makedirs(outdir)

    # numbering here isn't working . . through each group contained all related pairs . . .
    plotted_last = False
    plot_group = 1
    #for groups in finder.homologous_groups_mapped:
    for groups in finder_info['homologous_groups_mapped']:
        if plotted_last:
            plot_group += 1
            plotted_last = False
        
        for pair in groups:
            
            # first check if this pair includes any of the 
            # final selection of ambiguous regions (without
            # the short regions)
            plot_A = False
            plot_B = False
            sA, eA = pair['A']['start'], pair['A']['end']
            sB, eB = pair['B']['start'], pair['B']['end']
            for s,e in finder_info['ambiguous_ranges']:
                if (s < eA and e > sA):
                    plot_A = True
                if (s < eB and e > sB):
                    plot_B = True
            
            if plot_A and plot_B:
                plotted_last = True
                plotname = '%02d_%07d_%07d_vs_%07d_%07d' % (
                                plot_group, sA, eA, sB, eB)
                
                plot_output_path = _os.path.sep.join([outdir, _os.path.extsep.join([plotname, 'svg'])])
                print('Plotting {}'.format(plot_output_path))
                
                if not force and _os.path.exists(plot_output_path):
                    print('Found %s\nUse "force = True" to overwrite' % plot_output_path)
                else:
                    repeat_plotter = Repeats.Plotter(finder_info, plot_output_path)
                    repeat_plotter.doPlot(pair)

if __name__ == '__main__':
    main()
