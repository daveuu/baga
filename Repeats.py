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
Repeats module from the Bacterial and Archaeal Genome Analyzer (BAGA).

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
        e = 'This chromosome object contains no ORF information: Repeats finding not possible'
        assert len(genome.ORF_ranges) > 0, e
        
        locus_chrom_info = dict(self.genome.ORF_ranges.items())
        # same for rRNA but add suffix: '__rRNA' to prevent translation later
        locus_chrom_info.update((rRNA_id+'__rRNA',[s,e,strd,n]) for rRNA_id,(s,e,strd,n) in \
            self.genome.rRNA_ranges.items())
        self.locus_chrom_info = locus_chrom_info
        ## NB ORFs_ordered includes rRNA encoding regions!
        self.ORFs_ordered = sorted(locus_chrom_info, key = locus_chrom_info.get)

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
            for ORF,(s,e,strand,gene_name) in self.locus_chrom_info.items():
                if s < aln.reference_end and aln.reference_start < e:
                    # exclude self
                    if ORF != aln.query_name:
                        s_q, e_q, strand, gene_name = \
                                self.locus_chrom_info[aln.query_name]
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

        all_ORF_hits = set(a for b in self.ORF_hits.values() for a in b)

        self.ORFs_with_hits_ordered = sorted(all_ORF_hits, key = self.locus_chrom_info.get)
    def getAdjacentORFsWithHit(self, ORF, maxdist = 3, direction = 1, getall = False):
        '''return next ORF with a paralog and num ORFs to it within maxdist'''
        thisORFn_hits = self.ORFs_with_hits_ordered.index(ORF)
        thisORFn_all = self.ORFs_ordered.index(ORF)
        if direction == -1:
            maxdist -= 1
            s = thisORFn_all + maxdist * direction
            e = thisORFn_all
        else:
            s = thisORFn_all
            e = thisORFn_all + (maxdist) * direction

        ORFS_to_check = self.ORFs_ordered[s:e][::direction]
        next_hit_ORFs = {}
        for n in range(1,maxdist):
            if thisORFn_hits + n * direction >= len(self.ORFs_with_hits_ordered):
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
        '''
        Delineate contiguous groups of ORFs with homology elsewhere in chromosome.

        The homoblocks dict is built during the algorithm and indicates 
        ORF-with-homolog membership to it's homologous block.

        The homologous_groups dict is added to when the extension (merging) of the 
        current block ends.


        For each ORF in a chromosome-ordered list in which all have homology hits 
        to one or more other ORFs:
          I) if no block is currently being extended, initiate required data sturctures:
            i) "this_block" list of ORF IDs in currently extending block
            ii) "homoblocks" dict: values are contiguous blocks (lists) or ORFs 
               elsewhere in the chromosome with homology to the currently extending 
               block of ORFs; keys are the ORFs homologous to the current ORF in the
               current block
            iii) "extended" dict: values are True if that block was extended, False if
                 not to keep a record of continuously extended blocks; keys are the first
                 ORF in each of the homologous, contiguous blocks in "homoblocks"
          
          II) check whether the current ORF has a neighbour within max_dist that has a
              homolog else where (and is potentially part of a homologous block to be
              included in current block)
              Use defaults: maxdist = 3, in direction = 1 (same as ORF's 5' to 3')
          
          III) as long as the neighbour does not hit anything in current block (which would 
               be a tandem repeat) iterate through hits of the neighbour
                i) get ORF strand
                ii) for each hit of the focal ORF, collect their own hits with:
                  a) adjacency (within three), strand (same)
                iii) do any of the adjacents of the hit, hit the adjacent of the focal ORF?
                  a) if they do, extend their entry in homoblocks dict
          
          IV) finally record if any homoblocks were extended in "extended" dict. With reference
              to "extended" dict, if no extensions were done, store all continuously extended
              homologous in "block_to_homoblocks" dict and reset "this_block" to trigger (I)

        '''

        # some sanity checks
        f1 = False
        for n,ORF1 in enumerate(self.ORFs_with_hits_ordered[:-1]):
            for ORF2 in self.ORFs_with_hits_ordered[(n+1):]:
                if self.ORF_hits[ORF1] & self.ORF_hits[ORF2]:
                    # if there is an intersection, it should be a union
                    if self.ORF_hits[ORF1] != self.ORF_hits[ORF1] | self.ORF_hits[ORF2]:
                        f1 = True

        for n,ORF1 in enumerate(self.ORFs_with_hits_ordered[:-1]):
            for ORF2 in self.ORFs_with_hits_ordered[(n+1):]:
                if self.ORF_hits[ORF1] & self.ORF_hits[ORF2]:
                    self.ORF_hits[ORF1] = self.ORF_hits[ORF1] | self.ORF_hits[ORF2]
                    self.ORF_hits[ORF2] = self.ORF_hits[ORF1] | self.ORF_hits[ORF2]

        for n,ORF1 in enumerate(self.ORFs_with_hits_ordered[::-1][:-1]):
            for ORF2 in self.ORFs_with_hits_ordered[::-1][(n+1):]:
                if self.ORF_hits[ORF1] & self.ORF_hits[ORF2]:
                    self.ORF_hits[ORF1] = self.ORF_hits[ORF1] | self.ORF_hits[ORF2]
                    self.ORF_hits[ORF2] = self.ORF_hits[ORF1] | self.ORF_hits[ORF2]

        f2 = False
        for n,ORF1 in enumerate(self.ORFs_with_hits_ordered[:-1]):
            for ORF2 in self.ORFs_with_hits_ordered[(n+1):]:
                if self.ORF_hits[ORF1] & self.ORF_hits[ORF2]:
                    if self.ORF_hits[ORF1] != self.ORF_hits[ORF1] | self.ORF_hits[ORF2]:
                        print(ORF1,ORF2)
                        f2 = True

        assert not f2

        # below, we'll be building "homoblocks": contiguous blocks
        # of homologous loci
        # homoblocks is ORFn1b1 => (ORFn1b1,ORFn2b1,...ORFnib1)
        # homoblocks is ORFn1b2 => (ORFn1b2,ORFn2b1,...ORFnib2)
        # continues as long as at least one is extended
        # eventually those as long as "this_block" (usd below) are
        # retained in block_to_homoblocks


        def retain_homoblocks(this_block, homoblocks, extended):
            '''
            determine how much of each block to keep based on continuous extension
            
            BBBBBBBBBB <-block
            TTTTTTTTTT <-full length homo-block (always at least one of these)
            TTTTTTFFFF <-partial length homo-block
            TTTFFFFTTT <-partial length homo-block that continues after break
            
            All of the continuous pieces are collected. Sometimes new focal blocks
            are also returned to match start point with homoblocks that restart or 
            start late
            '''
            collected_homoblocks = sorted(homoblocks.values())
            keep_homoblocks = {}
            for this_homoblock in collected_homoblocks:
                this_cont_block = []
                this_cont_homoblock = []
                b = 0
                for n,added in enumerate(extended[this_homoblock[0]]):
                    if added:
                        this_cont_block += [this_block[n]]
                        this_cont_homoblock += [this_homoblock[b]]
                        # separate iterator because only homologs added nothing corresponds to Falses
                        b += 1
                    elif len(this_cont_block):
                        try:
                            keep_homoblocks[tuple(this_cont_block)] += [tuple(this_cont_homoblock)]
                        except KeyError:
                            keep_homoblocks[tuple(this_cont_block)] = [tuple(this_cont_homoblock)]
                        this_cont_block = []
                        this_cont_homoblock = []
                
            return(keep_homoblocks)

        check_any_exts = []
        check = False
        blocks = []
        this_block = []
        # need to know that at least one homologous block extended continuously
        extended = {}  # just count how many ORFs added
        block_to_homoblocks = {}
        for n,thisORF in enumerate(self.ORFs_with_hits_ordered):
            thisORFhits = sorted(set(self.ORF_hits[thisORF]) - set([thisORF]), key = self.locus_chrom_info.get)
            print('{} ({}/{}) has {} hits'.format(thisORF,n+1,len(self.ORFs_with_hits_ordered),len(thisORFhits)))
            if thisORF in ('PP_5SA__rRNA',):
                check = True
                #import ipdb; ipdb.set_trace()
            
            # (I)
            if not len(this_block):
                # new block
                # initiate a dict of extending history
                # at least one homologous block must be extended continuously
                # i.e., same length as this_block
                # record true/false for each block initiated at same point in time:
                # need to know contiguity of shorter alignments
                # first True is first ORF in self.ORF_hits[thisORF]
                extended = {firstORFhit:[True] for firstORFhit in thisORFhits}
                homoblocks = {firstORFhit:[firstORFhit] for firstORFhit in thisORFhits}
                # start with current ORF which must have hits
                # else current ORF was tested and added in previous iteration
                this_block = [thisORF]
                print('Initiated for: {} versus {}'.format(thisORF, '+'.join(sorted(self.ORF_hits[thisORF]))))
            
            for thisORFhit in thisORFhits:
                if thisORFhit not in homoblocks and thisORFhit not in extended:
                    homoblocks[thisORFhit] = [thisORFhit]
                    # need to front pad with Falses because should be same length as block
                    extended[thisORFhit] = [False]*(len(this_block)-1) + [True]
                    print('New block: {}'.format(thisORFhit))
            
            # use to check for extensions
            homoblocks_extended = []
            # (II)
            if n + 1 == len(self.ORFs_with_hits_ordered):
                # finish on last ORF
                # extending a block over the putative origin of replication not implemented
                # but separate blocks either side will be detected (and could be joined later)
                thisORFnearestWithHit, dist_to_next = False, False
            else:
                thisORFnearestWithHit, dist_to_next = self.getAdjacentORFsWithHit(thisORF)
                
            if thisORFnearestWithHit:
                thisORFnearestWithHit_hits = self.ORF_hits[thisORFnearestWithHit]
                if len(set(this_block) & thisORFnearestWithHit_hits):
                    # next ORF has a homolog in this block == tandem repeat
                    print('Next ORF is tandem repeat: do not attempt to extend')
                else:
                    print('{} is within max_dist of {} and has hits'.format(thisORFnearestWithHit,thisORF))
                    # found one nearby with a hit (thisORFnearestWithHit)
                    # get strand info for matching in any homologous blocks
                    thisORFnearestWithHit_strand = self.locus_chrom_info[thisORFnearestWithHit][2]
                    thisORF_strand = self.locus_chrom_info[thisORF][2]
                    if thisORFnearestWithHit_strand == thisORF_strand:
                        same_strand = True
                    else:
                        same_strand = False
                    # do thisORFnearestWithHit and thisORF both have paralogs within
                    # maxdist of each other?
                    # (and is strand/orientation conserved?)
                    # this is tackled from the direction of thisORF's hits and their
                    # adjacents; do those hit thisORFnearestWithHit?
                    # (III)
                    thisORFnearestWithHit_hits = sorted(thisORFnearestWithHit_hits, key = self.locus_chrom_info.get)
                    new_blocks = []
                    for thisORFhit in thisORFhits:
                        print('thisORFhit: {}'.format(thisORFhit))
                        # collect all ORFs within maxdist with homologs (not just the nearest)
                        thisORFhit_strand = self.locus_chrom_info[thisORFhit][2]
                        ## what about an inverted ORF within a contiguity here? <== is right to break blocks here? yes because will continue in a new block
                        if thisORF_strand == thisORFhit_strand:
                            # orientations matches: check in same direction relative to sequence in file
                            thisORFhit_nearestWithHits = self.getAdjacentORFsWithHit(thisORFhit, direction = 1, getall = True)
                        else:
                            # any homologous block will be on other strand: check in other direction
                            thisORFhit_nearestWithHits = self.getAdjacentORFsWithHit(thisORFhit, direction = -1, getall = True)
                        
                        # check whether neighbours hit any of hit's neighbours
                        for thisORFnearestWithHit_hit in thisORFnearestWithHit_hits:
                            # *** this is the pivotal test of extending continguous homologous block pairs
                            # all of the hits of the potential ORF for extension are tested in the same way
                            if thisORFnearestWithHit_hit in thisORFhit_nearestWithHits:
                                thisORFnearestWithHit_hit_strand = self.locus_chrom_info[thisORFnearestWithHit_hit][2]
                                print('Found potential homologues in neighbour')
                                if (thisORFnearestWithHit_hit_strand == thisORFhit_strand) == same_strand:
                                    # strand orientation between this and the next in this homologous block matches
                                    print('orientations match: {} , {} ; {} , {}'.format(
                                            thisORFnearestWithHit_hit_strand, thisORFhit_strand, 
                                            thisORFnearestWithHit_strand, thisORF_strand))
                                    
                                    # extending at least one of the homologous blocks
                                    dist_to_thisORF_nearestWithHit = thisORFhit_nearestWithHits[thisORFnearestWithHit_hit]
                                    A = 'Focus:    %s has %s at %s dist\n' % (
                                            thisORF, thisORFnearestWithHit, 
                                            dist_to_next)
                                    B = 'Paralogs: %s has %s at %s dist\n' % (
                                            thisORFhit, thisORFnearestWithHit_hit, 
                                            dist_to_thisORF_nearestWithHit)
                                    print(A+B)
                                    # can only extend (existing block) if started at beginning of this block
                                    if thisORFhit in homoblocks and thisORFnearestWithHit_hit not in homoblocks:
                                        print('{} in a homologous block that started with current block'.format(thisORFhit))
                                        # second conditional checks for tandems in homoblocks (and doesn't add if already in)
                                        # *** this is the extension of the other contiguous block
                                        # as we iterate along, homoblocks ORF ==> block
                                        # also needs to iterate all ORF keys
                                        # so they correspond to the "current" ORF
                                        homoblocks[thisORFhit] += [thisORFnearestWithHit_hit]
                                        assert thisORFnearestWithHit_hit not in homoblocks, \
                                                'thisORFnearestWithHit_hit already in homoblocks'
                                        homoblocks[thisORFnearestWithHit_hit] = homoblocks[thisORFhit]
                                        print('Added with: {} to {}'.format(
                                                thisORFnearestWithHit_hit, '+'.join(homoblocks[thisORFnearestWithHit_hit])))
                                        del homoblocks[thisORFhit]
                                        # keep a record of which homologous blocks merged in each iteration
                                        homo_block_firstORF = homoblocks[thisORFnearestWithHit_hit][0]
                                        homoblocks_extended += [homo_block_firstORF]
                                        # only add one at a time else tandem repeats get added incorrectly
                                        break
            
            for firstORF in extended:
                if firstORF in homoblocks_extended:
                    extended[firstORF] += [True]
                else:
                    extended[firstORF] += [False]
            
            # add to this block now to test for continuous merge
            # if OK iteration will just continue, else will store blocks and reset
            this_block += [thisORFnearestWithHit]
            
            # (IV)
            # check at least one homologous block has been continuously extended
            #continuous_merge = len([b for b,a in extended.items() if sum(a) == len(this_block)])
            # check for any merging, retain_homoblocks() will now extract continuous blocks
            any_extending = any([b[-1] for a,b in extended.items()])
            check_any_exts += [any_extending]
            #import ipdb; ipdb.set_trace()
            if not any_extending or not thisORFnearestWithHit:
                print('Break in merging <======')
                # not any_extending: assessed adjacent ORFs but no contiguous
                # (can this ever happen in new method? if there is 
                # thisORFnearestWithHit at least a non-contiguous merge must happen)
                # /\/\/\/\/\ check this
                # not thisORFnearestWithHit: no ORFs with homologs close enough to be considered contiguous
                # tandem_repeat: thisORFnearestWithHit hits something already in block.
                # Stop and store block but re-initiate before next iteration
                check = False
                # either no extension made within max_dist
                # or the only extension was in a homologous block that already broke contiguity
                print('merge info')
                print('extended:',sorted(extended.items()))
                print('homo_blocks:',sorted(homoblocks.items()))
                # blocks = [] need to record in this?
                print('this_block',this_block)
                # also save homologous contiguous blocks
                keep_homoblocks = retain_homoblocks(this_block, homoblocks, extended)
                print('=====> Storing:')
                print(keep_homoblocks,'\n\n\n')
                for b,h in keep_homoblocks.items():
                    assert b not in block_to_homoblocks
                    block_to_homoblocks[b] = h
                else:
                    # reset for new block
                    this_block = []

        # collect existing shorter homologous blocks
        # put into new block => homoblocks
        new_block_to_homoblocks = {}
        for block,homoblocks in sorted(block_to_homoblocks.items()):
            for homoblock in homoblocks:
                if len(homoblock) < len(block):
                    # homoblock shorter than current block, need to check and/or make new entry
                    # compile same length homoblocks for this shorter one in currrent entry
                    new_homoblocks = set([block[:len(homoblock)]])
                    for b in homoblocks:
                        if b != homoblock and len(b) >= len(homoblock):
                            new_homoblocks.add(b[:len(homoblock)])
                    
                    if homoblock in new_block_to_homoblocks:
                        new_block_to_homoblocks[homoblock].update(new_homoblocks)
                    else:
                        new_block_to_homoblocks[homoblock] = new_homoblocks


        # add full length blocks with same length homoblocks
        for block,homoblocks in sorted(block_to_homoblocks.items()):
            keeps = set()
            for homoblock in homoblocks:
                if len(homoblock) >= len(block):
                    assert len(homoblock) == len(block), 'block {} '\
                            'points to longer homologous block'.format(
                            '+'.join(block),'+'.join(homoblock))
                    keeps.add(homoblock)
            
            if block in new_block_to_homoblocks:
                new_block_to_homoblocks[block].update(keeps)
            else:
                new_block_to_homoblocks[block] = keeps

        check_ORF_hits = {}
        for block,homoblocks in new_block_to_homoblocks.items():
            for theseHomologs in zip(*([block] + list(homoblocks))):
                for H1 in theseHomologs:
                    for H2 in theseHomologs:
                        try:
                            check_ORF_hits[H1].add(H2)
                        except KeyError:
                            check_ORF_hits[H1] = set([H2])
                        try:
                            check_ORF_hits[H2].add(H1)
                        except KeyError:
                            check_ORF_hits[H2] = set([H1])

        assert sorted(check_ORF_hits) == sorted(self.ORF_hits), ''\
                'lost or gained some homologs?!'

        gained = set()
        lost = set()
        for thisORF in self.ORFs_with_hits_ordered:
            if check_ORF_hits[thisORF] != self.ORF_hits[thisORF]:
                gained.update(check_ORF_hits[thisORF] - self.ORF_hits[thisORF])
                lost.update(self.ORF_hits[thisORF] - check_ORF_hits[thisORF])
                print(self.ORF_hits[thisORF] - check_ORF_hits[thisORF])

        print('gained: {}'.format(len(gained)))
        print('lost: {}'.format(len(lost)))
        # gain or loss isn't acceptable!

        # blocks are collected in a 5'3' +ve strand direction even if they are on the reverse strand
        # this means these and any +ve strand homologous blocks are collected "backwards"
        # e.g. A, B, C on -ve with Z, Y, X on +ve,
        # but 5'3' order on -ve is really C, B, A
        # and 5'3' order on +ve is really X, Y, Z
        # but for alignment purposes, A-C is collected then reverse complemented
        # while X-Z is collected, hence:
        homologous_groups = set()
        for this,those in sorted(new_block_to_homoblocks.items()):
            # innermost sort ensures loci are in chromomosome order, 
            # not loci label order
            homologous_groups.add(tuple(sorted(map(
                    lambda x: tuple(sorted(x, key = self.locus_chrom_info.get)),
                    [this]+list(those)))))

        self.homologous_groups = sorted(homologous_groups)

        self.block_to_homoblocks = block_to_homoblocks

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


    def align_blocks(self, extend_len = 200,
                                max_extensions = 10,
                                min_pID = 0.95,
                                num_terminal_window_steps = 5):
        
        '''
        Align homologous blocks using the Needleman Wunch pairwise global alignment
        
        extension options:
            extend_len: initial extension length
            max_extensions: maximum extension iterations
            min_pID: min percent ID over . . .
            num_terminal_window_steps: . . . this number of moving window steps, 
        to continue extension
        '''
        

        def countGaps(seq_str, from_end = False):
            '''Count number of terminal gap symbols in a string'''
            if from_end:
                seq_str = seq_str[::-1]
            
            num_gaps = 0
            for char in seq_str:
                if char == '-':
                    num_gaps += 1
                else:
                    break
            return(num_gaps)

        def alignNW(seqrecords, exe_aligner = self.exe_aligner):
            '''Needleman-Wunsch alignment using seqalign
            
            seqrecords must be a list of two BioPython SeqRecord instances
            '''
            out_handle = _StringIO()
            _SeqIO.write(seqrecords, out_handle, 'fasta')
            cmd = [exe_aligner]
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
            # put aligned amino acids into a SeqRecord
            A = _SeqRecord(_Seq(A), id = 'A')
            B = _SeqRecord(_Seq(B), id = 'B')
            return([A,B])

        def alignNW_Nuc_as_AA(nuc_A_seq, nuc_B_seq):
            A_nuc, B_nuc = _SeqRecord(nuc_A_seq, id = 'A'), _SeqRecord(nuc_B_seq, id = 'B')
            A_aa, B_aa = _SeqRecord(nuc_A_seq.translate(), id = 'A'), _SeqRecord(nuc_B_seq.translate(), id = 'B')
            A_aa, B_aa = alignNW([A_aa, B_aa])
            # align nucleotides as codons to aligned amino acids
            alnd = self.Nuc2AA([A_nuc, B_nuc], [A_aa, B_aa])
            return(alnd['A'], alnd['B'])

        # for inter-ORF nucleotides:
        # alignNW(seqrecords)
        # for ORFs
        # alignNW_NUC_as_AA(nuc_A_seq, nuc_B_seq)

        def collectForAligning(repeated_loci_ranges, repeated_loci, genome_seq):
            repeated_locus_types = []
            repeated_locus_strings = []
            for i in range(len(repeated_loci_ranges)-1):
                start = repeated_loci_ranges[i]
                end = repeated_loci_ranges[i+1]
                repeated_locus_strings += [genome_seq[start:end]]
                # ORFs and inter-ORFs alternate
                if i % 2 == 0:
                    if repeated_loci[int(i/2)].endswith('__rRNA'):
                        repeated_locus_types += ['rRNA']
                    else:
                        repeated_locus_types += ['ORF']
                else:
                    repeated_locus_types += ['inter']
            
            return({'types':repeated_locus_types,'seqs':repeated_locus_strings})


        def reverseRange(start, end, sequence):
            '''given start and end ordinates and sequence, return equivalents on reverse strand'''
            return(len(sequence) - end, len(sequence) - start)

        # forward = _Seq('TGGCCAAAGACACCATCGTGCTGCGCGAATC')
        # reverse = forward.reverse_complement()
        # s, e = 12, 19
        # str(reverse[slice(*reverseRange(s,e,forward))]) == str(forward[s:e].reverse_complement())

        def collectLociRanges(repeated_loci, 
                              ORF_ranges, 
                              genome_sequence, 
                              strand):
            ORF_overlaps = []
            loci_ranges = []
            for i,locus in enumerate(repeated_loci[::strand]):
                start, end = ORF_ranges[locus][:2]
                ## ensure all coding regions are whole codons
                ## sometimes incomplete codons get through in
                ## annotated pseudogenes or errors
                ends_codon_diff = (end - start) % 3
                end = end - ends_codon_diff
                if strand == -1:
                    start, end = reverseRange(start, end, genome_sequence)
                
                if i > 0 and start < loci_ranges[-1]:
                    # this ORF start is prior to previous ORF end
                    # record overlaps for pair-specific correction
                    # record index of first one
                    ORF_overlaps += [i-1]
                
                loci_ranges += [start, end]
            
            return(loci_ranges, ORF_overlaps)


        def do_extension(genome_use_A, 
                         genome_use_B, 
                         extensionA_start, 
                         extensionB_start, 
                         extend_len,
                         max_extensions,
                         direction = 1):
            
            extensionA_end = int(extensionA_start)
            extensionB_end = int(extensionB_start)
            for extend_num in range(max_extensions):
                # collect next chunk
                extensionA_end += (extend_len * direction)
                # print(extensionA_start,extensionA_end,extend_len)
                # allow for extending at start
                ext_start_use, ext_end_use = sorted([extensionA_start,extensionA_end])
                SeqRecA = _SeqRecord(genome_use_A[ext_start_use:ext_end_use], id = 'A')
                extensionB_end += (extend_len * direction)
                # print(extensionB_start,extensionB_end,extend_len)
                ext_start_use, ext_end_use = sorted([extensionB_start,extensionB_end])
                SeqRecB = _SeqRecord(genome_use_B[ext_start_use:ext_end_use], id = 'B')
                # align next chunk
                SeqRecA_alnd,SeqRecB_alnd = alignNW([SeqRecA,SeqRecB])
                # print(str(SeqRecA_alnd.seq))
                # print(str(SeqRecB_alnd.seq))
                pIDs = self.get_percent_ID(SeqRecA_alnd.seq,SeqRecB_alnd.seq, window = 100, step = 20)
                # if end of extended bit is high ID, extend and test again until pID drops off
                mean_pID = sum([pID for pos,pID in pIDs][::direction][-num_terminal_window_steps:])/float(num_terminal_window_steps)
                print('mean pID: {:.0%}, min pID: {:.0%} after {:,} bp'.format(mean_pID,min_pID,(extend_num+1)*extend_len))
                if mean_pID < min_pID:
                    # extension divergent enough be at end of duplication
                    print('end of extension at ends of contiguities')
                    break
                
                if extend_num == max_extensions-1:
                    # should have broken by now: no more iterations
                    # add this aligned bit anyway
                    print('WARNING: extended %s x %s bp but regions not yet divergent enough' % (max_extensions, extend_len))
            
            return(str(SeqRecA_alnd.seq), extensionA_end, str(SeqRecB_alnd.seq), extensionB_end)


        ## this code block is to account for overlapping features
        ## (typically ORFs are non-overlapping)
        ## bases are only aligned once but at an overlap would
        ## be translated into two different AA sequences. To prevent
        ## ensure bases are only translated and aligned once, one
        ## ORF range is shortened.
        ## See notes below in main loop for more info

        def update_ORF_ends(
            overlap_nuc_len_A,
            preORFA_end, 
            postORFA_start,
            num_gaps_preA,
            num_gaps_postA,
            overlap_nuc_len_B,
            preORFB_end, 
            postORFB_start,
            num_gaps_preB,
            num_gaps_postB):
            '''Deal with several scenarios of ORF overlaps by updating starts, ends
            
            (i) overlap_nuc_len_A and overlap_nuc_len_B are identical in size and 
            there are no terminal gaps in num_gaps_preA, postA, preB, postB.
            
            (ii) overlap_nuc_len_A corresponds to number of gaps either num_gaps_preA
            of _postA.
            
            (iii) same as (ii) but for B
            '''
            # (i)
            if overlap_nuc_len_A > 0 and overlap_nuc_len_B > 0 and \
                    overlap_nuc_len_A == overlap_nuc_len_B:
                # increase post-ORFs this number of codons
                num_codons_diff = int(round(overlap_nuc_len_A / 3.0 + 0.5))
                new_postORFA_start = postORFA_start + (num_codons_diff*3)
                new_postORFB_start = postORFB_start + (num_codons_diff*3)
                print('attempted fix for same overlap in both blocks')
                #import pdb; pdb.set_trace()
                return([preORFA_end, new_postORFA_start], [preORFB_end, new_postORFB_start])
                # loci_ranges_updated_A += [preORFA_end, new_postORFA_start]
                # loci_ranges_updated_B += [preORFB_end, new_postORFB_start]
            
            # (ii) or (iii)
            # deal with overlaps independently in block A and B
            # NB overlaps in A are checked against terminal gaps in B
            # and vice versa
            overlap_info = []
            overlap_info += [('A', overlap_nuc_len_A, preORFA_end, postORFA_start, 
                    num_gaps_preB, num_gaps_postB)]
            overlap_info += [('B', overlap_nuc_len_B, preORFB_end, postORFB_start, 
                    num_gaps_preA, num_gaps_postA)]
            
            new_ranges = []
            for (block, overlap_nuc_len, preORF_end, postORF_start, num_gaps_pre, 
                    num_gaps_post) in overlap_info:
                new_preORF_end = preORF_end
                new_postORF_start = postORF_start
                if overlap_nuc_len > 0:
                    print('Overlap in block {} of {} bp'.format(block, overlap_nuc_len))
                    ## create new ranges creating an inter-ORF region
                    # first update according to any terminal gaps in codon alignment
                    # of other block for both the pre-ORF and post-ORF; if there is
                    # still an overlap, remove additional codon-multiple lengths until
                    # there isn't
                    if num_gaps_pre > 0:
                        print('Caused {} gaps in other block pre-ORF'\
                                ''.format(num_gaps_pre))
                        # late end in pre-ORF A
                        new_preORF_end -= num_gaps_pre
                    if num_gaps_post > 0:
                        print('Caused {} gaps in other block post-ORF'\
                                ''.format(num_gaps_post))
                        # early start in this block post-ORF: increment forward
                        # num_gaps_post is always codons
                        new_postORF_start += num_gaps_post
                    # test again for overlap
                    if not new_preORF_end < new_postORF_start:
                        # still overlapping: overlap is longer than terminal gaps
                        # in codon alignment
                        # remove required additional length from end of pre-ORF
                        # positive is overlap
                        remaining_overlap = new_preORF_end - new_postORF_start
                        # round up overlap to nearest codon in post ORF
                        num_codons_diff = int(round(remaining_overlap / 3.0 + 0.5))
                        new_preORF_end -= (num_codons_diff*3)
                        print('attempted fix for an overlap greater than gaps at end of codon alignment')
                        print('remaining overlap: {}; terminal gaps: {} pre, {} post'.format(
                                remaining_overlap, num_gaps_pre, num_gaps_post))
                        #import pdb; pdb.set_trace()
                    if not new_preORF_end < new_postORF_start:
                        # overlap is longer than terminal gaps in codon alignment
                        # remove required additional length from start of post-ORF
                        # positive is overlap
                        remaining_overlap = new_preORF_end - new_postORF_start
                        # round up overlap to nearest codon in post ORF
                        num_codons_diff = int(round(remaining_overlap / 3.0 + 0.5))
                        new_postORF_start += (num_codons_diff*3)
                        print('attempted fix for an overlap greater than gaps at end of codon alignment')
                        print('remaining overlap: {}; terminal gaps: {} pre, {} post'.format(
                                remaining_overlap, num_gaps_pre, num_gaps_post))
                        #import pdb; pdb.set_trace()
                
                new_ranges += [(new_preORF_end, new_postORF_start)]
            return(new_ranges)

        pID_window_size = 200

        genome_strands = {}
        genome_strands[1] = _Seq(self.genome.sequence.tostring())
        genome_strands[-1] = genome_strands[1].reverse_complement()

        homologous_groups_alnd = []
        for g_n,groups in enumerate(self.homologous_groups):
            print('\n\nHomologous Group {}'.format(g_n))
            # for each group, align pairwise combinations of homologous, contiguous blocks
            alignment_combos = {}
            for n,repeated_loci_A in enumerate(groups[:-1]):
                print('A: {}'.format(' - '.join(repeated_loci_A)))
                strandA = self.locus_chrom_info[repeated_loci_A[0]][2]
                # print(repeated_loci_A, strandA)
                genome_use_A = genome_strands[strandA]
                # collect ranges for these loci
                # locus|inter-locus|locus etc
                loci_ranges_use_A, ORF_overlaps_A = collectLociRanges(repeated_loci_A,
                                                                      self.locus_chrom_info, 
                                                                      genome_use_A,
                                                                      strandA)
                # print(loci_ranges_use_A, ORF_overlaps_A)
                for repeated_loci_B in groups[(n+1):]:
                    print('vs. B: {}'.format(' - '.join(repeated_loci_B)))
                    strandB = self.locus_chrom_info[repeated_loci_B[0]][2]
                    # print(repeated_loci_B, strandB)
                    genome_use_B = genome_strands[strandB]
                    # collect ranges for these loci
                    # locus|inter-locus|locus etc
                    
                    #### loci_ranges_use_B ==     [4486395, 4486716, 4486712, 4487048, 4487111, 4488647] <== proper codons
                    #### loci_ranges_updated_B == [4486395, 4486716, 4486716, 4487048, 4487111, 4488647] <== not proper codons so fails assertion . . .
                    
                    loci_ranges_use_B, ORF_overlaps_B = collectLociRanges(repeated_loci_B,
                                                                          self.locus_chrom_info, 
                                                                          genome_use_B,
                                                                          strandB)
                    # print(loci_ranges_use_B, ORF_overlaps_B)
                    if len(repeated_loci_B) > 1:
                        # more than one consecutive ORF means there might be overlaps
                        # (some nucleotides appearing in two different ORFs - but
                        # should only be aligned once for moving window percent identity)
                        loci_ranges_updated_A = [loci_ranges_use_A[0]]
                        loci_ranges_updated_B = [loci_ranges_use_B[0]]
                        for i in range(len(repeated_loci_A)-1):
                            # i is this ORF, i+1 is next
                            # pre- and post- is relative to the inter-ORF region
                            # (which doesn't really exist in the case of the ORF overlap)
                            preORFA_start, preORFA_end = loci_ranges_use_A[(i)*2:(i)*2+2]
                            preORFB_start, preORFB_end = loci_ranges_use_B[(i)*2:(i)*2+2]
                            postORFA_start, postORFA_end = loci_ranges_use_A[(i+1)*2:(i+1)*2+2]
                            postORFB_start, postORFB_end = loci_ranges_use_B[(i+1)*2:(i+1)*2+2]
                            if i not in ORF_overlaps_A and i not in ORF_overlaps_B:
                                # neither homologous block A nor block B overlap between these
                                # ORFs: nothing to adjust
                                # add current ORF ranges and move on
                                loci_ranges_updated_A += [preORFA_end, postORFA_start]  # [preORFA_end, new_postORFA_start]
                                loci_ranges_updated_B += [preORFB_end, postORFB_start]
                            else:
                                # get overlap lengths (positive means overlapping)
                                overlap_nuc_len_A = preORFA_end - postORFA_start
                                overlap_nuc_len_B = preORFB_end - postORFB_start
                                
                                ## this code block is to account for overlapping features
                                ## (typically ORFs are non-overlapping)
                                ## bases are only aligned once but at an overlap would
                                ## be translated into two different AA sequences. To prevent
                                ## ensure bases are only translated and aligned once, one
                                ## ORF range is shortened.
                                
                                ## if only overlap exists in only one block, and that overlap is caused by
                                ## an increase of decrease in an ORF in one block and not other,
                                ## terminal gaps in alignment between the ORF pair indicates which ORF should
                                ## be shortened: the one without the gaps
                                
                                ## Occasionally overlap affected regions are encountered which
                                ## are already off the end of one repeat so that the above
                                ## assumptions do not hold and the overlap is not fixed with
                                ## the below code: but if these are already at the end of a repeat
                                ## it doesn't really matter . . .
                                
                                # compare num AAs per homologous ORFs
                                # identify if length difference of prior or subsequent homologous ORFs are closest to overlap length
                                # align pre- and post- ORFs as original dimensions
                                # only align as codons if ORFs: could be encoded rRNA
                                if repeated_loci_A[i].endswith('__rRNA') or \
                                        repeated_loci_B[i].endswith('__rRNA'):
                                    print('========> found rRNA {}, {}'.format())
                                    preORFA_aln, preORFB_aln = alignNW_Nuc(
                                            genome_use_A[preORFA_start:preORFA_end], 
                                            genome_use_B[preORFB_start:preORFB_end])
                                else:
                                    preORFA_aln, preORFB_aln = alignNW_Nuc_as_AA(
                                            genome_use_A[preORFA_start:preORFA_end], 
                                            genome_use_B[preORFB_start:preORFB_end])
                                if repeated_loci_B[i+1].endswith('__rRNA') or \
                                        repeated_loci_B[i+1].endswith('__rRNA'):
                                    print('========> found rRNA {}, {}'.format())
                                    postORFA_aln, postORFB_aln = alignNW_Nuc(
                                            genome_use_A[postORFA_start:postORFA_end], 
                                            genome_use_B[postORFB_start:postORFB_end])
                                else:
                                    postORFA_aln, postORFB_aln = alignNW_Nuc_as_AA(
                                            genome_use_A[postORFA_start:postORFA_end], 
                                            genome_use_B[postORFB_start:postORFB_end])
                                
                                # return any gaps at end of each alignment to diagnose which ORF to
                                # shorten
                                num_gaps_preA = countGaps(preORFA_aln.seq, from_end = True)
                                num_gaps_preB = countGaps(preORFB_aln.seq, from_end = True)
                                num_gaps_postA = countGaps(postORFA_aln.seq, from_end = False)
                                num_gaps_postB = countGaps(postORFB_aln.seq, from_end = False)
                                
                                print('Assessing overlaps between ORFs {} and {} '\
                                        'in these homologous blocks'.format(i, i+1))
                                
                                # get non-overlapping pre-ORF end and post-ORF start
                                # maintain codon lengths
                                # first attempt to correct lengths by removing
                                # additional codons on one or other
                                
                                new_ORFA_ords, new_ORFB_ords = update_ORF_ends(
                                        overlap_nuc_len_A,
                                        preORFA_end, 
                                        postORFA_start,
                                        num_gaps_preA,
                                        num_gaps_postA,
                                        overlap_nuc_len_B,
                                        preORFB_end, 
                                        postORFB_start,
                                        num_gaps_preB,
                                        num_gaps_postB)
                                
                                loci_ranges_updated_A += new_ORFA_ords
                                loci_ranges_updated_B += new_ORFB_ords
                        
                        loci_ranges_updated_A += [postORFA_end]
                        loci_ranges_updated_B += [postORFB_end]
                        # check updated versions are complete and without overlaps
                        e1 = 'odd number of ORF-inter-ORF boundaries'
                        e2 = 'overlap correction failed'
                        assert len(loci_ranges_updated_A) == len(loci_ranges_use_A), e1
                        assert loci_ranges_updated_A == sorted(loci_ranges_updated_A), e2
                        assert len(loci_ranges_updated_B) == len(loci_ranges_use_B), e1
                        assert loci_ranges_updated_B == sorted(loci_ranges_updated_B), e2
                        # check all corrections yielded codons
                        assert all([(loci_ranges_updated_A[i+1]-loci_ranges_updated_A[i])%3==0 \
                                for i in range(0,len(loci_ranges_updated_A),2)]), 'not all A '\
                                'loci updated to codons lengths'
                        assert all([(loci_ranges_updated_B[i+1]-loci_ranges_updated_B[i])%3==0 \
                                for i in range(0,len(loci_ranges_updated_B),2)]), 'not all B '\
                                'loci updated to codons lengths'
                        
                    else:
                        loci_ranges_updated_A = list(loci_ranges_use_A)
                        loci_ranges_updated_B = list(loci_ranges_use_B)
                    
                    # all corrected ORF ranges should be codon multiple lengths
                    # but uncorrected (non-overlapping) might have been annotated erroneously
                    # or genuinely include partial codons . . . (pseudogenes?)
                    # so ensure all are codon length prior to translations
                    loci_ranges_updated_A_fixed = []
                    loci_ranges_updated_B_fixed = []
                    for i in range(len(loci_ranges_updated_A)/2):
                        s,e = loci_ranges_updated_A[(i*2):(i*2)+2]
                        ends_codon_diff = (e - s) % 3
                        loci_ranges_updated_A_fixed += [s, e - ends_codon_diff]
                        s,e = loci_ranges_updated_B[(i*2):(i*2)+2]
                        ends_codon_diff = (e - s) % 3
                        loci_ranges_updated_B_fixed += [s, e - ends_codon_diff]
                    
                    loci_ranges_updated_A = loci_ranges_updated_A_fixed
                    loci_ranges_updated_B = loci_ranges_updated_B_fixed
                    
                    ## now do the actual aligning
                    repeated_seqs2aln_A = collectForAligning(loci_ranges_updated_A, 
                            repeated_loci_A, genome_use_A)
                    repeated_seqs2aln_B = collectForAligning(loci_ranges_updated_B, 
                            repeated_loci_B, genome_use_B)
                    Aseq_all_alnd = []
                    Bseq_all_alnd = []
                    for i,(Aseq,Bseq) in enumerate(zip(repeated_seqs2aln_A['seqs'],
                                                       repeated_seqs2aln_B['seqs'])):
                        if repeated_seqs2aln_A['types'][i] == 'ORF':
                            Aseq_aln, Bseq_aln = alignNW_Nuc_as_AA(Aseq, Bseq)
                        else:
                            # inter-ORF or rRNA
                            Aseq_aln, Bseq_aln = alignNW([_SeqRecord(Aseq, id = 'A'), 
                                                          _SeqRecord(Bseq, id = 'B')])
                        
                        Aseq_all_alnd += [str(Aseq_aln.seq)]
                        Bseq_all_alnd += [str(Bseq_aln.seq)]
                    
                    ## now extend alignments at each end
                    if len(self.ORFs_in_tandem_repeats.intersection(
                            repeated_loci_A + repeated_loci_B)) > 0:
                        print('Not extending tandem repeats')
                    else:
                        ## extend at end
                        # check for near-100% identity at end and extend further as necessary
                        # but not if ORFs in tandem repeats because:
                        # alignments will be to other parts of repeat, not extending contiguities
                        # last one, extend to get to end of duplication
                        # get current (non-extended percent identity)
                        pIDs = self.get_percent_ID(Aseq_all_alnd[-1], Bseq_all_alnd[-1], 
                                window = pID_window_size, step = 20)
                        mean_pID = sum([pID for pos,pID in pIDs][::1][-num_terminal_window_steps:])/float(num_terminal_window_steps)
                        if mean_pID > min_pID:
                            # end of existing alignment not divergent enough to end it: extend
                            print('extending at end')
                            extensionA_start = int(loci_ranges_updated_A[-1])
                            extensionB_start = int(loci_ranges_updated_B[-1])
                            SeqRecA_alnd_seq, extensionA_end, SeqRecB_alnd_seq, extensionB_end = do_extension(genome_use_A, 
                                                                                                              genome_use_B, 
                                                                                                              extensionA_start, 
                                                                                                              extensionB_start, 
                                                                                                              extend_len,
                                                                                                              max_extensions,
                                                                                                              direction = 1)
                            
                            Aseq_all_alnd += [SeqRecA_alnd_seq]
                            Bseq_all_alnd += [SeqRecB_alnd_seq]
                            # update last position
                            loci_ranges_updated_A[-1] = extensionA_end
                            loci_ranges_updated_B[-1] = extensionB_end
                        else:
                            print('already divergent at end ({:.0%}), not extending'.format(mean_pID))
                        
                        ## extend at end
                        # need to reverse sequences
                        pIDs = self.get_percent_ID(Aseq_all_alnd[0], Bseq_all_alnd[0], 
                                window = pID_window_size, step = 20)
                        mean_pID = sum([pID for pos,pID in pIDs][::-1][-num_terminal_window_steps:])/float(num_terminal_window_steps)
                        if mean_pID > min_pID:
                            # end of existing alignment not divergent enough to end it: extend
                            print('extending at beginning')
                            extensionA_start = int(loci_ranges_updated_A[0])
                            extensionB_start = int(loci_ranges_updated_B[0])
                            SeqRecA_alnd_seq, extensionA_end, SeqRecB_alnd_seq, extensionB_end = do_extension(genome_use_A, 
                                                                                                              genome_use_B, 
                                                                                                              extensionA_start, 
                                                                                                              extensionB_start, 
                                                                                                              extend_len,
                                                                                                              max_extensions,
                                                                                                              direction = -1)
                            Aseq_all_alnd.insert(0, SeqRecA_alnd_seq)
                            Bseq_all_alnd.insert(0, SeqRecB_alnd_seq)
                            # update first position
                            loci_ranges_updated_A[0] = extensionA_end
                            loci_ranges_updated_B[0] = extensionB_end
                        else:
                            print('already divergent at start ({:.0%}), not extending'.format(mean_pID))
                    
                    ## now store for next stage
                    # keep sequences in 5-3 orientation as aligned
                    A = ''.join(Aseq_all_alnd)
                    B = ''.join(Bseq_all_alnd)
                    
                    # but store ordinates reversed for -ve strand
                    if strandA == 1:
                        delimitersA = loci_ranges_updated_A[0], loci_ranges_updated_A[-1]
                    else:
                        # from start last aligned ORF pair member (with extension) to end of first aligned pair member
                        delimitersA = reverseRange(loci_ranges_updated_A[0], loci_ranges_updated_A[-1], genome_use_A)
                    
                    if strandB == 1:
                        delimitersB = loci_ranges_updated_B[0], loci_ranges_updated_B[-1]
                    else:
                        # from start last aligned ORF pair member (with extension) to end of first aligned pair member
                        delimitersB = reverseRange(loci_ranges_updated_B[0], loci_ranges_updated_B[-1], genome_use_B)
                    
                    # delimiter relative position determine strand
                    alignment_combos[repeated_loci_A, repeated_loci_B] = ((delimitersA, A),(delimitersB, B))
            
            homologous_groups_alnd += [alignment_combos]

        self.homologous_groups_alnd = homologous_groups_alnd

    def get_percent_ID(self, A, B, window = 100, step = 20):
        pID_per_window = []
        for i in range(0, len(A)-window, step):
            Achunk = A[i:i+window]
            Bchunk = B[i:i+window]
            pID = sum([a == b for a,b in zip(Achunk,Bchunk)])/float(window)
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

        e =  'The start is after end in the aligned sequence range. '\
        'Reverse strand ranges should be be ascending as if on forward strand. '\
        'Reverse complement is then calculated. Start is {:,}, end is {:,}, '\
        'strand is {}'.format(
        aligned_seq['start'], aligned_seq['end'], strand)
        assert aligned_seq['start'] < aligned_seq['end'], e

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
                    chromchar = _Seq(
                            self.genome.sequence[alnd_chrom_pos0:alnd_chrom_pos0+1].\
                            tostring()).reverse_complement()[0]
                
                #print(char, chromchar)
                e = 'mismatch detected when assigning percent identity between '\
                'duplications to chromosome positions . . . :-(\nReference '\
                'chromosome used for plotting may not match that used for '\
                "aligning duplicate region pairs? (else there's a bug)\n"\
                'x at {}: {}={}\n{}'.format(aln_pos0, char, chromchar,aligned_seq['seq_str'])
                assert char == chromchar, e
                #print(strand,alnd_chrom_pos0,aln_pos0,char,chromchar,aligned_seq['start'],aligned_seq['end'])
                try:
                    # account for width of window
                    # within which percent ID calculated
                    pIDs_by_chrmpos0[alnd_chrom_pos0 + (window / 2 * strand)] = pIDs[aln_pos0]
                    #print('> at {}: {}={}'.format(aln_pos0, char, chromchar))
                except KeyError:
                    pass
            
                alnd_chrom_pos0 += strand   # this needs to decrements if strand == -1


        return(pIDs_by_chrmpos0)

    def map_alignments_to_chromosome(self):

        homologous_groups_mapped = []

        for g_n,group in enumerate(self.homologous_groups_alnd):
            # one pairwise: two homologs
            # three pairwise: three homologs
            these_pairwise = []
            for ORFs, ORFs_info in group.items():
                ORFsA, ORFsB = ORFs
                aligned_A = {}
                aligned_B = {}
                # start and end are base-0 chromosome sequence Python slices
                (aligned_A['start'], aligned_A['end']), aligned_A['seq_str'] = ORFs_info[0]
                (aligned_B['start'], aligned_B['end']), aligned_B['seq_str'] = ORFs_info[1]
                # strands of each region
                aligned_A['strand'] = self.locus_chrom_info[ORFsA[0]][2]
                aligned_B['strand'] = self.locus_chrom_info[ORFsB[0]][2]
                # percent identity over aligned region
                pIDs = dict(self.get_percent_ID(aligned_A['seq_str'], aligned_B['seq_str'], window = 100, step = 20))
                # map the pIDs to from alinmnet to chromosome
                # print('A: {} {} {}'.format(aligned_A['start'], aligned_A['end'], aligned_A['strand']))
                # print('B: {} {} {}'.format(aligned_B['start'], aligned_B['end'], aligned_B['strand']))
                #print('A: {}'.format(aligned_A['seq_str']))
                aligned_A['aln_pos0_2_chrm_pos0_pIDs'] = self.aln_pos0_2_chrm_pos0_pIDs(aligned_A, pIDs, window = 100)
                #print('B: {}'.format(aligned_B['seq_str']))
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

        if len(all_positions) > 0:
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
        else:
            print('No repeats found.')
            self.ambiguous_ranges = []




    def findRepeats(self, minimum_percent_identity = 0.95, 
                          minimum_repeat_length = 400,
                          max_extensions = 20,
                          exe_bwa = False, 
                          exe_samtools = False,
                          exe_exonerate = False,
                          local_repeats_path = ['repeats'], 
                          local_genomes_path = ['genome_sequences'], 
                          retain_bwa_output = True, 
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
            self.exe_bwa = \
                    _dependencies['bwa']['checker']['arguments']['path'][-1]
        else:
            self.exe_bwa = exe_bwa


        if not exe_samtools:
            self.exe_samtools = get_exe('samtools')
        elif exe_samtools == 'system':
            self.exe_samtools = \
                    _dependencies['samtools']['checker']['arguments']['path'][-1]
        else:
            self.exe_samtools = 'samtools'
            

        if not exe_aligner:
            self.exe_aligner = get_exe('seq-align')
        elif exe_bwa == 'system':
            self.exe_aligner = \
                    _dependencies['seq-align']['checker']['arguments']['path'][-1]
        else:
            self.exe_aligner = 'seq-align'


        local_repeats_path = _os.path.sep.join(local_repeats_path)
        local_genomes_path = _os.path.sep.join(local_genomes_path)

        genome_fna = '{}/{}.fna'.format(local_genomes_path, self.genome.id)
        genome_fna_ORFs = '{}/{}_ORFs.fna'.format(local_genomes_path,
                self.genome.id)

        try:
            _os.makedirs(local_genomes_path)
        except OSError:
            pass

        try:
            _os.makedirs(local_repeats_path)
        except OSError:
            pass


        print('Writing genome to FASTA')
        _SeqIO.write(_SeqRecord(_Seq(self.genome.sequence.tostring()), 
                id = self.genome.id), genome_fna, 'fasta')

        #### align ORFs to chromosome ####
        cmd = [self.exe_bwa, 'index', genome_fna]
        print('Called: {}'.format(' '.join(cmd)))
        try:
            _subprocess.call(cmd)
        except OSError:
            print('Problem running BWA at {}. Please use Dependencies module '\
                    'to install locally or check system path'.format(cmd[0]))


        # collect ORF sequences
        recs = []
        for ID,(s,e,strand,gene_name) in self.genome.ORF_ranges.items():
            if strand == 1:
                recs += [_SeqRecord(_Seq(self.genome.sequence[s:e].tostring()
                        ), id = ID)]
            else:
                recs += [_SeqRecord(_Seq(self.genome.sequence[s:e].tostring()
                        ).reverse_complement(), id = ID)]

        # same for rRNA but add suffix: '__rRNA' to prevent translation later
        for ID,(s,e,strand,gene_name) in self.genome.rRNA_ranges.items():
            if strand == 1:
                recs += [_SeqRecord(_Seq(self.genome.sequence[s:e].tostring()
                        ), id = ID+'__rRNA')]
            else:
                recs += [_SeqRecord(_Seq(self.genome.sequence[s:e].tostring()
                        ).reverse_complement(), id = ID+'__rRNA')]

        _SeqIO.write(recs, genome_fna_ORFs, 'fasta')

        ## relaxed alignments using BWA i.e. find modestly divergent ORFs
        cmd = [self.exe_bwa, 'mem',
            # Band width. Essentially, gaps longer than INT will not be found
            # ... also affected by the scoring matrix and the hit length [100]
            '-w', '100',
            # Off-diagonal X-dropoff (Z-dropoff). Stop extension when the
            # difference between the best and the current extension score is
            # ... avoids unnecessary extension ... [100]
            '-d', '100',
            # Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. Larger
            # value yields fewer seeds, which leads to faster alignment speed
            # but lower accuracy. [1.5]
            '-r', '1.5',
            # Discard a MEM if it has more than INT occurence in the genome.
            # This is an insensitive parameter. [10000]
            '-c', '10000',
            # Matching score. [1]
            '-A', '1',
            # Mismatch penalty. The sequence error rate is approximately:
            # {.75 * exp[-log(4) * B/A]}. [4]
            '-B', '2',
            # Gap open penalty. [6]
            '-O', '3',
            # Gap extension penalty. A gap of length k costs O + k*E
            # (i.e. -O is for opening a zero-length gap). [1]
            '-E', '1',
            # Clipping penalty. SW extension: best score reaching the end
            # of query is larger than best SW score minus the clipping
            # penalty, clipping not be applied. [5]
            '-L', '5',
            # Penalty for an unpaired read pair .....applied without -P?
            # set to zero here? [9]
            #-U 9
            # Assume the first input query file is interleaved paired-end FASTA/Q.
            #-p
            # Complete read group header line.
            #-R STR
            # Don't output alignment with score lower than INT. This option
            # only affects output. [30]
            '-T', '30',
            # Output all found alignments for single-end or unpaired paired-end
            # reads. These alignments will be flagged as secondary alignments.
            '-a',
            # Append append FASTA/Q comment to SAM output. This option can be
            # used to transfer read meta information (e.g. barcode) to the
            # SAM output.
            # -C
            # Use hard clipping 'H' in the SAM output. This option may
            # dramatically reduce the redundancy of output when mapping long
            # contig or BAC sequences.
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
        cmd = [self.exe_samtools, 'view', '-bS', _os.path.sep.join(
                [local_repeats_path, sam_name])]
        print('Called: {}\nPiping to {}'.format(' '.join(cmd), bam_name))
        with open(_os.path.sep.join([local_repeats_path, bam_name]), "wb") as out:
            _subprocess.call(cmd, stdout = out)

        # where does output end up (which path?)
        cmd = [self.exe_samtools, 'sort', 
                _os.path.sep.join([local_repeats_path, bam_name]), 
                _os.path.sep.join([local_repeats_path, sam_name[:-4]])]
        print('Called: {}'.format(' '.join(cmd)))
        _subprocess.call(cmd)

        bam_name_sorted = sam_name.replace('.sam', '.bam')
        cmd = [self.exe_samtools, 'index', 
                _os.path.sep.join([local_repeats_path, bam_name_sorted])]
        print('Called: {}'.format(' '.join(cmd)))
        _subprocess.call(cmd)

        if not retain_bwa_output:
            f = _os.path.sep.join([local_repeats_path, sam_name])
            print('Removing: {}'.format(f))
            _os.unlink(f)
            f = _os.path.sep.join([local_repeats_path, bam_name])
            print('Removing: {}'.format(f))
            _os.unlink(f)

        #### parse and analyse ORF-chromosome alignment ####
        print('Parsing and analysing the ORF-chromosome alignments . . .')
        self.parseBamForDoubleHits(_os.path.sep.join([local_repeats_path, 
                bam_name_sorted]))

        # there's something curious here . . . what and why is this happening?
        # something to do with one-way matches being made two-way (symmetrical)
        self.followHitsAndAdd()
        max_iterations = 10
        i = 0
        while len(self.ORF_hits) != len(set([a for b in self.ORF_hits.values() \
                for a in b])):
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
        # min_pID is used to find an initial set of repetitive and possibly 
        # divergent regions and should be lower than the eventual minimum
        # threshold: different to minimum_percent_identity which is used for final
        # selection of ambiguous regions below
        pre_min_pID = minimum_percent_identity * 0.85
        self.align_blocks(min_pID = pre_min_pID, max_extensions = max_extensions)
        self.map_alignments_to_chromosome()

        #### identify 98% (default) identical regions
        self.identify_ambiguous_regions(minimum_percent_identity = minimum_percent_identity)
        # repeats of unit length less than the insert size of the paired end 
        # sequenced fragments should be resolvable so omit them
        # (although tandem repeats of total length > insert size
        # may still be a problem . . .)
        self.merge_ambiguous_regions(minimum_repeat_length = minimum_repeat_length)

    def findRepeatsNucmer(self, minimum_percent_identity = 0.98, 
                          minimum_repeat_length = 400,
                          exe_nucmer = False, 
                          local_repeats_path = ['repeats'], 
                          local_genomes_path = ['genome_sequences'], 
                          force = False):
        '''
        Find repeats!
        Use the very fast nucmer method. This approach is not aware of protein coding regions.
        Initial assignment and alignment of homologous blocks requires >=95% nucleotide identity
        Final selection of ambiguous repeats for filtering selectes regions >=98% identity
        '''


        # collect paths to executables and data

        exe_nucmer = False

        def get_exe(name):
            from baga.Dependencies import dependencies as _dependencies
            exe = []
            exe += [_dependencies[name]['destination']]
            exe += _dependencies[name]['checker']['arguments']['path']
            exe = _os.path.sep.join(exe)
            return(exe)

        self.exe_nucmer = get_exe('mummer')

        local_genomes_path = _os.path.sep.join(local_genomes_path)

        genome_fna = '{}/{}.fna'.format(local_genomes_path, self.genome.id)

        try:
            _os.makedirs(local_genomes_path)
        except OSError:
            pass

        print('Writing genome to FASTA')
        _SeqIO.write(_SeqRecord(_Seq(self.genome.sequence.tostring()), 
                id = self.genome.id), genome_fna, 'fasta')


        # nucmer to self
        cmd = [self.exe_nucmer, '--maxmatch','--nosimplify','--prefix={}'\
                ''.format(self.genome.id), genome_fna, genome_fna]
        print('Called: {}'.format(' '.join(cmd)))
        try:
            _subprocess.call(cmd)
        except OSError:
            print('Problem running nucmer at {}. Please use Dependencies module '\
                    'to install locally or check system path'.format(cmd[0]))

        # collect coordinates
        exe = self.exe_nucmer.replace('nucmer','show-coords')
        cmd = [exe, '-r', '{}.delta'.format(self.genome.id)]
        coords = _subprocess.check_output(cmd)

        #print('these are the coords')
        #print(coords)
        #(minimum_percent_identity, minimum_repeat_length)

        # parse coords: simply retain ordinates of regions affected by repeats
        # for comparison with more detailed baga method
        all_positions = set()
        for n,line in enumerate(coords.split('\n')):
            #print(line)
            cells = line.split()
            if n > 3 and len(cells) > 5:
                s1,e1,s2,e2,l1,l2 = map(int,[cells[0], cells[1], cells[3], 
                        cells[4], cells[6], cells[7]])
                if s2 > e2:
                    # reverse strand
                    s2,e2 = e2,s2
                pid = float(cells[9])
                if minimum_percent_identity < (pid/100):
                    Apositions = set(range(s1-1,e1))
                    Bpositions = set(range(s2-1,e2))
                    A_nonself = Apositions - Bpositions
                    B_nonself = Bpositions - Apositions
                    x = False
                    # minimum_repeat_length should really be checked against
                    # contiguous lengths, not total bp which might now be non-
                    # contiguous after set subtractions above
                    #if len(A_nonself) >= minimum_repeat_length:
                    if len(A_nonself):
                        all_positions.update(A_nonself)
                        x = True
                    #if len(B_nonself) >= minimum_repeat_length:
                    if len(B_nonself):
                        all_positions.update(B_nonself)
                        x = True
                    if x:
                        pass
                        #print(line.rstrip())

        if len(all_positions) > 0:
            all_ambiguous_ranges = self.makeRanges(sorted(all_positions))
            initial_num_repeats = len(all_ambiguous_ranges)
            all_ambiguous_ranges = [(s,e) for s,e in all_ambiguous_ranges if \
                    e - s > minimum_repeat_length]
            s = sum([(e - s) for s,e in all_ambiguous_ranges])
            print('Dropped {} repeats less than {} basepairs; {} remain '\
                    'spanning {:,} basepairs'.format(initial_num_repeats, 
                    minimum_repeat_length, len(all_ambiguous_ranges), s))
            self.ambiguous_ranges = all_ambiguous_ranges
        else:
            print('No repeats found.')
            self.ambiguous_ranges = []

        self.homologous_groups_mapped = 'not implemented'


        #print(self.ambiguous_ranges)

        ### could adapt some of below . . .
        #### get contiguous blocks and their homologs (other blocks) ####
        #self.orderORFHits()
        #self.getHomologousContiguousBlocks()

        #### align contiguous homologous blocks ####
        #self.getORFsInTandemRepeats()
        # min_pID is used to find an initial set of repeatitive and possibly 
        # divergent regions and should be lower than the eventual minimum
        # threshold: different to minimum_percent_identity which is used for final
        # selection of ambiguous regions below
        #self.align_blocks(min_pID = 0.85, max_extensions = 15)
        #self.map_alignments_to_chromosome()

        #### identify 98% (default) identical regions
        #self.identify_ambiguous_regions(minimum_percent_identity = minimum_percent_identity)
        # repeats of unit length less than the insert size of the paired end 
        # sequenced fragments should be resolvable so omit them
        # (although tandem repeats of total length > insert size
        # may still be a problem . . .)
        #self.merge_ambiguous_regions(minimum_repeat_length = minimum_repeat_length)



    def compareRepeatRegions(self, minimum_repeat_length = 400):
        '''Compare repeat regions found by nucmer to a previous BAGA analysis'''

        from baga import Repeats

        print('Loading baga repeats analysis of genome %s' % self.genome.id)
        # genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
        # finder = Repeats.Finder(genome)
        filein = 'baga.Repeats.FinderInfo-{}.baga'.format(self.genome.id)
        baga_finder_info = Repeats.loadFinderInfo(filein)

        all_nucmer_positions = set()
        for s,e in self.ambiguous_ranges:
            all_nucmer_positions.update(range(s,e))

        all_baga_positions = set()
        for s,e in baga_finder_info['ambiguous_ranges']:
            all_baga_positions.update(range(s,e))

        print('nucmer positions: {}'.format(len(all_nucmer_positions)))
        print('baga positions: {}'.format(len(all_baga_positions)))

        print('nucmer not baga positions: {}'.format(len(all_nucmer_positions - all_baga_positions)))
        nucmer_not_baga = self.makeRanges(sorted(all_nucmer_positions - all_baga_positions))
        nucmer_not_baga = [(s,e) for s,e in nucmer_not_baga if e-s >= 400]
        print('nucmer not baga positions: {} (filtered at minimum repeat length {}bp)'\
                ''.format(sum((e-s) for s,e in nucmer_not_baga), minimum_repeat_length))

        print('baga not nucmer positions: {}'.format(len(all_baga_positions - all_nucmer_positions)))
        baga_not_nucmer = self.makeRanges(sorted(all_baga_positions - all_nucmer_positions))
        baga_not_nucmer = [(s,e) for s,e in baga_not_nucmer if e-s >= 400]
        print('baga not nucmer positions: {} (filtered at minimum repeat length {}bp)'\
                ''.format(sum((e-s) for s,e in baga_not_nucmer), minimum_repeat_length))


        print('these in nucmer not baga')
        for_csv = []
        for (s,e) in nucmer_not_baga:
            d = e-s
            print('{:,}-{:,} is {:,}bp'.format(s,e,d))
            for_csv += [('"nucmer_not_baga"',s,e,d)]

        print('these in baga not nucmer')
        for (s,e) in baga_not_nucmer:
            d = e-s
            print('{:,}-{:,} is {:,}bp'.format(s,e,d))
            for_csv += [('"baga_not_nucmer"',s,e,d)]

        for s,e in self.ambiguous_ranges:
            d = e-s
            for_csv += [('"nucmer"',s,e,d)]

        for s,e in baga_finder_info['ambiguous_ranges']:
            d = e-s
            for_csv += [('"baga"',s,e,d)]

        for_csv.sort(key = lambda x: x[1])
        for_csv = ['{1}{0}{2}{0}{3}{0}{4}'.format(',',n,s,e,d) for \
                n,s,e,d in for_csv]

        fout_name = 'baga.Repeats.nucmer_vs_baga.csv-{}.csv'.format(
                self.genome.id)
        with open(fout_name, 'w') as fout:
            fout.write('"detection method","start (bp)","end (bp)",'\
                    '"length (bp)"\n')
            fout.write('\n'.join(for_csv))


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
            add(self.genome.ORF_ranges, 'ORF_ranges')  ### to be replaced by "locus" version
            add(self.locus_chrom_info, 'locus_chrom_info')  ### to be updated in plotter
            add(self.genome.large_mobile_element_ranges, 'large_mobile_element_ranges')
            # repeats
            add(self.homologous_groups_mapped, 'homologous_groups_mapped')
            add(self.ambiguous_ranges, 'ambiguous_ranges')
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
        #for ID,(s,e,d,name) in self.finder_info['ORF_ranges'].items():
        for ID,(s,e,d,name) in self.finder_info['locus_chrom_info'].items():
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
        # establish lane plot area (always single for percent ID at duplicate regions)
        plot_spacing = 3
        lane_upper_y = lower_y - (lane_num + 1) * y_per_lane + plot_spacing
        lane_lower_y = lower_y - lane_num * y_per_lane - plot_spacing
        lane_plot_height = lane_lower_y - lane_upper_y

        # plot filled area for percent identity
        colour = ('80', '80', '80', '%')
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

        # add a darker line to aid visualisation through highlights
        linecolour = ('50', '50', '50', '%')
        commands = ['M %s %s' % (
                                    plotpositions_x[0] + plotstart_x, 
                                    lane_lower_y - lane_plot_height * plotdepths[0]
                                    )]

        for n,d in enumerate(plotdepths[1:]):
            # because everything is upside down in SVG minus goes up on page
            commands += ['L %s %s' % (
                                        plotpositions_x[n+1] + plotstart_x,
                                        lane_lower_y - lane_plot_height * d
                                        )]
        # don't end with 'z' because don't want path 'closed'

        plot_path = self.dwg.path(
                            d=' '.join(commands), stroke = self.rgb(*linecolour), 
                            stroke_linecap='round', stroke_width= '2', 
                            fill = 'none'
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
                                        plotstart_x - label_horizontal_offset * 4.5, 
                                        lane_lower_y - (lower_y - upper_y) * 1.4
                                    ), 
                                    fill='black', font_family=use_fontfamily, 
                                    text_anchor=textanc, 
                                    font_size = '{}pt'.format(int(font_size[:-2])*1.4), baseline_shift='-50%')

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
        colour = ('2', '25', '100', '%')
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
                                    d=' '.join(commands), stroke = 'none', #self.rgb(*colour), 
                                    stroke_linecap='round', stroke_width= '3', 
                                    fill = self.rgb(*colour), fill_rule='evenodd',
                                    fill_opacity = 0.45
                                    )
                
                self.dwg.add(plot_path)

        # first indicate regions with high identity in the current plot
        colour = ('100', '2', '36', '%')
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
                                    d=' '.join(commands), stroke = 'none', #self.rgb(*colour), 
                                    stroke_linecap='round', stroke_width= '3', 
                                    fill = self.rgb(*colour), fill_rule='evenodd',
                                    fill_opacity = 0.36
                                    )
                
                self.dwg.add(plot_path)


def loadFinderInfo(filein):
    finder_info = {}
    with _tarfile.open(filein, "r:gz") as tar:
        for member in tar:
            contents = _StringIO(tar.extractfile(member).read())
            if contents == u'not implemented':
                raise NotImplementedError('attempted to load Repeats Finder '\
                        'object without proper data content ("not implemented")')
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
def summariseRepeats(genome_name, force = False):
    '''Summarise repeats in a csv file'''

    # first load Repeats.Finder analysis results and related information

    from baga import Repeats

    filein = 'baga.Repeats.FinderInfo-{}.baga'.format(genome_name)
    finder_info = Repeats.loadFinderInfo(filein)


    # write to simple summary of regions for parsing
    # baga.Repeats.filter_regions-<genome>.csv
    # write to more detailed summary of regions indicating which repeats link to each other
    # baga.Repeats.details-<genome>.csv

    loaded = False
    try:
        import baga
        ambiguous_ranges = baga.bagaload('baga.Repeats.filter_regions-{}'.format(genome_name))
        loaded = True

    except IOError:
        print('Could not load baga.Repeats.filter_regions-{}'.format(genome_name))
        print('You probably need to find the repeats first with the --find option')

    if loaded:
        total = 0
        foutpath = 'baga.Repeats.filter_regions-{}.csv'.format(genome_name)
        with open(foutpath, 'w') as fout:
            fout.write('"length (bp)","start","end"\n')
            for s,e in ambiguous_ranges:
                print('{:,} bp from {:,} to {:,}'.format(e-s+1, s, e))
                total += e-s+1
                fout.write('{},{}\n'.format(s, e))
        
        print('\n{} repetitive, ambiguous regions spanning {:,} bp of chromosome\n'.format(len(ambiguous_ranges), total))
        print('Simple summary for parsing written to: {}'.format(foutpath))


    # numbering here isn't working . . through each group contained all related pairs . . .
    plotted_last = False
    plot_group = 1
    #for groups in finder.homologous_groups_mapped:
    foutpath = 'baga.Repeats.details-{}.csv'.format(genome_name)
    with open(foutpath, 'w') as fout:
        group_num = 0
        fout.write('"region A start","region A end","region B start","region B end"\n')
        for groups in finder_info['homologous_groups_mapped']:
            
            this_group = []
            
            for pair in groups:
                # first check if this pair includes any of the 
                # final selection of ambiguous regions (without
                # the short regions)
                report_A = False
                report_B = False
                sA, eA = pair['A']['start'], pair['A']['end']
                sB, eB = pair['B']['start'], pair['B']['end']
                for s,e in finder_info['ambiguous_ranges']:
                    if (s < eA and e > sA):
                        report_A = True
                    if (s < eB and e > sB):
                        report_B = True
                
                if report_A and report_B:
                    this_group += [((sA, eA),(sB, eB))]
            
            if len(this_group):
                group_num += 1
                lens = [a for b in this_group for a in b]
                lens = [a[1]-a[0] for a in lens]
                mean_len = int(sum(lens) / float(len(lens)))
                fout.write('Group {}, approx. {}bp:\n'.format(group_num,mean_len))
                for (sA, eA),(sB, eB) in this_group:
                    fout.write('{},{},{},{}\n'.format(sA, eA, sB, eB))
                fout.write('\n')
        
        print('More detailed summary written to: {}'.format(foutpath))

if __name__ == '__main__':
    main()
