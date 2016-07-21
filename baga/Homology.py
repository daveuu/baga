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
Homology module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains methods for a novel, fast and scalably algorithm and 
implementation for delineating families of homologous protein coding genes 
(open reading frames) among groups of genomes.
'''

# stdlib
from baga import _subprocess
from baga import _os
from baga import _pickle
from baga import _gzip
from baga import _tarfile
from baga import _array
from baga import _time
from baga import _glob
from baga import _StringIO
from baga import _BytesIO
from baga import _defaultdict
from baga import _re
from baga import _Counter

import uuid as _uuid
import mmap as _mmap
import random as _random
from itertools import islice as _islice
from concurrent.futures import ProcessPoolExecutor as _ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor as _ThreadPoolExecutor
from operator import itemgetter as _itemgetter

# external Python modules
from Bio import SeqIO as _SeqIO
from Bio.Seq import Seq as _Seq
from Bio.SeqRecord import SeqRecord as _SeqRecord
from roaringbitmap import RoaringBitmap as _rrbitmap
from roaringbitmap import MultiRoaringBitmap as _mrbitmap
import dendropy as _dendropy
import numpy as _np
import lmdb as _lmdb
import struct as _struct
from itertools import chain as _chain
import bcolz as _bcolz #;_bcolz.print_versions()
try:
    import numexpr as _ne
    have_numexpr = True
except ImportError:
    have_numexpr = False

# package specific
from baga import decide_max_processes as _decide_max_processes
from baga import report_time as _report_time
from baga import MetaSample as _MetaSample
from baga import PROGRESS
def main():
    pass
class Finder(_MetaSample):
    '''
    Methods to find groups of homologous protein coding genes from among genomes
    '''

    def __init__(self, analysis_name = 'homology', bagas = [], gffs = [], 
            gbks = [], inherit_from = False, **kwargs):
        module_name = __name__
        super(Finder, self).__init__(analysis_name, module_name, 
                inherit_from = inherit_from, **kwargs)
        # cannot save until have obtained data so set to False initially
        file_name = False
        if inherit_from != 'self':
            # insist on having something in these for instantiation?
            assert sum([len(bagas),len(gbks),len(gffs)]) >= 2, 'you must provide at '\
                    'least two genomes via bagas, gffs or gbks'
            self.bagas = bagas
            self.gbks = gbks
            self.gffs = gffs
            self.name = analysis_name
            self.first_round = True
            # in-memory results
            self.bm_selected_genomes = []
            self.bm_tested_genomes = []
            self.bm_clusters = []
            self.not_recruited_bm_ORFs = []
            self.bm_unassigned_ORFs_by_genome = {}
            # outputs
            self.folders = {}
            # mash sketch files
            self.folders['sketches'] = 'homology__{}__mash_sketches'.format(self.name)
            # temporary whole genome fastas for mash
            self.folders['DNAs'] = 'homology__{}__DNAs'.format(self.name)
            # all amino acids for CD-HIT reduction
            self.folders['AAs'] = 'homology__{}__AAs'.format(self.name)
            # bcolz folders for DNA, AA and RAA sequences in carrays
            self.folders['DNAs_bcolz'] = 'homology__{}__DNAs_bc'.format(self.name)
            self.folders['AAs_bcolz'] = 'homology__{}__AAs_bc'.format(self.name)
            self.folders['RAAs_bcolz'] = 'homology__{}__RAAs_bc'.format(self.name)
            # LMDB folder
            self.folders['metadata'] = 'homology__{}__metadata'.format(self.name)
            # output folder for results
            self.folders['results'] = 'homology__{}__results'.format(self.name)

    ## dictionary to convert codons into amino acids, an integer encoded, 
    ## compressed alphabet and a slightly reduced amino acid with higher
    ## specificity

    # this one is based on:
    # http://www.ingentaconnect.com/content/ben/cbio/2012/00000007/00000002/art00007
    # modified to use 10 slots:
    # [('G',), ('P',), ('I', 'V'), ('F', 'Y', 'W'), ('A', 'L', 'M'),
    # ('E', 'Q', 'R', 'K'), ('N', 'D'), ('H', 'S'), ('T',), ('C')]

    translate = {
    'GGT':('G',0),'GGC':('G',0),'GGA':('G',0),'GGG':('G',0),
    'CCT':('P',1),'CCC':('P',1),'CCA':('P',1),'CCG':('P',1),
    'ATC':('I',2),'ATT':('I',2),'ATA':('I',2),
    'GTT':('V',2),'GTC':('V',2),'GTA':('V',2),'GTG':('V',2),
    'TAT':('Y',3),'TAC':('Y',3),
    'TTT':('F',3),'TTC':('F',3),
    'TGG':('W',3),
    'GCT':('A',4),'GCC':('A',4),'GCA':('A',4),'GCG':('A',4),
    'TTA':('L',4),'TTG':('L',4),'CTT':('L',4),'CTC':('L',4),'CTA':('L',4),'CTG':('L',4),
    'ATG':('M',4),
    'GAA':('E',5),'GAG':('E',5),
    'CAA':('Q',5),'CAG':('Q',5),
    'CGT':('R',5),'CGC':('R',5),'CGA':('R',5),'CGG':('R',5),'AGA':('R',5),'AGG':('R',5),
    'AAA':('K',5),'AAG':('K',5),
    'AAT':('N',6),'AAC':('N',6),
    'GAT':('D',6),'GAC':('D',6),
    'CAT':('H',7),'CAC':('H',7),
    'TCT':('S',7),'TCC':('S',7),'TCA':('S',7),'TCG':('S',7),'AGT':('S',7),'AGC':('S',7),
    'ACT':('T',8),'ACC':('T',8),'ACA':('T',8),'ACG':('T',8),
    'TGT':('C',9),'TGC':('C',9),
    #'ATG':('-',-1),
    'TAA':('*',-2),'TGA':('*',-2),'TAG':('*',-2)}

    stop_codons = {
            '+':{b'TAA',b'TGA',b'TAG'},
            '-':{b'TTA',b'TCA',b'CTA'}}

    # for codon,(AA,C) in translate.items():
        # print(Seq(codon).translate(),AA,C)
        # assert str(Seq(codon).translate()) == AA

    # complement DNA including codes with redundancy
    complement = {
    'A':'T','T':'A','G':'C','C':'G',
    'R':'Y','Y':'R','K':'M','M':'K',
    'S':'S','W':'W'}

    def load_gff(self, gff_file):
        '''Extract genome sequence and ordinates from gff files'''
        prefasta = True
        ORFs = _defaultdict(dict)
        last_header = tuple()
        this_fasta = []
        DNAseqs = {}
        for line in open(gff_file):
            #if line.startswith('# Sequence Data:'):
            #    seqinfo = dict([b.split('=') for b in line.rstrip().split(';')])
            #    seqnames += [seqinfo['seqhdr'].strip('"')]
            if '##FASTA' in line:
                prefasta = False
            if prefasta and not line.startswith('#'):
                bits = line.rstrip().split('\t')
                otherbits = bits[8].split(';')
                # prodigal only?
                locus_ID = otherbits[0][3:]
                partial = otherbits[1].split('=')[1]
                if partial == '00':
                    seq_ID = bits[0]
                    s = int(bits[3])
                    e = int(bits[4])
                    strand = bits[6]
                    # need to omit stop codon while loading
                    if strand == '+':
                        ORFs[seq_ID][locus_ID] = s-1,e-3,strand
                    else:
                        ORFs[seq_ID][locus_ID] = s-1+3,e,strand
            else:
                if line.startswith('>'):
                    if len(this_fasta):
                        DNAseqs[last_header] = _np.fromiter(''.join(this_fasta), 
                                dtype='|S1')
                        this_fasta = []
                    last_header = tuple(line[1:].rstrip().split(' ',1))
                elif not line.startswith('#'):
                    this_fasta += [line.rstrip()]
        DNAseqs[last_header] = _np.fromiter(''.join(this_fasta),dtype='|S1')
        return(DNAseqs,dict(ORFs))

    def load_baga(self, baga_file):
        '''load a baga.CollectData.Genome file ready for baga.Homology analysis'''
        from baga import CollectData
        path_bits = baga_file.split(_os.path.sep)
        sample_id = baga_file.split('-')[-1][:-5]
        working_dir = _os.path.realpath(_os.path.curdir)
        if len(path_bits) > 1:
            genome_dir = _os.path.sep.join(path_bits[:-1])
            _os.chdir(genome_dir)
        
        newgenome = CollectData.Genome(sample_id, inherit_from = 'self')
        if len(path_bits) > 1:
            _os.chdir(working_dir)
        
        ## need to update baga to include things like 
        ## chromosome_ID,chromosome_name per chromosome/plasmid/contig?
        DNAseqs = {}
        for n,(seq_ID,seqarray) in enumerate(sorted(newgenome.sequence.items())):
            DNAseqs[(seq_ID,'mol_{}'.format(n))] = _np.frombuffer(seqarray, 
                    dtype='S1')[::4]
        
        ORFs = {}
        for seq_ID,info in newgenome.annotations.items():
            these_ORFs = {}
            for ORF_ID,(s,e,strand,name) in info[0].items():
                # need to omit stop codon from bagas
                if strand == -1:
                    these_ORFs[ORF_ID] = s+3,e,'-'  # less efficient changing bagas here?
                else:
                    these_ORFs[ORF_ID] = s,e-3,'+'
            ORFs[seq_ID] = these_ORFs
        
        return(DNAseqs,ORFs)

    def load_gff_old(self, gff_file):
        '''Extract genome sequence and ordinates from gff files'''
        prefasta = True
        ORFs = _defaultdict(dict)
        fasta = []
        seqnames = []
        for line in open(gff_file):
            if line.startswith('# Sequence Data:'):
                seqinfo = dict([b.split('=') for b in line.rstrip().split(';')])
                seqnames += [seqinfo['seqhdr'].strip('"')]
            if '##FASTA' in line:
                prefasta = False
            if prefasta and not line.startswith('#'):
                bits = line.rstrip().split('\t')
                otherbits = bits[8].split(';')
                # prodigal only?
                ID = otherbits[0][3:]
                partial = otherbits[1].split('=')[1]
                if partial == '00':
                    mol = bits[0]
                    s = int(bits[3])
                    e = int(bits[4])
                    strand = bits[6]
                    ORFs[mol][ID] = s-1,e,strand
            else:
                fasta += [line]
        
        data = _StringIO()
        data.writelines(fasta)
        data.seek(0)
        seqs_DNA_arrays = {g.id:_array('u',g.seq) for g in \
                _SeqIO.parse(data,'fasta')}
        return(seqs_DNA_arrays,dict(ORFs),seqnames)

    def load_gbk(self, gbk_file):
        ORFs = _defaultdict(dict)
        seqs_DNA_arrays = {}
        seqnames = []
        for g in SeqIO.parse(gbk_file,'genbank'):
            seqs_DNA_arrays[g.id] = _array('u',g.seq)
            for f in g.features:
                if f.type == 'source':
                    seqnames += [f.qualifiers['organism'][0]]
                if f.type == 'CDS':
                    s,e = f.location.nofuzzy_start,f.location.nofuzzy_end
                    strand = f.strand
                    ORFs[g.id][f.qualifiers['locus_tag'][0]] = s,e,strand
        
        return(seqs_DNA_arrays,dict(ORFs),seqnames)

    def load_baga_old(self, baga_file):
        '''load a baga.CollectData.Genome file ready for baga.Homology analysis'''
        from baga import CollectData
        path_bits = baga_file.split(_os.path.sep)
        sample_id = baga_file.split('-')[-1][:-5]
        working_dir = _os.path.realpath(_os.path.curdir)
        if len(path_bits) > 1:
            genome_dir = _os.path.sep.join(path_bits[:-1])
            _os.chdir(genome_dir)
        
        newgenome = CollectData.Genome(sample_id, inherit_from = 'self')
        if len(path_bits) > 1:
            _os.chdir(working_dir)
        
        seqs_DNA_arrays = newgenome.sequence
        ORFs = {}
        for seq_ID,info in newgenome.annotations.items():
            these_ORFs = {}
            for ORF_ID,(s,e,strand,name) in info[0].items():
                # need to omit stop codon from bagas
                if strand == -1:
                    these_ORFs[ORF_ID] = s+3,e,'-'  # less efficient changing bagas here?
                else:
                    these_ORFs[ORF_ID] = s,e-3,'+'
            ORFs[seq_ID] = these_ORFs
        
        return(seqs_DNA_arrays,ORFs,[newgenome.organism])

    def saveArrayDict(self, ArrayDict, filename):
        with _tarfile.open(filename, "w:gz", compresslevel = 3) as tar:
            for arr_name, arr in ArrayDict.items():
                ioob = _BytesIO(arr.tobytes())
                ioob.seek(0, _os.SEEK_END)
                length = ioob.tell()
                ioob.seek(0)
                thisone = _tarfile.TarInfo(name = arr_name)
                thisone.size = length
                tar.addfile(tarinfo = thisone, fileobj = ioob)

    def loadArrayDict(self, filename, datatype = 'H'):
        ArrayDict = {}
        with _tarfile.open(filename, "r:gz") as tar:
            for member in tar:
                contents = _BytesIO(tar.extractfile(member).read())
                ArrayDict[member.name] = _array(datatype, contents.getvalue())
        return(ArrayDict)

    def _make_folders(self):
        for folder,path in self.folders.items():
            try:
                _os.makedirs(path)
            except FileExistsError:
                self.logger.debug('Folder of {} already exists: {}'.format(
                        folder, path))
                pass


    def verify_files_restored(self, checkfiles = [('sketches','msh')]):
        missing = []
        for folder,ext in checkfiles:
            for genome_num in self.ORFs_per_genome:
                try:
                    filename = '{}/{:0{}d}.{}'.format(self.folders[folder], 
                            genome_num, self.padding, ext)
                    open(filename)
                except FileNotFoundError:  # PermissionDeniedError? IOError as e; e.errornum ?
                    missing += [filename]
        
        return(missing)


    def _transfer_ORFs2fasta(self, ORF_IDs, sourcefile, destinationfile):
        '''Write a new fasta given an iterable of sequence IDs and a fasta'''
        # collect file positions of all ORFs
        ORF_ranges = {}
        ID = 'delme'
        s = 0
        e = 0
        with open(sourcefile) as fasta:
            i = 0
            for line in fasta:
                if line.startswith('>'):
                    e = i
                    ORF_ranges[ID] = s,e
                    i += len(line)
                    s = i
                    ID = line[1:].rstrip()
                else:
                    i += len(line)
        ORF_ranges[ID] = s,i
        # read via memory map and write them to new file
        with open(sourcefile, "r+b") as fin, open(destinationfile, 'w') as fout:
            # memory-map the file, size 0 means whole file
            mm = _mmap.mmap(fin.fileno(), 0)
            # read content via slice notation
            for ORF in ORF_IDs:
                s,e = ORF_ranges[ORF]
                fout.write('>{}\n{}'.format(ORF,mm[s:e].decode('utf-8')))

    ### write ORF family sequences to fasta files
    def _write_clusters(self, cluster_IDs, cluster_members, cl_locus_IDs, 
            cl_ranges, sequence_array, ORF_padding, file_ext = 'fasta', 
            chars_per_line = 50):
        '''
        Write a fasta file per cluster from sequence carrays
        
        Parameters
        ----------
        cluster_IDs     : list
            Names (IDs) of each cluster
        cluster_members : list
            lists of tuples for ORF members
        cl_locus_IDs    : list
            lists of strings of input locus IDs
        cl_ranges       : list
            lists of tuple of on-disk bcolz carray slices
        sequence_array  : list
            the on-disk bcolz carray to take sequence from
        file_ext        : str
            Fasta file extension e.g., 'faa' or 'fna'
        '''
        for cID, ORF_IDs, locus_IDs, ranges in zip(cluster_IDs, 
                cluster_members, cl_locus_IDs, cl_ranges):
            filename = '{}/{}.{}'.format(self.folders['results'], cID, 
                    file_ext)
            seqs = [_np.unicode(sequence_array[s:e], 'utf-8') for s,e in ranges]
            self.logger.debug('Writing {} sequences to {}'.format(len(seqs), filename))
            records = []
            for (g,o), locus_ID, seq in \
                    zip(ORF_IDs, locus_IDs, seqs):
                ORF_ID = '{:0{}d}.{:0{}d}'.format(g, self.padding, 
                        o, ORF_padding)
                num_lines = round(len(seq)/chars_per_line+0.5)
                seqstr = '\n'.join([seq[n*50:(n+1)*50] for n in \
                        range(num_lines)])
                records += ['>{} {}\n{}'.format(ORF_ID, locus_ID, 
                        seqstr)]
            
            open(filename, 'w').write('\n'.join(records))

    def rolling_window(self, array, k):
        shape = array.shape[:-1] + (array.shape[-1] - k + 1, k)
        strides = array.strides + (array.strides[-1],)
        return(_np.lib.stride_tricks.as_strided(array, shape = shape, 
                strides = strides))

    def call_orfs_prodigal(self, seqs_DNA_arrays):
        from baga.Dependencies import dependencies
        d = dependencies['prodigal']['default_source']
        exe_prodigal = _os.path.join(dependencies['prodigal'][d]['exe']['main'])
        filename = '{}/{}.fna'.format(self.folders['DNAs'], _uuid.uuid4())
        while _os.path.exists(filename):
            # if parallel prepareORFs()
            # (remote chance of race condition: try?)
            filename = '{}/{}.fna'.format(self.folders['DNAs'], _uuid.uuid4())
        chars_per_line = 50
        with open(filename, 'w') as fout:
            # add longest first
            for seqID,seq in sorted(seqs_DNA_arrays.items(), 
                    key = lambda x:len(x[1]), reverse = True):
                num_lines = int(len(seq)/chars_per_line+0.5)
                seqstr = '\n'.join([seq[n*50:(n+1)*50].tounicode() for n in \
                        range(num_lines)])
                fout.write('>{}\n{}\n'.format(seqID, seqstr))
        
        gff_file = filename[-4:]+'.gff'
        cmd = [exe_prodigal, '-i', filename, '-c', '-f', 'gff', '-o', gff_file]
        self.logger.log(PROGRESS, 'Calling ORFs with prodigal: "{}"'.format(
                ' '.join(cmd)))
        proc_returncode,stdout = self._launch_external(cmd, 'prodigal', 
                self.log_folder, main_logger = self.logger)
        ORFs = _defaultdict(dict)
        for line in open(gff_file):
            if not line.startswith('#'):
                bits = line.rstrip().split('\t')
                otherbits = bits[8].split(';')
                # prodigal only?
                ID = otherbits[0][3:]
                partial = otherbits[1].split('=')[1]
                if partial == '00':
                    mol = bits[0]
                    s = int(bits[3])
                    e = int(bits[4])
                    strand = bits[6]
                    ORFs[mol][ID] = s-1,e,strand
        
        _os.unlink(filename)
        _os.unlink(gff_file)
        
        return(dict(ORFs))

    def _union_find(self, data):
        '''Create Disjoint Data Structure from list of lists'''
        parents = {}
        def find(i):
            j = parents.get(i, i)
            if j == i:
                return i
            k = find(j)
            if k != j:
                parents[i] = k
            return k
        for l in filter(None, data):
            parents.update(dict.fromkeys(map(find, l), find(l[0])))
        merged = {}
        for k, v in parents.items():
            merged.setdefault(find(v), []).append(k)
        return list(merged.values())

    # Sorts child nodes ascending by the number of children each child node has.
    # Sort by edge lengths to better capture phylogenetic diversity? <===
    def _ladderize_diversity(self, tree, ascending=True):
        """
        Sorts child nodes in ascending (if ``ascending`` is `False`) or
        descending (if ``ascending`` is `False`) order in terms of edge
        lengths descending from each child node.
        """
        node_desc_counts = {}
        for nd in tree.postorder_node_iter():
            if len(nd._child_nodes) == 0:
                if nd.edge_length is None:
                    node_desc_counts[nd] = 0
                else:
                    node_desc_counts[nd] = nd.edge_length
            else:
                if nd.edge_length is not None:
                    total = nd.edge_length
                else:
                    total = 0
                for child in nd._child_nodes:
                    total += node_desc_counts[child]
                    if child.edge_length is not None:
                        total += child.edge_length
                node_desc_counts[nd] = total
                nd._child_nodes.sort(key=lambda n: node_desc_counts[n], 
                        reverse=not ascending)

    def _collect_cluster_info(self, clusters):
        '''
        Query LMDB ORF metadata on a per cluster basis
        
        Return an array of broad cluster membership, and minimum and
        maximum ORF lengths in nucleotides with a single ORF as a 
        representative.
        
        Parameters
        ----------
        clusters : list
            List of lists of ORFs
        '''
        # info for parsing LMDB metadata database
        ORFs_s_k = self.ORFs_struct_key
        ORFs_s_v = self.ORFs_struct_value
        AA_s_sl = self.ORFs_struct_slices['AA']
        AA_struct = ORFs_s_v[AA_s_sl]
        AA_b_sl = self.ORFs_bytes_slices['AA']
        BC_s_sl = self.ORFs_struct_slices['broad_cluster_num']
        BC_struct = ORFs_s_v[BC_s_sl]
        BC_b_sl = self.ORFs_bytes_slices['broad_cluster_num']
        # store cluster info and a representative for this group
        cluster_info_by_ORFrep = _np.empty((len(clusters),5), dtype = _np.uint32)
        # keep a note of ranges for when eventually collecting RAAs from
        # contiguous carrays
        cluster_carray_ranges = {}
        # could have a single context for lmdb_env outside loops?
        # query LMDB for ORF lens and IDs
        with _lmdb.Environment(self.folders['metadata'], readonly = True, 
                create = False, max_dbs = 3) as lmdb_env:
            ORFs_metadata = lmdb_env.open_db(b'ORFs')
            with lmdb_env.begin(db = ORFs_metadata) as txn:
                for c,cluster in enumerate(clusters):
                    ORFs = sorted(cluster)
                    raw = [txn.get(_struct.pack(ORFs_s_k, g, o)) \
                            for g,o in ORFs]
                    ranges = {ORF:_struct.unpack(AA_struct, data[AA_b_sl]) \
                            for ORF,data in zip(ORFs,raw)}
                    broad_clusters = {ORF:_struct.unpack(BC_struct, data[BC_b_sl])[0] \
                            for ORF,data in zip(ORFs,raw)}
                    cluster_carray_ranges.update(ranges)
                    lengths = [e-s for s,e in ranges.values()]
                    # arbitrary selection of a representative ORF
                    thisORF = ORFs[0]
                    # assert len(set(broad_clusters.values())) == 1, \
                            # 'broad clusters incongruent with a best match group cluster!'
                    # genome, ORF, broad, min length, max length
                    cluster_info_by_ORFrep[c,:] = thisORF[0], thisORF[1], \
                            broad_clusters[thisORF], min(lengths), \
                            max(lengths)
        
        return cluster_info_by_ORFrep, cluster_carray_ranges

    def _collect_cluster_ranges(self, clusters):
        '''
        Query LMDB database per cluster for sequence ranges and loci ID
        
        Input must be list of lists of ORFs, info returned in lists of same 
        order. Return lists of ranges in AA and DNA carrays and input loci.
        
        Parameters
        ----------
        clusters : list
            List of lists of ORFs
        '''
        # info for parsing LMDB metadata database
        ORFs_s_k = self.ORFs_struct_key
        ORFs_s_v = self.ORFs_struct_value
        AA_s_sl = self.ORFs_struct_slices['AA']
        AA_struct = ORFs_s_v[AA_s_sl]
        AA_b_sl = self.ORFs_bytes_slices['AA']
        DNA_s_sl = self.ORFs_struct_slices['DNA']
        DNA_struct = ORFs_s_v[DNA_s_sl]
        DNA_b_sl = self.ORFs_bytes_slices['DNA']
        # last in bytes string is variable length is locus ID
        # get last stop among slices
        locus_ID_start = max(Slice.stop for Slice in self.ORFs_bytes_slices.values())
        # keep a note of ranges for collecting AAs and/or DNAs from
        # contiguous carrays
        DNAranges = []
        AAranges = []
        locus_IDs = []
        # could have a single context for lmdb_env outside loops?
        # query LMDB for ORF lens and IDs
        with _lmdb.Environment(self.folders['metadata'], readonly = True, 
                create = False, max_dbs = 3) as lmdb_env:
            ORFs_metadata = lmdb_env.open_db(b'ORFs')
            with lmdb_env.begin(db = ORFs_metadata) as txn:
                for c,cluster in enumerate(clusters):
                    raw = [txn.get(_struct.pack(ORFs_s_k, g, o)) \
                            for g,o in cluster]
                    AAranges += [[_struct.unpack(AA_struct, data[AA_b_sl]) \
                            for data in raw]]
                    DNAranges += [[_struct.unpack(DNA_struct, data[DNA_b_sl]) \
                            for data in raw]]
                    locus_IDs += [[data[locus_ID_start:].decode('utf-8') \
                            for data in raw]]
        
        return DNAranges, AAranges, locus_IDs


    def _mrb_jaccard(self, mrb, As, Bs, nprocs = 1):
        '''multiroaringbitmap jaccard distances to numpy array with threading'''
        if nprocs == 1:
            results = _np.frombuffer(
                    mrb.jaccard_dist(As,Bs), 
                    dtype = _np.float64).astype(_np.float16)
        elif nprocs > 1:
            prechunksize = int((len(As)+nprocs)/nprocs)
            # slightly more efficient to create separate arrays above?
            As_list = (As[p*prechunksize:(p+1)*prechunksize] \
                    for p in range(nprocs))
            Bs_list = (Bs[p*prechunksize:(p+1)*prechunksize] \
                    for p in range(nprocs))
            # jaccard_dist() has nogil so threading works
            # allows efficient use of memory and cores
            with _ThreadPoolExecutor(max_workers=nprocs) as executor:
                results = _np.concatenate([
                        _np.frombuffer(a, dtype = _np.float64).astype(
                        _np.float16) for a in executor.map(mrb.jaccard_dist, 
                        As_list, Bs_list, chunksize=1)])
        return(results)

    def _collect_tophits(self, from_hits):
        tops = {}           # shortest A to B hit (possibly draws)
        nontops = {}        # other A to B hits but not shortest (ordered and could also be draws)
        zeros = set()       # those shortest A to B hits that are zero distance (no chance of inparalogs)
        for A,Bs in from_hits.items():
            # get increasing unique distances
            incr_distances = sorted(set(Bs.values()))
            closest = incr_distances[0]
            # collect all ORFs with top score (could be a draw for top match)
            tops[A] = sorted([B_ORF for B_ORF,dist in Bs.items() if dist == closest])
            if closest == 0:
                zeros.add(A)
            # collect other ORFs matching increasing scores in order
            these_nontops = []
            for this_dist in incr_distances[1:]:
                these_nontops += [sorted([B_ORF for B_ORF,dist in Bs.items() \
                        if dist == this_dist])]
            if len(these_nontops):
                nontops[A] = these_nontops
        return(tops,zeros,nontops)

    def _build_mrbs(self, all_ORFs_order, batch1_ranges, batch2_ranges, 
            RAAs_carray):
        '''
        Work out each genome range to collect that spans required ORFs in carray.
        
        Whole ORF collection may span many genomes. Instead of slicing the 
        potentially disk-based carray for each ORF, first organise ORFs by genome, 
        then extract the region of each genome required, then slice each ORF from 
        within that in-memory numpy array. This reduces disk IO (one search per 
        genome) without much RAM cost (one partial genome loaded at a time).
        
        Parameters
        ----------
        ORF_ranges  : dict
            ORF ID tuple with tuple of range in the carray of RAAs
        '''
        # for converting array of 0-9 integers into single base-10 integer
        k_multiplier = _np.fromiter([10**n for n in range(self.k-1,-1,-1)], 
                dtype = self.kmer_collections_dtype)
        ORF_ranges = dict(batch1_ranges.items())
        ORF_ranges.update(batch2_ranges.items())
        # find start and end in RAA carray on a per genome basis for IO efficiency
        max_ordinate = max((e for ORF,(s,e) in ORF_ranges.items()))
        all_genomes = sorted(set((g for g,o in ORF_ranges)))
        genome_ranges = {genome:[max_ordinate,0] for genome in all_genomes}
        for (genome,ORF),(s,e) in ORF_ranges.items():
            if genome_ranges[genome][0] > s:
                genome_ranges[genome][0] = s
            if genome_ranges[genome][1] < e:
                genome_ranges[genome][1] = e
        
        # collect bitmpas in a dictionary by ID then build multibitmap
        # in required ORF order
        rb_kmers = {}
        for genome,(genome_s,genome_e) in genome_ranges.items():
            # minimize slices to carray by slicing once per genome
            genome_contiguous_RAAs = RAAs_carray[genome_s:genome_e]
            # then slice ORFs in memory for k-mer counting
            # adjust per ORF ranges to within this slice for this genome
            # (could vectorise this in a numpy array for all genomes?)
            these_ranges = sorted([ORFinfo for ORFinfo in ORF_ranges.items() \
                    if ORFinfo[0][0] == genome])
            new_ranges = (((ORF,(AA_s-genome_s,AA_e-genome_s)) \
                    for ORF,(AA_s,AA_e) in these_ranges))
            # extract k-mers by striding numpy array, straight to roaringbitmap
            # put in dict by key for ordering later
            rb_kmers.update((ORF,_rrbitmap((self.rolling_window(
                    genome_contiguous_RAAs[s:e], self.k)*k_multiplier)\
                    .sum(1, dtype = self.kmer_collections_dtype))) for ORF,(s,e) in new_ranges)
        
        # sorted in order corresponding to pairwise comparisons
        rb_kmers = [rb_kmers[ORF] for ORF in sorted(rb_kmers)]
        mrb = _mrbitmap(rb_kmers)
        
        return mrb

    def _compare_ORF_pairs(self, mrb, pairwise_test_inds, pairwise_test_IDs, 
            jacc_dist_max, nprocs = 1):
        '''do jaccard distances retaining only hits'''
        self.logger.debug('Compiling pairs for calculations')
        self.logger.debug('Calculating {:,} jaccard distances in bitmap'\
                ''.format(len(pairwise_test_inds[0])))
        b1_to_b2_dists = self._mrb_jaccard(mrb, pairwise_test_inds[0], 
                pairwise_test_inds[1], nprocs = 1)
        # not all v all so always fewer distances than ORFs in mrb
        self.logger.debug('Parsing {:,} jaccard distances for a pair of genome '\
                'groups'.format(len(pairwise_test_inds[0])))
        
        # retain only "hits" below threshold distance (would need calibrating
        # to AA pID distance)
        if nprocs == 1 or not have_numexpr:
            hits = b1_to_b2_dists < jacc_dist_max
        else:
            # numexpr quicker on 2 or more cores
            hits = _ne.evaluate('b1_to_b2_dists < jacc_dist_max')
        
        # results retaining only hits
        b1_to_b2_dists = b1_to_b2_dists[hits]
        
        # used in _parse_distances_per_genome() to get genome specific entries
        # These indicate which items in batch 1 and batch 2 correspond to
        # which genomes also contains. (make a single array?)
        # makes a 2 by x array from list of 2-tuples
        pairwise_hit_IDs = []
        pairwise_hit_IDs += [_np.array(pairwise_test_IDs[0], dtype = _np.uint32)[hits]]
        pairwise_hit_IDs += [_np.array(pairwise_test_IDs[1], dtype = _np.uint32)[hits]]
        
        return b1_to_b2_dists, pairwise_hit_IDs


    def _parse_distances(self, b1_to_b2_dists, pairwise_hit_IDs):
        '''
        Extract per ORF distances to other batch of ORFs from lists
        
        From a list of distances with lists of pairs create dictionaries of 
        distances. Genome of ORF ignored: only input batches matter.
        
        Parameters
        ----------
        b1_to_b2_dists : numpy.array
            distance between pairs of ORFs
        pairwise_hit_IDs : list
            of 1d numpy.array ORF IDs corresponding to distances
        
        removed:
        b2_ORFs_in_mrb : numpy.array
            indexes in bitmap of ORFs of second member of pair
        ORF_IDs_in_mrb : numpy.array
            ORF IDs in same order as in bitmap
        '''
        # record the hits and their scores
        b1_hits = {}
        b2_hits = {}
        if pairwise_hit_IDs[0].shape[0]:
            # only proceed if there were any hits
            # iterate through each (unique) batch 1 ORF cluster rep
            # probably a numpy quicker way?
            ORFs_in_b1 = _np.vstack({tuple(ORF) for ORF in pairwise_hit_IDs[0]})
            for b1_cluster_rep in ORFs_in_b1:
                # get scores for this b1 cluster rep against all b2 cluster reps
                # indexes of distance featuring this ORF
                these_hit_inds = _np.all(pairwise_hit_IDs[0]==b1_cluster_rep, axis=1)
                # scores featuring this ORF
                hit_scores = b1_to_b2_dists[these_hit_inds]
                hit_IDs = pairwise_hit_IDs[1][these_hit_inds]
                these_hits = {tuple(ID):s for ID,s in zip(hit_IDs,hit_scores)}
                b1_hits[tuple(b1_cluster_rep)] = these_hits
            # iterate through each (unique) batch 2 ORF cluster rep
            ORFs_in_b2 = _np.vstack({tuple(ORF) for ORF in pairwise_hit_IDs[1]})
            for b2_cluster_rep in ORFs_in_b2:
                # get scores for this b1 cluster rep against all b2 cluster reps
                these_hit_inds = _np.all(pairwise_hit_IDs[1]==b2_cluster_rep, axis=1)
                hit_scores = b1_to_b2_dists[these_hit_inds]
                hit_IDs = pairwise_hit_IDs[0][these_hit_inds]
                these_hits = {tuple(ID):s for ID,s in zip(hit_IDs,hit_scores)}
                b2_hits[tuple(b2_cluster_rep)] = these_hits
        else:
            self.logger.debug('no hits to process in _parse_distances()')
        return b1_hits, b2_hits


    ## remove if not needed
    # def _parse_distances_per_genome(self, b1_to_b2_dists, pairwise_tests, 
                # b1_ORFs_in_mrb, genome1, b2_ORFs_in_mrb, genome2, 
                # ORF_IDs_in_mrb):
        # '''
        # Extract per ORF distances to other genome from lists
        
        # From a list of distances with lists of pairs and genomes of interest
        # create dictionaries of distances.
        
        # Parameters
        # ----------
        # b1_to_b2_dists
        # pairwise_tests
        # b1_ORFs_in_mrb
        # genome1
        # b2_ORFs_in_mrb
        # genome2
        # ORF_IDs_in_mrb
        # '''
        # # ORF_IDs_in_mrb
        # ORF2BM = {}
        # dists_these_genomes = (pairwise_tests[0][:,0]==genome1) & (pairwise_tests[1][:,0]==genome2)
        # if not any(dists_these_genomes):
            # return {}, {}, ORF2BM
        # # get distances for this genome combination
        # these_results = b1_to_b2_dists[dists_these_genomes]
        # # get indexes in multibitmap input for this genome combination
        # these_b1_inds = b1_ORFs_in_mrb[dists_these_genomes]
        # these_b2_inds = b2_ORFs_in_mrb[dists_these_genomes]
        # # iterate through each (unique) batch 1 ORF cluster rep
        # # record the hits and their scores
        # b1_hits = {}
        # for b1_cluster_rep in _np.unique(these_b1_inds):
            # # get scores for this b1 cluster rep against all b2 cluster reps
            # these_hits = these_b1_inds==b1_cluster_rep
            # hit_scores = these_results[these_hits]
            # hit_indexes = these_b2_inds[these_hits]
            # hit_IDs = ORF_IDs_in_mrb[hit_indexes,:]
            # these_hits = {tuple(ID):s for ID,s in zip(hit_IDs,hit_scores)}
            # b1_hits[tuple(ORF_IDs_in_mrb[b1_cluster_rep])] = these_hits
            # # retain relevent ID to index mappings for additional jaccard dist calcs
            # ORF2BM.update(((tuple(ID),int(i)) for ID,i in zip(hit_IDs,hit_indexes)))
        # # iterate through each (unique) batch 2 ORF cluster rep
        # # record the hits and their scores
        # b2_hits = {}
        # for b2_cluster_rep in _np.unique(these_b2_inds):
            # # get scores for this b1 cluster rep against all b2 cluster reps
            # these_hits = these_b2_inds==b2_cluster_rep
            # hit_scores = these_results[these_hits]
            # hit_indexes = these_b1_inds[these_hits]
            # hit_IDs = ORF_IDs_in_mrb[hit_indexes,:]
            # these_hits = {tuple(ID):s for ID,s in zip(hit_IDs,hit_scores)}
            # b2_hits[tuple(ORF_IDs_in_mrb[b2_cluster_rep])] = these_hits
            # # retain relevent ID to index mappings for additional jaccard dist calcs
            # ORF2BM.update(((tuple(ID),int(i)) for ID,i in zip(hit_IDs,hit_indexes)))
        # return(b1_hits, b2_hits, ORF2BM)

    def _parse_hits(self, 
            b1_tops, b1_zeros, b1_nontops, b1_hits, 
            b2_tops, b2_zeros, b2_nontops, b2_hits):
        '''
        parse best hits.
        
        these are orthologs plus maybe inparalogs from "tops"
        they might have non-top hit inparalogs for recruitment: check nontops
            may include i) draws identical paralogs at RAA resolution (may have DNA differences)
                             ==> attempt resolution by position later
                       ii) separate orthologous groups with out-paralogous between group relationships
                      iii) single copy families or with inparalogs resolved from among "nontops"
        '''
        these_recip_best_hits = \
                set((tuple((b2s[0],b1)) for b1,b2s in b1_tops.items())) & \
                set((tuple((b2,b1s[0])) for b2,b1s in b2_tops.items()))
        # add to cumulative list
        these_recip_best_hits = [list(pair) for pair in these_recip_best_hits]
        # for below, record which ORFs collected in this call to genome_pw_best_matches <=== currently not returned but maybe needed in Map(genome_pw)
        collected_ORFs = set([a for b in these_recip_best_hits for a in b])
        ### parse non-best hits (make a function to share with genome_pw)
        # some might have reciprocal best hits elsewhere
        # i.e., ensure these potential inparalogs are not known orthologs
        # if the best hit was zero distance, these can only be out-paralogs
        # (singletons, at least in this context)
        # Also merge to single dict
        all_zeros = b1_zeros | b2_zeros
        nontops_retained = {}
        singletons = set()
        for A,nontopBs in _chain(b1_nontops.items(),b2_nontops.items()):
            keeps = []
            for theseBs in nontopBs:
                # ensure ORFs are not best hits elsewhere
                notbests = sorted(set(theseBs) - collected_ORFs)
                if len(notbests):
                    keeps += [notbests]
            if len(keeps):
                if A in all_zeros:
                    # if no distance between orthologs, keep this as an outparalog
                    # i.e., family of one at this stage because no best hits
                    # elsewhere against this other genome
                    # (or two if identical distance i.e., near-identical ORFs)
                    singletons.update((tuple(keep) for keep in keeps))
                else:
                    # else it could be an inparalog
                    # keep to compare paralog distance with ortholog distance
                    nontops_retained[A] = keeps
        
        # add to cumulative list
        if len(singletons):
            self.logger.debug('Added {} outparalogs to near-identical orthologs'.format(
                    len(singletons)))
            these_recip_best_hits += (list(singleton) for singleton in singletons)
        
        ### collect ORFs for testing in or out paralogy after
        # now see if dist to paralogs is shorter than to ortholog
        # if yes: inparalog
        # if no: outparalog
        # if close: phylogenetic test in wider context? <== bootstrap?
        these_para_pairs = []
        for query,hits in nontops_retained.items():
            for n,hit in enumerate(hits):
                # [0] here allows for draws: all have same score so treat the same
                # use this_hit as key but retain hit for all drawn ORFs
                this_hit = hit[0]
                # don't know if query was in genome A or B so try both
                try:
                    top_hit = b1_tops[query][0]
                    top_hit_dist = b1_hits[query][b1_tops[query][0]]
                    this_hit_dist = b1_hits[query][this_hit]
                except KeyError:
                    top_hit = b2_tops[query][0]
                    top_hit_dist = b2_hits[query][b2_tops[query][0]]
                    this_hit_dist = b2_hits[query][this_hit]
                # retain:
                #    query (ortholog 1), 
                #    top_hit (ortholog 2), 
                #    hit (potential inparalog(s) to ortholog 2, include all but test one)
                #    query to top_hit distance <== greater than this for outparalog, less for in-
                ### use ORF2BM here to collect mrb indexes?
                these_para_pairs += [(query,top_hit,hit,top_hit_dist)]
        
        return these_recip_best_hits, these_para_pairs

    def _test_inparalogy(self, mrb, para_pairs, ORF2BM):
        # prepare indices of pairs to measure jaccard distances between
        recip_best_hits = []
        test_b1s = _array('L')
        test_b2s = _array('L')
        ortho_dists = _np.empty((len(para_pairs),), _np.float16)
        for i,(query,top_hit,other_hit,top_hit_dist) in enumerate(para_pairs):
            # not all ORFs per genome retained after CD-HIT so need to look up position
            # and add offset for start of that genome in mrb
            # otherhit is one or more at same distance to ortholog to be tested
            # for in/out paralogy
            # testAs.append(self.bm_group_genome_ranges[top_hit[0]][0] + \
                    # self.bm_group_genome_ORFids[top_hit[0]].index(top_hit))
            # testBs.append(self.bm_group_genome_ranges[top_hit[0]][0] + \
                    # self.bm_group_genome_ORFids[other_hit[0][0]].index(other_hit[0]))
            test_b1s.append(ORF2BM[top_hit])
            test_b2s.append(ORF2BM[other_hit[0]])
            ortho_dists[i] = top_hit_dist
        
        para_pairs = _np.array(para_pairs, dtype = _np.object)
        para_dists = self._mrb_jaccard(mrb, test_b1s, test_b2s, nprocs = 1)
        # possible but unlikely to be same distance
        # if so, should go into same cluster and test at ORF synteny level
        inparalogs = para_dists <= ortho_dists
        
        ### collect inparalogs
        for i,(query,top_hit,other_hit,top_hit_dist) in enumerate(para_pairs[inparalogs]):
            self.logger.debug('for query {0} top hit was {1} at {2:.03f} while this '\
                    'potential inparalog {3} is only {4:.03f} from the top '\
                    'hit ortholog {1}: it is inparalogous to {3}'.format(str(query), 
                    str(top_hit), top_hit_dist, str(other_hit[0]), para_dists[inparalogs][i]))
            
            recip_best_hits += [[top_hit, inparalog] for inparalog in other_hit]
        
        ### collect the outparalogs as singletons in this context
        for i,(query,top_hit,other_hit,top_hit_dist) in enumerate(para_pairs[~inparalogs]):
            recip_best_hits += [other_hit]
        return recip_best_hits

    def processORFs(self, call_orfs = False, do_checks = True, force = False, 
            fail_prop_limit = 0.2, min_ORFs_per_genome = 100, store_to_disk = True, 
            assign_narrow_clusters = True, unclustered_clearence_freq = 15,
            ):
        '''
        Per genome, 
            a) read in all DNA seqs (e.g. genbank, gff, exclude partials) convert 
                    straight to convert to unicode array
            b) write to an fna per genome for mash and sketch (combine later)
            d) convert to AA char and numeric reduced AA alphabet tuples using 
                    single dict and unzip to numeric np arrays and python arrays 
                    pass on KeyError for Ns
            e) write each AA sequence to giant FASTA (try excluding new lines for 
                    now - marginally quicker if CD-Hit doesnt mind)
            f) write a dictionary of reduced AA alphabet numeric np arrays to disk

        Parameters
        ----------
        fail_prop_limit : float
            Above this proportion of failed ORFs, abandon whole genome. Default: 0.2
        min_ORFs_per_genome : int
            Minimum number of ORFs annotated in input, else abandon whole 
            genome. Useful for quality control. Default 100
        unclustered_clearence_freq : int
            how often (number of genomes) single member clusters should be cleared
            from the in-memory store of narrow clusters. Default 15

        store_to_disk : True  <== update here and function for option to store all in memory including carrays
            if True, the LMDB key-value store of ORF metadata will be updated. Use False
            to store ORF broad cluster affiliation in a dictionary in memory. If the whole
            analysis can fit in memory, it is faster not to store_to_disk. To allow larger
            analyses than memory permits, set store_to_disk as True.
        '''

        # max threads is slightly different to cpus
        # . . can probably use more
        # max_processes = _decide_max_processes( max_cpus )

        self._make_folders()

        from baga.Dependencies import dependencies
        d = dependencies['mash']['default_source']
        exe_mash = _os.path.join(dependencies['mash'][d]['exe']['main'])

        genomes_to_process = []
        for baga in self.bagas:
            genomes_to_process += [(baga,self.load_baga)]

        for gff in self.gffs:
            genomes_to_process += [(gff,self.load_gff)]

        for gbk in self.gbks:
            genomes_to_process += [(gff,self.load_gbk)]

        # then store locations of saved, processed files

        padding = len(str(len(self.bagas)+len(self.gffs)+len(self.gbks)))
        # for storage with genome metadata
        genome_num2organism = {}
        # for easy access to total ORFs per genome and generation of ORF IDs
        # because they are simply integer ranges. Total here is base-1 while
        # IDs start at zero
        ORFs_per_genome = {}

        ## threaded?
        # from multiprocessing import Pool
        # def process_a_genome(x):
            # # load data, convert to RAA etc
            # return x*x
        # with Pool(5) as p:
            # print(p.map(process_a_genome, genomes_to_process))
            # # could include sleep(n) at beginning to stagger IO ?
            # http://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments/5443941#5443941
            # might have to fiddle with the logging . . . at least make a note of PID or process num
            # only have progress at beginning and end, debug elsewhere.
            # even something clever like observing whether input file ahead is open for reading or not and wait for it to be closed
            # os.readlink('/proc/<pid>/0-n') or just touch then unlink a lock file while files are open . . .

        def perform_checks( fwd_slices, fwd_strand_mask_A, fwd_strand_mask_B, 
                            rev_slices, rev_strand_mask_A, rev_strand_mask_B, 
                            fwd_slices_A, fwd_slices_B, 
                            rev_slices_A, rev_slices_B):
            # check are codon shaped
            e = 'coding region mask is not codon multiple'
            assert fwd_strand_mask_A.sum()%3 == 0, e
            assert fwd_strand_mask_B.sum()%3 == 0, e
            assert rev_strand_mask_A.sum()%3 == 0, e
            assert rev_strand_mask_B.sum()%3 == 0, e
            # check are same length as input
            assert sum([e-s for ORF,(s,e) in fwd_slices]) == \
                    fwd_strand_mask_A.sum() + \
                    fwd_strand_mask_B.sum(), 'forward strand coding region '\
                    'input != output'
            assert sum([e-s for ORF,(s,e) in rev_slices]) == \
                    rev_strand_mask_A.sum() + \
                    rev_strand_mask_B.sum(), 'reverse strand coding region '\
                    'input != output'
            # check are non-overlapping
            e = 'overlaps remain'
            if fwd_slices_A:
                assert max(_Counter([range(s,e) \
                        for ORF,(s,e),fstrand in fwd_slices_A]).values()) == 1, e
            if fwd_slices_B:
                assert max(_Counter([range(s,e) \
                        for ORF,(s,e),fstrand in fwd_slices_B]).values()) == 1, e
            if rev_slices_A:
                assert max(_Counter([range(s,e) \
                        for ORF,(s,e),fstrand in rev_slices_A]).values()) == 1, e
            if rev_slices_B:
                assert max(_Counter([range(s,e) \
                        for ORF,(s,e),fstrand in rev_slices_B]).values()) == 1, e

        if assign_narrow_clusters:
            # make streaming narrow clusters at genome parse time
            narrows = _defaultdict(list)

        if store_to_disk:
            ### store ORF, chromsome and genome metadata in LMDB ###
            # https://lmdb.readthedocs.org/en/release/
            # create lmdb environment: on-disk file that can contain more
            # than one key-value DB
            # use some risky but reasonable optimisations
            lmdb_env_options = {
                    # 2GB (actually suitable for in-memory on even small machines!)
                    # implemented using LMDB to guarantee future scalability to huge
                    # Later work how many bacterial genomes worth of metadata fits in
                    # 2GB - can be unnecessarily big for safety.
                    'map_size':1024*1024*2000, 
                    # put content in subdir
                    'subdir':True,
                    'readonly':False, 
                    # metasync=False is optimisation,
                    # last thing before a crash is lost: fine here
                    'metasync':False, 
                    'sync':False, 
                    'map_async':False, 
                    # file creation mode . . permission bitmask?
                    'mode':493, 
                    'create':True, 
                    'readahead':True, 
                    # writemap=True another slightly risky optimisation: fine here 
                    # because write once, read many. "Note that sync=False, 
                    # writemap=True leaves the system with no hint for when to 
                    # write transactions to disk, unless sync() is called"
                    # <== do this and sync() after each genome.
                    'writemap':True, 
                    'meminit':True, 
                    'max_readers':16, 
                    # max_dbs=4 Max DBs available. 0 assumes env used as a single database.
                    # Setting to 4 for main + ORFs, chromsomes, genome, additional metadata
                    'max_dbs':4, 
                    'max_spare_txns':1, 
                    'lock':True}
            
            self.lmdb_env_options = lmdb_env_options
            try:
                _os.unlink(_os.path.sep.join([self.folders['metadata'], 'data.mdb']))
            except FileNotFoundError:
                pass
            try:
                _os.unlink(_os.path.sep.join([self.folders['metadata'], 'lock.mdb']))
            except FileNotFoundError:
                pass
            lmdb_env = _lmdb.Environment(self.folders['metadata'], **lmdb_env_options)
            # primary keys are ORFs, genomes and chromosomes for corresponding sub databases 
            # of metadata. 
            # LMDB sorts key in lexicographical order but integers packed into bytes do not
            # sort into base-10 integer order, although zero does come first.
            # Genomes start at one but ORFs start at zero so the first ORF in a genomes can
            # be reliably found by (genome_num, 0).
            # The point of all this is for simple (and fast) iterating of the cursor through
            # the keys to collect e.g., a genome of ORFs at a time by seeking to first and
            # iterating for known number of ORFs
            ORFs_metadata = lmdb_env.open_db(b'ORFs')
            genomes_metadata = lmdb_env.open_db(b'genomes')
            chromosomes_metadata = lmdb_env.open_db(b'chromosomes')
        else:
            ### store ORF, chromsome and genome metadata in dict in memory ###
            genomes_metadata = {}
            chromosomes_metadata = {}
            ORFs_metadata = {}


        ### structs:
        # ORF metadata struct

        ### key: genome_num, ORF_num ###
        self.ORFs_struct_key = 'LH'
        # up to 65k ORFs per genome using H
        # L for >65k genomes (up to >4 billion)

        ### value: ###
        # for interpreting metadata in LMDB
        # (or dicts in memory if store_to_disk = False)

        self.ORFs_struct_value = 'LLLLLLLL?'
        # up to 4 billion total nucleotides via DNA_e as L
        # up to 4 billion broad clusters as L (lot's of rare ORFs?)
        # up to 4 billion chromosomes via chrm_num as L
        # (drafts have 100s contigs leading to >65k cumulative IDs
        # so chrm_num cannot be H)
        ## NB order of struct chunks has to decrease with size? bug?
        ## i.e., an L on end does not work

        # these slices must correspond to the widths listed in the struct module docs
        # https://docs.python.org/3/library/struct.html
        # and the format in ORF_value_struct
        struct_collect_slices = [
                # DNA_s, DNA_e,
                ('DNA',(0,2)),
                # AA_s, AA_e,
                ('AA',(2,4)),
                # chrm_s, chrm_e,
                ('chrm',(4,6)),
                # broad_cluster_num,
                ('broad_cluster_num',(6,7)),
                # chromosome_num,
                ('chrm_num',(7,8)),
                # f_strand,
                ('strand',(8,9)),
                ]

        # generate corresponding slices in structs of key and value for ORFs info
        struct_slices = {}
        bytes_slices = {}
        last_bytes_end = 0
        for n,(s,e) in struct_collect_slices:
            struct_slices[n] = slice(s, e)
            this_bytes_end = last_bytes_end + \
                    _struct.calcsize(self.ORFs_struct_value[struct_slices[n]])
            bytes_slices[n] = slice(last_bytes_end, this_bytes_end)
            last_bytes_end = this_bytes_end


        # initially zero, starting from 1 for actual clusters which are
        # calculated later after which this database is updated
        unassigned_broad_num = 0

        # better renamed self.ORFs_struct_slices as self.ORFs_struct_key_sl ?
        self.ORFs_struct_slices = struct_slices
        # better renamed self.ORFs_bytes_slices as self.ORFs_struct_value_sl ?
        self.ORFs_bytes_slices = bytes_slices


        # confirmation and demo of packing to bytes for LMDB
        unpacked = 41032143, 41032143, 13677381, 13677381, 7347, 8460, 132156, 1, False
        assert len(_struct.pack(self.ORFs_struct_value, *unpacked)) == \
                _struct.calcsize(self.ORFs_struct_value)
        packed = _struct.pack(self.ORFs_struct_value, *unpacked)
        assert unpacked[:2] == _struct.unpack(
                self.ORFs_struct_value[self.ORFs_struct_slices['DNA']], 
                packed[self.ORFs_bytes_slices['DNA']])
        assert unpacked[2:4] == _struct.unpack(
                self.ORFs_struct_value[self.ORFs_struct_slices['AA']], 
                packed[self.ORFs_bytes_slices['AA']])
        assert unpacked[4:6] == _struct.unpack(
                self.ORFs_struct_value[self.ORFs_struct_slices['chrm']], 
                packed[self.ORFs_bytes_slices['chrm']])
        assert unpacked[6] == _struct.unpack(
                self.ORFs_struct_value[self.ORFs_struct_slices['broad_cluster_num']], 
                packed[self.ORFs_bytes_slices['broad_cluster_num']])[0]
        assert unpacked[7] == _struct.unpack(
                self.ORFs_struct_value[self.ORFs_struct_slices['chrm_num']], 
                packed[self.ORFs_bytes_slices['chrm_num']])[0]
        assert unpacked[8] == _struct.unpack(
                self.ORFs_struct_value[self.ORFs_struct_slices['strand']], 
                packed[self.ORFs_bytes_slices['strand']])[0]


        ### replicon and genome metadate keys and values ###

        # stored as bytes strings in LMDB or in memory in dicts
        # chromosome_num: no reset between genomes, IDs integers get high
        self.chrms_key_struct = 'L'
        # ORF_num to ORF_num: ORF ID range in chromosome from this genome
        self.chrms_val_struct = 'HH'
        # allow > 65k genomes
        self.gnoms_key_struct = 'L'
        # chromosome_num to chromosome_num in this genome, total ORFs
        self.gnoms_val_struct = 'LLH'
        # (longest first for struct)


        ## set up options for sequence storage in bcolz.carrays
        ## these not conserved on re-open? reverts to blosclz . .
        ## tuen the compression options here . . DNA and AA have different information content density
        DNAs_cparams = _bcolz.cparams(clevel=2,cname='zlib')
        AAs_cparams = _bcolz.cparams(clevel=3,cname='lz4')
        RAAs_cparams = _bcolz.cparams(clevel=2,cname='zlib')
        if store_to_disk:
            DNAs_carray_folder = self.folders['DNAs_bcolz']
            AAs_carray_folder = self.folders['AAs_bcolz']
            RAAs_carray_folder = self.folders['RAAs_bcolz']
        else:
            DNAs_carray_folder = None
            AAs_carray_folder = None
            RAAs_carray_folder = None


        rejected_genomes = []
        num_included_ORFs = 0
        included_genome_files = []
        try:
            filename = '{}/all.faa'.format(self.folders['AAs'])
            allAAs = open(filename,'w')
            self.logger.debug('Will write all amino acid ORF sequences to "{}"'\
                    ''.format(filename))
        except IOError:
            self.logger.error('Folders must be created first: make_folders()')

        # for processing and quality control of strings of actual nucleotides and
        # amino acids
        unambiguous_nucs = {b'A', b'C', b'G', b'T'}

        translate = dict([(bytes(k,encoding='utf-8'),v) for k,v in self.translate.items()])

        # for ordinates in main sequence DB
        # only for very first entries at DB creation
        DNA_s = 0
        AA_s = 0
        DNA_e = 0
        AA_e = 0

        cum_chromosome_num = 0
        # genomes start at 1; this will be iterated once before first one
        genome_num = 0

        for n,(genome_file,load) in enumerate(genomes_to_process):
            # encode genome names efficiently
            # could make dict to proper name if included anywhere
            self.logger.log(PROGRESS, 'Parsing {} of {}: {}'.format(n+1, 
                    len(genomes_to_process), genome_file))
            # load using baga, gff or gbk function as determined above
            # <== allow for multifasta i.e., no ORFs+prodigal below (what about ordinates)
            # update other loads
            # stop codons should be omitted in load function
            DNAseqs,ORFs_info = load(genome_file)
            ## need to test if input correct for call_orfs <===
            if call_orfs:
                ORFs = self.call_orfs_prodigal(
                        DNAseqs.view('S{}'.format(len(DNAseqs)))[0])
            
            # check total ORF counts
            total_ORFs = 0
            for seqID,ORF_info in ORFs_info.items():
                total_ORFs += len(ORF_info)
                if len(ORF_info) == 0:
                    self.logger.warning('No ORFs found in {0} from {1}.'.format(
                            seqID, genome_file))
            if total_ORFs == 0:
                self.logger.warning('No ORFs found at all in {}. Excluding.'\
                        ''.format(genome_file))
                rejected_genomes += [genome_file+': no ORFs']
                continue
            elif total_ORFs < min_ORFs_per_genome:
                self.logger.warning('Too few ORFs ({}) found in {}. Need at '\
                        'least {}. Excluding.'.format(total_ORFs, genome_file, 
                        min_ORFs_per_genome))
                rejected_genomes += [genome_file+': too few ORFs ({} < {}, '\
                        'min. required)'.format(total_ORFs, 
                        min_ORFs_per_genome)]
                continue
            
            # first do a QC on whole genome
            # for speed, do a relatively quick per ORF check and reject here,
            # instead of detailed parse+analysis only to find low quality after
            bad_ORFs = {}
            dump_replicons = []
            for (chromosome_ID,chromosome_name),nucs in sorted(DNAseqs.items()):
                ORF_info = ORFs_info[chromosome_ID]
                ## ORF pre-filter for
                ##     non-codon multiples
                ##     ambiguous characters
                ##     stop codons
                # bad annotations can include ORFs off the ends of contigs
                # so filter them out before main list comprehension
                ORF_info_use = [(ORF,(s,e,strd)) for ORF,(s,e,strd) in \
                        ORF_info.items() if e < len(nucs)]
                bad_ORFs[chromosome_ID] = [ORF for ORF,(s,e,strd) in \
                        ORF_info_use if \
                        (e-s)%3!=0 or \
                        set(nucs[s:e]) != unambiguous_nucs or \
                        self.stop_codons[strd] & set(_np.array(nucs[s:e]).view('S3'))]
                # then add the ORFs off the end of the contig to omit
                bad_ORFs[chromosome_ID] += [ORF for ORF,(s,e,strd) in \
                        ORF_info.items() if e >= len(nucs)]
                
                if set(bad_ORFs[chromosome_ID]) == set(ORF_info):
                    # if all ORFs bad or if there are no ORFs called, 
                    # will dump this chromosome from further analysis
                    dump_replicons += [(chromosome_ID,chromosome_name)]
                    if len(bad_ORFs[chromosome_ID]):
                        self.logger.warning('All {} open reading frames failed '\
                        'QC because of ambiguous characters or partial terminal codons '\
                        'in {} from {}: {}'.format(len(ORF_info), 
                        chromosome_ID, genome_file, 
                        ', '.join(sorted(bad_ORFs[chromosome_ID]))))
                else:
                    num_fails = len(bad_ORFs[chromosome_ID])
                    if num_fails:
                        self.logger.warning('{} of {} open reading frames failed '\
                        'QC because of ambiguous characters, partial terminal codons '\
                        'or stop codons in {} from {}: {}'.format(num_fails, 
                        len(ORF_info), chromosome_ID, genome_file, 
                        ', '.join(sorted(bad_ORFs[chromosome_ID]))))
                        # not writing fails, just reporting their locus IDs
                
            prop_fails = sum((len(ORFs) for chromsome_ID,ORFs in \
                    bad_ORFs.items())) / total_ORFs
            if prop_fails > fail_prop_limit:
                self.logger.warning('Excluding {} because {:.00%} of ORFs \
                        failed QC (limit is <{:.00%})'.format(genome_file, 
                        prop_fails, fail_prop_limit))
                rejected_genomes += [genome_file+': too many low quality '\
                        'ORFs ({:.00%} > {:.00%} limit)'.format(prop_fails, 
                        fail_prop_limit)]
                continue
            
            # delete remaining chromosomes without any adequate quality ORFs
            # cannot get this far without any chromosomes surviving this cull
            # i.e., prop_fails <= fail_prop_limit
            for (chromosome_ID,chromosome_name) in dump_replicons:
                del DNAseqs[(chromosome_ID,chromosome_name)]
            
            # there should be some chromsomes (contigs) left
            assert DNAseqs
            
            # genome passed
            genome_num += 1
            if 'mol_' in [seqID[:4] for _,seqID in DNAseqs]:
                # organism_name unknown because molecule names unknown
                organism_name = genome_file
            else:
                # get proper organism names from load() functions: this won't work for contigs
                organism_name = '; '.join([c[1] for c in sorted(DNAseqs, 
                        key = lambda x: len(DNAseqs[x]))])
            assembly_acc_num = 'unknown'
            genome_num2organism[genome_num] = organism_name
            self.logger.info('Processing genome {} of {} input genomes: {}'.format(
                    genome_num, len(genomes_to_process), 
                    genome_num2organism[genome_num]))
            # reset per genome, not between chromosomes
            # ORFs IDs start at zero (not 1)
            # ORF_num = 0 delme
            cum_ORF_num = 0
            genome_DNAs = []
            genome_AAs = []
            genome_RAAs = []
            these_chromosomes_metadata = {}
            for chromosome_num,((chromosome_ID,chromosome_name),nucs) in \
                    enumerate(sorted(DNAseqs.items()), cum_chromosome_num):
                
                ### parse sequence data: slice out ORFs and encode to AA and RAA ###
                ORF_info = ORFs_info[chromosome_ID]
                ##==> can handle double overlapping ORFs but not 3-way (prophage?)
                ## get lists of slices per strand (omit stop codons)
                fwd_slices = [(ORF,(s,e-3)) for ORF,(s,e,strd) in ORF_info.items() \
                        if strd != '-' and ORF not in bad_ORFs[chromosome_ID]]
                rev_slices = [(ORF,(s+3,e)) for ORF,(s,e,strd) in ORF_info.items() \
                        if strd == '-' and ORF not in bad_ORFs[chromosome_ID]]
                fwd_slices.sort(key = lambda x: x[1][0])
                rev_slices.sort(key = lambda x: x[1][0])
                ## separate overlapping
                if len(fwd_slices) == 1:
                    # unless there's only one ORF here (can happen on small contigs)
                    fwd_slices_A = [(ORF,(s,e),True) for ORF,(s,e) in fwd_slices]
                    fwd_slices_B = []
                else:
                    fwd_slices_A = []
                    fwd_slices_B = []
                    overlapped = False
                    for (ORF1,(s1,e1)),(ORF2,(s2,e2)) in \
                            zip(fwd_slices[:-1], fwd_slices[1:]):
                        if overlapped:
                            overlapped = False
                            continue
                        fwd_slices_A += [(ORF1,(s1,e1),True)]
                        if e1 > s2:
                            fwd_slices_B += [(ORF2,(s2,e2),True)]
                            overlapped = True
                    if not overlapped and len(fwd_slices) > 1:
                        fwd_slices_A += [(ORF2,(s2,e2),True)]
                
                if len(rev_slices) == 1:
                    rev_slices_A = [(ORF,(s,e),True) for ORF,(s,e) in rev_slices]
                    rev_slices_B = []
                else:
                    rev_slices_A = []
                    rev_slices_B = []
                    overlapped = False
                    for (ORF1,(s1,e1)),(ORF2,(s2,e2)) in \
                            zip(rev_slices[:-1],rev_slices[1:]):
                        if overlapped:
                            overlapped = False
                            continue
                        rev_slices_A += [(ORF1,(s1,e1),False)]
                        if e1 > s2:
                            rev_slices_B += [(ORF2,(s2,e2),False)]
                            overlapped = True
                    if not overlapped and len(rev_slices) > 1:
                        rev_slices_A += [(ORF2,(s2,e2),False)]
                
                ## get masks of ORFs on each strand
                fwd_strand_mask_A = _np.zeros(nucs.shape,dtype=_np.bool)
                fwd_strand_mask_B = _np.zeros(nucs.shape,dtype=_np.bool)
                [fwd_strand_mask_A.__setitem__(((slice(s,e)),),_np.True_) \
                        for ORF,(s,e),fstrand in fwd_slices_A]
                [fwd_strand_mask_B.__setitem__(((slice(s,e)),),_np.True_) \
                        for ORF,(s,e),fstrand in fwd_slices_B]
                rev_strand_mask_A = _np.zeros(nucs.shape,dtype=_np.bool)
                rev_strand_mask_B = _np.zeros(nucs.shape,dtype=_np.bool)
                [rev_strand_mask_A.__setitem__(((slice(s,e)),),_np.True_) \
                        for ORF,(s,e),fstrand in rev_slices_A]
                [rev_strand_mask_B.__setitem__(((slice(s,e)),),_np.True_) \
                        for ORF,(s,e),fstrand in rev_slices_B]
                # copy strand for reverse
                rev_nucs = _np.array(nucs)
                # complement whole strand (no need to mask coding regions)
                rev_nucs[(nucs==b'A')] = b'T'
                rev_nucs[(nucs==b'T')] = b'A'
                rev_nucs[(nucs==b'C')] = b'G'
                rev_nucs[(nucs==b'G')] = b'C'
                # make a new reversed version of the chromosome
                rev_nucs = _np.array(rev_nucs[::-1])
                # and reverse the ORF masks
                rev_strand_mask_A = rev_strand_mask_A[::-1]
                rev_strand_mask_B = rev_strand_mask_B[::-1]
                # reverse the slices
                rev_slices_A = [(ORF,(nucs.shape[0]-e,nucs.shape[0]-s),fstrand) \
                        for ORF,(s,e),fstrand in rev_slices_A]
                rev_slices_B = [(ORF,(nucs.shape[0]-e,nucs.shape[0]-s),fstrand) \
                        for ORF,(s,e),fstrand in rev_slices_B]
                
                if do_checks:
                    # check against the reversed masks
                    # and sequence
                    try:
                        perform_checks(
                                fwd_slices, fwd_strand_mask_A, fwd_strand_mask_B, 
                                rev_slices, rev_strand_mask_A, rev_strand_mask_B, 
                                fwd_slices_A, fwd_slices_B, 
                                rev_slices_A, rev_slices_B)
                    except AssertionError:
                        self.logger.error('Processing of open reading frames failed '\
                                'checks for this genome')
                
                ## try codon replacement across whole chromosome at once
                # get continuous codon views of coding regions from whole chromosome
                fwd_nucs_ORFs_A = nucs[fwd_strand_mask_A].view('|S3')
                fwd_nucs_ORFs_B = nucs[fwd_strand_mask_B].view('|S3')
                rev_nucs_ORFs_A = rev_nucs[rev_strand_mask_A].view('|S3')
                rev_nucs_ORFs_B = rev_nucs[rev_strand_mask_B].view('|S3')
                
                # translate each codon into continuous array of AAs and RAAs
                # (quicker than boolean arrays as done for complementation above
                # (because '|S3' not optimized?), could be some way of sticking
                # to '|S1' and indexing start of codons to index AAs? stride 
                # _np.where() in 3s + AND on offsets for positions 2 and 3?)
                # S1 (bytes in Python3) much more compact than utf-8
                # and numpy.unicode() much faster than ''.join() at other end.
                # Cannot iterate generator twice as separate list comprehensions
                # would require so put in list and . . .
                # use compound datatype to collect both AA letter and RAA integer
                # this necessitates safe = True for bcolz below - something to do 
                # with striding? is a bit slower so speed gain of single dict looks
                # might be eroded or cancelled because of that.
                dt = _np.dtype([('AA', '|S1'), ('RAA', _np.uint8)])
                AA_RAA_ORFs = _np.fromiter(_chain.from_iterable([
                        map(translate.get,fwd_nucs_ORFs_A),
                        map(translate.get,fwd_nucs_ORFs_B),
                        map(translate.get,rev_nucs_ORFs_A),
                        map(translate.get,rev_nucs_ORFs_B)]), dtype = dt)
                # combine all coding region nucs into continuous array
                DNAs = _np.concatenate((fwd_nucs_ORFs_A.view('|S1'), 
                        fwd_nucs_ORFs_B.view('|S1'), rev_nucs_ORFs_A.view('|S1'), 
                        rev_nucs_ORFs_B.view('|S1')), axis=0)
                
                # sanity check <== pass
                # rev_start = int((fwd_strand_mask_A.sum() + fwd_strand_mask_B.sum()) / 3)
                # n = 0
                # s,e = rev_slices_A[-(n+1)][1]
                # l = int((e-s)/3)
                # from_baga = AA_RAA_ORFs['AA'][rev_start:rev_start+l].tobytes().decode('utf-8')
                # from_file = str(_Seq(nucs[slice(
                        # *info[rev_slices_A[-(n+1)][0]][:2])].tobytes().decode('utf-8'))\
                        # .reverse_complement().translate())[:-1]
                # from_baga,from_file
                
                # save for later
                genome_DNAs += [DNAs]
                genome_AAs += [AA_RAA_ORFs['AA']]
                genome_RAAs += [AA_RAA_ORFs['RAA']]
                
                ### store ORF metadata for this chromosome ###
                
                these_ORFs_metadata = {}
                # first the forward strand ORFs with slices as they will
                # appear in carrays, then reverse which need to be reversed
                # because masks are in strand order . . .
                chrm_slices = _chain.from_iterable(
                        [fwd_slices_A, fwd_slices_B, 
                        rev_slices_A[::-1], rev_slices_B[::-1]])
                
                # keep track of slices in current chromosome array
                # with these_AAs_s:these_AAs_e which are reset for each
                # chromosome, DNA_s:DNA_e grow across all genomes parsed
                # and will be used for recalling sequence on demand 
                # (optionally on-disk carray of all sequence) while
                # chrm_s:chrm_e refer to chromosome sequence being parsed
                these_AAs_s = 0
                these_AAs_e = 0
                for ORF_num,(locus_ID,(chrm_s,chrm_e),f_strand) in \
                        enumerate(chrm_slices, cum_ORF_num):
                    key = _struct.pack(self.ORFs_struct_key, genome_num, ORF_num)
                    ORF_len = chrm_e - chrm_s
                    AA_len = int(ORF_len/3)
                    DNA_e += ORF_len
                    AA_e += AA_len
                    these_AAs_e += AA_len
                    # sanity check ========> rev_slices_A immediately fail <===========
                    # from_baga = AA_RAA_ORFs['AA'][these_AAs_s:these_AAs_e].tobytes().decode('utf-8')
                    # strnd = info[locus_ID][2]
                    # if strnd == '+':
                        # from_file = str(_Seq(nucs[slice(*info[locus_ID][:2])].tobytes().decode('utf-8')).translate())[:-1]
                    # else:
                        # from_file = str(_Seq(
                                # nucs[slice(*info[locus_ID][:2])].tobytes().decode('utf-8')
                                # ).reverse_complement().translate())[:-1]
                    # assert from_baga == from_file
                    
                    # fixed width packed integer values at beginning
                    # variable width locus ID at end
                    value = _struct.pack(self.ORFs_struct_value, DNA_s, DNA_e, 
                            AA_s, AA_e, chrm_s, chrm_e, unassigned_broad_num, 
                            chromosome_num, f_strand) + locus_ID.encode('utf-8')
                    
                    these_ORFs_metadata[key] = value
                    
                    if assign_narrow_clusters:
                        # narrow clusters are matching reduced amino acid sequences
                        # allows common SNPs - synonymous and nonsynonymous within RAA
                        narrows[AA_RAA_ORFs['RAA'][these_AAs_s:these_AAs_e].tobytes()] \
                                += [(genome_num, ORF_num)]
                    
                    ### write to full AA database
                    allAAs.write('>{:0{}d}__{:05d}\n{}\n'.format(genome_num, 
                            padding, ORF_num, _np.unicode(
                            AA_RAA_ORFs['AA'][these_AAs_s:these_AAs_e].tobytes(), 
                            'utf-8')))
                    
                    DNA_s += ORF_len
                    AA_s += AA_len
                    these_AAs_s += AA_len
                    num_included_ORFs += 1
                
                if store_to_disk:
                    # in LMDB if store_to_disk = True else dict in memory
                    with lmdb_env.begin(db = ORFs_metadata, write=True) as transaction:
                        # transaction.putmulti(these_ORFs_metadata.items()) not available
                        for key, value in these_ORFs_metadata.items():
                            transaction.put(key, value)
                        # need this, because metasync=False, sync=False
                        # i.e., manually sending data to disk because this use case is
                        # simple and predictable.
                        lmdb_env.sync()
                    # docs imply force=True required when created with sync=False
                    # but not available . . .
                else:
                    # keep in one big dict in memory
                    ORFs_metadata.update(these_ORFs_metadata)
                
                # put metametadata per chromosome in a dict
                RefSeqID = chromosome_ID
                GenBankID = 'unknown'  # <===== update "load" functions to collect organism ID etc <== handle contigs properly also
                # key == chromosome_num which is unique across whole dataset
                # value starts with ORF_num_start, ORF_num_end i.e., range(ORF_num_start, ORF_num_end)
                # for generating IDs and calculating total by subtraction
                # L because cumulative unique chromsome IDs can go over >65k when
                # using draft genomes i.e., lots of contigs
                key = _struct.pack(self.chrms_key_struct, chromosome_num)
                value = _struct.pack(self.chrms_val_struct, cum_ORF_num, ORF_num) + \
                        '|'.join([GenBankID, RefSeqID, chromosome_name]).encode("utf-8")
                these_chromosomes_metadata[key] = value
                if these_ORFs_metadata:
                    # only iterate if there were any ORFs in this sequence
                    # can happen with short contigs . . .
                    cum_ORF_num = ORF_num + 1
            
            ######## store number of ORFs
            
            ### per genome, write genome and chromosome level metadata and sequences to disk ###
            # metametadata per genome: organism_name, assembly_acc_num, input filename
            genome_key = _struct.pack(self.gnoms_key_struct, genome_num)
            # key = genome_num which is unique across whole dataset
            # value starts with total ORFs, chromsome_num_start, chromsome_num_end 4 bytes each
            # i.e., range(chromsome_num_start, chromsome_num_end)
            # for generating IDs for lookups and calculating total by subtraction
            genome_value = _struct.pack(
                    self.gnoms_val_struct, cum_chromosome_num, chromosome_num,
                    ORF_num+1) + '|'.join([organism_name, assembly_acc_num, 
                    genome_file]).encode("utf-8")
            # iterate cumulative chromsome number so next genome starts
            # counting chromosomes at the next integer for unique IDs
            cum_chromosome_num = chromosome_num + 1
            
            if store_to_disk:
                with lmdb_env.begin(db = genomes_metadata, write = True) as transaction:
                    transaction.put(genome_key,genome_value)
                    lmdb_env.sync()
                # chromosome_num to GenBankID|RefSeqID|contigname
                # empty as required (e.g. if GenBank only or contig etc)
                with lmdb_env.begin(db = chromosomes_metadata, write = True) as transaction:
                    for key,value in these_chromosomes_metadata.items():
                        transaction.put(key, value)
                    lmdb_env.sync()
            else:
                genomes_metadata[genome_key] = genome_value
                chromosomes_metadata.update(these_chromosomes_metadata)
            
            if genome_num == 1:
                # first round, add all pieces, for expected_len use +50% * num_genomes
                # safe = False because we always add 1d numpy arrays so no additional 
                # coersion needed
                DNAs_exp_len = int(sum([DNAs.shape[0] \
                        for DNAs in genome_DNAs]) * len(genomes_to_process) * 1.5)
                DNAs_carray = _bcolz.carray(genome_DNAs[0], 
                        rootdir = DNAs_carray_folder, cparams = DNAs_cparams, 
                        expectedlen = DNAs_exp_len, mode = 'w', safe = False)
                for DNAs in genome_DNAs[1:]:
                    DNAs_carray.append(DNAs)
                
                # actually, safe = True for RAAs and AAs because something about array 
                # construction makes carray corruption (misaligned?) without. 
                # Would be good to work out why sometime because safe = True
                # is a bit slower. Probably compound data type used to do single
                # dict lookup when translating
                AAs_exp_len = int(sum([AAs.shape[0] \
                        for AAs in genome_AAs]) * len(genomes_to_process) * 1.5)
                AAs_carray = _bcolz.carray(genome_AAs[0], 
                        rootdir = AAs_carray_folder, cparams = AAs_cparams, 
                        expectedlen = AAs_exp_len, mode = 'w', safe = True)
                for AAs in genome_AAs[1:]:
                    AAs_carray.append(AAs)
                
                RAAs_carray = _bcolz.carray(genome_RAAs[0], 
                        rootdir = RAAs_carray_folder, cparams = RAAs_cparams, 
                        expectedlen = AAs_exp_len, mode = 'w', safe = True)
                for RAAs in genome_RAAs[1:]:
                    RAAs_carray.append(RAAs)
            
            else:
                # might be more efficient to keep open between genome loops?
                for DNAs in genome_DNAs:
                    DNAs_carray.append(DNAs)
                for AAs in genome_AAs:
                    AAs_carray.append(AAs)
                for RAAs in genome_RAAs:
                    RAAs_carray.append(RAAs)
            
            DNAs_carray.flush()
            AAs_carray.flush()
            RAAs_carray.flush()
            
            self.logger.info('{} ({}) passed QC with {:.00%} of {} ORFs '\
                    'omitted (limit is <{:.00%}). Retained {}'.format(
                    genome_file, genome_num, prop_fails, total_ORFs, 
                    fail_prop_limit, ORF_num+1))
            
            # useful for set ops on ORFs because IDs are integer range
            # start at 0 not 1 range(ORF_num) returns all IDs for sets
            ORFs_per_genome[genome_num] = ORF_num + 1
            
            ## do minhash mash sketch of genome sequence
            sketch_out = '{}/{:0{}d}.msh'.format(self.folders['sketches'], 
                    genome_num, padding)
            if not _os.path.exists(sketch_out) or force:
                # save genome sequence to temporary fasta for mash sketch
                fna_out = '{}/{:0{}d}.fna'.format(self.folders['DNAs'], 
                        genome_num, padding)
                lines = ['>{:0{}d}\n'.format(genome_num, padding)]
                # join chromosomes to one contiguous sequence for minhash calcs
                lines += [''.join((g.tobytes().decode('utf-8') 
                        for g in DNAseqs.values()))]
                cmd = [exe_mash,'sketch','-o',sketch_out, fna_out]
                self.logger.log(PROGRESS, 'Doing Mash minhash sketch of DNA '\
                        'sequences for {} ({}): "{}"'.format(genome_file, 
                        genome_num, ' '.join(cmd)))
                open(fna_out,'w').write('\n'.join(lines))
                # try: here?
                proc_returncode,stdout = self._launch_external(cmd, 'mash', 
                        self.log_folder, main_logger = self.logger)
                _os.unlink(fna_out)
            else:
                self.logger.log(PROGRESS, 'Found mash sketch for {} ({}) at '\
                        '{}, use --force to generate a new one'.format(
                        genome_file, genome_num, sketch_out))
            
            if (genome_num+1) % unclustered_clearence_freq == 0:
                # every e.g. 10 genomes, clear single copy (yet to cluster) sequences
                # to save memory. Most sequence redundancy will be in the core genome
                # and will be encountered after only a few genomes.
                before = len(narrows)
                narrows = _defaultdict(list,{seq:ORFs for seq,ORFs in narrows.items() \
                        if len(ORFs) > 1})
                after = len(narrows)
                self.logger.debug('Removed {} unclustered ORFs from narrow clustering '\
                        'leaving {} clusters'.format(before-after, after))
                
                # could also write broad and narrow cluster data to disk at these intervals
                # for continuing an interrupted analysis

        allAAs.close()


        ## should put this in LMDB too?
        ## same time as above? as per ORF ID?
        ## would lower memory footprint . . .
        if assign_narrow_clusters:
            #### retaining all narrow clustered ORFs
            self.narrow_clusters = [tuple(narrow) for narrow in narrows.values() if \
                    len(narrow) > 1]
            #### retaining only single copy narrow clusters . . 
            # single_copy_narrows = []
            # ORFs_clustered = 0
            # # ensure only single copy ORFs are clustered so that double copy ORFs can be
            # # properly assessed for in- or out-paralogy by distance measures
            # notkept = 0
            # for seq,ORFs in narrows.items():
                # g2o = {g:(g,o) for g,o in ORFs}
                # g_freq = _Counter(g for g,o in ORFs)
                # keep = tuple(g2o[g] for (g,f) in g_freq.items() if f == 1)
                # notkeep = tuple((g,f) for (g,f) in g_freq.items() if f > 1)
                # if len(keep) > 1:
                    # # if after removing multi-copy ORFs more than one ORFs remain in
                    # # the cluster, save it (but not the integer encoded sequence)
                    # single_copy_narrows += [keep]
                    # ORFs_clustered += len(keep)
                # if notkeep:
                    # self.logger.debug('Omitting from narrows . . .')
                    # for g,f in notkeep:
                        # self.logger.debug('{} from genomes {}'.format(f,g))
                        # notkept += f
            
            # self.logger.info('{:,} similar single copy ORFs assigned to {:,} clusters'\
                    # ''.format(ORFs_clustered, len(single_copy_narrows)))
            # self.logger.info('{:,} ORFs removed from narrows because multicopy'\
                    # ''.format(notkept))
            # self.narrow_clusters = single_copy_narrows

        # genomes and ORFs encoded in contiguous integer range from zero
        # len(self.ORFs_per_genome) is total genomes
        # ORFs_per_genome is count from 1, genome and ORF IDs are from 1
        assert set(range(1, len(ORFs_per_genome)+1)) == set(ORFs_per_genome), ''\
                'Genome IDs should be a continuous range of integers from 1'
        self.ORFs_per_genome = ORFs_per_genome
        # for writing per-genome files
        self.padding = padding

        DNAs_carray.flush()
        AAs_carray.flush()
        RAAs_carray.flush()

        if store_to_disk:
            # on disk, can be re-opened on demand
            ## causes
            ## SyntaxError: can not delete variable 'RAAs_carray' referenced in nested scope
            ## sometimes - garbage collection etc can probably handle this anyway
            # del DNAs_carray
            # del AAs_carray
            # del RAAs_carray
            # close after last genome (could do between genomes but .flush() and .sync()
            # enough to be safe?)
            lmdb_env.close()
        else:
            # store state in Finder instance instead
            self.sequences = {}
            self.sequences['DNA'] = DNAs_carray
            self.sequences['AA'] = AAs_carray
            self.sequences['RAA'] = RAAs_carray
            self.genomes_metadata = genomes_metadata
            self.chromosomes_metadata = chromosomes_metadata
            self.ORFs_metadata = ORFs_metadata

        self.logger.info('{} genomes processed. {} failed QC leaving {} for '\
                'analysis of {:,} ORFs.'.format(len(genomes_to_process), 
                len(rejected_genomes), len(genomes_to_process) - len(rejected_genomes), 
                num_included_ORFs))

        if rejected_genomes:
            self.logger.info('These input genomes failed QC: ' + \
                    '; '.join(rejected_genomes))


    def make_genome_tree(self, retain_individual_sketches = False):
        '''Make a minhash-distance neighbor-joining genome tree

        Infer a tree-like phylogeny of genomes using neighbor-joining from Mash 
        minhash distances of concatentated DNA sequences.

        Parameters
        ----------
        retain_individual_sketches : bool
            Do not delete each genome Mash sketch generated by processORFs() 
            after pasting them to single file
        '''

        # max threads is slightly different to cpus
        # . . can probably use more
        # max_processes = _decide_max_processes( max_cpus )

        from baga.Dependencies import dependencies
        d = dependencies['mash']['default_source']
        exe_mash = _os.path.join(dependencies['mash'][d]['exe']['main'])

        genome_padding = len(str(len(self.ORFs_per_genome)))

        # write list of input mashes
        list_of_sketches = '{}/for_pasting.txt'.format(self.folders['sketches'])
        with open(list_of_sketches,'w') as fout:
            # genome IDs start at 1
            for genome in range(1, len(self.ORFs_per_genome)+1):
                fout.write('{}/{:0{}d}.msh\n'.format(self.folders['sketches'], 
                        genome, self.padding))

        allsketches = '{}/allseqs'.format(self.folders['sketches'])
        try:
            _os.unlink(allsketches+'.msh')
        except FileNotFoundError:
            pass
        cmd = [exe_mash, 'paste', allsketches, '-l', list_of_sketches]
        self.logger.log(PROGRESS,'Combining Mash sketches: "{}"'.format(
                ' '.join(cmd)))
        proc_returncode,stdout = self._launch_external(cmd, 'mash', 
                self.log_folder, main_logger = self.logger)

        _os.unlink('{}/for_pasting.txt'.format(self.folders['sketches']))

        if not retain_individual_sketches:
            self.logger.debug('Deleting the {} individual Mash sketches after '\
                    'pasting them to single file'.format(self.num_genomes))
            # genome IDs start at 1
            for genome in range(1, len(self.ORFs_per_genome)+1):
                _os.unlink('{}/{:0{}d}.msh'.format(self.folders['sketches'], 
                        genome, self.padding))

        allsketches += '.msh'
        cmd = [exe_mash, 'dist', allsketches, allsketches]
        self.logger.log(PROGRESS,'Calculating MinHash Mash distances: "{}"'\
                ''.format(' '.join(cmd)))
        proc_returncode,dists = self._launch_external(cmd, 
                'mash', self.log_folder, main_logger = self.logger,
                return_stdout = True)

        self.logger.log(PROGRESS,'Calculating NJ tree from MinHash Mash distances')
        # convert to a phylip matrix
        phylip_dists_dict = {}
        for line in dists.decode('utf-8').split('\n'):
            if len(line):
                cells = line.split('\t')
                A = cells[0].split('/')[-1].split('.')[0]
                B = cells[1].split('/')[-1].split('.')[0]
                try:
                    phylip_dists_dict[A][B] = float(cells[2])
                except KeyError:
                    phylip_dists_dict[A] = {B:float(cells[2])}
                try:
                    phylip_dists_dict[B][A] = float(cells[2])
                except KeyError:
                    phylip_dists_dict[B] = {A:float(cells[2])}

        phylip_dists = ['   {}'.format(len(phylip_dists_dict))]
        for A,Bs in sorted(phylip_dists_dict.items()):
            this_row = [A]
            for B,dist in sorted(Bs.items()):
                this_row += [str(dist)]
            phylip_dists += ['\t'.join(this_row)]

        #### for very big trees: rapidNJ (else DendroPy, below)
        # dists_file = '{}/genome_dists.phy'.format(self.folders['sketches'])
        # open(dists_file,'w').write('\n'.join(phylip_dists))
        # d = dependencies['rapidNJ']['default_source']
        # exe_rapidNJ = _os.path.join(dependencies['rapidNJ'][d]['exe']['main'])
        # try:
            # o = _subprocess.check_output([exe_rapidNJ, dists_file,
                    # '--input-format', 'pd'])
        # except _subprocess.CalledProcessError as e:
            # print(e)
        # print(o)
        # t = _dendropy.Tree.get_from_string(o.decode('utf-8'),'newick')

        # could output phylip distances for convenience . . .

        for_dendropy = ['\t'.join(['']+sorted(phylip_dists_dict))]
        for_dendropy += phylip_dists[1:]

        pdm = _StringIO('\n'.join(for_dendropy))
        pdm.seek(0)
        pdm = _dendropy.PhylogeneticDistanceMatrix.from_csv(
                src=pdm, delimiter="\t")

        tree = pdm.nj_tree()
        tree.reroot_at_midpoint()
        tree_outfile = '{}/mash_tree.newick'.format(self.folders['sketches'])
        self.logger.log(PROGRESS, 'Wrote Mash minhash genome tree to {}'.format(tree_outfile))
        tree.write_to_path(tree_outfile, 'newick')

        self.genome_tree = tree
        self.mash_distances = phylip_dists_dict

    def phylo_select_genomes(self, max_genomes_per_group = 10, 
            genome_pairs_shared = 1):
        '''Make phylogenetically optimised genome selections

        From the ladderized tip order along the genome tree, alternately
        select a genome for each group. This is part of the "map" algorithm 
        in the context of the MapReduce paradigm.

        Parameters
        ----------
        max_genomes_per_group : int
            Preferred genome group size for pairwise comparisons including those
            shared between groups.
        genome_pairs_shared   : int
            Number of genome pairs shared between groups. Allows merging of genomes
            at the "Reduce" stage. More shared genomes means more pairwise 
            comparisons which scales poorly.
        '''


        # "Phylogenetic MapReduce for scalable gene family delineation across millions of bacterial genomes"
        # (2 million on 16 cores . . .)

        ## 1) divide and conquer is required for scalability

        ## 2) phylogeny-aware: divide even PD selection with "neighbour"-wise overlapping pairs for spanning.
        ## Strikes a balance between spanning rare, local families by best match (highly xenologous ORFs less so)
        ## and maximising gain at union-find stage
        ## benefits from all sub-groups spanning deepest split in mash tree i.e., same context
        # objective: delineate groups for separate best match analyses that
            # i) encompass common, core families: high phylogenetic diversity
            # ii) encompass rare, accessory families: low phylogenetic diversity
            # iii) minimises number of genomes in each analysis <= to solve scaling problem
            # iv) overlap <= to allow joining of families by transience of orthology or inparalogy without unnecessary analysis

        # several optima to aim at: including simply group size versus number of groups
        # even sized groups are better than skewed sized groups
            # skew includes larger groups that are disporportionately slower to analyse: scaling
        # small groups size is better than large
            # larger groups are disporportionately slower to analyse: scaling
        # large groups include more direct comparisons, less merging by overlap
            # merging by overlap only works for families present in common genomes and means more direct analyses afterwards

        ## main dilemma is how big should the overlap be?
        ## Bigger overlap:
            # larger, slower PW comparisons
            # more successful merges at reduce 1 unionfind
        ## small, rare, patchy families always most difficult to capture by unionfind
        ## more like not possible: will always need second round of PW by family reps so don't worry too much: stick with 2

        assert max_genomes_per_group > genome_pairs_shared*2, "requiring {} "\
                "shared genomes but max. group size {} will not leave any "\
                "unique genomes per group. max_genomes_per_group > "\
                "genomes_shared is required.".format(max_genomes_per_group, 
                genomes_shared)

        use_genomes_per_group = max_genomes_per_group - genome_pairs_shared*2
        num_groups = round((len(self.ORFs_per_genome) / use_genomes_per_group) + 0.5)

        # reset any previously selected genomes
        self.bm_selected_genomes = []
        # reset any previously tested genomes
        self.bm_tested_genomes = []

        self.logger.debug('Distributing {} genomes into {} groups for '\
                'separate pair-wise analysis'.format(len(self.ORFs_per_genome), 
                num_groups))

        # mash genome tree
        if not self.genome_tree.is_rooted:
            # actually important to be midpoint rooted
            # arbitrary rooting could be bad for phylogenetic assumptions
            # accurate rooting probably better than midpoint approximation
            # should have been midpoint rerooted in .make_genome_tree()
            self.genome_tree.reroot_at_midpoint()

        ## 1) ladderize by diversity
        # Similar to conventional ladderize for plotting which sorts by
        # number of child nodes but this sorts by cumulative edge length
        # i.e., diversity.
        # Change in similarity between adjacent genomes increases approximately
        # gradually allowing organisation into groups of similar phylogenetic
        # diversity with little computational cost.

        self._ladderize_diversity(self.genome_tree, ascending=False)
        #print(self.genome_tree.as_ascii_plot(plot_metric = 'length'))

        # only select from genomes with some narrow cluster representatives
        # some very similar genomes won't be represented at all
        ladderized_order = [n.taxon for n in self.genome_tree.leaf_nodes()]
        assert len(ladderized_order) == len(self.ORFs_per_genome)
        #print(ladderized_order)

        ## 2) put genomes into groups by selecting alternative genomes along 
        ## phylogeny representation
        # argument is genomes_per_groupm not total groups, to ensure each task 
        # is a managable size regardless of total number of genomes

        groups = []
        for g in range(num_groups):
            groups += [ladderized_order[g::num_groups]]

        #groupsizes = _Counter((len(group) for group in groups))
        assert sum(map(len,groups)) == len(self.ORFs_per_genome)

        # ensure all groups have at least two of deepest clades for
        # proper out-paralog context
        # ==> "balance" shared genomes here instead of always adding same one?
        root_node = self.genome_tree.mrca(taxa=self.genome_tree.taxon_namespace)
        deepest_bitmasks = [node.bipartition.split_bitmask for node in root_node.child_nodes()]
        txn = self.genome_tree.taxon_namespace
        for group in groups:
            this_bitmask = txn.get_taxa_bitmask(taxa=group)
            span_children = [True if this_bitmask & deep_bitmask > 0 else False \
                    for deep_bitmask in deepest_bitmasks]
            if sum(span_children) < 2:
                for n,child in enumerate(deepest_bitmasks):
                    if not span_children[n]:
                        spanning_taxon = txn.bitmask_taxa_list(deepest_bitmasks[n])[0]
                        self.logger.debug('To enforce spanning of deepest split, '\
                                'adding {} to group: {}'.format(spanning_taxon.label, 
                                ','.join(t.label for t in group)))
                        group += [spanning_taxon]
                        break


        ## 3) copy a number of genome pairs from each group to one other to establish
        ## ovelaps i.e., each pair is a shared genome from both ladderised "neighbours"
        # genome_pairs_shared
        self.logger.debug('Distributing {} pair(s) of shared genomes between '\
                'groups'.format(genome_pairs_shared))

        groupsizes = ['There were ']
        for size,freq in _Counter((len(group) for group in groups)).most_common():
            groupsizes += ['{} groups of size {}'.format(freq,size),'; ']
        groupsizes = ''.join(groupsizes[:-1])
        self.logger.debug(groupsizes)

        orig_groups = list([list(group) for group in groups])
        for s in range(genome_pairs_shared):
            for n,(group,nextgroup) in enumerate(zip(orig_groups[:-1],orig_groups[1:])):
                # this group gets a genome from the end of other group
                groups[n] += [nextgroup[-(1+s)]]
                # gives from beginning of self to other group
                groups[n+1] += [group[s]]
            # last group gets end of first group
            groups[-1] += [orig_groups[0][-(1+s)]]
            # last group gives beginning to first group
            groups[0] += [orig_groups[-1][s]]

        # occasionally, (2) has already added a genome that would be added here
        # so ensure unique by using set()
        groupsizes = ['After genome sharing there are ']
        for size,freq in _Counter((len(set(group)) for group in groups)).most_common():
            groupsizes += ['{} groups of size {}'.format(freq,size),'; ']
        groupsizes = ''.join(groupsizes[:-1])
        self.logger.debug(groupsizes)

        self.bm_selected_genomes = []
        for group in groups:
            self.bm_selected_genomes += [_array('H', set(int(taxon.label) for taxon in group))]

        all_selected = set(a for b in self.bm_selected_genomes for a in b)
        assert set(self.ORFs_per_genome) == all_selected, \
                'Not all genomes selected for analysis :-( . . this is a bug'

    def process_narrow_clusters(self):
        '''Select among narrow clusters a representative ORF

        Bias representatives towards genomes featured in more best match analysis groups
        i.e., those linking transitive homology between analyses.

        ==> not sure this actually matters?
        '''


        if not self.bm_selected_genomes:
            self.logger.warning('phylo_select_genomes() must be run to populate '\
                    '"bm_selected_genomes" list')
            return

        self.logger.info('Selecting narrow ORF cluster representatives')

        # for very large scaling, these can be put into LMDBs
        # total ORFs in narrow clusters will depend on number and
        # relatedness of input genomes: more closely related, more and
        # larger narrow clusters. However, very large 2-tuple:list(2-tuples)
        # should fit in memory.
        cluster_expansions = {}
        cluster_noreps = {}
        cluster_rep_bygenome = {}

        # Select a representative per narrow cluster prioritising genomes 
        # shared between groups. This will skew genome
        # choice so more gene families have an ORF in two or more genome
        # groups and can be merged in the Reduce 1 stage.

        # always prefer higher frequency, more shared genome
        genome_rep_counts = _Counter((genome for cluster in \
                self.bm_selected_genomes for genome in cluster))
        freqs = sorted(set(genome_rep_counts.values()), reverse = True)
        shared_genomes = []
        for thisfreq in freqs:
            shared_genomes += [[genome for genome,freq in \
                    genome_rep_counts.items() if freq == thisfreq]]

        # all genomes should be present in this dictionary
        # even if set is empty i.e., no ORFs clustered
        # it is used to determine which ORFs should be included
        # in comparisons when extracting kmers per ORF later
        clustered_bygenome = {genome:set() for genome in self.ORFs_per_genome}

        c = 0
        for cluster in self.narrow_clusters:
            assert len(cluster) > 1
            cluster_represented = False
            # select a single ORF per genome represented in this cluster
            # (should be single copy only . . .)
            genomes2reps = dict(((genome_num,(genome_num,ORF_num)) \
                    for genome_num,ORF_num in cluster))
            # currently not considered: RAA jaccard is limit of resolution
            # during initial clustering - paralogous draws are resolved later
            # assert len(genomes2reps) == len(cluster), 'multi-copy narrow cluster!!'
            # decide which genome should have the cluster representative
            # first consider one of those at the highest frequency
            for prioritised in shared_genomes:
                options = sorted(set(prioritised) & set(genomes2reps))
                if options:
                    genome = _random.choice(options)
                    representative_ORF = genomes2reps[genome]
                    cluster_expansions[representative_ORF] = cluster
                    c += 1
                    cluster_represented = True
                    # collect representatives
                    try:
                        cluster_rep_bygenome[genome] += [representative_ORF]
                    except KeyError:
                        cluster_rep_bygenome[genome] = [representative_ORF]
                    # record which are represented, 
                    # excluding representatives and singles
                    for genome2,ORF in cluster:
                        if representative_ORF != (genome2,ORF):
                            clustered_bygenome[genome2].add((genome2,ORF))
                    break
            # a representative should always have been selected
            # so assert this is the case
            assert representative_ORF in cluster
            assert cluster_represented, 'a narrow cluster is not '\
                    'represented by any genome . . . this is a bug :-('


        self.logger.info('{:,} narrow clusters are represented'\
                ''.format(len(cluster_expansions)))

        total_represented = len([a for b in cluster_expansions.values() for a in b]) - c
        self.logger.info('{:,} ORFs will be indirectly represented by narrow clusters '\
                'after database reduction'.format(total_represented))

        total_remaining = sum((num_ORFs - len(clustered_bygenome[g]) \
                for g,num_ORFs in self.ORFs_per_genome.items()))
        self.logger.info('{:,} ORFs will be directly analysed after database reduction'\
                ''.format(total_remaining))

        # some genome-to-genome best matches already done within narrow clusters
        reps_bm_done = {}
        # some high-similarity groups already found where all genomes present within a
        # a narrow cluster. (small risk of multiple, mutually exclusive out-paralogous groups?)
        # most probably orthologs, extra copies probably inparalogs for resolution by DNA and/or
        # chromosome position, but possibly out-paralogs which could be resolved by DNA jaccards . . .
        bm_complete = []
        complete_reps = []
        for rep,mems in cluster_expansions.items():
            if set([genome for genome,ORF_num in mems]) == set(self.ORFs_per_genome):
                # every genome represented in this CD-HIT high identity cluster
                # close to zero diversity (less than narrow cluster limit 
                # setting) across all genomes in this gene family
                bm_complete += [mems]
                # must also exclude rep from further clustering (all best matches
                # and possibly some additional - already found)
                # so record as clustered (represented already added above)
                # "clustered_bygenome" ultimatley determines which ORFs are
                # included n later best match clustering
                clustered_bygenome[rep[0]].add(rep)
                complete_reps += [rep]
            else:
                # make a note of these unambiguous best matches
                # i.e. which genomes has this representative already been matched to
                # within a narrow cluster to be expanded at the end
                hits_to_others = set([genome for genome,ORF_num in mems \
                        if (genome,ORF_num) != rep])
                if len(hits_to_others):
                    reps_bm_done[rep] = hits_to_others

        # will be added at end so can be removed here
        for genome,ORF_num in complete_reps:
            del cluster_expansions[genome,ORF_num]
            try:
                cluster_rep_bygenome[genome].remove((genome,ORF_num))
            except KeyError:
                # genome might be fully clustered i.e., absent in representatives because
                # always represented by other genomes because very similar
                pass

        if len(bm_complete):
            self.logger.info('{:,} narrow clusters have members spanning all '\
                    'genomes (one per genome) and are assumed to be orthologous. '\
                    'There is a small risk these are a pair of separate '\
                    'orthologous families, in mutually exclusive genomes, that '\
                    'are out-paralogous to each other.'.format(
                    len(bm_complete)))


        # dictionary of ORFs representing CD-HIT clusters

        self.reduced_db = {}
        self.reduced_db['expansions'] = cluster_expansions
        # same as expansions, organised by genomes
        # ORFs not added to CD-HIT clusters + 
        # ORFs removed from CD-HIT clusters because of in/out-paralog ambiguity +
        # these == ORFs for best match analysis
        self.reduced_db['reps_by_genome'] = cluster_rep_bygenome
        self.reduced_db['clustered_by_genome'] = clustered_bygenome
        # these are used during best match tests to check whether a
        # comparison needs to be done
        self.reduced_db['bm_genomes_within_narrows'] = reps_bm_done
        # these are added at the end at expansion time
        self.reduced_db['bm_clusters_within_narrows'] = bm_complete

        # work out which genomes have no ORFs for analysis after narrow cluster
        # reduction for very similar genomes . . .

        # genomes with ORFs representing other ORFs
        genomes_in_representing = set([genome_num for genome_num,ORF_num in \
                self.reduced_db['expansions']])

        # genomes with ORFs represented by ORFs from self or others
        genomes_in_represented = set()
        for cluster in self.reduced_db['expansions'].values():
            for genome_num,ORF_num in cluster:
                genomes_in_represented.add(genome_num)

        # collect genomes with ORFs neither in nor representing a narrow cluster
        # i.e., ORFs not involved in narrow clustering at all
        # then compare with representing ORFs in each genome
        genomes_w_nonrepresenting = set()
        # iterate through genome IDs
        # genome IDs start at 1
        for get_genome_num in self.ORFs_per_genome:
            if get_genome_num not in self.reduced_db['reps_by_genome']:
                # not ORFs are representatives at all
                genomes_w_nonrepresenting.add(get_genome_num)
            else:
                # generate ORF IDs and see if any (how many?) are not representatives
                # use knowledge of total ORFs in this genome with range()
                these_ORFs = set(tuple([genome_num,ORF_num]) \
                        for ORF_num in range(self.ORFs_per_genome[get_genome_num]))
                if len(these_ORFs - set(self.reduced_db['reps_by_genome'][get_genome_num])):
                    genomes_w_nonrepresenting.add(get_genome_num)
                else:
                    e = "failed to find genome {} in DB".format(get_genome_num)
                    print(e)
                    #self.logger.error(e)
                    # raise LookupError

        # genomes with all ORFs represented in clusters
        # can happen with two very similar genomes
        self.genomes_not_representing = sorted(
                genomes_in_represented - \
                genomes_in_representing - \
                genomes_w_nonrepresenting)

        if len(self.genomes_not_representing):
            self.logger.warn('{} genome(s), ID(s): {}, are entirely represented '\
                    'by others after narrow cluster database reduction. Some genomes '\
                    'are identical or nearly so.'.format(
                    len(self.genomes_not_representing), 
                    ', '.join(self.genomes_not_representing)))
            
            # update genomes selected for best match analysis
            new_bm_selected_genomes = []
            left_over_single_genomes = []
            for i,genome_group in enumerate(self.bm_selected_genomes):
                genome_group = [genome for genome in genome_group if \
                        genome not in self.genomes_not_representing]
                if len(genome_group) > 1:
                    new_bm_selected_genomes += [_array('H', genome_group)]
                elif len(genome_group) == 1:
                    left_over_single_genomes += [genome_group[0]]
            
            # redistribute any single genomes left over
            while left_over_single_genomes:
                for group in new_bm_selected_genomes:
                    group += [left_over_single_genomes.pop()]
                    if not left_over_single_genomes:
                        break
            
            if len(new_bm_selected_genomes) < self.bm_selected_genomes:
                self.logger.debug('There are {} fewer genome groups for best matching '\
                        'because some very high genome similarity')
            
            groupsizes = ['After processing narrow clusters there are ']
            for size,freq in _Counter((len(set(group)) for group in groups)).most_common():
                groupsizes += ['{} groups of size {}'.format(freq,size),'; ']
            groupsizes = ''.join(groupsizes[:-1])
            self.logger.debug(groupsizes)
            
            self.bm_selected_genomes = new_bm_selected_genomes


    def AA_DB_broad_clustering(self, iterations = [(0.80,5,0.6), (0.68,4,0.6), (0.56,4,0.6)], 
            store_to_disk = True, max_cpus = -1, max_memory = 8):
        '''Create broad clusters of distant sequences using CD-HIT

        Parameters
        ----------
        iterations : list
            tuples of (float,int,float) representing (identity, word length, 
            length difference) for subsequent iterations. Identity should decrease,
            Word length should decrease with identity as recommended in the CD-Hit 
            documentation. Length difference can remain constant.
        store_to_disk : True
            if True, the LMDB key-value store of ORF metadata will be updated. Use False
            to store ORF broad cluster affiliation in a dictionary in memory. If the whole
            analysis can fit in memory, it is faster not to store_to_disk. To allow larger
            analyses than memory permits, set store_to_disk as True.
        max_cpus : int
            Maximum CPUs or threads to use
        max_memory : int
            Maximum memory in GB to use
        '''

        # create input DB from .reduced_db['reps_by_genome'] ORFs
        # this is remaining AA sequences after narrow clusters collapsed
        # i.e., narrow cluster representatives and those unaffected by 
        # narrow clustering (representing themselves only)
        sourcefile = '{}/all.faa'.format(self.folders['AAs'])
        outfile = '{}/reduced.faa'.format(self.folders['AAs'])

        ORF_IDs = []
        for genome,num_ORFs in self.ORFs_per_genome.items():
            # generate all ORF IDs for this genome
            remaining_ORFs = set(tuple([genome,ORF_num]) \
                    for ORF_num in range(num_ORFs))
            # subtract ORFs within narrow clusters
            try:
                remaining_ORFs.difference_update(self.reduced_db['clustered_by_genome'][genome])
            except KeyError:
                self.logger.debug('genome {} has no ORFs in narrow clusters . . .'.format(genome))
            # restore narrow cluster representatives. This leaves ORFs left
            # after narrow clustering for subsequent analysis.
            try:
                remaining_ORFs.update(self.reduced_db['reps_by_genome'][genome])
            except KeyError:
                self.logger.debug('genome {} has no ORFs representing narrow clusters . . .'.format(genome))
            ORF_IDs += [('{:0{}d}__{:05}'.format(genome,self.padding,ORF) \
                    for genome,ORF in remaining_ORFs)]

        ORF_IDs = _chain(*ORF_IDs)
        # fast way of transferring a subset of fasta file entries to another file by ID
        # uses memory map of initial fasta and writes appropriate slices to new file
        self._transfer_ORFs2fasta(ORF_IDs, sourcefile, outfile)

        from baga.Dependencies import dependencies
        d = dependencies['cd-hit']['default_source']
        exe_cdhit = _os.path.join(dependencies['cd-hit'][d]['exe']['main'])

        max_processes = _decide_max_processes( max_cpus )
        max_memory = max_memory*1000

        # the maximum distance for broad CD-hit clustering needs to cluster 
        # all potential orthologs and inparalogs. Insufficent clustering here 
        # guarantees false negative orthologs later

        broad_clusters = []
        first = True
        for identity, word_length, length in iterations:
            if first:
                infilename = '{}/reduced.faa'.format(self.folders['AAs'])
                outfile = '{}/broad_{}'.format(self.folders['AAs'],identity)
                first = False
            else:
                infilename = outfile
                outfile = '{}/broad_{}'.format(self.folders['AAs'],identity)
            cmd = [exe_cdhit,
                    '-i',infilename,
                    '-o',outfile,
                    '-c',str(identity),
                    '-aL',str(length),
                    '-n',str(word_length),
                    '-M',str(max_memory),
                    '-d','0', # omit description in .clstr output file
                    '-T',str(max_processes)]
            
            self.logger.log(PROGRESS, 'Running CD-Hit broad clustering on all '\
                    'input ORFs in {}.faa, output to {}.faa ({:.01%} identity, '\
                    '{:.01%} length)'.format(infilename, outfile, identity, length))
            proc_returncode,stdout = self._launch_external(cmd, 'cd-hit', 
                    self.log_folder, main_logger = self.logger)
            self.logger.info('Parsing CD-Hit clusters . . .')
            these_clusters = []
            this_cluster = []
            for line in open('{}.clstr'.format(outfile)):
                if line.startswith('>'):
                    if len(this_cluster):
                        these_clusters += [this_cluster]
                        this_cluster = []
                else:
                    this_cluster += [line.split('>')[1].split('...')[0]]
            
            # add last cluster
            if len(this_cluster):
                these_clusters += [this_cluster]
            
            # merge with smaller clusters from previous round to maintain full
            # membership throughout iterations
            broad_clusters = self._union_find(broad_clusters + these_clusters)

        ORF2broad_cluster = {}
        for n,cluster in enumerate(broad_clusters):
            for member in cluster:
                ORF2broad_cluster[tuple(map(int,member.split('__')))] = n

        self.logger.info('There were {} final broad clusters'.format(
                len(ORF2broad_cluster)))


        ## update ORF broad cluster affiliation in metadata store ##
        # slice in bytes string value
        broad_cluster_num_sl = self.ORFs_bytes_slices['broad_cluster_num']
        # encoding letter for struct: e.g. H for uint32
        broad_cluster_num_enc = \
                self.ORFs_struct_value[self.ORFs_struct_slices['broad_cluster_num']]
        if store_to_disk:
            # update metadata in LMDB with ORF broad cluster affiliation
            # is updating entire database in one transaction OK?
            # or should this be done genome by genome?
            # ==> might be better to create a new database here then delete the old
            # one because LMDB never overwrites data: original DB will just grow . .
            lmdb_env = _lmdb.Environment(self.folders['metadata'], 
                    **self.lmdb_env_options)
            ORFs_metadata = lmdb_env.open_db(b'ORFs')
            with lmdb_env.begin(db = ORFs_metadata, write = True) as transaction:
                cursor = transaction.cursor()
                for k_bytes,v_bytes in cursor:
                    ORF = _struct.unpack(self.ORFs_struct_key, k_bytes)
                    try:
                        broad_cluster_num = ORF2broad_cluster[ORF]
                        # assign new values via a mutable bytes array
                        new_v_bytes = bytearray(v_bytes)
                        new_v_bytes[broad_cluster_num_sl] = \
                                _struct.pack(broad_cluster_num_enc, broad_cluster_num)
                        transaction.put(k_bytes, bytes(new_v_bytes))
                    except KeyError:
                        pass
            lmdb_env.close()
        else:
            # update self.ORFs_metadata with ORF broad cluster affiliation
            for ORF,broad_cluster_num in ORF2broad_cluster.items():
                k_bytes = _struct.pack(self.ORFs_struct_key, *ORF)
                new_v_bytes = bytearray(self.ORFs_metadata[k_bytes])
                new_v_bytes[broad_cluster_num_sl] = \
                        _struct.pack(broad_cluster_num_enc, broad_cluster_num)
                self.ORFs_metadata[k_bytes] = new_v_bytes
            
            # alternatively just use a separate dict directly for broad clusters
            # (use slightly more memory - cannot resume later)
            # single dict or key-value store for all ORF metadata should save on
            # lookups later, but updating it here requires lookups anyway . . .
            # . . . but dict and LMDB lookups relatively cheap?
            # self.broad_clusters = ORF2broad_cluster

        del ORF2broad_cluster


    def extract_RAA_kmers(self, k = 5, nprocs = 1):
        '''Extract k-mers in reduced amino acid alphabet ORF sequences

        For the genomes that are about to be pair-wise compared by best-match, 
        load sequences and extract k-mers. Frequencies of each k-mers are not 
        counted, each k-mer is simply added to a set (roaringbitmap) per ORF.

        Parameters
        ----------
        k : int
            k-mer (word) length, up to 19.
        '''

        self.logger.info('Extracting k-mers in reduced amino acid alphabet ORF '\
                'sequences for round {} of best match testing'.format(
                len(self.bm_tested_genomes)+1))

        if k <= 2:
            dtype = _np.uint8  # (0 to 255)
        elif 2 < k <= 4:
            dtype = _np.uint16 # (0 to 65535)
        elif 4 < k <= 9:
            dtype = _np.uint32 # (0 to 4294967295)
        elif 9 < k <= 19:
            dtype = _np.uint64 # (0 to 18446744073709551615)
            assert dtype != _np.uint64, 'RoaringBitmaps do not support k > 9. '\
                    'Would have to use conventional Set objects but that is not '\
                    'currently implemented'

        self.logger.debug('Selected data type "{}" for k = {}'.format(dtype, k))
        # for converting list of 0-9 integers into single integer
        k_multiplier = _np.fromiter([10**n for n in range(k-1,-1,-1)], dtype = dtype)


        ### attempts at getting this picklable for multiprocessing k-mer extraction <==
        # cannot even pickle this dummy function for multiprocessing?!
        def extract_kmers(genome_num, ORF_ranges, job_num, jobs_total, path_to_array):
            return('a','b','c','d')
        # tried to remove all references to self and other objects that might contain
        # nested functions but even dummy function above cannot be pickled . . .
        def extract_kmers_parallel(genome_num, ORF_ranges, job_num, jobs_total, path_to_array):
            # select just those RAAs required from carray
            # extract k-mers and put into roaringbitmaps
            ORF_ranges.sort(key = lambda x: x[1][0])
            # read off disk via bcolz carray: single slice for whole genome
            # (includes some not needed but saves 1000s slices)
            genome_s,genome_e = ORF_ranges[0][1][0],ORF_ranges[-1][1][1]
            with _bcolz.open(path_to_array, mode = 'r') as RAAs_carray:
                genome_contiguous_RAAs = RAAs_carray[genome_s:genome_e]
            # adjust per ORF ranges to within this slice
            new_ranges = [(AA_s-genome_s,AA_e-genome_s) for ORF_ID,(AA_s,AA_e) in ORF_ranges]
            # extract k-mers by striding numpy array, straight to roaringbitmap
            # rb_kmers = [_rrbitmap((rolling_window(genome_contiguous_RAAs[s:e], k)*m).sum(1, dtype = dtype)) 
                    # for s,e in new_ranges]
            rb_kmers = 'dummy'
            # collect full lengths: number k-mers will be less than sequence length
            # i.e., repeat insensitive . . good or bad?
            ORF_lens = [(AA_e - AA_s) for ORF_ID,(AA_s,AA_e) in this_genome]
            ORFids = [ORF_ID for ORF_ID,AA_se in this_genome]
            return(rb_kmers, ORF_lens, ORFids, genome_num)


        ## used
        def extract_kmers(genome_num, job_num, jobs_total):
            self.logger.log(PROGRESS, 'Collecting info for {}-mer extraction of genome {} ({} of {})'\
                    ''.format(k, genome_num, job_num, jobs_total))
            # collect ranges from LMDB
            total_ORFs = self.ORFs_per_genome[genome_num]
            genome_start = _struct.pack(ORFs_s_k, genome_num, 0)
            ORF_ords = []
            ORFids = []
            ORF_broad_clusters = []
            # make a set of key bytes we are *not* interested in
            # (handled by narrow clusters)
            narrow_clustered_ORFs = set([_struct.pack(ORFs_s_k, *ORF) for ORF in \
                    self.reduced_db['clustered_by_genome'][genome_num]])
            with lmdb_env.begin(db = ORFs_metadata) as txn:
                cursor = txn.cursor()
                if cursor.set_key(genome_start):
                    # seek to start of first genome: collect all,
                    # iterating probably cheap especially on big DB.
                    # cf fetch per ORF versus building small dict then
                    # querying it
                    itr = _islice(iter(cursor), total_ORFs)
                    for k_bytes,v_bytes in itr:
                        if k_bytes not in narrow_clustered_ORFs:
                            # only record ORF if it's a representative of a narrow cluster
                            # or not affected by narrow clusters
                            # (might be quicker to grab all data first, retain some, then unpack?)
                            # (or to get each key for required ORFs: make LMDB do more lookup work?)
                            genome_num, ORF_num = _struct.unpack(ORFs_s_k, k_bytes)
                            AA_s,AA_e = _struct.unpack(ORFs_s_v[AA_s_sl], 
                                    v_bytes[AA_b_sl])
                            ORF_ords += [(AA_s,AA_e)]
                            ORFids += [(genome_num,ORF_num)]
                            ORF_broad_clusters += [_struct.unpack(ORFs_s_v[BC_s_sl], 
                                    v_bytes[BC_b_sl])[0]]
                else:
                    e = "failed to find genome {} in ORF DB".format(genome_num)
                    self.logger.error(e)
                    raise LookupError(e)
            
            self.logger.debug('Extracting {}-mers for {} ORFs in genome {} ({} of {})'\
                    ''.format(k, len(ORFids), genome_num, job_num, jobs_total))
            # select just those RAAs required from carray
            # extract k-mers and put into roaringbitmaps
            # read off disk via bcolz carray: single slice for whole genome
            # (includes some not needed but saves 1000s slices)
            genome_s = min([s for s,e in ORF_ords])
            genome_e = max([e for s,e in ORF_ords])
            genome_contiguous_RAAs = RAAs_carray[genome_s:genome_e]
            # adjust per ORF ranges to within this slice
            new_ranges = [(AA_s-genome_s,AA_e-genome_s) for (AA_s,AA_e) in ORF_ords]
            # extract k-mers by striding numpy array, straight to roaringbitmap
            rb_kmers = [_rrbitmap((rolling_window(genome_contiguous_RAAs[s:e], 
                    k)*k_multiplier).sum(1, dtype = dtype)) for s,e in new_ranges]
            # collect full lengths: number k-mers will be less than sequence length
            # i.e., repeat insensitive . . good or bad?
            ORF_lens = [(AA_e - AA_s) for AA_s,AA_e in ORF_ords]
            return(rb_kmers, ORF_lens, ORFids, ORF_broad_clusters, genome_num)

        these_genomes = self.bm_selected_genomes[0]

        # define rolling_window at "top-level" for pickling to multiprocessing
        # also defined as self.rolling_window() need a trick for this
        def rolling_window(array, k):
            shape = array.shape[:-1] + (array.shape[-1] - k + 1, k)
            strides = array.strides + (array.strides[-1],)
            return(_np.lib.stride_tricks.as_strided(array, shape = shape, 
                    strides = strides))

        def lmdb_env_begin(lmdb_env):
            return(lmdb_env.begin())

        # open RAA DB
        RAAs_carray = _bcolz.open(self.folders['RAAs_bcolz'], mode = 'r')
        # info for parsing LMDB metadata database
        ORFs_s_k = self.ORFs_struct_key
        ORFs_s_v = self.ORFs_struct_value
        AA_s_sl = self.ORFs_struct_slices['AA']
        AA_b_sl = self.ORFs_bytes_slices['AA']
        BC_s_sl = self.ORFs_struct_slices['broad_cluster_num']
        BC_b_sl = self.ORFs_bytes_slices['broad_cluster_num']

        i = 0
        # open existing lmdb
        with _lmdb.Environment(self.folders['metadata'], readonly = True, 
                create = False, max_dbs = 4) as lmdb_env:
            # used in extract_kmers() below
            ORFs_metadata = lmdb_env.open_db(b'ORFs')
            # save ORF k-mers to list of roaringbitmaps for multiroaringbitmap
            bm_group_rb_kmers = []
            # with corresponding lengths
            # ORFs to be compared (minus CD-HIT represented)
            # total ORFs involved in this analysis are those not consumed in the
            # narrow clusters for the genomes in this group
            num_ORFs_for_comparisons = 0
            for genome_num in these_genomes:
                if genome_num in self.genomes_not_representing:
                    # some genomes entirely removed by narrow clusters
                    # (v. similar to another genome)
                    continue
                ORFs_in_this_genome = self.ORFs_per_genome[genome_num] - \
                        len(self.reduced_db['clustered_by_genome'][genome_num])
                num_ORFs_for_comparisons += ORFs_in_this_genome
            
            bm_group_ORF_lens = _np.empty((num_ORFs_for_comparisons,), dtype=_np.uint16)
            # ORF broad cluster affiliations
            bm_group_ORF_broad_clusters = \
                    _np.empty((num_ORFs_for_comparisons,), dtype=_np.uint32)
            # record genome ranges within array and mrb
            bm_group_genome_ranges = {}
            # record ORF IDs corresponding to genome ranges within array and mrb
            #bm_group_genome_ORFids = {}
            bm_group_genome_ORFids = _np.empty((num_ORFs_for_comparisons,2), \
                    dtype=_np.uint32)
            # a function here would need:
            # self.bm_selected_genomes, num_ORFs, lmdb_env
            if nprocs == 1:
                for n,genome_num in enumerate(these_genomes):
                    if genome_num in self.genomes_not_representing:
                        # some genomes entirely removed by narrow clusters
                        # (v. similar to another genome)
                        continue
                    # also get broad cluster from LMDB
                    rb_kmers, ORF_lens, ORFids, ORF_broad_clusters, genome_num_done = \
                            extract_kmers(genome_num, job_num = n+1, 
                            jobs_total = len(these_genomes))
                    bm_group_rb_kmers += rb_kmers
                    bm_group_ORF_lens[i:i+len(ORF_lens)] = ORF_lens
                    bm_group_ORF_broad_clusters[i:i+len(ORF_lens)] = ORF_broad_clusters
                    bm_group_genome_ranges[genome_num] = i, i+len(ORF_lens)
                    ## makes more sense to keep IDs in same ordered list as bitmap etc
                    ## not in dict of ID lists paired with dict of genome ranges? (as above)
                    bm_group_genome_ORFids[i:i+len(ORF_lens)] = ORFids
                    i += len(ORF_lens)
                    #bm_group_genome_ORFids[genome_num] = ORFids
            elif nprocs > 1:
                raise NotImplementedError("this is not working: cannot make a picklable extract_kmers()")
                # this would need updating to serial version if used
                use_these_genomes = [genome_num \
                        for genome_num in these_genomes if \
                        genome_num not in self.genomes_not_representing]
                all_ORF_ranges = []
                for n,genome_num in enumerate(these_genomes):
                    if genome_num in self.genomes_not_representing:
                        # some genomes entirely removed by CD-HIT
                        # this code would not be required if CD-HIT before PDA
                        continue
                    total_ORFs = self.ORFs_per_genome[genome_num]
                    use_ORFs = set(self.reduced_db['reps_by_genome'][genome_num])
                    # collect ranges from LMDB
                    genome_start = _struct.pack(ORFs_s_k, genome_num, 0)
                    this_genome = []
                    with lmdb_env.begin() as txn:
                        cursor = txn.cursor()
                        if cursor.set_key(genome_start):
                            # seek to start of first genome: collect all,
                            # iterating probably cheap especially on big DB.
                            # cf fetch per ORF versus building small dict then
                            # querying it
                            itr = _islice(iter(cursor), total_ORFs)
                            for k_bytes,v_bytes in itr:
                                genome_num,ORF_num = _struct.unpack(ORFs_s_k, k_bytes)
                                if (genome_num,ORF_num) in use_ORFs:
                                    # only record info if it's an ORF we are interested in
                                    AA_s,AA_e = _struct.unpack(ORFs_s_v[AA_s_sl], 
                                            v_bytes[AA_b_sl])
                                    this_genome += [((genome_num,ORF_num), (AA_s,AA_e))]
                        else:
                            e = "failed to find genome {} in ORF DB".format(genome_num)
                            self.logger.error(e)
                            raise LookupError(e)
                    all_ORF_ranges += [this_genome]
                
                #lmdb_envs = [lmdb_env]*jobs_total
                #RAAs_carrays = [RAAs_carray]*jobs_total
                job_nums = list(range(1,jobs_total+1))
                # record the order they come out in, could vary for concurrancy
                genome_order = []
                # use processes not threads because of GIL
                with _ProcessPoolExecutor(max_workers=nprocs) as executor:
                    for rb_kmers, ORF_lens, ORFids, genome_num in executor.map(
                            extract_kmers_parallel, 
                            use_these_genomes,
                            all_ORF_ranges,
                            job_nums, 
                            [jobs_total]*jobs_total,
                            [self.folders['RAAs_bcolz']]*jobs_total):
                        bm_group_rb_kmers += rb_kmers
                        bm_group_ORF_lens[i:i+len(ORF_lens)] = ORF_lens
                        bm_group_genome_ranges[genome_num] = i, i+len(ORF_lens)
                        i += len(ORF_lens)
                        bm_group_genome_ORFids[genome_num] = ORFids
                        genome_order += [genome_num]
                # update genome order if it changed
                these_genomes = genome_order
                ######## test with genome_order <== just a copy if in serial

        # create multiroaringbitmap per genome best match group
        bm_genome_group_mrb_kmers = _mrbitmap(bm_group_rb_kmers)
        del bm_group_rb_kmers

        # close
        RAAs_carray.free_cachemem()
        # del RAAs_carray # scoping issues for deleting; could be just Python 2

        self.bm_genome_group_mrb_kmers = bm_genome_group_mrb_kmers
        self.bm_group_genome_ORFids = bm_group_genome_ORFids
        self.bm_group_ORF_lens = bm_group_ORF_lens
        self.bm_group_ORF_broad_clusters = bm_group_ORF_broad_clusters
        self.bm_group_genome_ranges = bm_group_genome_ranges
        self.kmer_collections_dtype = dtype
        self.k = k


    def genome_pw_best_matches(self, lower_len_lim = 0.85, 
            upper_len_lim = 1.15, jacc_dist_max = 0.8, 
            nprocs = 1):
        '''Calculate shared k-mer ORF to ORF distances between genomes

        Parameters
        ----------
        lower_len_lim : float
            proportional length difference to allow sequence comparison
        upper_len_lim : float
            proportional length difference to allow sequence comparison
        '''


        # Selecting orthologs by chromosome position later.
        # Distances below RAA detection threshold considered not significantly
        # different with respect to stochastic mutation+selection processes i.e.,
        # not enough evidence to declare an out-paralog: select ortholog by
        # chromosome position.

        ## values around 0.99 might appear by chance depending on length of each sequence
        ## ==> could be calced or just measured empirically?

        self.logger.info('Calculating ORF distances between {} genomes, round {}'\
                ''.format(len(self.bm_group_genome_ranges), 
                len(self.bm_tested_genomes)+1))

        self.logger.log(PROGRESS, 'Preparing jaccard distance pairs for round '\
                '{}'.format(len(self.bm_tested_genomes)+1))

        genome_ranges_sorted = sorted(self.bm_group_genome_ranges.items())

        # pairwise bi-directional genome comparisons
        total_comparisons = 0
        for n,(genome1,(s1,e1)) in enumerate(genome_ranges_sorted[:-1]):
            for genome2,(s2,e2) in genome_ranges_sorted[n+1:]:
                total_comparisons += (e1-s1)*(e2-s2)

        total_ORFs = sum([(e-s) for g,(s,e) in genome_ranges_sorted])

        # bcolz quite slow here, could try snappy at low ratio if memory becomes an issue?
        #all_A_lens = _bcolz.carray(_np.zeros((full_lens,),dtype=_np.uint32), expectedlen=full_lens)
        #all_B_lens = _bcolz.carray(_np.zeros((full_lens,),dtype=_np.uint32), expectedlen=full_lens)

        # uint32 up to 4,294,967,295 ORFs <== this would need uint64 first
        if total_ORFs > 18446744073709551615:
            self.logger.error('Too many ORFs among genomes ({:,})! Compare fewer '\
                    'genomes per batch . . .'.format(total_ORFs))
            raise NotImplementedError

        # prepare arrays for pairwise indices for all input ORFs
        A_indexes = _array('L')
        B_indexes = _array('L')

        # slower python array extension
        # but much more efficient for memory
        genomeAs = _array('L')
        genomeBs = _array('L')
        for n,(genomeA,(s1,e1)) in enumerate(genome_ranges_sorted[:-1]):
            lowers = self.bm_group_ORF_lens[s1:e1]*lower_len_lim
            uppers = self.bm_group_ORF_lens[s1:e1]*upper_len_lim
            b_clusters = self.bm_group_ORF_broad_clusters[s1:e1]
            for genomeB,(s2,e2) in genome_ranges_sorted[n+1:]:
                # for this genome pair, collect ORF pairs for comparison
                # restrict to being within lower-upper limit
                for i,lower,upper,b_cluster in zip(range(s1,e1),lowers,uppers,b_clusters):
                    Bs_before = len(B_indexes)
                    B_indexes.extend(_np.where(
                            # length and broad cluster restrictions for ORF-to-ORF tests
                            # 0 is unassigned broad clusters for short ORFs
                            (lower < self.bm_group_ORF_lens[s2:e2]) & \
                            (self.bm_group_ORF_lens[s2:e2] < upper) & \
                            ((self.bm_group_ORF_broad_clusters[s2:e2] == b_cluster) | \
                                    (self.bm_group_ORF_broad_clusters[s2:e2] == 0))
                            )[0] + s2)
                    Bs_after = len(B_indexes)
                    num_pairs = Bs_after-Bs_before
                    A_indexes.extend([i]*num_pairs)
                    genomeAs.extend([int(genomeA)]*num_pairs)
                    genomeBs.extend([int(genomeB)]*num_pairs)

        # some ORFs not even considered for pair-wise tests
        # too diverse or differing too much in length
        all_indexes = set(A_indexes)
        all_indexes.update(B_indexes)
        not_tested_indexes = set(range(self.bm_group_genome_ORFids.shape[0])) - \
                all_indexes

        # calculate jaccard distances
        self.logger.log(PROGRESS, 'Calculating {:,} jaccard distance pairs for round '\
                '{}'.format(len(genomeAs), len(self.bm_tested_genomes)+1))
        self.logger.debug('Would have been {:,} pairs for all ORF lengths and without '\
                'broad clustering for round {}'.format(total_comparisons, 
                len(self.bm_tested_genomes)+1))

        # use function here? (see reduce 2)
        if nprocs == 1:
            results = _np.frombuffer(
                    self.bm_genome_group_mrb_kmers.jaccard_dist(A_indexes,B_indexes), 
                    dtype = _np.float64).astype(_np.float16)
        elif nprocs > 1:
            prechunksize = int((len(A_indexes)+nprocs)/nprocs)
            # slightly more efficient to create separate arrays above?
            A_indexes_list = (A_indexes[p*prechunksize:(p+1)*prechunksize] for p in range(nprocs))
            B_indexes_list = (B_indexes[p*prechunksize:(p+1)*prechunksize] for p in range(nprocs))
            # jaccard_dist() has nogil so threading works
            # allows efficient use of memory and cores
            with _ThreadPoolExecutor(max_workers=nprocs) as executor:
                results = _np.concatenate([_np.frombuffer(a,dtype = _np.float64).astype(_np.float16) \
                        for a in executor.map(self.bm_genome_group_mrb_kmers.jaccard_dist, 
                        A_indexes_list, B_indexes_list, chunksize=1)])

        self.logger.log(PROGRESS, 'Parsing {:,} jaccard distances for round '\
                '{}'.format(len(genomeAs), len(self.bm_tested_genomes)+1))

        # retain only "hits" below threshold distance (would need calibrating
        # to AA pID distance)
        if nprocs == 1 or not have_numexpr:
            hits = results < jacc_dist_max
        else:
            # numexpr quicker on 2 or more cores
            hits = _ne.evaluate('results < jacc_dist_max')

        # (uint32 for ORF indexes of best match allows millions of inputs and
        # matches 'L' used above)
        # these arrays contain ORF indicies in original multi-roaringbitmap
        # list-like input. _np.frombuffer mismatched types causes length mismatch.
        # array.array('L',[]) == _np.uint64
        # array.array('H',[]) == _np.uint16
        # buffer interface is very fast and implemented for Python arrays
        # _np.frombuffer(array.array) >>>> _np.fromiter(array.array)

        A_indexes = _np.frombuffer(A_indexes, 
                dtype = _np.uint64)[hits].astype(_np.uint32)
        B_indexes = _np.frombuffer(B_indexes, 
                dtype = _np.uint64)[hits].astype(_np.uint32)
        # these indicate which items in As and Bs correspond to which genomes
        # allow for millions of genomes with uint32 :-)
        genomeAs = _np.frombuffer(genomeAs, 
                dtype = _np.uint64)[hits].astype(_np.uint32)
        genomeBs = _np.frombuffer(genomeBs, 
                dtype = _np.uint64)[hits].astype(_np.uint32)

        results = results[hits]

        # record indexes of ORFs that never hit another
        miss_indexes = set(range(self.bm_group_genome_ORFids.shape[0]))
        miss_indexes -= set(A_indexes)
        miss_indexes -= set(B_indexes)

        ### parse results from 1-d array and delineate out-paralogous groups
        ### (orthologs, inparalogs and xenologs within) per genome pair

        # these relationships may change with variations of input genomes
        # because of the non-transitive character of orthology

        def collect_tophits(from_hits):
            tops = []           # shortest A to B hit (possibly draws)
            nontops = {}        # other A to B hits but not shortest (ordered and could also be draws)
            zeros = set()       # those shortest A to B hits that are zero distance (no chance of inparalogs)
            for A_ORF,Bs in from_hits.items():
                # get increasing unique distances
                incr_distances = sorted(set(Bs.values()))
                closest = incr_distances[0]
                # collect all ORFs with top score (could be a draw for top match)
                tops += sorted([(A_ORF,B_ORF) for B_ORF,dist in Bs.items() if dist == closest])
                if closest == 0:
                    zeros.add(A_ORF)
                # collect other ORFs matching increasing scores in order
                these_nontops = []
                for this_dist in incr_distances[1:]:
                    these_nontops += [sorted([B_ORF for B_ORF,dist in Bs.items() if dist == this_dist])]
                if len(these_nontops):
                    nontops[A_ORF] = these_nontops
            return(tops,zeros,nontops)

        genome_sorted = [g for g,r in genome_ranges_sorted]
        # list of lists for union find (orthologous and inparalogous)
        recip_best_hits = []
        # ID to index in multiroaringbitmap for collecting cluster union and intersections
        ORF2BM = {}
        # record these for preparing mrb retest
        zero_draws = {}
        # for collecting k-mers
        all_draw_ORFs = set()
        # collect info to see if dist to paralogs is shorter than to ortholog
        # if yes: inparalog
        # if no: outparalog
        # if close: phylogenetic test in wider context?
        tops = {}
        nontops = {}
        para_pairs = []

        ## somewhere around here might want to assess draws wrt bootstrapped distances . . .
        ## e.g. when is a DNA distance significantly greater than zero?

        ## this section meanders a bit too much between Python dicts and NumPy arrays
        ## one day might benefit from sticking to arrays . . .

        for n,genomeA in enumerate(genome_sorted[:-1]):
            #start_i_A = self.bm_group_genome_ranges[genomeA][0]
            genomeA_ORFids = self.bm_group_genome_ORFids[
                    slice(*self.bm_group_genome_ranges[genomeA])]
            these_tops = {}
            these_nontops = {}
            for genomeB in genome_sorted[n+1:]:
                #start_i_B = self.bm_group_genome_ranges[genomeB][0]
                #bm_group_genome_ORFids_a[bm_group_genome_ORFids_a[:,0]==genomeB]
                #genomeB_ORFids = self.bm_group_genome_ORFids[genomeB]
                genomeB_ORFids = self.bm_group_genome_ORFids[
                        slice(*self.bm_group_genome_ranges[genomeB])]
                # get hit indexes for this AvB genome combination
                this_genome_pair = (genomeAs==int(genomeA)) & (genomeBs==int(genomeB))
                # get distances for this genome combination
                these_results = results[this_genome_pair]
                # get indexes in multibitmap input for this genome combination
                these_A_indexes = A_indexes[this_genome_pair]
                these_B_indexes = B_indexes[this_genome_pair]
                # iterate through each (unique) A ORF
                # record the hits and their scores
                A_hits = {}
                for A_index in _np.unique(these_A_indexes):
                    # get scores for this A against Bs
                    these_scores = these_results[(these_A_indexes==A_index)]
                    these_indexes = these_B_indexes[(these_A_indexes==A_index)]
                    these_IDs = self.bm_group_genome_ORFids[these_indexes]
                    #these_IDs = [genomeB_ORFids[B - start_i_B] for B in these_indexes]
                    these_hits = {tuple(ID):s for ID,s in zip(these_IDs,these_scores)}
                    #A_hits[genomeA_ORFids[A_index - start_i_A]] = these_hits
                    A_hits[tuple(self.bm_group_genome_ORFids[A_index])] = these_hits
                    # retain relevent ID to index mappings for additional jaccard dist calcs
                    for ID,i in zip(these_IDs,these_indexes):
                        # numpy.uint16 to int for indexing multiroaringbitmap
                        ORF2BM[tuple(ID)] = int(i)
                # iterate through each (unique) B ORF
                # record the hits and their scores
                B_hits = {}
                for B_index in _np.unique(these_B_indexes):
                    # get scores for this B against As
                    these_scores = these_results[(these_B_indexes==B_index)]
                    these_indexes = these_A_indexes[(these_B_indexes==B_index)]
                    these_IDs = self.bm_group_genome_ORFids[these_indexes]
                    #these_IDs = [genomeA_ORFids[A - start_i_A] for A in these_indexes]
                    these_hits = {tuple(ID):s for ID,s in zip(these_IDs,these_scores)}
                    #B_hits[genomeB_ORFids[B_index - start_i_B]] = these_hits
                    B_hits[tuple(self.bm_group_genome_ORFids[B_index])] = these_hits
                    # retain relevent ID to index mappings for additional jaccard dist calcs
                    for ID,i in zip(these_IDs,these_indexes):
                        ORF2BM[tuple(ID)] = int(i)
                A_tops,A_zeros,A_nontops = collect_tophits(A_hits)
                B_tops,B_zeros,B_nontops = collect_tophits(B_hits)
                # A_tops and B_tops: orthologous
                # (multiple entries for RAA identical distances)
                # A_zeros and B_zeros: RAA identical orthologs
                # (refers to all entries in A_tops and B_tops)
                all_zeros = A_zeros | B_zeros
                # from A_nontops: potential inparalogs in B
                # from B_nontops: potential inparalogs in A
                
                ## parse best hits
                # these are orthologs plus maybe paralogs from "tops"
                # they might have non-top hit inparalogs for recruitment: check nontops
                ## may include i) draws identical paralogs at RAA resolution (may have DNA differences)
                ##                     ==> attempt resolution by position later
                ##            ii) separate orthologous groups with out-paralogous between group relationships
                ##           iii) single copy families or with inparalogs resolved from among "nontops"
                these_recip_best_hits = \
                        set((tuple((A,B)) for A,B in A_tops)) & \
                        set((tuple((A,B)) for B,A in B_tops))
                # record which ORFs collected as best matches here for exclusion from 
                # assessment of inparalogs below
                collected_ORFs = set([a for b in these_recip_best_hits for a in b])
                #if checkthese & set(these_A_indexes) or checkthese & set(these_B_indexes):
                # if checkthese & set(ORF2BM):
                    # import pdb; pdb.set_trace()
                
                ### parse non-best hits
                # some might have reciprocal best hits elsewhere
                # i.e., ensure these potential inparalogs are not known orthologs
                # if the best hit was zero distance, these can only be out-paralogs
                # (singletons, at least in this context)
                # Also merge to single dict
                nontops_retained = {}
                singletons = set()
                for I,nontopJs in _chain(A_nontops.items(),B_nontops.items()):
                    keeps = []
                    for theseJs in nontopJs:
                        # ensure ORFs are not best hits elsewhere
                        notbests = sorted(set(theseJs) - collected_ORFs)
                        if len(notbests):
                            keeps += [notbests]
                    if len(keeps):
                        if I in all_zeros:
                            # if no distance between orthologs, keep this as an outparalog
                            # i.e., family of one at this stage because no best hits
                            # elsewhere against this other genome
                            # (or two if identical distance i.e., near-identical ORFs)
                            singletons.update((tuple(keep) for keep in keeps))
                        else:
                            # else it could be an inparalog
                            # keep to compare paralog distance with ortholog distance
                            # could be a draw so keep as list
                            nontops_retained[I] = keeps
                
                # add to cumulative list
                recip_best_hits += (list(pair) for pair in these_recip_best_hits)
                if len(singletons):
                    self.logger.debug('Added {} outparalogs to near-identical orthologs'.format(
                            len(singletons)))
                    # retested in Reduce 2 against ORFs in other groups
                    recip_best_hits += (list(singleton) for singleton in singletons)
                
                ### collect ORFs for testing in or out paralogy after
                # now see if dist to paralogs is shorter than to ortholog
                # if yes: inparalog
                # if no: outparalog
                # if close: phylogenetic test in wider context? <== bootstrap?
                for query,hits in nontops_retained.items():
                    for n,hit in enumerate(hits):
                        # [0] here allows for draws: all have same score so treat the same
                        # use this_hit as key but retain hit for all drawn ORFs
                        this_hit = hit[0]
                        # don't know if query was in genome A or B so try both
                        try:
                            #top_hit = A_tops[query][0]
                            top_hit = [B_ORF for A_ORF,B_ORF in A_tops if A_ORF == query][0]
                            top_hit_dist = A_hits[query][top_hit]
                            this_hit_dist = A_hits[query][this_hit]
                        except IndexError:
                            #top_hit = B_tops[query][0]
                            top_hit = [A_ORF for B_ORF,A_ORF in B_tops if B_ORF == query][0]
                            top_hit_dist = B_hits[query][top_hit]
                            this_hit_dist = B_hits[query][this_hit]
                        # retain:
                        #    query (ortholog 1), 
                        #    top_hit (ortholog 2), 
                        #    hit (potential inparalog(s) to ortholog 2, include all but test one)
                        #    query to top_hit distance <== greater than this for outparalog, less for in-
                        ### use ORF2BM here to collect mrb indexes?
                        para_pairs += [(query,top_hit,hit,top_hit_dist)]
                        
                        # best to avoid logging too deep within loops
                        # self.logger.debug('for query {0} top hit was {1} at {2:.03f} while this '\
                        #        'potential inparalog {3} was {4:.03f}'.format(str(query), 
                        #        str(top_hit), top_hit_dist, str(this_hit), this_hit_dist))

        ### test for inparalogs: inparalogs should be closer than orthologs
        # prepare indices of pairs to measure jaccard distances between
        testA_inparas = _array('L')
        testB_inparas = _array('L')
        ortho_dists = _np.empty((len(para_pairs),), _np.float16)
        for i,(query,top_hit,other_hit,top_hit_dist) in enumerate(para_pairs):
            # not all ORFs per genome retained after CD-HIT so need to look up position
            # and add offset for start of that genome in mrb
            # otherhit is one or more at same distance to ortholog to be tested
            # for in/out paralogy
            testA_inparas.append(ORF2BM[top_hit])
            testB_inparas.append(ORF2BM[other_hit[0]])
            ortho_dists[i] = top_hit_dist

        para_pairs = _np.array(para_pairs, dtype = _np.object)
        para_dists = _np.frombuffer(self.bm_genome_group_mrb_kmers.jaccard_dist(
                testA_inparas, testB_inparas), dtype = _np.float64).astype(_np.float16)

        # possible but unlikely to be same distance
        # if so, should go into same cluster and test at ORF synteny level
        inparalogs = para_dists <= ortho_dists

        ### collect inparalogs
        for i,(query,top_hit,other_hit,top_hit_dist) in enumerate(para_pairs[inparalogs]):
            self.logger.debug('for query {0} top hit was {1} at {2:.03f} while this '\
                    'potential inparalog {3} is only {4:.03f} from the top '\
                    'hit ortholog {1}: it is inparalogous to {3}'.format(str(query), 
                    str(top_hit), top_hit_dist, str(other_hit[0]), para_dists[inparalogs][i]))
            
            recip_best_hits += [[top_hit, inparalog] for inparalog in other_hit]

        ### collect the outparalogs as singletons in this context
        for i,(query,top_hit,other_hit,top_hit_dist) in enumerate(para_pairs[~inparalogs]):
            recip_best_hits += [other_hit]

        # for cluster in recip_best_hits:
            # if not isinstance(cluster[0][0],int):
                # print('C',cluster)
                # import pdb; pdb.set_trace()


        self.logger.log(PROGRESS, 'Creating clusters of homologous ORFs for round '\
                '{}'.format(len(self.bm_tested_genomes)+1))

        ###############
        # some ORFs will have been included in PW but not hit anything
        # others will have been excluded from PW by lengths and/or broad clusters
        # ensure they are included at this stage for between group comparisons
        # for genome,ORFs in self.bm_group_genome_ORFids.items():
            # # narrow clustered already excluded
            # recip_best_hits += ([ORF] for ORF in ORFs)

        # add those not tested as singletons here
        recip_best_hits += [[tuple(self.bm_group_genome_ORFids[i])] for i \
                in not_tested_indexes]

        # add those tested but non-hits as singletons here
        recip_best_hits += [[tuple(self.bm_group_genome_ORFids[i])] for i \
                in miss_indexes]

        clusters = self._union_find(recip_best_hits)

        try:
            self.bm_clusters += [clusters]
        except AttributeError:
            self.bm_clusters = [clusters]


        # record which have been prepared for testing and which remain
        self.bm_tested_genomes += [self.bm_selected_genomes.pop(0)]

    def reduce_clusters(self, lower_len_lim = 0.85, 
            upper_len_lim = 1.15, jacc_dist_max = 0.8, 
            nprocs = 1):
        '''Reduce clusters from each (mapped) genome group.

        Combine by disjoint-set/union-find the clusters created in each analysis 
        of the separate genome groups exploiting the overlap of shared genomes 
        between groups. Then address clusters which have not been tested together 
        because they were in separate genome groups but no common genome contained
        an ortholog by performing pangenome versus pangenome best match analyses.

        Parameters
        ----------
        lower_len_lim : float
            proportional length difference to allow sequence comparison
        upper_len_lim : float
            proportional length difference to allow sequence comparison
        '''


        self.logger.info('Merging groups from separate genome group analyses')
        # this exploits the genome membership overlap between initial groups of genomes
        # ideally, all homologous groups are formed here. Exceptions will be where
        # "linking" genomes have lost an ORF in which case those two clusters need to be
        # tested against each other at this pangenome versus pangenome comparison here
        merged_clusters = self._union_find(
                (a for b in self.bm_clusters for a in b))

        ###### no duplicates after merging best match clusters
        assert len([a for b in merged_clusters for a in b]) == len(set([a for b in merged_clusters for a in b]))

        self.logger.info('Assessing inter-genome group cluster relationships '\
                'for further merging')

        thegroups = [set(group) for group in self.bm_tested_genomes]
        # open RAA DB
        RAAs_carray = _bcolz.open(self.folders['RAAs_bcolz'], mode = 'r')

        # non-overlapping group vs group i.e., pangenome v pangenome of genome groups

        # singletons after first round that could be inparalogous but never tested
        # because no hits within initial group. Can be tested here where both hits are
        # in same genome. 
        putative_para_singletons = set([cluster[0] for cluster in merged_clusters if len(cluster) == 1])

        recip_best_hits = []

        for n,group1 in enumerate(thegroups[:-1]):
            # ===> build all of group 1 bitmaps here
            # then select by length and broad cluster as needed
            for m,group2 in enumerate(thegroups[n+1:],n+1):
                # ===> build all of group 2 bitmaps here
                # then select by length and broad cluster as needed
                # and store in memory or on disk for when this group 2 becomes group 1
                # check effect on memory usage
                self.logger.debug('Collecting required ORFs for comparing '\
                        'genome groups: {} and {}'.format(n+1, m+1))
                batch1 = []
                batch2 = []
                # select appropriate clusters for this group pair
                # any overlap with one but not other includes within-group
                # plus others but exlcudes those spanning both groups which
                # includes core and near-core clusters
                for cluster in merged_clusters:
                    genomes = set(genome for genome,_ in cluster)
                    group1_overlap = genomes & group1
                    group2_overlap = genomes & group2
                    if group1_overlap and not group2_overlap:
                        batch1 += [cluster]
                    elif group2_overlap and not group1_overlap:
                        batch2 += [cluster]
                # collect cluster info and a representative for each cluster
                # in this genome group.
                # (genome, ORF, broad cluster, min, max ORF lengths)
                # Also keep a note of ranges for when eventually collecting
                # RAAs from contiguous carrays
                batch1_cluster_info, batch1_ranges = self._collect_cluster_info(batch1)
                self.logger.debug('batch1_cluster_info.shape: {} by {}'.format(*batch1_cluster_info.shape))
                # adjust min and max cluster sequence lengths for between
                # cluster comparisons in which only those with an overlap of 
                # min-max lengths plus original PW comparison length range are
                # compared (numpy recasting so not using *=)
                batch1_cluster_info[:,3] = batch1_cluster_info[:,3]*lower_len_lim
                batch1_cluster_info[:,4] = batch1_cluster_info[:,4]*upper_len_lim
                # same for batch 2
                batch2_cluster_info, batch2_ranges = self._collect_cluster_info(batch2)
                self.logger.debug('batch2_cluster_info.shape: {} by {}'.format(*batch2_cluster_info.shape))
                batch2_cluster_info[:,3] = batch2_cluster_info[:,3]*lower_len_lim
                batch2_cluster_info[:,4] = batch2_cluster_info[:,4]*upper_len_lim
                # now we have batch1_cluster_info vs batch2_cluster_info
                # containing ORF_ID of cluster representative (0,1),
                # broad_cluster (2), minlen (3), maxlen (4).
                b1_ORFs_ind_for_tests = _array('L')
                b2_ORFs_ind_for_tests = _array('L')
                b1_ORFs_IDs_for_tests = []
                b2_ORFs_IDs_for_tests = []
                # among ORFs for this genome group comparison
                # collect broad clusters
                broad_clusters = _np.unique(_np.concatenate(
                        [batch1_cluster_info[:,2], batch2_cluster_info[:,2]]))
                self.logger.debug('{} broad clusters among ORFs for comparison'\
                        ''.format(len(broad_clusters)))
                # compile a list of ORFs for testing here
                ORFs_for_testing1 = [tuple(ORF) for ORF in batch1_cluster_info[:,:2]]
                ORFs_for_testing2 = [tuple(ORF) for ORF in batch2_cluster_info[:,:2]]
                # ORFs will be in this order in the bitmap
                # (is sorting actually necessary here if order is stored?)
                all_ORFs_order = {ORF:n for n,ORF in \
                        enumerate(sorted(ORFs_for_testing1 + ORFs_for_testing2))}
                # assign bitmap indexes here while building pairwise comparisons
                batch1_ORF_indices = _np.fromiter((all_ORFs_order[ORF] 
                        for ORF in ORFs_for_testing1), _np.uint32)
                batch2_ORF_indices = _np.fromiter((all_ORFs_order[ORF] \
                        for ORF in ORFs_for_testing2), _np.uint32)
                batch1_cluster_info = _np.column_stack(
                        (batch1_cluster_info, batch1_ORF_indices))
                batch2_cluster_info = _np.column_stack(
                        (batch2_cluster_info, batch2_ORF_indices))
                for broad_cluster in broad_clusters:
                    # collect ORFs in each broad cluster for this genome group pair
                    these = batch1_cluster_info[:,2] == broad_cluster
                    g1_this_brdcls = batch1_cluster_info[these,]
                    these = batch2_cluster_info[:,2] == broad_cluster
                    g2_this_brdcls = batch2_cluster_info[these,]
                    if g1_this_brdcls.shape[0] and g2_this_brdcls.shape[0]:
                        # at least some ORFs shared: if lengths close enough,
                        # store for PW comparison
                        for genome_num, ORF_num, bc, min_len, max_len, bitmap_i in \
                                g1_this_brdcls:
                            keeps = g2_this_brdcls[\
                                    (g2_this_brdcls[:,3] < max_len) & \
                                    (g2_this_brdcls[:,4] > min_len)]
                            if keeps.shape[0]:
                                # same batch 1 cluster rep: ORF ID and its index in multibitmap
                                b1_ORFs_ind_for_tests.extend((bitmap_i,)*keeps.shape[0])
                                b1_ORFs_IDs_for_tests += \
                                        [tuple([genome_num, ORF_num])]*keeps.shape[0]
                                # against a selection of batch 2 cluster reps
                                b2_ORFs_ind_for_tests.extend(keeps[:,-1])
                                b2_ORFs_IDs_for_tests += \
                                        [tuple(ORF_id) for ORF_id in keeps[:,:2]]
                
                pairwise_test_IDs = b1_ORFs_IDs_for_tests,b2_ORFs_IDs_for_tests
                pairwise_test_inds = b1_ORFs_ind_for_tests,b2_ORFs_ind_for_tests
                # among the many ORFs excluded from comparisons here are those
                # in other clusters that span into respective groups (would
                # have been already merged), ORFs in a different "broad 
                # cluster" i.e., those we know are non-hits with very different
                # sequence, those with very different lengths. By process of 
                # elimination we have mostly hits left.
                self.logger.debug('Collected {:,} ORF pairs for comparison between '\
                        'groups: {} ({} clusters) and {} ({} clusters)'.format(
                        len(pairwise_test_IDs[0]), n+1, batch1_cluster_info.shape[0], 
                        m+1, batch2_cluster_info.shape[0]))
                
                # we now know which ORFs can be compared between these genome groups
                # (appropriate genomes, broad clusters, length restrictions)
                # build multiroaring bitmap and convert paired IDs to indices
                # in bitmap
                self.logger.debug('Building bitmaps . . .')
                mrb = self._build_mrbs(all_ORFs_order, batch1_ranges, batch2_ranges, RAAs_carray)
                # then apply the pairwise tests using this multiroaring bitmap and retain only hits
                self.logger.debug('Calculating distances . . .')
                # all PW ORF indices and IDs in, hit IDs out with distances
                b1_to_b2_dists, pairwise_hit_IDs = \
                        self._compare_ORF_pairs(mrb, pairwise_test_inds, pairwise_test_IDs, 
                        jacc_dist_max = jacc_dist_max, nprocs = nprocs)
                
                self.logger.debug('Parsing distances . . .')
                # collect info to see if dist to paralogs is shorter than to ortholog
                # if yes: inparalog
                # if no: outparalog
                # if close: phylogenetic test in wider context?
                b1_hits, b2_hits = self._parse_distances(b1_to_b2_dists, pairwise_hit_IDs)
                b1_tops,b1_zeros,b1_nontops = self._collect_tophits(b1_hits)
                b2_tops,b2_zeros,b2_nontops = self._collect_tophits(b2_hits)
                # b1_tops and b2_tops: orthologous
                # b1_zeros and b2_zeros: RAA identical orthologs
                # from A_nontops: potential inparalogs in B
                # from B_nontops: potential inparalogs in A
                
                ## parse hits
                these_recip_best_hits, para_pairs = self._parse_hits(
                        b1_tops, b1_zeros, b1_nontops, b1_hits, 
                        b2_tops, b2_zeros, b2_nontops, b2_hits)
                recip_best_hits += these_recip_best_hits
                ### test for inparalogs: closer than orthologs
                # this is still valid in this pangenome-v-pangenome context because 
                # retains outparalogs as singletons but excluded here because
                # already in
                para_pairs_test = []
                for query,top_hit,other_hit,top_hit_dist in para_pairs:
                    if top_hit[0] == other_hit[0][0]:
                        # only inparalogs if in same genome
                        if top_hit in putative_para_singletons or \
                                any([hit in putative_para_singletons \
                                        for hit in other_hit]):
                            # only need testing here if singletons after
                            # within group analyses
                            para_pairs_test += [(query, top_hit, other_hit, 
                                    top_hit_dist)]
                
                if len(para_pairs_test):
                    recip_best_hits += (pair for pair in self._test_inparalogy(
                            mrb, para_pairs_test, all_ORFs_order) if len(pair) > 1)

        self.reduced_clusters = self._union_find(
                _chain(merged_clusters, recip_best_hits))


    def expand_clusters(self, sort_by_num_genomes = True, 
            sort_by_num_ORFs = False):
        '''
        Expand families containing representatives of CD-HIT clusters

        For very large analyses sorting might be too time consuming but 
        defaults to families spanning the most genomes reported first.

        Parameters
        ----------
        sort_by_num_genomes : bool
            put ORF families in order of total genomes spanned
        sort_by_num_ORFs    : bool
            put ORF families in order of total ORFs contained (ignored if 
            sort_by_num_genomes is True)
        '''

        ## expand from narrow clusters (matching RAA sequences)
        # generate arbitrary family IDs here: index in list
        self.logger.info('Expanding narrow clusters')

        full_clusters = []
        for cluster in self.reduced_clusters:
            full_cluster = []
            for ORF in cluster:
                try:
                    full_cluster += self.reduced_db['expansions'][ORF]
                except KeyError:
                    full_cluster += [ORF]
            
            full_cluster.sort()
            full_clusters += [full_cluster]

        # add the one-per genome complete families with matching RAA sequences
        for mems in self.reduced_db['bm_clusters_within_narrows']:
            full_clusters += [sorted(mems)]

        if sort_by_num_genomes:
            # by number genomes
            full_clusters.sort(key=lambda c: len(set(g for g,_ in c)), 
                    reverse = True)
        elif sort_by_num_ORFs:
            full_clusters.sort(key=len, reverse = True)

        self.full_clusters = full_clusters


    def generate_results(self, generate_fastas_AA = False, 
            generate_fastas_DNA = False, chars_per_line = 50, 
            families_per_write_batch = 100):
        '''Generate cluster results

        For now a simple fasta-like output:

        >C01
        G1__ORF_01 G2__ORF_02 G3__ORF_04
        >C02
        G1__ORF_02 G2__ORF_05 G3__ORF_01

        This uses internal IDs for now, will expand later
        '''

        self.logger.info('Writing output for analysis {} to {}'.format(
                self.name, self.folders['results']))

        ### collect genome metadata: name, input file, number ORFs
        ### write to text file <results folder>/genome_encoding.baga
        ## TODO: include assembly_acc_num <== (needs fixing at other end)
        # for unpacking bytes of total ORFs and genome ID number
        # (could change if total genomes very high? Depend on code at other end)
        # determined near the end of .processORFs() when genome and chromosome 
        # metadata stored

        # these are set in Finder.processORFs()
        # self.chrms_key_struct = 'L'
        # ORF_num to ORF_num: ORF ID range in chromosome from this genome
        # self.chrms_val_struct = 'HH'
        # allow > 65k genomes
        # self.gnoms_key_struct = 'L'
        # total ORFs, chromosome_num to chromosome_num in this genome
        # self.gnoms_val_struct = 'LLH'
        ORF_per_gnom_struct = self.gnoms_val_struct[2:3]
        s = _struct.calcsize(self.gnoms_val_struct[:2])
        e = _struct.calcsize(self.gnoms_val_struct[:3])
        ORFs_per_gnom_sl = slice(s,e)
        gnom_str_data_sl = slice(_struct.calcsize(self.gnoms_val_struct), None)

        filename = '{}/genome_encoding.baga'.format(self.folders['results'])
        max_ORFs = 0
        with _lmdb.Environment(self.folders['metadata'], readonly = True, 
                create = False, max_dbs = 4) as lmdb_env, open(filename, 'w') as fout:
            fout.write('"Genome ID"\t"Total ORFs"\t"Source file"\t"Organism name"\n')
            genomes_metadata = lmdb_env.open_db(b'genomes')
            with lmdb_env.begin() as txn:
                csr = txn.cursor(genomes_metadata)
                # LMDB keys are sorted
                for g,info in csr.iternext(keys = True, values = True):
                    genome_num = _struct.unpack(self.gnoms_key_struct, g)[0]
                    num_ORFs = _struct.unpack(ORF_per_gnom_struct, 
                            info[ORFs_per_gnom_sl])[0]
                    if num_ORFs > max_ORFs:
                        max_ORFs = num_ORFs
                    organism_name, assembly_acc_num , genome_file = \
                            info[gnom_str_data_sl].decode("utf-8").split('|')
                    fout.write('"{:0{}d}"\t{}\t"{}"\t"{}"\n'.format(genome_num, 
                            self.padding, num_ORFs, genome_file, organism_name))

        cluster_padding = len(str(len(self.full_clusters)))
        ORF_padding = len(str(max_ORFs))

        ### write results: clusters ID and ORF membership to text file
        ### <results folder>/families.baga
        filename = '{}/families.baga'.format(self.folders['results'])
        with open(filename, 'w') as fout:
            for cID,members in enumerate(self.full_clusters):
                members_str = ['{:0{}d}.{:0{}d}'.format(g, self.padding, 
                        o, ORF_padding) for g,o in sorted(members)]
                fout.write('>{:0{}d}\n{}\n'.format(cID, cluster_padding, 
                        ' '.join(members_str)))


        ## probably quicker to load all ranges for AA from LMDB into dict or
        ## nested lists because random access to in-memory probably quicker
        ## then write all of AAs, then do same for DNAs.
        ## currently lots of opening and closing and random access of LMDB

        if generate_fastas_AA or generate_fastas_DNA:
            num_batches = round((len(self.full_clusters)/families_per_write_batch)+0.5)
            AAs_carray = _bcolz.open(self.folders['AAs_bcolz'], mode = 'r')
            DNAs_carray = _bcolz.open(self.folders['DNAs_bcolz'], mode = 'r')
            for i in range(num_batches):
                s,e = i*families_per_write_batch,(i+1)*families_per_write_batch
                cluster_members = self.full_clusters[s:e]
                cluster_IDs = ['{:0{}d}'.format(c, cluster_padding) \
                        for c in range(s,e)]
                # get ORF info for this batch of clusters
                # cl_locus_IDs are external IDs and may not indicate genome
                cl_ranges_DNA, cl_ranges_AA, cl_locus_IDs = \
                        self._collect_cluster_ranges(cluster_members)
                if generate_fastas_AA:
                    self.logger.log(PROGRESS, 'Writing ORF families protein sequences in batch '\
                            '{} of {}'.format(i+1, num_batches))
                    self._write_clusters(cluster_IDs, cluster_members, cl_locus_IDs, 
                            cl_ranges_AA, AAs_carray, ORF_padding, 'faa', chars_per_line)
                if generate_fastas_DNA:
                    self.logger.log(PROGRESS, 'Writing ORF families nucleotide sequences in batch '\
                            '{} of {}'.format(i+1, num_batches))
                    self._write_clusters(cluster_IDs, cluster_members, cl_locus_IDs, 
                            cl_ranges_DNA, DNAs_carray, ORF_padding, 'fna', chars_per_line)

        if generate_fastas_AA:
            AAs_carray.free_cachemem(); del AAs_carray

        if generate_fastas_DNA:
            DNAs_carray.free_cachemem(); del DNAs_carray




    def do_de_novo(self, max_genomes_per_group = 10, genome_pairs_shared = 2, 
            bm_lower_len_lim = 0.8, bm_upper_len_lim = 1.25, k = 5, 
            jacc_dist_max = 0.65, max_cpus = -1, max_memory = 8, 
            retain_individual_sketches = False, performchecks = True):
        '''Perform a complete inference of homologous groups of protein coding genes'''

        # make a minhash tree using mash distances and dendropy NJ
        self.make_genome_tree(
                retain_individual_sketches = retain_individual_sketches)
        # could report a tree at this point?
        # print(fam_finder_S.mash_NJ_tree.as_ascii_plot(plot_metric='length'))
        # select groups of genomes for separate analyses
        self.phylo_select_genomes(max_genomes_per_group, genome_pairs_shared)

        self.process_narrow_clusters()

        # could add assertions: see do_de_novo_old()

        ## need to decide on when/how to decide on CPUs/processes etc
        nprocs = _decide_max_processes(max_cpus)

        self.AA_DB_broad_clustering(
                iterations = [(0.80,5,0.6), (0.68,4,0.6), (0.56,4,0.6)], 
                store_to_disk = True, max_cpus = 1, max_memory = 8)

        while len(self.bm_selected_genomes):
            # collects for genomes in fam_finder.bm_selected_genomes[0]
            # could be quicker by putting ranges in two column table and
            # vectorise k-mer counting?
            # >1 procs/cpus/threads/cores not yet implemented
            self.extract_RAA_kmers(k, nprocs = 1)
            self.genome_pw_best_matches(lower_len_lim = bm_lower_len_lim, 
                    upper_len_lim = bm_upper_len_lim, jacc_dist_max = jacc_dist_max, 
                    nprocs = nprocs)

        # Reduce as in MapReduce. Includes some more ORF comparisons
        self.reduce_clusters(lower_len_lim = bm_lower_len_lim, 
                upper_len_lim = bm_upper_len_lim, jacc_dist_max = jacc_dist_max, 
                nprocs = nprocs)

        # expand narrow, high identity clusters
        self.expand_clusters()

if __name__ == '__main__':
    main()
