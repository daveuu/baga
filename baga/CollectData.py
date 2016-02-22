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
CollectData module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains functions to collect reference genomes and query read sets 
from a the Internet via a URL or from a local path.
'''

### Import Dependencies for this module ###

# use underscore prefix to keep these private

# stdlib
from ftplib import FTP as _FTP
from time import sleep as _sleep
import urllib2 as _urllib2
from glob import glob as _glob

# not sure if there is any advantage in going via the package __init__.py for 
# widely used imports
from baga import _cPickle
from baga import _gzip
from baga import _os
from baga import _sys
from baga import _re
from baga import _tarfile
from baga import _json
from baga import _StringIO
from baga import _array
from baga import _time
from baga import _md5

from baga import report_time as _report_time
from baga import MetaSample as _MetaSample
from baga import PROGRESS
from baga import PY3 as _PY3

# external Python modules
from Bio import Entrez as _Entrez
from Bio import SeqIO as _SeqIO
def main():
    pass
class Genome_old:
    '''
    Collect a reference genome sequence from the Internet or a local path
    '''
    def __init__(self,  accession = False, 
                        user_email = False, 
                        local_path = False, 
                        format = False):
        '''
        Option 1: provide "accession" number and "user_email":
            
        Download an annotated genome sequence from NCBI Nucleotide 'nuccore' 
        database. Requires accession number of entry e.g. 'FM209186' for 
        Pseudomonas aeruginosa LESB58 and your email address so NCBI can contact 
        you in the event of a problem. Uses BioPython's Bio.Entrez module.

        Option 2: provide "local_path" and "format" as either "genbank" or "baga".

        If format = "genbank": Load an annotated genome sequence from a local 
        genbank file. For example, one downloaded previously from 
        ftp://ftp.ncbi.nih.gov/genomes/Bacteria/

        If format = "baga": Reload a genome that was previously processed with 
        either of the two options above for use with the Bacterial and Archaeal 
        Genome Analyser.

        Option 3:

        Just provide "accession" to autoload a locally saved genome.
        '''
        
        # sanity checks
        e = ["Provide accession OR local_path and format"]  #,accession, user_email, local_path, format
        e += ["Provide accession and user_email to download version from NCBI via Entrez"]
        e += ["accession = {}, user_email = {}, local_path = {}, format = {}".format(accession, user_email, local_path, format)]
        assert ((accession) or (local_path and format)), '\n'.join(e)
        

        def extractLoci(seq_record):
            '''
            Extract some ORF information from genbank record 
            and put into convenient dictionaries
            '''
            
            ORF_ranges = {}
            large_mobile_element_ranges = {}
            
            GI_prophage = _re.compile('[Ii]sland|[Pp]hage|GI')
            for f in seq_record.features:
                if f.type == 'CDS':
                    # until AA and detailed annotations of ORFs needed, don't save other qualifiers
                    try:
                        thisgene = f.qualifiers['gene'][0]
                    except KeyError:
                        thisgene = ''
                    
                    ORF_ranges[f.qualifiers['locus_tag'][0]] = (f.location.start.position,f.location.end.position,f.location.strand, thisgene)
                    
                if f.type == 'misc_feature':
                    try:
                        # not all misc_features have a "note"
                        feature_note = f.qualifiers['note'][0]
                    except KeyError:
                        continue
                    
                    if _re.search(GI_prophage, feature_note) and f.location.nofuzzy_end - f.location.nofuzzy_start > 10000:
                        large_mobile_element_ranges[f.qualifiers['note'][0]] = f.location.nofuzzy_start, f.location.nofuzzy_end
            
            # filter ORFs within ORFs (artifacts? PLES_21351 and PLES_21361 in LESB58)
            ORF_ranges_sorted = sorted(ORF_ranges.items(), key = lambda x: x[1][0])
            inner_ORFs = set()
            for n, (ORF1, (s1, e1, strnd1, genename1)) in enumerate(ORF_ranges_sorted[:-1]):
                for ORF2, (s2, e2, strnd2, genename2) in ORF_ranges_sorted[(n+1):]:
                    if s1 < s2 and e2 < e1:
                        print('%s is within %s; dumping former' % (ORF2,ORF1))
                        inner_ORFs.add(ORF2)
                    if s2 < s1 and e1 < e2:
                        print('%s is within %s; dumping former' % (ORF1,ORF2))
                        inner_ORFs.add(ORF2)
            
            for ORF in inner_ORFs:
                del ORF_ranges[ORF]
            
            return(ORF_ranges, large_mobile_element_ranges)

        def DL(url, verbose = True):
            req = _urllib2.urlopen(url)
            CHUNK = 16 * 1024 * 32
            data = _StringIO()
            c = 0
            for chunk in iter(lambda: req.read(CHUNK), ''):
                c += CHUNK
                if verbose:
                    print("{:,} bytes".format(c))
                data.write(chunk)
            
            if verbose:
                print('Download complete . . .')
            data.seek(0)
            return(data)

        def getFromEntrez(search_id, user_email):
            '''
            download a genome sequence given a search ID
            
            search_id is recommended to be a refseq or genbank accession number
            or other unambiguous ID that will return a single result
            '''
            
            from Bio.Entrez.Parser import ValidationError as _ValidationError
            
            if '.' in search_id:
                search_id_unversioned,requested_ver = search_id.split('.')
            else:
                search_id_unversioned,requested_ver = search_id,None
            
            if '_' in search_id_unversioned:
                search_id_is_refseq = True
            else:
                search_id_is_refseq = False
            
            _Entrez.email = user_email
            handle = _Entrez.esearch(db = "assembly", term = search_id_unversioned)
            result = _Entrez.read(handle)
            if len(result['IdList']) != 1:
                print('WARNING: Your search ID: "{}" returned {} assembly results '\
                'from ncbi.nlm.nih.gov/assembly but a single result is required.'.format(
                        search_id, len(result['IdList'])))
                raise LookupError
            Assembly_ID = result['IdList'][0]
            handle = _Entrez.esummary(db = "assembly", id = Assembly_ID)
            # some ways of handling unexpected content from NCBI
            try:
                raw = _Entrez.read(handle, validate=True)
                info = raw['DocumentSummarySet']['DocumentSummary'][0]
            except _ValidationError as e:
                print('WARNING: The information about this genome returned by NCBI Entrez failed validation (ValidationError):\n{}'.format(e))
                print('Trying without validation . . .')
                handle = _Entrez.esummary(db = "assembly", id = Assembly_ID)
                info = _Entrez.read(handle, validate=False)['DocumentSummarySet']['DocumentSummary'][0]
            
            print('Found: {} ({})'.format(info['Organism'],info['AssemblyStatus']))
            
            # collect download links
            try:
                genbank_ftp = _re.findall(
                        '<FtpPath type="GenBank">([^<]+)</FtpPath>', 
                        info['Meta'])[0]
                print('Found Genbank link:\n{}'.format(genbank_ftp))
            except IndexError:
                genbank_ftp = False
                print('GenBank link not found')
            
            try:
                refseq_ftp = _re.findall(
                        '<FtpPath type="RefSeq">([^<]+)</FtpPath>', 
                        info['Meta'])[0]
                print('Found RefSeq link:\n{}'.format(refseq_ftp))
            except IndexError:
                refseq_ftp = False
                print('RefSeq link not found')
            
            e = 'Failed to retrieve FTP download links from MetaData:\n{}'.format(info['Meta'])
            assert genbank_ftp or refseq_ftp, e
            
            if refseq_ftp:
                use_link = refseq_ftp
            elif genbank_ftp:
                use_link = genbank_ftp
            
            # collect accessions and versions
            refseq_ass_acc = info['AssemblyAccession']
            e = 'No RefSeq assembly found for {}. You can double check at http://www.ncbi.nlm.nih.gov/assembly'.format(search_id)
            assert refseq_ass_acc[:3] == 'GCF'
            genbank2refseq = {}
            genbank2version = {}
            refseq2genbank = {}
            refseq2version = {}
            data = DL('ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/{}.assembly.txt'.format(
                    refseq_ass_acc), verbose = False)
            ID_info = data.readlines()
            
            for line in ID_info:
                if line[0] != '#' and len(line) > 0:
                    cells = line.split('\t')
                    genbank_acc, gb_ver = cells[4].split('.')
                    refseq_acc, rs_ver = cells[6].split('.')
                    genbank2refseq[genbank_acc] = refseq_acc
                    genbank2version[genbank_acc] = gb_ver
                    refseq2genbank[refseq_acc] = genbank_acc
                    refseq2version[refseq_acc] = rs_ver
            
            if search_id_is_refseq:
                use_name = search_id_unversioned + '.' + refseq2version[search_id_unversioned]
                if requested_ver is None:
                    print('Found version {} of RefSeq accession {}'.format(
                            refseq2version[search_id_unversioned], search_id_unversioned))
                elif requested_ver != refseq2version[search_id_unversioned]:
                    print('RefSeq accession {} version {} was requested, '\
                    'but version {} is the current version and will be used instead'.format(
                            search_id_unversioned, requested_ver, 
                            refseq2version[search_id_unversioned]))
            else:
                use_refseq = genbank2refseq[search_id_unversioned]
                print('Will use RefSeq accession {} (latest version {}) which '\
                'corresponds to provided GenBank accession {}'.format(
                        use_refseq, refseq2version[use_refseq], search_id_unversioned))
                
                use_name = use_refseq + '.' + refseq2version[use_refseq]
            
            ### could collect other replicons in this genome . . .
            if len(refseq2version) > 1:
                print('(this is 1 of {} replicons in this genome)'.format(len(refseq2version)))
            else:
                print('(this is the only replicon in this genome)')
            
            # download checksums
            data = DL(use_link + '/md5checksums.txt', verbose = False)
            checksum = [l.split('  ./') for l in data.readlines() if '_genomic.gbff.gz' in l][0][0]
            # download sequences and annotations
            use_link += '/' + use_link.split('/')[-1] + '_genomic.gbff.gz'
            print('Downloading from:\n{}'.format(use_link))
            data = DL(use_link, verbose = True)
            hasher = _md5()
            buff = data.read(65536)
            while len(buff) > 0:
                hasher.update(buff)
                buff = data.read(65536)
            
            e = '. . . checksum fail!'
            assert hasher.hexdigest() == checksum, e
            print('. . . checksum {} passed!'.format(checksum))
            data.seek(0)
            archive = _gzip.GzipFile(mode="rb", fileobj = data)
            records = list(_SeqIO.parse(archive, 'genbank'))
            for seq_record in records:
                if use_name == seq_record.id:
                    self.ORF_ranges, self.large_mobile_element_ranges = extractLoci(seq_record)
                    self.sequence = _array('c', seq_record.seq)
                    self.id = seq_record.id

        def getFromEntrezNucleotide(accession, user_email):
            print("WARNING: NCBI's Entrez query system used here can be unreliable "\
                    "for data download. If the download does not start (if you "\
                    "don't see '524,288 bytes ...') within a few seconds, press "\
                    "ctrl-c and issue the same command again (up-arrow, enter, "\
                    "usually works)\n")
            
            _Entrez.email = user_email
            handle = _Entrez.efetch(db = "nuccore", rettype = "gb", retmode = "text", 
                    id = accession)
            try:
                records = list(_SeqIO.parse(handle, 'genbank'))
                #seq_record = _SeqIO.read(handle, "genbank")
            except ValueError as error_message:
                print("There was a problem with the genome (accession: {}) downloaded "\
                        "from NCBI via Entrez: {}. Retry because Entrez can be "\
                        "unreliable, or try loading from a .gbk file downloaded "\
                        "manually from e.g., ftp://ftp.ncbi.nih.gov/genomes/Bacteria/"\
                        "".format(accession, error_message))
            handle.close()
            # self.ORF_ranges, self.large_mobile_element_ranges = extractLoci(seq_record)
            # self.sequence = _array('c', seq_record.seq)
            # self.id = seq_record.id
            for seq_record in records:
                if accession in seq_record.id:
                    self.ORF_ranges, self.large_mobile_element_ranges = extractLoci(seq_record)
                    self.sequence = _array('c', seq_record.seq)
                    self.id = seq_record.id

        def loadFrombaga(local_path):
            with _tarfile.open(local_path, "r:gz") as tar:
                for member in tar:
                    contents = _StringIO(tar.extractfile(member).read())
                    try:
                        # either json serialised conventional objects
                        contents = _json.loads(contents.getvalue())
                    except ValueError:
                        # or longer python array.array objects
                        contents = _array('c', contents.getvalue())
                    
                    setattr(self, member.name, contents)

        def loadFromGBK(local_path):
            seq_record = list(_SeqIO.parse(local_path, "genbank"))[0]
            self.ORF_ranges, self.large_mobile_element_ranges = extractLoci(seq_record)
            self.sequence = _array('c', seq_record.seq)
            self.id = seq_record.id
        
        if accession and user_email:
            try:
                getFromEntrez(accession, user_email)
            except LookupError:
                print('Falling back to download from www.ncbi.nlm.nih.gov/nuccore '\
                        'database via NCBI Entrez.')
                getFromEntrezNucleotide(accession, user_email)
            
        elif accession and not user_email:
            success = False
            try:
                tryfilename = 'baga.CollectData.Genome-{}.baga'.format(accession)
                loadFrombaga(tryfilename)
                success = True
            except IOError:
                pass
            
            if not success:
                try:
                    tryfilename = '{}.gbk'.format(accession)
                    loadFromGBK(tryfilename)
                    success = True
                except IOError:
                    pass
            
            if success:
                print('Successfully loaded a genome from: {}'.format(tryfilename))
                print('If you would like to download the NCBI version, provide an email address for use with Entrez.')
            else:
                print('Could not find any local versions of genome with accession number: {}'.format(accession))
                print('If you would like to download the NCBI version, provide an email address for use with Entrez.')
            
        elif local_path and format == 'genbank':
            loadFromGBK(local_path)
            
        elif local_path and format == 'baga':
            loadFrombaga(local_path)
            
        elif local_path and format:
            print('Format "{}" not supported (try "baga" or "genbank")'.format(format))
            
        elif local_path:
            print('Please specify format: "baga" or "genbank")')
    def saveLocal(self, name = False):
        '''
        Save a reference genome to a local compressed baga file. This saves 
        Internet bandwidth if downloading from NCBI and time if loading a 
        genbank file.
        'filename' can exclude extension: .baga will be added
        A .baga file is mostly Python dictionaries in JSON strings and
        array.array objects in a tar.gz format.
        '''
        
        if name:
            fileout = 'baga.CollectData.Genome-{}.baga'.format(name)
        else:
            fileout = 'baga.CollectData.Genome-{}.baga'.format(self.id)
        
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
                elif isinstance(att, dict) or isinstance(att, str):
                    # ensure only dicts or strings for genome objects but shouldn't be anything else anyway
                    io = _StringIO()
                    _json.dump(att, io)
                    io.seek(0, _os.SEEK_END)
                    length = io.tell()
                    io.seek(0)
                    thisone = _tarfile.TarInfo(name = att_name)
                    thisone.size = length
                    tar.addfile(tarinfo = thisone, fileobj = io)
class Genome(_MetaSample):
    '''
    Collect one or more chromosome sequences from the Internet or a local path
    '''

    def __init__(self, sample_name = 'getgenomes', inherit_from = False, 
            **kwargs):
        module_name = __name__
        super(Genome, self).__init__(sample_name, module_name, 
                inherit_from = inherit_from, **kwargs)
        # self.file_name is updated depending on source, usually to an accession
        # cannot save until have obtained data so set to False initially
        file_name = False
        

    def queryEntrezAssembly(self, search_term, user_email, 
            exclude_draft = False, exclude_complete = False):
        '''
        Obtain download links and IDs for NCBI Assembly genome sequences

        One or more search terms can be provided. Assembly accession IDs are the 
        most precise search term: information to retrieve all chromosomes or 
        contigs will be obtained per genome.

        To download a complete genome consisting of a single chromosome, a genbank
        or refseq accession can be used. If additional chromosomes or plasmids
        are found in the genome corresponding to the provided sequence accession
        they will also be downloaded to the same file unless single_sequences is set 
        to True.

        Taxonomic names can also be provided, but beware of ambiguous search terms
        because download details for all assemblies will be obtained.

        Your email address is required so NCBI can contact you in the event of a 
        problem. If you search 'Bacteria', you might get an email from NCBI for 
        using too much bandwidth (try to avoid doing this).

        To do the actual downloading call .downloadNext() for each item.

        Uses BioPython's Bio.Entrez module.
        '''
        from Bio.Entrez.Parser import ValidationError as _ValidationError

        _Entrez.email = user_email
        use_db = "assembly"

        # if search term is string, use directly
        # if search term is integer, treat as taxonomy ID and convert to
        # name if lower than genus or collect genera within ID if search
        # on each genus

        if exclude_draft and exclude_complete:
            raise LookupError('Excluding both draft and complete genomes will '\
                    'return nothing!')

        try:
            handle = _Entrez.esearch(db = use_db, term = search_term, retmax = 100000)
            result = _Entrez.read(handle)
            handle.close()
        except RuntimeError as e:
            error = 'There might be a problem with your search term ({}): '\
                    'BioPython returned a RuntimeError ({})'\
                    ''.format(search_term, e)
            self.logger.error(error)
            raise LookupError(error)

        self.logger.info('{} assembly results for "{}"'.format(len(result['IdList']), 
                search_term))
        if len(result['IdList']) == 0:
            error = 'Your search ID: "{}" returned {} assembly results from '\
                    'ncbi.nlm.nih.gov/assembly. Try getFromEntrezNucleotide()'.format(
                    search_id, len(result['IdList']))
            self.logger.error(error)
            raise LookupError(error)

        self.logger.debug('Collecting summary for {}: {}'.format(len(result['IdList'])
                ,','.join(result['IdList'])))

        # limited by length of URI created
        perfetch = 100
        do_validation = True
        # some ways of handling unexpected content from NCBI
        try:
            handle = _Entrez.esummary(db = use_db, 
                    id = ','.join(result['IdList'][:perfetch]))
            raw = _Entrez.read(handle, validate = do_validation)
            handle.close()
            summaries = raw['DocumentSummarySet']['DocumentSummary']
        except _ValidationError as e:
            self.logger.warning('WARNING: The information about this genome '\
                    'returned by NCBI Entrez failed validation '\
                    '(ValidationError): {}'.format(e))
            self.logger.warning('Trying without validation . . .')
            do_validation = False
            handle = _Entrez.esummary(db = use_db, 
                    id = ','.join(result['IdList'][:perfetch]))
            raw = _Entrez.read(handle, validate = do_validation)
            handle.close()
            summaries = raw['DocumentSummarySet']['DocumentSummary']

        # get the remainder with or without validation
        for i in range(1,int(len(result['IdList'])/perfetch)+1):
            self.logger.log(PROGRESS,'Fetching summaries for items {}-{}'\
                    ''.format(i*perfetch,i*perfetch+perfetch))
            handle = _Entrez.esummary(db = use_db, 
                    id = ','.join(result['IdList'][i*perfetch:i*perfetch+perfetch]))
            result2 = _Entrez.read(handle, validate = do_validation)
            handle.close()
            summaries += result2['DocumentSummarySet']['DocumentSummary']

        ass_acc_problems = {}
        assemblies_info = {}

        for summary in summaries:
            self.logger.info('Found: {} ({}, uid:{})'.format(summary['Organism'],
                    summary['AssemblyStatus'], summary.attributes['uid']))
            
            # "Complete Genome" or "Scaffold"
            # self.logger.debug('Assembly status is "{}"'.format(
                    # summary['AssemblyStatus']))
            # self.logger.debug('Taxonomy ID is: {}'.format(summary['Taxid']))
            # self.logger.debug('Species Tax. ID is: {}'.format(
                    # summary['SpeciesTaxid']))
            if summary['AssemblyStatus'] == "Scaffold" and exclude_draft:
                self.logger.debug('Skipping uid:{} because "{}"'.format(
                        summary.attributes['uid'], summary['AssemblyStatus']))
                continue
            if summary['AssemblyStatus'] == "Complete Genome" and exclude_complete:
                self.logger.debug('Skipping uid:{} because "{}"'.format(
                        summary.attributes['uid'], summary['AssemblyStatus']))
                continue
            
            info = {}
            # collect download links for sequences
            refseq_ftps = _re.findall(
                    '<FtpPath type="RefSeq">([^<]+)</FtpPath>', 
                    summary['Meta'])
            if len(refseq_ftps) > 1:
                raise NotImplementedError('More than one RefSeq FTP URL '\
                        'found for uid:{}, {}'.format(summary.attributes['uid'], 
                        ', '.join(refseq_ftps)))
            elif len(refseq_ftps) == 1 and len(refseq_ftps[0]) > 0:
                refseq_ftp = refseq_ftps[0]
                self.logger.info('Found RefSeq link: {}'.format(refseq_ftp))
            else:
                refseq_ftp = False
                self.logger.info('RefSeq link not found')
            genbank_ftps = _re.findall(
                    '<FtpPath type="GenBank">([^<]+)</FtpPath>', 
                    summary['Meta'])
            if len(genbank_ftps) > 1:
                raise NotImplementedError('More than one GenBank FTP URL '\
                        'found for uid:{}, {}'.format(summary.attributes['uid'], 
                        ', '.join(genbank_ftps)))
            elif len(genbank_ftps) == 1 and len(genbank_ftps[0]) > 0:
                genbank_ftp = genbank_ftps[0]
                self.logger.info('Found Genbank link: {}'.format(genbank_ftp))
            else:
                genbank_ftp = False
                self.logger.info('GenBank link not found')
            
            if not(genbank_ftp or refseq_ftp):
                e = 'Failed to retrieve FTP download links for either RefSeq nor '\
                        'GenBank assemlies from MetaData:\n{}'.format(summary['Meta'])
                self.logger.error(e)
                ass_acc_problems[summary.attributes['uid']] = summary['Organism']
                continue
            elif refseq_ftp:
                use_link = refseq_ftp
            elif genbank_ftp:
                use_link = genbank_ftp
            
            # collect some stats
            def getStat(name, meta):
                patt = '<Stat category="{}" sequence_tag="all">([^<]+)</Stat>'.format(name)
                n = _re.findall(patt, meta)[0]
                return(int(n))
            
            if summary['AssemblyStatus'] == "Complete Genome":
                for name in ["replicon_count","non_chromosome_replicon_count"]:
                    try:
                        info[name] = getStat(name, summary['Meta'])
                    except (IndexError, TypeError):
                        self.logger.warning("Couldn't find {} in MetaData: {}"\
                                "".format(name, summary['Meta']))
            elif summary['AssemblyStatus'] == "Scaffold" or \
                    summary['AssemblyStatus'] == "Contig":
                for name in ["contig_count", "contig_n50", "contig_l50", 
                        "total_length", "ungapped_length"]:
                    try:
                        info[name] = getStat(name, summary['Meta'])
                    except (IndexError, TypeError):
                        self.logger.warning("Couldn't find {} in MetaData: {}"\
                                "".format(name, summary['Meta']))
            
            # collect accessions and check whether RefSeq or GenBank
            ass_acc = summary['AssemblyAccession']
            if ass_acc[:3] == 'GCF':
                info['collection'] = 'RefSeq'
            else:
                self.logger.warning('No RefSeq assembly found for uid:{}, '\
                        'accession:{}. You can double check by searching at '\
                        'http://www.ncbi.nlm.nih.gov/assembly'.format(
                        summary.attributes['uid'], ass_acc))
                if ass_acc[:3] == 'GCA':
                    info['collection'] = 'GenBank'
                else:
                    self.logger.warning('No GenBank assembly found for uid:{}, '\
                            'accession:{}. You can double check by searching at '\
                            'http://www.ncbi.nlm.nih.gov/assembly'.format(
                            summary.attributes['uid'],ass_acc))
                    self.logger.warning('Not collecting anything for {}: does not '\
                            'start with GCF nor GCA'.format(ass_acc))
                    ass_acc_problems[summary.attributes['uid']] = summary['Organism']
                    continue
            
            info['url'] = use_link
            info['AssemblyAccession'] = ass_acc
            info['AssemblyStatus'] = summary['AssemblyStatus']
            info['Taxid'] = summary['Taxid']
            info['SpeciesTaxid'] = summary['SpeciesTaxid']
            info['Organism'] = summary['Organism']
            # summary['Synonym']['Genbank']
            # try:
                # if summary['Synonym']['Similarity'] == 'identical':
                    # self.logger.warning('GenBank and RefSeq assemblies '\
                            # 'non-identical for {}, {}'.format(info['Synonym']['Genbank'],
                            # info['Synonym']['RefSeq']))
            # except KeyError:
                # self.logger.debug('No "Synonym" entry for "{}" in Assembly DB'\
                        # ''.format(search_id))
            
            assemblies_info[summary.attributes['uid']] = info

        if len(ass_acc_problems):
            self.logger.warning('For {} results, details could not be retrieved '\
                    'and are ignored ({} were successful)'.format(
                    len(ass_acc_problems), len(assemblies_info)))

        self.assemblies_info = assemblies_info
        self.assemblies_problems = ass_acc_problems


    def downloadNextAssembly(self, include_complete = True, 
            include_scaffolds = True, include_contigs = True,
            refseq_only = False, force = False, retain_gbk = False,
            use_genome_identifier = False, path = '.'):
        '''
        Download the next genome in a sorted dict generated by queryEntrezAssembly()

        This function acts similar to a generator (an may one day be implemented 
        as such). A genome is downloaded via FTP from the NCBI Assembly database. 
        The GenBank file is decompressed and parsed in memory, feature ranges 
        and names are stored and the sequence is converted from a Python string 
        to a Python unicode array. This provides large improvements in speed and 
        memory efficiency. Sequences, feature ranges (a tuple of ORFs and large 
        feature dicts) and full name of each replicon are saved in a dictionaries 
        as attributes of the genome object (.sequence, .annotations, .names).

        The intended use is to call the `saveLocal()` method after fetching a 
        genome and repeating both calls until all genomes have been retrieved. 
        On each call, the current genome information is removed from the .info 
        dict attribute. When all genomes have been processed and the .info dict 
        attribute is empty, StopIteration is raised.


        Parameters
        ----------
        include_complete : bool
            download 'Complete Genome' assemblies (set false to exclude)
        include_scaffolds : bool
            download 'Scaffold' assemblies (contigs joined with NNNs)
        include_contigs : bool
            download 'Contig' assemblies
        refseq_only : bool
            only download if RefSeq version available, not GenBank
        force : bool
            if a version is already saved to disk, repeat download, else skip
        retain_gbk : bool
            save complete GenBank file to disk as well as optimised 
            internal version
        use_genome_identifier : str, optional
            use this to genome identifier to build filename instead of 
            assembly accession


        Returns
        -------
        genome_identifier : str or None
            If a genome was downloaded, the identifier is returned. This is 
            usually the assembly accession. `file_name` attribute is also set 
            for the .saveLocal method. If no download was made, None is returned
            and `file_name` attribute is also set to False (bool).


        Raises
        ------
        StopIteration
            If no genomes in the `assemblies_info` attribute dict.
        IOError
            If md5checksums.txt file from NCBI Assemblies FTP server cannot 
            be parsed. The NCBI Assembly uid is included in the error.
        ValueError
            If md5 checksum for downloaded genbank file not match. The NCBI 
            Assembly uid is included in the error.
        _URLError
            If a download fails
        '''


        if len(self.assemblies_info) == 0:
            raise(StopIteration)
            # see https://docs.python.org/3/glossary.html#term-iterator
            # below should be .__next__() ?

        # collect next
        uid,info = sorted(self.assemblies_info.items())[0]
        del self.assemblies_info[uid]
        # this could be over-ridden
        if use_genome_identifier:
            genome_identifier = use_genome_identifier
        else:
            genome_identifier = info['AssemblyAccession']

        if refseq_only and info['collection'] == 'GenBank':
            self.logger.info('Not downloading GenBank only assembly for uid:{}, '\
                    'accession:{}, organism:{}'.format(uid, 
                    info['AssemblyAccession'], info['Organism']))
            self.file_name = False
            return(genome_identifier)
        elif any([
                (not include_complete and info['AssemblyStatus'] == 'Complete Genome'),
                (not include_scaffolds and info['AssemblyStatus'] == 'Scaffold'),
                (not include_contigs and info['AssemblyStatus'] == 'Contig')]):
            self.logger.info('Not downloading because AssemblyStatus is "{}" for '\
                    'uid:{}, accession:{}, organism:{}'.format(info['AssemblyStatus'], 
                    uid, info['AssemblyAccession'], info['Organism']))
            self.file_name = False
            return(genome_identifier)
        else:
            use_filename = '{}.{}-{}.baga'.format(__name__, type(self).__name__, 
                    genome_identifier)
            use_path = _os.path.sep.join([path,use_filename])
            if _os.path.exists(use_path) and \
                    _os.path.getsize(use_path) > 0 and not force:
                self.logger.info('File for {} already exists and force = False: '\
                        'not downloading. Use force = True to download again'\
                        ''.format(use_path))
                self.file_name = False
                return(genome_identifier)
            # download checksums
            self.logger.info('Downloading checksum from: {}'.format(
                    info['url'] + '/md5checksums.txt'))
            data = self._DL(info['url'] + '/md5checksums.txt', verbose = False)
            try:
                checksum = [l.decode('utf-8').split('  ./') for l in \
                        data.readlines() if '_genomic.gbff.gz' in \
                        l.decode('utf-8')][0][0]
            except ValueError as e:
                error = 'Problem parsing checksum for {}: {}'\
                        ''.format(uid, data.read().decode('utf-8'))
                self.logger.error(error)
                # return uid with error message as part of `args` attribute
                raise IOError(error, uid)
            
            # download sequences and annotations
            use_link = info['url'] + '/' + info['url'].split('/')[-1] + \
                    '_genomic.gbff.gz'
            self.logger.info('Downloading sequence data from:\n{}'.format(use_link))
            data = self._DL(use_link, verbose = True)
            hasher = _md5()
            buff = data.read(65536)
            while len(buff) > 0:
                hasher.update(buff)
                buff = data.read(65536)
            
            # raise some other error here? ValueError?
            if hasher.hexdigest() != checksum:
                error = 'checksum fail!'
                self.logger.error(error)
                self.file_name = False
                raise ValueError(error, uid)
            
            self.logger.log(PROGRESS, '. . . checksum {} passed!'.format(checksum))
            
            data.seek(0)
            if _PY3:
                gbk_content = _TextIOWrapper(_gzip.GzipFile(mode="rb", fileobj = data), 
                        encoding = 'utf-8')
            else:
                gbk_content = _gzip.GzipFile(mode="rb", fileobj = data)
                        #encoding = 'utf-8')
            self.loadFromGBK(gbk_content)
            if retain_gbk:
                gbk_content.seek(0)
                gbk_filename = use_link.split('/')[-1].replace('.gz','')
                self.logger.info('Also writing {}'.format(gbk_filename))
                try:
                    self.logger.info('Writing genbank file to {}.'\
                            ''.format(gbk_filename))
                    open(gbk_filename, 'w').write(gbk_content.read())
                except OSError:
                    self.logger.warning('Seems like writing {} failed. '\
                            'Continuing with any others for download.'\
                            ''.format(gbk_filename))
            
            self.file_name = use_filename
            self.logger.debug('will use filename: {}'.format(self.file_name))
            self.sample_name = genome_identifier
            self.source = use_link
            self.organism = info['Organism']
            return(genome_identifier)


    def downloadFromList(self, accs_urls_dict, force = False, retain_gbk = False,
            user_email = False, path = '.'):
        '''
        Download genomes from a list of accession numbers optionally with URLs

        These are usually from a previous baga "CollectData --genomes" run and 
        is designed to simplify reproducibity. The same query is likely to return
        a larger selection if repeated in the future.

        Problems with downloads are logged as warnings but attempts on other 
        downloads are continued.

        Parameters
        ----------
        accs_urls_dict : dict
            keys are assumed to be accession numbers, values are corresponding 
            URLs for download. The URLs are optional and values can be zero 
            length strings or None in which case Entrez is queried with the 
            accessions to get the URL for each. Keys can be any search term that
            return a single result and will be used for "sample_name" and in the
            file name.

        Returns
        -------
        accs_downloaded : list
            A list of the accessions successfully downloaded.
        '''


        # if any accessions are without download FTPs, we'll need an email address
        missing = [a for a,u in accs_urls_dict.items() if len(u) == 0 or u is None]
        if len(missing):
            if not user_email:
                error = 'No email address provided for Entrez but some accessions are '\
                        'missing download URLs: {}'.format(', '.join(missing))
                self.logger.error(error)
                raise ValueError(error)
            else:
                _Entrez.email = user_email
                use_db = "assembly"

        downloaded = []
        for accession,url in accs_urls_dict.items():
            use_filename = '{}.{}-{}.baga'.format(__name__, type(self).__name__, 
                    accession)
            use_path = _os.path.sep.join([path,use_filename])
            if _os.path.exists(use_path) and \
                    _os.path.getsize(use_path) > 0 and not force:
                self.logger.info('File for {} already exists and force = False: '\
                        'not downloading. Use force = True to download again'\
                        ''.format(use_path))
                downloaded += [accession]
                continue
            
            self.logger.info('Attempting to download {} from a list.'\
                    ''.format(accession))
            
            if url is None or len(url) == 0:
                # no download url: query Entrez with accession to get it
                self.logger.info('No URL provided, querying NCBI Assembly via '\
                        'Entrez with "{}"'.format(accession))
                try:
                    handle = _Entrez.esearch(db = use_db, term = accession, 
                            retmax = 100)
                    result = _Entrez.read(handle)
                    handle.close()
                except RuntimeError as e:
                    error = 'Skipping {0} - there might be a problem with your '\
                            'search term ({0}): BioPython returned a RuntimeError '\
                            '({1})'.format(accession, e)
                    self.logger.warning(error)
                    continue
                
                if len(result['IdList']) != 1:
                    error = 'Skipping - your search ID: "{}" returned {} '\
                            'assembly results from ncbi.nlm.nih.gov/assembly, '\
                            'but one is required. Using full accessions should '\
                            'return a single result'.format(accession, 
                            len(result['IdList']))
                    self.logger.warning(error)
                    continue
                else:
                    uid = result['IdList'][0]
            
                # some ways of handling unexpected content from NCBI
                try:
                    handle = _Entrez.esummary(db = use_db, id = uid)
                    raw = _Entrez.read(handle, validate = True)
                    handle.close()
                    summary = raw['DocumentSummarySet']['DocumentSummary'][0]
                except _ValidationError as e:
                    self.logger.warning('WARNING: The information about this genome '\
                            '(uid:{}) returned by NCBI Entrez failed validation '\
                            '(ValidationError): {}'.format(uid,e))
                    self.logger.warning('Trying without validation . . .')
                    handle = _Entrez.esummary(db = use_db, id = uid)
                    raw = _Entrez.read(handle, validate = False)
                    handle.close()
                    summary = raw['DocumentSummarySet']['DocumentSummary'][0]
                
                self.logger.info('Found: {} ({}, uid:{})'.format(summary['Organism'],
                        summary['AssemblyStatus'], summary.attributes['uid']))
                
                info = {}
                # collect download links for sequences
                refseq_ftps = _re.findall(
                        '<FtpPath type="RefSeq">([^<]+)</FtpPath>', 
                        summary['Meta'])
                if len(refseq_ftps) > 1:
                    raise NotImplementedError('More than one RefSeq FTP URL '\
                            'found for uid:{}, {}'.format(summary.attributes['uid'], 
                            ', '.join(refseq_ftps)))
                elif len(refseq_ftps) == 1 and len(refseq_ftps[0]) > 0:
                    refseq_ftp = refseq_ftps[0]
                    self.logger.info('Found RefSeq link: {}'.format(refseq_ftp))
                else:
                    refseq_ftp = False
                    self.logger.info('RefSeq link not found')
                genbank_ftps = _re.findall(
                        '<FtpPath type="GenBank">([^<]+)</FtpPath>', 
                        summary['Meta'])
                if len(genbank_ftps) > 1:
                    raise NotImplementedError('More than one GenBank FTP URL '\
                            'found for uid:{}, {}'.format(summary.attributes['uid'], 
                            ', '.join(genbank_ftps)))
                elif len(genbank_ftps) == 1 and len(genbank_ftps[0]) > 0:
                    genbank_ftp = genbank_ftps[0]
                    self.logger.info('Found Genbank link: {}'.format(genbank_ftp))
                else:
                    genbank_ftp = False
                    self.logger.info('GenBank link not found')
                
                if not(genbank_ftp or refseq_ftp):
                    e = 'Skipping - failed to retrieve FTP download links for either '\
                            'RefSeq nor GenBank assemlies from MetaData:\n{}'.format(
                            summary['Meta'])
                    self.logger.warning(e)
                    continue
                elif refseq_ftp:
                    use_link = refseq_ftp
                elif genbank_ftp:
                    use_link = genbank_ftp
                
                url = use_link + '/' + use_link.split('/')[-1] + \
                    '_genomic.gbff.gz'
            
            # have download URL
            # download checksums
            md5_url = '/'.join(url.split('/')[:-1]) + '/md5checksums.txt'
            self.logger.info('Downloading checksum from: {}'.format(md5_url))
            try:
                data = self._DL(md5_url, verbose = False)
            except _URLError:
                self.logger.warning('Skipping {} - failed to download '\
                        'MD5 checksums from: {}'.format(accession, md5_url))
                continue
            try:
                checksum = [l.decode('utf-8').split('  ./') for l in \
                        data.readlines() if '_genomic.gbff.gz' in \
                        l.decode('utf-8')][0][0]
            except ValueError as e:
                error = 'Skipping {} - problem parsing checksum for {}: {}'\
                        ''.format(uid, data.read().decode('utf-8'))
                self.logger.warning(error)
                continue
            
            # download sequences and annotations
            self.logger.info('Downloading sequence data from: {}'.format(url))
            try:
                data = self._DL(url, verbose = True)
            except _URLError:
                self.logger.warning('Skipping {} - failed to download '\
                        'sequence data from: {}'.format(accession, md5_url))
                continue
            hasher = _md5()
            buff = data.read(65536)
            while len(buff) > 0:
                hasher.update(buff)
                buff = data.read(65536)
            
            # raise some other error here? ValueError?
            if hasher.hexdigest() != checksum:
                error = 'checksum fail!'
                self.logger.error(error)
                self.file_name = False
                raise ValueError(error, uid)
            
            self.logger.log(PROGRESS, '. . . checksum {} passed!'.format(checksum))
            data.seek(0)
            gbk_content = _TextIOWrapper(_gzip.GzipFile(mode="rb", fileobj = data), 
                    encoding = 'utf-8')
            self.loadFromGBK(gbk_content)
            if retain_gbk:
                gbk_content.seek(0)
                gbk_filename = url.split('/')[-1].replace('.gz','')
                self.logger.info('Also writing {}'.format(gbk_filename))
                try:
                    self.logger.info('Writing genbank file to {}.'\
                            ''.format(gbk_filename))
                    open(gbk_filename, 'w').write(gbk_content.read())
                    if path != '.':
                        new_path = _os.path.sep([path,gbk_filename])
                        self.logger.log(PROGRESS, 'Moving to new path: {}'.format(new_path))
                        _os.rename(use_filename, new_path)
                except OSError:
                    self.logger.warning('Seems like writing {} failed. '\
                            'Continuing with any others for download.'\
                            ''.format(gbk_filename))
            self.file_name = use_filename
            self.logger.debug('will use filename: {}'.format(self.file_name))
            self.sample_name = accession
            self.source = url
            #self.organism = info['Organism']   ## how to get this from list? gbk should contain under: SOURCE and/or ORGANISM
            self.saveLocal(exclude = ["assemblies_info", "assemblies_problems"])
            if path != '.':
                # avoid having full path in self.file_name so saved
                # data is not tied to one location
                new_path = _os.path.sep.join([path,use_filename])
                self.logger.log(PROGRESS, 'Moving to new path: {}'.format(new_path))
                _os.rename(use_filename, new_path)
            downloaded += [accession]

        return(downloaded)


    def getFromEntrezNucleotide(self, search_term, user_email, 
            use_genome_identifier = False):
        '''
        download a genome sequence given an accession search ID

        Download an annotated genome sequence from NCBI Assembly
        database. Requires RefSeq or GenBank accession number of entry e.g. 
        'FM209186' for Pseudomonas aeruginosa LESB58 and your email address 
        so NCBI can contact you in the event of a problem. Uses BioPython's 
        Bio.Entrez module.
        '''
        _Entrez.email = user_email

        use_db = 'nuccore'
        handle = _Entrez.esearch(db = use_db, term = search_term, retmax = 100000)
        result = _Entrez.read(handle)
        handle.close()
        self.logger.debug('Search term "{}" returned {} IDs from {}'.format(
                search_term, len(result['IdList']), use_db))
        handle = _Entrez.esummary(db = use_db, 
                id = ','.join(result['IdList']))

        do_validation = True
        raw = _Entrez.read(handle, validate = do_validation)
        handle.close()
        refseq_chromosomes = []
        genbank_chromosomes = []
        for summary in raw:
            if 'plasmid' in summary['Title'] or 'genome' in summary['Title']:
                if summary['Caption'][:3] == 'NC_':
                    refseq_chromosomes += [summary]
                elif summary['Caption'][:2] == 'CP':
                    genbank_chromosomes += [summary]

        if len(refseq_chromosomes):
            use_entries = refseq_chromosomes
            self.logger.debug('Search term "{}" returned {} RefSeq replicon sequences from {}'.format(
                search_term, len(refseq_chromosomes), use_db))
            collection = 'RefSeq'
        elif len(genbank_chromosomes):
            use_entries = genbank_chromosomes
            self.logger.debug('Search term "{}" returned {} GenBank replicon sequences from {}'.format(
                search_term, len(refseq_chromosomes), use_db))
            collection = 'GenBank'
        else:
            raise LookupError('Failed to find any RefSeq or GenBank accession '\
                    'numbers for "{}" in nuccore'.format(search_term))

        # first get annotations
        ids = ','.join([summary['Id'] for summary in use_entries])
        handle = _Entrez.efetch(db = use_db, id = ids, rettype = 'gb', 
                retmode = "text")
        data = _StringIO()
        data.write(handle.read())
        handle.close()
        data.seek(0)
        loci = {}
        try:
            for rec in _SeqIO.parse(data,'genbank'):
                accession = rec.id
                loci[rec.id] = self._extractLoci(rec)
        except ValueError as error_message:
            error = "There was a problem with the genome annotaion as genbank "\
                    "(accession: {}) downloaded from NCBI via Entrez: {}. Retry "\
                    "because Entrez can be unreliable, or try loading from a "\
                    ".gbk file downloaded manually from e.g., ftp://ftp.ncbi.nih.gov/"\
                    "genomes/Bacteria/".format(accession, error_message)
            self.logger.error(error)
            raise LookupError(error)

        # then get sequences (sometimes gb is annotations only)
        handle = _Entrez.efetch(db = use_db, id = ids, rettype = 'fasta')
        data = _StringIO()
        data.write(handle.read())
        handle.close()
        data.seek(0)
        seqs = {}
        names = {}
        try:
            for rec in _SeqIO.parse(data,'fasta'):
                accession = rec.id.split('|')[3]
                seqs[accession] = _array('u', str(rec.seq))
                names[accession] = rec.description
        except ValueError as error_message:
            error = "There was a problem with the genome sequence as fasta "\
                    "(accession: {}) downloaded from NCBI via Entrez: {}. Retry "\
                    "because Entrez can be unreliable, or try loading from a "\
                    ".gbk file downloaded manually from e.g., ftp://ftp.ncbi.nih.gov/"\
                    "genomes/Bacteria/".format(accession, error_message)
            self.logger.error(error)
            raise LookupError(error)

        if sorted(loci) != sorted(seqs):
            raise LookupError('Failed to collect matching annotations and sequences '\
                    'from nuccore for "{}"'.format(search_term))

        self.sequence = seqs
        self.annotations = loci
        self.names = names

        ### how many sequences will be saved?
        ### all that match search term?
        ### that allows singles for specific accessions and groups for specific organisms
        ### single result: accession which will probably match search term
        ### multiple results . . sanitised search term?
        ### and optionally --use_my_id <> from CLI ?

        if use_genome_identifier:
            genome_identifier = use_genome_identifier
        elif len(seqs) == 1:
            genome_identifier = sorted(seqs.keys())[0]
        else:
            # could use a pathname sanitiser here?
            genome_identifier = search_term.replace(' ','_')

        self.sample_name = genome_identifier

        self.file_name = '{}.{}-{}.baga'.format(__name__, type(self).__name__, 
                genome_identifier)

        return(genome_identifier)

    def loadFromGBK(self, data):
        '''
        Load an annotated genome sequence from a genbank file-like object.

        The attributes: "sequence", "annotations" and "names" are populated 
        which are dicts with an entry for each replicon or contig. The 
        attributes "file_name", "source" and "sample_name" are not populated
        but are required by most other methods and functions so should be 
        defined.

        GenBank files can be manually downloaded from e.g.,
        ftp://ftp.ncbi.nih.gov/genomes/Bacteria/

        Parameters
        ----------
        data: str or file
            The file-like object can be a local path, a file handle or a 
            file-like object such as generated by the io or StringIO modules
            because BioPython deals with that internally.
        '''
        records = _SeqIO.parse(data, 'genbank')
        seqs = {}
        loci = {}
        names = {}
        c = 0
        #could catch UnicodeDecodeError and maybe others for failed BioPython parsing?
        if _PY3:
            arraytype = 'u'
        else:
            arraytype = 'c'

        for seq_record in records:
            ORF_ranges, large_mobile_element_ranges, ordinate_offset = \
                    self._extractLoci(seq_record)
            if ordinate_offset:
                use_seq = str(seq_record.seq[-ordinate_offset:]) + \
                        str(seq_record.seq[:-ordinate_offset])
                seqs[seq_record.id] = _array(arraytype, use_seq)
            else:
                seqs[seq_record.id] = _array(arraytype, str(seq_record.seq))
            loci[seq_record.id] = ORF_ranges, large_mobile_element_ranges
            names[seq_record.id] = seq_record.description
            self.logger.debug('parsed: {}: {}'.format(
                    seq_record.id, seq_record.description))
            c += 1

        self.logger.log(PROGRESS, '{} records in genbank data {}'.format(c,data))
        self.sequence = seqs
        self.annotations = loci
        self.names = names

    def loadFrombaga(self, local_path):
        '''
        Reload a genome that was previously processed with a CollectData.Genome 
        object.
        '''
        with _tarfile.open(local_path, "r:gz") as tar:
            for member in tar:
                contents = _StringIO(tar.extractfile(member).read())
                try:
                    # either json serialised conventional objects
                    contents = _json.loads(contents.getvalue())
                except ValueError:
                    # or longer python array.array objects
                    contents = _array('c', contents.getvalue())
                
                setattr(self, member.name, contents)

    def _extractLoci(self, seq_record):
        '''
        Extract some ORF information from genbank record and put into convenient 
        dictionaries. Checks for features spanning start/end of circular 
        sequences and applies an offset so the start < end and start == 0.
        
        Returns
        -------
        ORF_ranges : dict
        large_mobile_element_ranges : dict
        ordinate_offset : int
            the offset integer to deal with compound features spanning 
            start/end (0 if not applied)
        '''
        ORF_ranges = {}
        large_mobile_element_ranges = {}
        ordinate_offset = 0
        GI_prophage = _re.compile('[Ii]sland|[Pp]hage|GI')
        for f in seq_record.features:
            if f.type == 'CDS':
                # until AA and detailed annotations of ORFs needed, don't save
                # other qualifiers
                try:
                    thisgene = f.qualifiers['gene'][0]
                except KeyError:
                    thisgene = ''
                
                if len(f.location.parts) == 2:
                    if f.location.parts[0].start.position == 0 and \
                            f.location.parts[1].end.position == len(seq_record.seq):
                        # this feature spans 0 on chromosome sequence in file
                        # will use an offset to simplify feature access
                        # by removing overlap
                        s = f.location.parts[1].start.position
                        e = f.location.parts[1].end.position
                        ordinate_offset = e - s
                        s = 0
                        e = f.location.parts[0].end.position + ordinate_offset
                        self.logger.info('Compound feature ({}) across start '\
                                'of sequence {} detected. All ordinates will '\
                                'include an offset'.format(
                                f.qualifiers['locus_tag'][0], seq_record.id))
                    else:
                        self.logger.warning('{} is a compound feature and not '\
                                'supported: ignoring'.format(
                                f.qualifiers['locus_tag'][0]))
                        continue
                else:
                    s = f.location.start.position + ordinate_offset
                    e = f.location.end.position + ordinate_offset
                ORF_ranges[f.qualifiers['locus_tag'][0]] = (s, e, 
                        f.location.strand, thisgene)
                
            if f.type == 'misc_feature':
                try:
                    # not all misc_features have a "note"
                    feature_note = f.qualifiers['note'][0]
                except KeyError:
                    continue
                
                if _re.search(GI_prophage, feature_note) and \
                        f.location.nofuzzy_end - f.location.nofuzzy_start > 10000:
                    large_mobile_element_ranges[f.qualifiers['note'][0]] = (
                            f.location.nofuzzy_start + ordinate_offset, 
                            f.location.nofuzzy_end + ordinate_offset)
        
        # filter ORFs within ORFs (artifacts? PLES_21351 and PLES_21361 in LESB58)
        ORF_ranges_sorted = sorted(ORF_ranges.items(), key = lambda x: x[1][0])
        inner_ORFs = set()
        for n, (ORF1, (s1, e1, strnd1, genename1)) in enumerate(ORF_ranges_sorted[:-1]):
            for ORF2, (s2, e2, strnd2, genename2) in ORF_ranges_sorted[(n+1):]:
                if s1 < s2 and e2 < e1:
                    self.logger.log(PROGRESS, '{} ({}, {}-{}) is within {} '\
                            '({}, {}-{}); dumping former'.format(ORF2, genename2, 
                            s2, e2, ORF1, genename1, s1, e1))
                    inner_ORFs.add(ORF2)
                if s2 < s1 and e1 < e2:
                    self.logger.log(PROGRESS, '{} ({}, {}-{}) is within {} '\
                            '({}, {}-{}); dumping former'.format(ORF1, genename1, 
                            s1, e1, ORF2, genename2, s2, e2))
                    inner_ORFs.add(ORF1)
        
        for ORF in inner_ORFs:
            del ORF_ranges[ORF]
        
        self.logger.log(PROGRESS, '{} ORFs, {} large features parsed in {}'\
                ''.format(len(ORF_ranges), len(large_mobile_element_ranges), 
                seq_record.id))
        return(ORF_ranges, large_mobile_element_ranges, ordinate_offset)

    def _DL(self, url, verbose = True):
        CHUNK = 16 * 1024 * 32
        if _PY3:
            req = _request.urlopen(url)
            data = _BytesIO()
        else:
            req = _urllib2.urlopen(url)
            data = _StringIO()
        c = 0
        for chunk in iter(lambda: req.read(CHUNK), ''):
            if c == 0:
                self.logger.info('Download started: {}'.format(url))
            c += CHUNK
            if verbose:
                self.logger.log(PROGRESS,"{:,} bytes".format(c))
            data.write(chunk)
            if len(chunk) == 0:
                break
        
        if verbose:
            self.logger.info('Download complete . . .')
        data.seek(0)
        return(data)
class Reads:
    '''
    Download reads from your local read archive
    '''
    def __init__(self):
        pass

    def getFromENA(self, run_acc_list, 
                         ftp_server_url = 'ftp.sra.ebi.ac.uk', 
                         local_reads_path = ['reads']):
        '''
        Given a list of 'run' accession numbers for paired end short read analyses, 
        download the read files from the European Nucleotide Archive.

        If using a mirror server, supply an alternative for 'ftp_server_url'.

        'local_reads_path' can be a path string or list or folder names.
        '''
        if isinstance(local_reads_path, list):
            local_reads_path = _os.path.sep.join(local_reads_path)

        if not _os.path.exists(local_reads_path):
            _os.makedirs(local_reads_path)

        print('Logging in to %s' % ftp_server_url)
        ftp = _FTP(ftp_server_url)
        # anonymous login
        print(ftp.login())

        def check_connection(ftp):
            try:
                print('FTP: %s' % ftp.voidcmd("NOOP"))
                # http://docs.python.org/2/library/ftplib.html
                return(True)
            except IOError as e:
                print('Seems to be a problem with the connection to FTP server:')
                print('I/O error({0}): {1}'.format(e.errno, e.strerror) )
                return(False)

        def calc_checksum(filepath):
            hasher = _md5()
            handle = open(filepath, 'rb')
            buff = handle.read(65536)
            while len(buff) > 0:
                hasher.update(buff)
                buff = handle.read(65536)
            
            return(hasher.hexdigest())

        downloaded_read_files = {}

        start_time = _time.time()
        failed = []
        for cnum,run_acc in enumerate(run_acc_list):
            
            query_url_base = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query='
            success = False
            tries = 0
            max_tries = 5
            while not success:
                rest_req = '"run_accession=%s"&result=read_run&fields=fastq_ftp,fastq_md5&display=report' % run_acc
                print('Sending query to ENA:\n%s' % rest_req)
                result = _urllib2.urlopen(query_url_base + rest_req).read()
                print('ENA accession numbers query result:\n%s' % result)
                if result.count('ERR') == 7:
                    success = True
                else:
                    print('Query result from ENA was unexpected on attempt %s of %s' % (tries, max_tries))
                    _time.sleep(0.5)
                    tries += 1
                    if tries == max_tries:
                        print('Attempt %s failed. Try again later and if problem persists, report bug.' % tries)
                        failed += [run_acc]
                        break
                        #_sys.exit(1)
            
            if not success:
                continue
            
            md5s = result.split('\n')[-2].split('\t')[-1][:-1].split(';')
            
            ENA_paths = result.split('\n')[-2].split('\t')[-2][:-1].split(';')
            
            ENA_reads_pair_paths = {}
            ENA_reads_pair_paths[1] = ENA_paths[0].replace(ftp_server_url, '')
            ENA_reads_pair_paths[2] = ENA_paths[1].replace(ftp_server_url, '')
            
            local_reads_pair_paths = {}
            local_reads_pair_paths[1] = local_reads_path + \
                                        _os.path.sep + \
                                        ENA_reads_pair_paths[1].split('/')[-1]
            local_reads_pair_paths[2] = local_reads_path + \
                                        _os.path.sep + \
                                        ENA_reads_pair_paths[2].split('/')[-1]
            
            downloaded_read_files[run_acc] = {}
            
            for f in (1,2):
                # ensure connection is still open
                while not check_connection(ftp):
                    _sleep(0.5)
                    print('Attempting to re-establish connection . . .')
                    ftp = _FTP(ftp_server_url)
                    # anonymous login
                    print(ftp.login())
                    pass
                
                expected_checksum = md5s[f - 1]
                
                exists = _os.path.exists(local_reads_pair_paths[f])
                if exists:
                    print('File %s for %s exists locally: %s' % (f, run_acc, local_reads_pair_paths[f]))
                    actual_checksum = calc_checksum(local_reads_pair_paths[f])
                    if actual_checksum == expected_checksum:
                        print('File checksum matches: %s. Skipping download' % (expected_checksum))
                        downloaded_read_files[run_acc][f] = local_reads_pair_paths[f]
                        continue
                    else:
                        print('Checksum mismatch')
                
                print('Downloading via %s: %s' % (ftp_server_url, ENA_reads_pair_paths[f]))
                res = ftp.retrbinary('RETR %s' % ENA_reads_pair_paths[f], 
                                     open(local_reads_pair_paths[f], 'wb').write)
                print('FTP: %s' % res)
                
                print('Calculating checksum . . .')
                actual_checksum = calc_checksum(local_reads_pair_paths[f])
                
                if actual_checksum == expected_checksum:
                    print('File checksum matches: %s.' % (expected_checksum))
                    downloaded_read_files[run_acc][f] = local_reads_pair_paths[f]
                else:
                    print('Checksum mismatch for: %s')
            
            if len(run_acc_list) > 1:
                # report durations, time left etc
                _report_time(start_time, cnum, len(run_acc_list))

        if len(failed) > 0:
            print('WARNING: some accession numbers did not return a result from ENA')
            print('Try searching http://www.ebi.ac.uk/ena in a web-browser for:')
            print(', '.join(failed))

        self.read_files = downloaded_read_files

    def getFromPath(self, path_to_fastq):
        '''
        Given a path to pairs of fastq short read files, parse them ready for analysis 
        with the Bacteria and Archaea Genome (BAG) Analyser.
        '''    

        use_files = []
        if isinstance(path_to_fastq, str):
            use_paths = [path_to_fastq]
        else:
            use_paths = path_to_fastq

        for path in use_paths:
            if _os.path.isdir(path):
                print('Checking in {}'.format(path))
                # supplied with path to folder - need to check contents
                path1 = _os.path.sep.join([path, '*.fastq'])
                file_list = _glob(path1)
                path2 = _os.path.sep.join([path, '*.fq'])
                file_list += _glob(path2)
                file_list.sort()
                
                path3 = _os.path.sep.join([path, '*.fastq.gz'])
                file_list_gz = _glob(path3)
                path4 = _os.path.sep.join([path, '*.fq.gz'])
                file_list_gz += _glob(path4)
                file_list_gz.sort()
                
                if len(file_list) == 0 and len(file_list_gz) == 0:
                    print('WARNING: did not find any files at {}, {}, {}, nor {}'.format(path1, path2, path3, path4))
                    
                elif len(file_list) == 0 and len(file_list_gz) > 0:
                    print('Found {} total gzipped fastq files'.format(len(file_list_gz)))
                    use_files += file_list_gz
                    
                elif len(file_list) > 0 and len(file_list_gz) == 0:
                    print('Found {} total uncompressed fastq files'.format(len(file_list)))
                    use_files += file_list
                    
                else:
                    print('Found compressed and uncompressed fastq files.\nUsing {} gzipped files'.format(len(file_list_gz)))
                    # could select from a combination without doubling up . . .
                    # preference for uncompressed:
                    # use_files = sorted(list(set(file_list_gz) - set([f+'.gz' for f in file_list])) + file_list)
                    use_files += file_list_gz
            else:
                try:
                    test = open(path, 'r')
                    test.close()
                    # part of a list of reads or shell expansion
                    use_files += [path]
                except IOError:
                    print('WARNING: did not find any files at {}'.format(path))

        use_files.sort()


        # check potential non-pair files
        keep_use_files = []
        for f in use_files:
            if 'singletons' in f:
                print('Assuming {} is not part of a pair: ignoring'.format(f))
                continue
            
            keep_use_files += [f]

        use_files = keep_use_files

        # check filenames for inclusion of known baga downstream files
        keep_use_files = []
        for f in use_files:
            this_suffix = ''
            for suffix in ('_subsmp','_adpt','_qual')[::-1]:
                this_suffix = suffix + this_suffix
                for f2 in use_files:
                    if f2 != f:
                        if this_suffix in f2 and f2.replace(this_suffix,'') == f:
                            error = 'ERROR: {} appears to be a file from a previous baga run that included {}. Try being more specific with the supplied path expansion to read VCFs (i.e., without baga suffixes allowed, e.g. "reads/*_[12].*"), or remove files generated in previous analyses'.format(f2, f)
                            _sys.exit(error)
            
            keep_use_files += [f]

        use_files = keep_use_files

        if len(use_files) == 0:
            print('Error: could not find any files at {}'.format(', '.join(path_to_fastq)))
            print('Please check paths and try again . . .')
            _sys.exit(1)

        if len(use_files) % 2 != 0:
            print('Please supply an even number of paired files. Found {}:\n{}'.format(len(use_files), '\n'.join(use_files)))
            _sys.exit(1)

        error_explanation = 'Problem parsing read files: ensure pairs are numbered '\
        '1 and 2\n'\
        'BAGA looks for a "1" or "2" labelling in read pair filenames and takes '\
        'the last digit in the filename (excluding the set number if present e.g., '\
        '_001.fastq).\n E.g. *R1.fastq.gz and *R2.fastq.gz would be OK, 1_thesereads1'\
        '.fastq.gz and 2_thesereads1.fastq.gz would not. (Leading digits OK for sample '\
        'numbering: 1_* 2_* 3_* etc but must each have 1 or 2 elsewhere in file '\
        'name)\n . . else please report as bug'

        # Illumina filename scheme:
        # <sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz
        # http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm

        # match pairs
        filepairs = {}
        for path in use_files:
            path_bits = path.split(_os.path.sep)
            filename_ext = path_bits[-1]
            # for now dump set number (if present? not always present?)
            # really need to deal with multiple sets and all likely versions of CASAVA filename schemes
            # _<set number (0-padded to 3 digits>.f
            use_filename_ext = _re.sub('(_[0-9]{3})(\.[fF])', r'\2', filename_ext)
            filename, ext = _re.findall('(.+)(\.fastq\.gz|\.fastq|\.fq\.gz|\.fq)$', use_filename_ext)[0]
            ones_and_twos = list(_re.finditer('[12]', filename))
            assert len(ones_and_twos) > 0, '{}. Problem filename: {}'.format(
                                                                        error_explanation, 
                                                                        filename)
            # make name for each pair that is consistant parts of file name
            # joining with space caused problems when incorporating into a filename downstream
            # and joining with underscore risks introducing double underscore which would cause splitting on __ later to fail
            s,e = ones_and_twos[-1].span()
            pairmember = ones_and_twos[-1].group()
            # omit the 1 or 2 relevent to pairing from the name
            part1,part2 = filename[:s],filename[e:]
            if len(part1) and len(part2):
                pairname = '-'.join([part1,part2])
            elif len(part1) and not len(part2):
                pairname = part1
            else:
                pairname = part2
            for known_suffix in ['.fastq.gz','.fq.gz','.fastq','.fq']:
                thismatch = _re.findall('('+known_suffix+')$', pairname)
                if thismatch:
                    pairnamenew = _re.sub('('+thismatch[0]+')$', '', pairname)
                    #print('Removed {} from {} == {}'.format(thismatch, pairname, pairnamenew))
                    pairname = pairnamenew.rstrip(' ')
                    continue
            
            # store with keys 1 or 2
            try:
                filepairs[pairname][int(pairmember)] = path
            except KeyError:
                filepairs[pairname] = {int(pairmember): path}

        # check pairs are accessible
        checked_read_files = {}
        for pairname,files in filepairs.items():
            assert len(files) == 2, '{}. Problem filename(s): {}'.format(
                    error_explanation, ', '.join(files.values()))
            
            print('Collected pair "{}": {} and {}'.format(
                    pairname, files[1], files[2]))
            
            try:
                if _os.path.getsize(files[1]) == 0:
                    print('File access fail (empty file): {}'.format(files[1]))
                    _sys.exit(1)
            except OSError:
                print('File access fail: {}'.format(files[1]))
                _sys.exit(1)
            
            try:
                if _os.path.getsize(files[2]) == 0:
                    print('File access fail (empty file): {}'.format(files[2]))
                    _sys.exit(1)
            except OSError:
                print('File access fail: {}'.format(files[2]))
                _sys.exit(1)
            
            checked_read_files[pairname] = files

        print('Total read pairs: {}'.format(len(checked_read_files)))

        self.read_files = checked_read_files

    def saveLocal(self, name):
        '''
        Save a downloaded read info to a local compressed pickle file.
        'name' can exclude extension: .baga will be added
        '''
        fileout = 'baga.CollectData.Reads-%s.baga' % name
        print('Saving to %s' % fileout)
        _cPickle.dump(self, _gzip.open(fileout, 'wb'))

if __name__ == '__main__':
    main()
