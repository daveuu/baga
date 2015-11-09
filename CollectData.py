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
CollectData module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains functions to collect reference genomes and query read sets 
from a the Internet via a URL or from a local path.
'''

### Import Dependencies for this module ###

# use underscore prefix to keep these private

# stdlib
from ftplib import FTP as _FTP
from time import sleep as _sleep
from hashlib import md5 as _md5
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

from baga import report_time as _report_time

# external Python modules
from Bio import Entrez as _Entrez
from Bio import SeqIO as _SeqIO
def main():
    pass
class Genome:
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
            e = 'Your search ID: "{}" returned {} assembly results '\
            'from ncbi.nlm.nih.gov/assembly but a single result is required.'.format(
                    search_id, len(result['IdList']))
            assert len(result['IdList']) == 1, e
            Assembly_ID = result['IdList'][0]
            handle = _Entrez.esummary(db = "assembly", id = Assembly_ID)
            info = _Entrez.read(handle)['DocumentSummarySet']['DocumentSummary'][0]
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
            refseq_ass_acc = info[u'AssemblyAccession']
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
                            search_id_unversioned, refseq2version[search_id_unversioned]))
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
            seq_record = _SeqIO.read(local_path, "genbank")
            self.ORF_ranges, self.large_mobile_element_ranges = extractLoci(seq_record)
            self.sequence = _array('c', seq_record.seq)
            self.id = seq_record.id
        
        
        
        
        
        
        if accession and user_email:
            getFromEntrez(accession, user_email)
            
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


        # match pairs
        filepairs = {}
        for path in use_files:
            path_bits = path.split(_os.path.sep)
            filename = path_bits[-1]
            # {1} prevents splitting at 1 or 2 first character
            bits = _re.split('([^0-9]{1}[12][^0-9])', filename)
            e = 'Problem parsing read files: ensure pairs are numbered 1 and 2\n'
            e += 'BAGA looks for 1 and 2 labelling in read pair filenames, excluding at the start and 1 or 2 among other digits.\n'
            e += 'E.g. *R1.fastq.gz and *R2.fastq.gz would be OK, 1_thesereads.fastq.gz and 2_thesereads.fastq.gz would not.\n'
            e += '(leading digits OK for sample numbering: 1_* 2_* 3_* etc but must each have 1 or 2 elsewhere in file name)\n'
            e += ' . . else please report as bug. Problem filename: {}'.format(filename)
            assert len(bits) == 3, e
            known_suffixes = ['.fastq.gz','.fq.gz','.fastq','.fq']
            # make name for each pair that is consistant parts of file name
            # joining with space caused problems when incorporating into a filename downstream
            # and joining with underscore risks introducing double underscore which would cause splitting on __ later to fail
            pairname = '-'.join([bits[0],bits[2]])
            for known_suffix in known_suffixes:
                thismatch = _re.findall('('+known_suffix+')$', pairname)
                if thismatch:
                    pairnamenew = _re.sub('('+thismatch[0]+')$', '', pairname)
                    #print('Removed {} from {} == {}'.format(thismatch, pairname, pairnamenew))
                    pairname = pairnamenew.rstrip(' ')
                    continue
            
            # store with keys 1 or 2
            try:
                filepairs[pairname][int(bits[1][1])] = path
            except KeyError:
                filepairs[pairname] = {int(bits[1][1]): path}

        # check pairs are accessible
        checked_read_files = {}
        for pairname,files in filepairs.items():
            print('Collected pair: {} and {}'.format(files[1], files[2]))
            
            try:
                if _os.path.getsize(files[1]) == 0:
                    print('File access fail: {}'.format(files[1]))
                    _sys.exit(1)
            except OSError:
                print('File access fail: {}'.format(files[1]))
                _sys.exit(1)
            
            try:
                if _os.path.getsize(files[2]) == 0:
                    print('File access fail: {}'.format(files[2]))
                    _sys.exit(1)
            except OSError:
                print('File access fail: {}'.format(files[2]))
                _sys.exit(1)
            
            checked_read_files[pairname] = files


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
