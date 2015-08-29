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
CollectData module from the Bacterial and Archaeal Genome (BAG) Analyser.

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
        
        def getFromEntrez(accession, user_email):
            _Entrez.email = user_email
            handle = _Entrez.efetch(db = "nuccore", rettype = "gb", retmode = "text", id = accession)
            seq_record = _SeqIO.read(handle, "genbank")
            handle.close()
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

        if len(path_to_fastq) == 1:
            # supplied with path to folder - need to check contents
            file_list = _glob(_os.path.sep.join([path_to_fastq[0], '*.fastq']))
            file_list += _glob(_os.path.sep.join([path_to_fastq[0], '*.fq']))
            file_list.sort()
            
            file_list_gz = _glob(_os.path.sep.join([path_to_fastq[0], '*.fastq.gz']))
            file_list_gz += _glob(_os.path.sep.join([path_to_fastq[0], '*.fq.gz']))
            file_list_gz.sort()
            
            if len(file_list) == 0 and len(file_list_gz) == 0:
                print('Error: did not find any files at {} nor {}'.format(file_list, file_list_gz))
                print('Please check paths and try again . . .')
                _sys.exit(1)
                
            elif len(file_list) == 0 and len(file_list_gz) > 0:
                print('Found {} total gzipped fastq files'.format(len(file_list_gz)))
                use_files = file_list_gz
                
            elif len(file_list) > 0 and len(file_list_gz) == 0:
                print('Found {} total uncompressed fastq files'.format(len(file_list)))
                use_files = file_list
                
            else:
                print('Found compressed and uncompressed fastq files.\n\
            Using {} gzipped files'.format(len(file_list_gz)))
                # could select from a combination without doubling up . . .
                # preference for uncompressed:
                # use_files = sorted(list(set(file_list_gz) - set([f+'.gz' for f in file_list])) + file_list)
                use_files = file_list_gz
        else:
            # supplied with list of reads or shell expansion
            use_files = sorted(path_to_fastq)

        if len(use_files) % 2 != 0:
            print('Please supply an even number of paired files')
            _sys.exit(1)

        # match pairs
        checked_read_files = {}
        for n,f in enumerate(use_files[::2]):
            #print(f,use_files[n*2 + 1])
            p1, p2 = f,use_files[n*2 + 1]
            try:
                if _os.path.getsize(p1) == 0:
                    print('File access fail: {}'.format(p1))
                    sys.exit(1)
            except OSError:
                print('File access fail: {}'.format(p1))
                sys.exit(1)
            
            try:
                if _os.path.getsize(p2) == 0:
                    print('File access fail: {}'.format(p2))
                    _sys.exit(1)
            except OSError:
                print('File access fail: {}'.format(p2))
                _sys.exit(1)
            
            p1_f = p1.split(_os.path.sep)[-1]
            p2_f = p2.split(_os.path.sep)[-1]
            for a in range(len(p1)):
                if p1_f[a] != p2_f[a]:
                    break
            
            pairname = p1_f[:a]
            print('Collected pair: {} and {}'.format(p1_f, p2_f))
            checked_read_files[pairname] = {1: p1, 2: p2}

        # no check on fails currently
        # if len(failed) > 0:
            # print('These failed:')
            # print(', '.join(failed))

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
