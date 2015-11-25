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
Dependencies module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains functions to check whether required progams are available.
It also contains functions to get some of the dependencies.
'''


# stdlib
from baga import _subprocess
from baga import _os
import sys as _sys
import urllib2 as _urllib2
import cStringIO as _cStringIO
import zipfile as _zipfile
import tarfile as _tarfile
import hashlib as _hashlib
import shutil as _shutil

# put local modules at start to over-ride local versions
_sys.path.insert(0,_os.path.sep.join([_os.path.abspath(_os.path.curdir),'local_packages']))
# need to fix these hacks:
_sys.path.insert(0,_os.path.sep.join([_os.path.abspath(_os.path.curdir),'local_packages/biopython-1.65']))
_sys.path.insert(0,_os.path.sep.join([_os.path.abspath(_os.path.curdir),'local_packages/svgwrite-1.1.6']))
_sys.path.insert(0,_os.path.sep.join([_os.path.abspath(_os.path.curdir),'local_packages/pyparsing-2.0.3']))
_sys.path.insert(0,_os.path.sep.join([_os.path.abspath(_os.path.curdir),'local_packages/DendroPy-4.0.3']))



def get_git(name, description, source, url, commit, checksum, destination, preparation, checker):
    '''
    Get a dependency from git
    '''
    
    if _os.path.realpath(_os.path.curdir) != destination:
        try:
            _os.chdir(destination)
        except OSError:
            _os.makedirs(destination)
            _os.chdir(destination)
    
    try:
        # clear any previous verions
        _shutil.rmtree(url.split('/')[-1].replace('.git',''))
    except OSError:
        pass
    
    git_server = url.replace('https://','').replace('http://','').split('/')[0]
    print('Downloading {} via git from {} . . .'.format(name, git_server))
    _subprocess.call(['git', 'clone', url])
    _os.chdir(url.split('/')[-1].replace('.git',''))
    _subprocess.call(['git', 'checkout', commit])
    # if repo uses git submodules, those will be set to the correct revisions for this commit
    # else will do nothing
    _subprocess.call(['git', 'submodule', 'update', '--init'])
    
    working_dir = _os.path.sep.join([destination,url.split('/')[-1].replace('.git','')])
    
    if preparation is not None:
        for do_this in preparation:
            if isinstance(do_this['arguments'], dict):
                do_this['function'](**do_this['arguments'])
            else:
                do_this['function'](*do_this['arguments'])
            
            # restore position in path if a prepare changed it
            if working_dir != _os.path.realpath(_os.path.curdir):
                _os.chdir(working_dir)
    
    _os.chdir(_os.path.pardir)
    _os.chdir(_os.path.pardir)

def get_download(name, description, source, url, commit, checksum, destination, preparation, checker):
    '''
    Download and unpack a dependancy
    '''
    ## 
    try:
        _os.chdir(destination)
    except OSError:
        _os.makedirs(destination)
        _os.chdir(destination)
    
    if checksum:
        hasher_algorithm = checksum.split('=')[0]
        if hasher_algorithm == 'md5':
            hasher = _hashlib.md5()
        elif hasher_algorithm == 'sha1':
            hasher = _hashlib.sha1()
        elif hasher_algorithm == 'sha224':
            hasher = _hashlib.sha224()
        elif hasher_algorithm == 'sha256':
            hasher = _hashlib.sha256()
        elif hasher_algorithm == 'sha384':
            hasher = _hashlib.sha384()
        elif hasher_algorithm == 'sha512':
            hasher = _hashlib.sha512()
        else:
            print("{} checksums not implemented in Python's hashlib!".format(hasher_algorithm))
    
    print('Downloading: %s' % url)
    req = _urllib2.urlopen(url)
    CHUNK = 16 * 1024 * 16
    data = _cStringIO.StringIO()
    c = 0
    for chunk in iter(lambda: req.read(CHUNK), ''):
        c += CHUNK
        print("{:,} bytes".format(c))
        data.write(chunk)
    
    print('Download complete . . .')
    data.seek(0)
    
    if checksum:
        buff = data.read(65536)
        while len(buff) > 0:
            hasher.update(buff)
            buff = data.read(65536)
        
        e = '. . . checksum fail!'
        assert hasher.hexdigest() == checksum.split('=')[1], e
        print('. . . checksum passed!')
        data.seek(0)
    
    if url[-6:] == 'tar.gz':
        archive = _tarfile.open(mode="r:gz", fileobj = data)
    elif url[-7:] == 'tar.bz2':
        archive = _tarfile.open(mode="r:bz2", fileobj = data)
    elif url[-4:] == '.zip':
        archive = _zipfile.ZipFile(data)
    
    if destination == 'local_packages':
        # extract as a pypi python package
        release = url.split('/')[-1][:-7]
        print('Extracting {} to {}'.format(release, _os.path.sep.join([destination,name])))
        c = 0
        nostrip = {'pysam'}
        if name in nostrip:
            try:
                _shutil.rmtree(archive.getnames()[0])
            except OSError:
                pass
            #_shutil.rmtree(_os.path.sep.join([destination, archive.getnames()[0]]))
            # some python modules should not be stripped . . more complex install
            for member in archive.getmembers():
                if member.isreg():
                    archive.extract(member)
                    print(member.name)
                    c += 1
        else:
            # others don't need additional compilation
            check_path1 = '{}/{}'.format(release,name)
            for member in archive.getmembers():
                if member.isreg() and check_path1 in member.name:
                    member.name = _os.path.sep.join(member.name.split(_os.path.sep)[1:])
                    archive.extract(member)
                    c += 1
        print('Extracted {} files'.format(c))
    
    else:
        # extract as a generic external program
        archive.extractall()
    
    if preparation:
        for do_this in preparation:
            if 'just_packages' in do_this['arguments']:
                # this is the only thing that differentiates this prepare()
                # from others that need some chdir <== this should be improved
                # see dep dict
                _os.chdir(_os.path.pardir)
                do_this['function'](*do_this['arguments']['package_list'])
            else:
                extracted_base_dir = archive.getnames()[0].split(_os.path.sep)[0]
                _os.chdir(extracted_base_dir)
                do_this['function'](**do_this['arguments'])
                _os.chdir(_os.path.pardir)
                _os.chdir(_os.path.pardir)
    else:
        _os.chdir(_os.path.pardir)

def get_other_packages(packages):
    '''Python package dependencies for python packages'''
    from baga import Dependencies
    for package in packages:
        get_download(**Dependencies.dependencies[package])

def get_GATK():
    print("Visit https://www.broadinstitute.org/gatk/download/auth?package=GATK and download after reviewing and agreeing license (and signing up to the forum). Place jar file to external_programs/GenomeAnalysisTK/GenomeAnalysisTK.jar")

def pysam_install():
    _subprocess.call([_sys.executable, 'setup.py', 'build'])
    for f in _os.listdir('build/lib.linux-x86_64-2.7/pysam'):
        if f[-3:] == '.so':
            #print(f)
            _shutil.copy('build/lib.linux-x86_64-2.7/pysam/{}'.format(f), 'pysam/{}'.format(f))
    
    try:
        _shutil.rmtree('../pysam')
    except OSError:
        pass
    
    _os.rename('pysam', '../pysam')


def discosnp_install():
    print(_os.listdir('.'))
    _subprocess.call(['bash','./compile_discoSnp++.sh'])
    # ensure VCF generation works if system default python is version 3
    # VCF code is python 2 only
    fixed = open('run_VCF_creator.sh').read().replace('python ','python2 ')
    open('run_VCF_creator.sh','w').write(fixed)


def prep_python_install(extras):
    print('Installing via setup.py . . .')
    # could try here?
    _subprocess.call([_sys.executable, 'setup.py', 'install'] + extras)


def prep_simple_make(path = False, configure = False, alt_command = False):
    if path:
        _os.chdir(_os.path.sep.join(path))
    
    if configure:
        _subprocess.call(['./configure'])
    
    if alt_command:
        _subprocess.call([alt_command])
    else:
        _subprocess.call(['make'])
    
    if path:
        _os.chdir(_os.path.sep.join([_os.path.pardir]*len(path)))

def check_no_error(path = '', 
                   system = False, 
                   java_commands = False, 
                   extra = False, 
                   extra_pre = False,
                   success_returncode = 1):
    '''
    check a binary executable can be called.
    
    A failure is considered as OSError, not a non-zero exit status.
    Can check in default baga location or a specified path or the system path 
    or java.
    Additional commands and be appended with 'extra' list or pre-pended with 
    'extra-pre' list e.g., specifying 'python2' for when python3 defaults but 
    python >3.3 is installed but not supported.
    '''
    if java_commands:
        p = _subprocess.Popen(['java'] + java_commands, stdout=_subprocess.PIPE, stderr=_subprocess.STDOUT)
        output = p.stdout.read()
        if 'Error:' in output:
            print(output)
            return(False)
        else:
            print(output)
            return(True)
    else:
        if system:
            cmd = path[-1:]
        elif path[0][0] == _os.path.sep or path[0] == '':
            # absolute path to elsewhere provided
            cmd = [_os.path.sep.join(path)]
        else:
            cmd = [_os.path.sep.join(['external_programs'] + path)]
        
        if extra:
            cmd += extra
        
        if extra_pre:
            cmd = extra_pre + cmd
        
        try:
            o = _subprocess.check_output(cmd)
            o = '\n'.join(['external> '+l for l in o.split('\n')])+'\n'
            print(o) # eventually put in the debugging log <== too much noise
            # some programs output a help screen if commands not given
            # . . . they are installed and available
            return(True)
            
        except _subprocess.CalledProcessError as e:
            if e.returncode == success_returncode:
                # some programs that exit return non-zero code if commands not given
                # usually this is 1, but occasionally it is 2 (e.g. cutadapt)
                # but a fail on 2 is e.g. python calling a non-existent python program
                # should be set in the dependencies dictionary if not 1 (default)
                # . . . but they are installed and available
                return(True)
                
            else:
                return(False)
            
        except OSError as e:
            print('{}: {}'.format(cmd[0], e))
            return(False)



def check_python_package(package, maj_version = False, system = False, package_path = 'local_packages'):
    '''Check if a package can be imported from the locally installed collection
    
    Print the results to stderr so they can be collected separately if calling 
    baga from the outside. Calling baga from baga is required when a package 
    has just been installed because a running python instance does not seem to 
    be able to detect newly installed packages.
    '''
    if system:
        # try system installed <== currently not implemented in CLI nor tested here
        x = _sys.path.pop(0)
        try:
            globals()['_'+package] = __import__(package)
            # this doesn't seem like a good idea if testing other system packages after this one . . .
            _sys.path.insert(0,'local_packages')
            
        except ImportError:
            print('Failed to import {} from system install'.format(package))
            _sys.path.insert(0,'local_packages')
            return(False)
    else:
        try:
            # import as with initial underscore to avoid name collisions
            # maybe should use importlib here
            globals()['_'+package] = __import__(package)
            
        except ImportError as e:
            print('Failed to import {} from system {}'.format(package, package_path))
            print(e)
            return(False)
        
    if hasattr(globals()['_'+package], '__version__'):
        v = getattr(globals()['_'+package], '__version__')
    elif hasattr(globals()['_'+package], 'version'):
        v = getattr(globals()['_'+package], 'version')
    else:
        v = 'unknown'
    
    if maj_version:
        if isinstance(v, str):
            try:
                this_maj_version = int(v.split('.')[0].split('_')[0])
            except ValueError:
                result = 'Could not check verions of this package: {} (either '\
                'package does not report version in a conventional way (unknown) '\
                'or this is a baga bug! Please raise an issue at '\
                'github.com/daveuu/baga.)'.format(v)
                _sys.stderr.write(result+'\n')
                #print(result)
        elif isinstance(v, tuple):
            this_maj_version = v[0]
        else:
            result = 'Could not check verions of this package: {} (either package '\
            'does not report version in a conventional way (unknown) or this is a '\
            'baga bug! Please raise an issue at github.com/daveuu/baga.'.format(v)
            _sys.stderr.write(result+'\n')
            #print(result)
        if maj_version == this_maj_version:
            result = 'Successfully imported {} required version {} from {}'.format(package, 
                                            v, getattr(globals()['_'+package], '__file__'))
            _sys.stderr.write(result+'\n')
            #print(result)
            return(True)
        elif maj_version > this_maj_version or maj_version < this_maj_version:
            result = 'Successfully imported {} at version {} from {}\n'.format(package, 
                                            v, getattr(globals()['_'+package], '__file__'))
            result += 'baga is not tested with this version. Version {} is required '\
            'which can be installed locally with the --get option'.format(maj_version)
            _sys.stderr.write(result+'\n')
            return(False)
    else:
        print('Successfully imported {} version {} from {}'.format(package, v, getattr(globals()['_'+package], '__file__')))
        print("baga's Dependencies module doesn't specify a version for this package so we'll assume it's OK")
        result = 'Successfully imported {} version {} from {}'.format(package, 
                                            v, getattr(globals()['_'+package], '__file__'))
        result += "baga's Dependencies module doesn't specify a version for this package "\
        "so we'll assume it's OK"
        _sys.stderr.write(result+'\n')
        return(True)



def check_git():
    try:
        version = _subprocess.check_output(['git','--version'])
    except OSError:
        print('Dependency error: git not found.')
        print('Please install git version 2.4 or more recent to continue')

    version = version.rstrip().split(' ')[-1].split('.')
    version_A, version_B, version_C = map(int, version)
    if not (version_A >= 2 and version_B >= 4):
        print('Dependency warning: git version < 2.4')
        print('''Proceeding but if git-related issues encountered,
        try updating git to version 2.4 or more recent''')
    else:
        print('git %s found: ok' % '.'.join(version))
# a method for local import that is explicit but fails for dendropy 4:
def local_import(dependency, module = False):
    '''Can provide a named module in package'''
    import imp as _imp
    package_path = [dependency['destination'], dependency['name']]
    get_module = dependency['name']
    package_path = _os.path.sep.join(package_path)
    _use_module = _imp.load_source(get_module, package_path)
    
    if module:
        _use_module = getattr(_use_module, module)
    
    return(_use_module)



def main():

    # get
    if dependencies[name]['source'] == 'git':
        get_git(**dependencies[name])
    elif dependencies[name]['source'] == 'download':
        get_download(**dependencies[name])

    # check (once checked, is loaded and cannot change between local and system)
    if dependencies[name]['destination'] == 'local_packages':
        checker = dependencies[name]['checker']['function']
        checker_args = dependencies[name]['checker']['arguments']
        checker(dependencies[name]['name'], **checker_args)
        #checker(**checker_args)
        #checker(dependencies[name]['name'], system = True)
    elif dependencies[name]['destination'] == 'external_programs':
        checker = dependencies[name]['checker']['function']
        checker_args = dependencies[name]['checker']['arguments']
        checker(**checker_args)
        #checker_args['system'] = True
        #checker(**checker_args)
dependencies = {}

# https://download.gnome.org/sources/librsvg/2.40/librsvg-2.40.9.tar.xz
# librsvg-2.40.9.sha256sum

# these need to be absolute for installing dependencies
destination_programs = _os.path.realpath('external_programs')
destination_packages = _os.path.realpath('local_packages')

dependencies['sickle'] = {
    'name': 'sickle',
    'description': 'per position quality score trimming of short reads',
    'source': 'git',
    'url': 'https://github.com/najoshi/sickle',
    'commit': 'f3d6ae3460e94ed8162d314b74fb470c27624432',
    'checksum': None,
    'destination': destination_programs,
    'preparation': [{'function': prep_simple_make,
                'arguments': {'path': False, 'configure': False}}],
    'checker': {'function': check_no_error,
                'arguments': {'path': ['sickle', 'sickle']}}
    }


dependencies['cutadapt'] = {
    'name': 'cutadapt',
    'description': 'remove library preparation artifacts from short reads',
    'source': 'git',
    'url': 'https://github.com/marcelm/cutadapt',
    'commit': '84bd634f674e09bb0cc61a5852ff3c87e0c8dea7',
    'checksum': None,
    'destination': destination_programs,
    # this one is awkward - could only install it to user's folder
    'preparation': [{'function': prep_python_install,
                'arguments': {'extras': ['--user']}}],
    'checker': {'function': check_no_error,
                'arguments': {'path': [_os.path.expanduser('~'), '.local','bin', 'cutadapt'],
                              # usually this is not required: 0 return code
                              # else is 1, the default success for non-zero
                              # but occassionally is 2 as in this case
                              'success_returncode':2}}
    }


dependencies['svgwrite'] = {
    'name': 'svgwrite',
    'description': 'generate plots in SVG files',
    'source': 'download',
    'url': 'https://pypi.python.org/packages/source/s/svgwrite/svgwrite-1.1.6.tar.gz',
    'commit': None,
    'checksum': 'md5=0d54ccf5584dd1f98fe22b7ac172bef1',
    'destination': destination_packages,
    'preparation': [{'function': get_other_packages,
                'arguments': {'package_list':[['pyparsing']], # really an argument list of length one which is another list
                              # this is the only thing that differentiates this prepare() from others that need some chdir in get_download()
                              # this should be improved
                              'just_packages':True}}],
    'checker': {'function': check_python_package,
                'arguments': {'maj_version':1}}
    }

dependencies['pyparsing'] = {
    'name': 'pyparsing',
    'description': 'fancy parsing for svgwrite!',
    'source': 'download',
    'url': 'https://pypi.python.org/packages/source/p/pyparsing/pyparsing-2.0.3.tar.gz',
    'commit': None,
    'checksum': 'md5=0fe479be09fc2cf005f753d3acc35939',
    'destination': destination_packages,
    'preparation': None,
    'checker': {'function': check_python_package,
                'arguments': {'maj_version':2}}
    }

dependencies['dendropy'] = {
    'name': 'dendropy',
    'description': 'phylogenetic computation',
    'source': 'download',
    'url': 'https://pypi.python.org/packages/source/D/DendroPy/DendroPy-4.0.3.tar.gz',
    'commit': None,
    'checksum': 'md5=aecd45bb1c5a04a35a812e764001d7ba',
    'destination': destination_packages,
    'preparation': None,
    'checker': {'function': check_python_package,
                'arguments': {'maj_version':4}}
    }

dependencies['bwa'] = {
    'name': 'bwa',
    'description': 'Burrows-Wheeler Aligner',
    'source': 'git',
    'url': 'https://github.com/lh3/bwa.git',
    # git rev-parse --short 981252192272bb5c5381af84ba68207ea8b24816 == 9812521
    'commit': '981252192272bb5c5381af84ba68207ea8b24816',
    'checksum': None,
    'destination': destination_programs,
    'preparation': [{'function': prep_simple_make,
                'arguments': {'path': False, 'configure': False}}],
    'checker': {'function': check_no_error,
                'arguments': {'path': ['bwa', 'bwa']}}
    }

dependencies['htslib'] = {
    'name': 'htslib',
    'description': 'htslib for samtools etc.',
    'source': 'git',
    'url': 'https://github.com/samtools/htslib',
    'commit': '1f45b1944cc5a9d9132328252c3a99aa259f4df6',
    'checksum': None,
    'destination': destination_programs,
    'preparation': [{'function': prep_simple_make,
                'arguments': {'path': False, 'configure': False}}],
    'checker': {'function': None,
                'arguments': {'path': None}}
    }

dependencies['samtools'] = {
    'name': 'samtools',
    'description': 'samtools for processing SAM and BAM files',
    'source': 'git',
    'url': 'https://github.com/samtools/samtools',
    # git rev-parse --short b0525f3e74836004f522c96a26bc7e660e032de4 == b0525f3
    'commit': 'b0525f3e74836004f522c96a26bc7e660e032de4',
    'checksum': None,
    'destination': destination_programs,
    'preparation': [
                {'function': get_git,
                 'arguments': dependencies['htslib']},
                              # having 'just_packages' key is the only thing that differentiates this 
                              # prepare() from others that need some chdir in get_download()
                              # this should be improved
                              #'just_packages':True}},
                {'function': prep_simple_make,
                 'arguments': {'path': False, 'configure': False}}
    ],
    'checker': {'function': check_no_error,
                'arguments': {'path': ['samtools', 'samtools']}}
    }

dependencies['biopython'] = {
    'name': 'Bio',
    'description': 'BioPython for reading, writing and downloading sequence data',
    'source': 'download',
    'url': 'https://pypi.python.org/packages/source/b/biopython/biopython-1.65.tar.gz',
    'commit': None,
    'checksum': 'md5=274a00e5629a135e84a8a7dfc0389935',
    'destination': destination_packages,
    'preparation': None,
    'checker': {'function': check_python_package,
                'arguments': {'maj_version':1}}
    }

dependencies['picard'] = {
    'name': 'picard',
    'description': 'manipulation of SAM and BAM files',
    'source': 'download',
    'url': 'https://github.com/broadinstitute/picard/releases/download/1.135/picard-tools-1.135.zip',
    'commit': None,
    'checksum': None,
    'destination': destination_programs,
    'preparation': None,
    'checker': {'function': check_no_error,
                'arguments': {'java_commands': ['-Xmx5g', '-jar', _os.path.sep.join(['external_programs','picard-tools-1.135','picard.jar']), 'MarkDuplicates', '--version']}}
    }

dependencies['seq-align'] = {
    'name': 'seq-align',
    'description': 'global pairwise sequence alignment',
    'source': 'git',
    'url': 'https://github.com/noporpoise/seq-align',
    # git rev-parse --short 05893831281ac26cfe299777416ef28f29f908f1 == 0589383
    'commit': '2472ddd774c9d01aed0bbe34b2b64724f6f2bc07',
    'checksum': None,
    'destination': destination_programs,
    'preparation': [{'function': prep_simple_make,
                'arguments': {'path': False, 'configure': False}}],
    'checker': {'function': check_no_error,
                'arguments':{'path': ['seq-align', 'bin', 'needleman_wunsch']}}
    }

dependencies['pysam'] = {
    'name': 'pysam',
    'description': 'for processing SAM and BAM files',
    'source': 'download',
    'url': 'https://pypi.python.org/packages/source/p/pysam/pysam-0.8.3.tar.gz',
    'commit': None,
    'checksum': 'md5=b1ae2a8ec3c6d20be30b2bc1aa995bbf',
    'destination': destination_packages,
    'preparation': [{'function': pysam_install,
                'arguments': {}}],
    #'preparation': None,
    'checker': {'function': check_python_package,
                'arguments': {}}
    }

dependencies['clonalframeml'] = {
    'name': 'clonalframeml',
    'description': 'Recombination Detection',
    'source': 'git',
    'url': 'https://github.com/xavierdidelot/clonalframeml',
    # git rev-parse --short 2d793a36f5160e311a6d15f18bf5cf89360d5e70 == 2d793a3
    'commit': '2d793a36f5160e311a6d15f18bf5cf89360d5e70',
    'checksum': None,
    'destination': destination_programs,
    'preparation': [{'function': prep_simple_make,
                'arguments': {'path': ['src'], 'configure': False}}],
    'checker': {'function': check_no_error,
                'arguments': {'path': ['clonalframeml', 'src', 'ClonalFrameML']}}
    }

dependencies['phyml'] = {
    'name': 'phyml',
    'description': 'Maximum Likelihood phylogeny inference and parameter estimation',
    'source': 'git',
    'url': 'https://github.com/stephaneguindon/phyml',
    # git rev-parse --short e71c55362c32f9fca07b9c2f2288dc05699fed41 == e71c553
    'commit': 'e71c55362c32f9fca07b9c2f2288dc05699fed41',
    'checksum': None,
    'destination': destination_programs,
    'preparation': [{'function': prep_simple_make,
                'arguments': {'path': [], 'configure': True, 'alt_command': './confphy'}}],
    'checker': {'function': check_no_error,
                'arguments': {'path': ['phyml', 'src', 'phyml'], 'extra': ['--version']}}
    }

# retired
# dependencies['exonerate'] = {
    # 'name': 'exonerate',
    # 'description': 'sequence alignment',
    # 'source': 'git',
    # 'url': 'https://github.com/nathanweeks/exonerate',
    # 'commit': 'v2.4.0',
    # 'checksum': None,
    # 'destination': 'external_programs',
    # 'preparation': {'function': prep_simple_make,
                # 'arguments': {'path': False, 'configure': True}},
    # 'checker': {'function': check_no_error,
                # 'arguments':{'path': ['exonerate-2.2.0-x86_64', 'bin', 'exonerate']}}
    # }


dependencies['spades'] = {
    'name': 'spades',
    'description': 'short read de novo assembler',
    'source': 'download',
    'url': 'http://spades.bioinf.spbau.ru/release3.6.1/SPAdes-3.6.1-Linux.tar.gz',
    'commit': None,
    'checksum': None,
    'destination': destination_programs,
    'preparation': None,
    'checker': {'function': check_no_error,
                'arguments':{'path': ['SPAdes-3.6.1-Linux', 'bin', 'spades.py'],
                             'extra_pre': ['python2']}}
    }

dependencies['discosnp'] = {
    'name': 'discosnp',
    'description': 'short read de novo variant caller',
    'source': 'download',
    'url': 'http://gatb-tools.gforge.inria.fr/versions/src/DiscoSNP++-2.2.1-Source.tar.gz',
    'commit': None,
    'checksum': None,
    'destination': destination_programs,
    'preparation': [{'function': discosnp_install,
                     'arguments': {}}],
    'checker': {'function': check_no_error,
                'arguments':{'path': ['DiscoSNP++-2.2.1-Source', 'run_discoSnp++.sh']}}
    }
dependencies_notes = {}

dependencies_notes['sickle'] = {
''
    }


dependencies_notes['cutadapt'] = {
''
    }

dependencies_notes['svgwrite'] = {
''
    }

dependencies_notes['dendropy'] = {
''
    }

dependencies_notes['bwa'] = {
''
    }

dependencies_notes['htslib'] = {
''
    }

dependencies_notes['samtools'] = {
''
    }

dependencies_notes['biopython'] = {
''
    }

dependencies_notes['picard'] = {
''
    }

dependencies_notes['exonerate'] = {
''
    }

dependencies_notes['exonerate'] = {
''
    }

dependencies_notes['seq-align'] = {
''
    }

dependencies_notes['pysam'] = {
''
    }

dependencies_notes['clonalframeml'] = {
''
    }

dependencies_notes['phyml'] = {
'Compilation from Git repository (i.e., from github.com) requires libtool among other things. Sometimes compute clusters are lacking a few essentials . . . Pre-compiled executables are available from http://www.atgc-montpellier.fr/phyml/binaries.php but require registration. After downloading, place binary in '
    }

dependencies_notes['clonalframeml'] = {
'Requires Python, supports: v2.4, 2.5, 2.6, 2.7, 3.2 and 3.3 (not 3.4 or 3.5)'
    }


dependencies_by_task = {}
dependencies_by_task['CollectData'] = [
'biopython',
]

dependencies_by_task['PrepareReads'] = [
'sickle',
'cutadapt',
'biopython',
]

dependencies_by_task['AlignReads'] = [
'bwa',
'samtools',
'picard',
# 'GATK' checked separately when path specified to GATK
]

dependencies_by_task['Repeats'] = [
'bwa',
'samtools',
'biopython',
'seq-align',
'svgwrite',
'pysam',
]

dependencies_by_task['Structure'] = [
'biopython',
'pysam',
'svgwrite',
'spades',
]

dependencies_by_task['CallVariants'] = [
# 'GATK' checked separately when path specified to GATK
]

dependencies_by_task['ApplyFilters'] = [

]

dependencies_by_task['ComparativeAnalysis'] = [
'phyml',
'clonalframeml',
'dendropy',
'biopython',
'svgwrite',
]


if __name__ == '__main__':
    main()
