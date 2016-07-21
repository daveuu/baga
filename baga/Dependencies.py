#! /usr/bin/env python3
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
Dependencies module from the Bacterial and Archaeal Genome Analyzer (BAGA).

This module contains functions to check whether required progams are available.
It also contains functions to get some of the dependencies.
'''

from baga import PY3 as _PY3

# stdlib
if _PY3:
    from baga import _BytesIO
    venv_dir = 'venv'
else:
    from baga import _StringIO
    venv_dir = 'venv2'

# stdlib
from baga import _subprocess
from baga import _os
from baga import _glob
from baga import _sys
from baga import _datetime
from baga import _logging
from baga import _urlopen

import zipfile as _zipfile
import tarfile as _tarfile
import hashlib as _hashlib
import shutil as _shutil
import stat as _stat

from baga import PROGRESS

# these need to be absolute for installing dependencies
destination_programs = _os.path.realpath('external_programs')
# local clones of git repos done manually
destination_dvcs_packages = _os.path.realpath('local_dvcs_package_repos')
# via pip even from a git repo
destination_packages = _os.path.realpath(_os.path.join(venv_dir,'lib',
        'python{}.{}'.format(_sys.version_info.major,_sys.version_info.minor),
        'site-packages'))
# get working folder
main_working_dir = _os.path.abspath(_os.curdir)
venv_path = _os.path.realpath(venv_dir)

# put venv modules at start to over-ride system versions
_sys.path.insert(0,destination_packages)




# could add non-zero exit status checks from ./configure etc with useful feedback?

# also add a handler for main i.e., to screen at DEBUG level? <==

def _configure_logging(log_folder, external_program, filename_suffix = False):
    if not filename_suffix:
        filename_suffix = 'for_{}'.format(external_program)
    time_start = _datetime.now()
    time_stamp = time_start.strftime("%d_%H-%M-%S-%f")
    filename = _os.path.join(log_folder, time_stamp+'__'+filename_suffix+'.log')
    this_ext_file_handler = _logging.FileHandler(filename, mode='w', encoding='utf-8')
    this_ext_file_handler.setLevel(_logging.DEBUG)
    this_module = __name__.split('.')[-1]
    logger_name = '{}_{}'.format(this_module,external_program)
    logger_for_external = _logging.getLogger(logger_name)
    logger_for_external.setLevel(_logging.DEBUG)
    logger_for_external.addHandler(this_ext_file_handler)
    # main logger name (which could be Task-level) will always match module 
    # name here (Dependencies)
    main_logger = _logging.getLogger(this_module)
    main_logger = _logging.LoggerAdapter(main_logger, {'task': this_module})
    main_logger.info('Writing output for {} to {}'.format(external_program, 
            filename))
    return(main_logger, logger_for_external)


# this is a lot like MetaSample.launch_external()
def _call_a_logged_process(cmd, logger, return_returncode = False):
    '''call a process via subprocess module and log outputs
    
    
    Logging for both BAGA Task and this module to timestamped files
    and terminal, depending on requested verbosity if called from the
    commandline interface.
    
    Parameters
    ----------
    cmd : list
        Command tokens for subprocess.Popen()
    logger : logger object
        With handlers etc for logging to various files
    return_returncode : bool
        Return the return code, or True for zero, False for non-zero
    '''
    this_module = __name__.split('.')[-1]
    main_logger = _logging.getLogger(this_module)
    main_logger = _logging.LoggerAdapter(main_logger, {'task': this_module})
    starting = 'Will call: {}'.format(' '.join(cmd))
    logger.log(PROGRESS, starting)
    main_logger.info(starting)
    try:
        proc = _subprocess.Popen(cmd, stdout = _subprocess.PIPE, 
                stderr = _subprocess.PIPE)
        # wait for completion
        stdout, stderr = proc.communicate()
        
        logger.log(PROGRESS, 'STDOUT for {} ({}) was:'.format(
                ' '.join(cmd), proc.pid))
        for line in stdout.decode('utf-8').split('\n'):
            logger.log(PROGRESS, '| '+line)
        
        logger.log(PROGRESS, 'STDERR for {} ({}) was:'.format(
                ' '.join(cmd), proc.pid))
        for line in stderr.decode('utf-8').split('\n'):
            logger.log(PROGRESS, '| '+line)
        
        logger.debug('Return code was {}'.format(proc.returncode))
        logger.debug('Call successful')
        if return_returncode:
            return(proc.returncode)
        else:
            return(True)
    except _subprocess.CalledProcessError as e:
        logger.warning('CalledProcessError for {}, return code: {}.'\
                ''.format(external_program, e.returncode))
        logger.debug(e)
        if return_returncode:
            return(e.returncode)
        else:
            return(False)
    except OSError as e:
        logger.warning('OSError for {}:'.format(cmd[0]))
        logger.debug(e)
        if return_returncode:
            # Typically FileNotFoundError which would be 2 if calling python
            return(2)
        else:
            return(False)


def checkmake_venv(venv_path = './venv', system_site_packages = False, with_pip = True):
    '''Check if a venv is present, create if not
    
    
    site-packages and local packages cannot be mixed i.e., 
    system_site_packages = True is not recommended for good version 
    control
    
    Parameters
    ----------
    venv_path : string
        Path to venv
    system_site_packages : bool
        To be passed to venv.create()
    with_pip : bool
        To be passed to venv.create()
    '''
    this_module = __name__.split('.')[-1]
    main_logger = _logging.getLogger(this_module)
    main_logger = _logging.LoggerAdapter(main_logger, {'task': this_module})
    if not _os.path.exists(_os.path.join(venv_path,'pyvenv.cfg')):
        main_logger.info('Building venv for Python packages')
        # from venv import EnvBuilder
        # eb = EnvBuilder(system_site_packages = system_site_packages, 
                # with_pip = with_pip)
        # eb.create(main_working_dir+_os.path.sep+'venv')
        import venv
        venv.create(venv_path, system_site_packages = system_site_packages, 
                with_pip = with_pip)
    else:
        main_logger.log(PROGRESS, 'Found venv for Python packages at "{}"'\
                ''.format(venv_path))

### specific get functions (called by get() below)


def get_pypi(pypi_name, name, log_folder = '.', pypi_version = False, **kwargs):
    '''python3.x pip built in venv'''
    checkmake_venv(venv_path)
    main_logger, logger_for_external = _configure_logging(log_folder, name, 
            filename_suffix = '{}_pip_get_from_pypi'.format(name))
    
    # **kwargs not getting dictionary supplied . . 
    # have to list all kwargs explicitly which is unexpected
    # try:
        # use_pypi_name = pypi_name+'=='+pypi_version
        # main_logger.info('Found a pypi version: {}'.format(use_pypi_name))
    # except NameError:
        # # if no version provided, install latest
        # use_pypi_name = pypi_name
        # main_logger.info('No pypi specified for {}'.format(use_pypi_name))
    if pypi_version:
        use_pypi_name = pypi_name+'=='+pypi_version
        main_logger.info('Found a pypi version: {}'.format(use_pypi_name))
    else:
        # if no version provided, install latest
        use_pypi_name = pypi_name
        main_logger.info('No pypi version specified for {}'.format(use_pypi_name))
    
    pip_exe = _os.path.join(venv_path,'bin','pip')
    _call_a_logged_process([pip_exe, 'install', use_pypi_name], logger_for_external)


def get_git(name, git, destination, preparation, log_folder = '.', 
        commit = False, **kwargs):
    '''
    Get a dependency from git (e.g. github.com or bitbucket.org)
    
    Preparation functions run after install (before use)
    '''
    main_logger, logger_for_external = _configure_logging(log_folder, name, 
            filename_suffix = '{}_git_clone'.format(name))
    
    # depency dict keys need reviewing
    # use of **kwargs means keys have to be match arguments exactly
    # key=git,value=repo_url forces argument name to git not url
    url = git
    
    if _os.path.realpath(_os.path.curdir) != destination:
        try:
            _os.chdir(destination)
        except OSError:
            _os.makedirs(destination)
            logger_for_external.debug('Created folder: {}'.format(destination))
            _os.chdir(destination)
    
    repo_name = url.split('/')[-1].replace('.git','')
    remove = True
    if _os.path.exists(repo_name):
        _os.chdir(repo_name)
        this_origin = _subprocess.check_output(['git', 'remote', 'get-url', 
                'origin']).decode('utf-8').rstrip()
        _os.chdir(_os.path.pardir)
        if this_origin == url:
            logger_for_external.info('repository in {} matches the name and '\
                    'origin of the requested: {}'.format(_os.path.realpath('.'), 
                    url))
            remove = False
        else:
            logger_for_external.info('repository in {} has a different origin '\
            '({}) to that requested: {}'.format(_os.path.realpath('.'), 
            this_origin, url))
    
    if remove:
        try:
            # clear any previous verions
            # wouldn't this be better updating from remote then checking out a specific version?
            _shutil.rmtree(repo_name)
            logger_for_external.debug('Removed previous git repository')
        except OSError:
            logger_for_external.warning('Failed to remove previous git repository!')
            pass
        git_server = url.replace('https://','').replace('http://','').split('/')[0]
        logger_for_external.log(PROGRESS,'Downloading {} via git from {} . . .'\
                ''.format(name, git_server))
        cmd = ['git', 'clone', url]
        call_successful = _call_a_logged_process(cmd, logger_for_external)
        if not call_successful:
            # raise a more general exception for the CLI?
            pass
    
    _os.chdir(repo_name)
    new_version = True
    last_commit = _subprocess.check_output(['git', 'show', '--name-status']
            ).decode('utf-8').split('\n')[0]
    logger_for_external.log(PROGRESS,'Last commit is: {}'\
            ''.format(last_commit))
    if not commit:
        logger_for_external.log(PROGRESS,'No commit requested')
    elif commit in last_commit:
        logger_for_external.log(PROGRESS, 'Already at requested commit: {}'\
                ''.format(commit))
        new_version = False
    else:
        logger_for_external.log(PROGRESS, 'Checking out requested commit: {}'\
                ''.format(commit))
        _call_a_logged_process(['git', 'checkout', commit], logger_for_external)
    
    # if repo uses git submodules, those will be set to the correct revisions 
    # for this commit else will do nothing
    _call_a_logged_process(['git', 'submodule', 'update', '--init'], 
            logger_for_external)
    
    working_dir = _os.path.sep.join([destination,repo_name])
    
    if preparation is not None:
        main_logger.info('Running preparation scripts for {} . . .'.format(name))
        for do_this in preparation:
            main_logger.log(PROGRESS, '{}()'.format(do_this['function'].__name__))
            ### if not new_version: do not run "make clean" else do <== to implement
            if 'package_list' in do_this['arguments']:
                do_this['function'](*do_this['arguments']['package_list'], 
                        logger = logger_for_external)
            elif isinstance(do_this['arguments'], dict):
                do_this['function'](logger = logger_for_external, 
                        **do_this['arguments'])
            else:
                do_this['function'](*do_this['arguments'], 
                        logger = logger_for_external)
            
            # restore position in path if a prepare changed it
            if working_dir != _os.path.realpath(_os.path.curdir):
                _os.chdir(working_dir)
    
    _os.chdir(_os.path.pardir)
    _os.chdir(_os.path.pardir)

def get_download(name, url, destination, preparation, checksum, 
        log_folder = '.', **kwargs):
    '''
    Download and unpack a dependancy
    '''
    # set up logging
    # set up handler for this subprocess
    # send to a unique time-stamped log file
    # also add a handler to the main log
    main_logger, logger_for_external = _configure_logging(log_folder, name, 
            filename_suffix = '{}_get_download'.format(name))
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
            e = "{} checksums not implemented in Python's hashlib!"\
                    "".format(hasher_algorithm)
            main_logger.error(e)
            raise NotImplementedError(e)
    
    main_logger.info('Downloading: %s' % url)
    req = _urlopen(url)
    CHUNK = 16 * 1024 * 64
    if _PY3:
        data = _BytesIO()
    else:
        data = _StringIO()
    c = 0
    for chunk in iter(lambda: req.read(CHUNK), ''):
        c += CHUNK
        main_logger.debug("{:,} bytes".format(c))
        data.write(chunk)
        if len(chunk) == 0:
            break
    
    main_logger.info('Download complete . . .')
    data.seek(0)
    
    if checksum:
        buff = data.read(65536)
        while len(buff) > 0:
            hasher.update(buff)
            buff = data.read(65536)
        
        e = '. . . checksum fail!'
        if hasher.hexdigest() == checksum.split('=')[1]:
            main_logger.error(e)
            raise Exception(e)
        else:
            main_logger.log(PROGRESS,'. . . checksum passed!')
        data.seek(0)
    
    if url[-6:] == 'tar.gz' or url[-3:] == 'tgz':
        archive = _tarfile.open(mode="r:gz", fileobj = data)
    elif url[-7:] == 'tar.bz2':
        archive = _tarfile.open(mode="r:bz2", fileobj = data)
    elif url[-4:] == '.zip':
        archive = _zipfile.ZipFile(data)
    else:
        raise Exception('did not recognise archive type for : {}'.format(url))
    
    if destination == destination_packages:
        # extract as a pypi python package
        release = url.split('/')[-1][:-7]
        main_logger.info('Extracting {} to {}'.format(release, 
                _os.path.sep.join([destination,name])))
        c = 0
        # this is obsolete with virtualenv
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
            # py3 vs py2: .getmembers() and .isreg() gone?
            # may be tar.gz specific not zip?
            for member in archive.getmembers():
                if member.isreg() and check_path1 in member.name:
                    member.name = _os.path.sep.join(member.name.split(
                            _os.path.sep)[1:])
                    archive.extract(member)
                    c += 1
        main_logger.log(PROGRESS, 'Extracted {} files'.format(c))
    
    else:
        # extract as a generic external program
        try:
            archive.extractall()
        except PermissionError as e:
            main_logger.log(PROGRESS, 'Attempting to overwrite existing version . . .')
            for f in archive.namelist()[::-1]:
                if _os.path.exists(f):
                    st = _os.stat(f)
                    _os.chmod(f, st.st_mode | _stat.S_IXUSR | _stat.S_IWUSR)
                    print(f,st.st_mode)
                    if _os.path.isfile(f):
                        _os.unlink(f)
                    else:
                        _os.rmdir(f)
            
            archive.extractall()
    
    if preparation:
        working_dir = _os.path.realpath(_os.path.curdir)
        main_logger.info('Running preparation scripts for {} . . .'.format(name))
        for do_this in preparation:
            if 'just_packages' in do_this['arguments']:
                # this is the only thing that differentiates this prepare()
                # from others that need some chdir <== this should be improved
                # see dep dict
                _os.chdir(_os.path.pardir)
                do_this['function'](*do_this['arguments']['package_list'], 
                        logger = logger_for_external)
            elif 'package_list' in do_this['arguments']:
                _os.chdir(_os.path.pardir)
                ### is this an alternative to above:
                ### check on latest dependencies dict format
                do_this['function'](*do_this['arguments']['package_list'], 
                        logger = logger_for_external)
            else:
                try:
                    extracted_base_dir = archive.namelist()[0].split(_os.path.sep)[0]
                except AttributeError:
                    extracted_base_dir = archive.getnames()[0].split(_os.path.sep)[0]
                _os.chdir(extracted_base_dir)
                if isinstance(do_this['arguments'], dict):
                    do_this['function'](logger = logger_for_external, **do_this['arguments'])
                else:
                    do_this['function'](*do_this['arguments'], logger = logger_for_external)
            
            # restore position in path (allows for unknown changes by a prepare)
            if working_dir != _os.path.realpath(_os.path.curdir):
                _os.chdir(working_dir)


def get_repo(name, url, destination, preparation, repo_exe = 'git', commit = False, 
        log_folder = '.', **kwargs):
    '''
    Get a dependency from a git or hg repository
    
    e.g. for git: github.com or bitbucket.org
    e.g. for mercurial (hg): bitbucket.org
    
    Preparation functions run after install (before use)
    '''
    # set up logging
    # set up handler for this subprocess
    # send to a unique time-stamped log file
    # also add a handler to the main log
    external_program = name
    time_start = _datetime.now()
    time_stamp = time_start.strftime("%d_%H-%M-%S-%f")
    # ("get_from_git" == clone and make and any other preparation)
    logfilename = '{}__{}_get_from_repo_{}.log'.format(time_stamp, external_program, repo_exe)
    filename = _os.path.join(log_folder, logfilename)
    this_ext_file_handler = _logging.FileHandler(filename, mode='w', encoding='utf-8')
    this_ext_file_handler.setLevel(_logging.DEBUG)
    logger_name = 'Dependencies_{}'.format(external_program)
    logger_for_external = _logging.getLogger(logger_name)
    logger_for_external.setLevel(_logging.DEBUG)
    logger_for_external.addHandler(this_ext_file_handler)
    # main logger name (which could be Task-level) will always match module 
    # name here (Dependencies)
    main_logger = _logging.getLogger('Dependencies')
    main_logger = _logging.LoggerAdapter(main_logger, {'task': 'Dependencies'})
    main_logger.info('Writing output for {} {} cloning and preparation to {}'\
            ''.format(external_program, repo_exe, filename))
    
    if _os.path.realpath(_os.path.curdir) != destination:
        try:
            _os.chdir(destination)
        except OSError:
            _os.makedirs(destination)
            logger_for_external.debug('Created folder: {}'.format(destination))
            _os.chdir(destination)
    
    try:
        # clear any previous verions
        # wouldn't this be better updating from remote then checking out a specific version?
        _shutil.rmtree(url.split('/')[-1].replace('.git','').replace('.hg',''))
        logger_for_external.debug('Removed previous {} repository'.format(repo_exe))
    except OSError:
        pass
    
    server = url.replace('https://','').replace('http://','').split('/')[0]
    logger_for_external.log(PROGRESS,'Downloading {} via {} from {} . . .'.format(name, repo_exe, server))
    cmd = [repo_exe, 'clone', url]
    call_successful = _call_a_logged_process(cmd, logger_name)
    if not call_successful:
        # raise a more general exception for the CLI?
        pass
    
    repo_name = url.split('/')[-1].replace('.git','').replace('.hg','')
    _os.chdir(repo_name)
    if commit:
        if repo_exe == 'git':
            _call_a_logged_process(['git', 'checkout', commit], logger_name)
        else:
            _call_a_logged_process(['hg', 'update', commit], logger_name)
    
    if repo_exe == 'git':
        # if repo uses git submodules, those will be set to the correct revisions 
        # for this commit else will do nothing
        _call_a_logged_process(['git', 'submodule', 'update', '--init'], logger_name)
        # apparently subrepo update in mercurial is 'on demand' i.e., automatic?
    
    working_dir = _os.path.sep.join([
            destination,url.split('/')[-1].replace('.git','').replace('.hg','')])
    
    if preparation is not None:
        main_logger.info('Running preparation scripts for {} . . .'.format(name))
        for do_this in preparation:
            main_logger.log(PROGRESS, '{}()'.format(do_this['function'].__name__))
            if 'package_list' in do_this['arguments']:
                do_this['function'](*do_this['arguments']['package_list'], logger_name = logger_name)
            elif isinstance(do_this['arguments'], dict):
                do_this['function'](logger_name = logger_name, **do_this['arguments'])
            else:
                do_this['function'](*do_this['arguments'], logger_name = logger_name)
            
            # restore position in path if a prepare changed it
            if working_dir != _os.path.realpath(_os.path.curdir):
                _os.chdir(working_dir)
    
    _os.chdir(_os.path.pardir)
    _os.chdir(_os.path.pardir)

def get_git_pip(name, 
                description, 
                source, 
                url, 
                destination, 
                preparation, 
                checker, 
                log_folder = '.',
                commit = False, 
                #pip = False, 
                #pypi_name = None, 
                #checksum = None
                ):
    '''
    Get a dependency from git using pip: preparation functions run before pip
    
    pip should handle dependencies where they are not handled or are non-Python, 
    use preparation to install them.
    '''
    
    # set up logging
    # set up handler for this subprocess
    # send to a unique time-stamped log file
    # also add a handler to the main log
    # set up logging
    # set up handler for this subprocess
    # send to a unique time-stamped log file
    # also add a handler to the main log
    main_logger, logger_for_external = _configure_logging(log_folder, name, 
            filename_suffix = '{}_get_git_pip'.format(name))
    
    working_dir = _os.path.abspath(_os.curdir)
    if preparation is not None:
        main_logger.info('Running {} preparation script(s)'\
                ''.format(len(preparation)))
        for do_this in preparation:
            # this will install pre-install dependencies for this git
            # python package
            if 'package_list' in do_this['arguments']:
                do_this['function'](*do_this['arguments']['package_list'], log_folder = log_folder)
            elif isinstance(do_this['arguments'], dict):
                do_this['function'](log_folder = log_folder, **do_this['arguments'])
            else:
                do_this['function'](*do_this['arguments'], log_folder = log_folder)
            
            # restore position in path if a prepare changed it
            if working_dir != _os.path.realpath(_os.path.curdir):
                _os.chdir(working_dir)
    
    checkmake_venv()
    
    cmd = ['{}/venv/bin/pip'.format(main_working_dir),
           'install',
           #'--editable', 
           # puts to: venv/src/biom-format/biom/__init__.py
           # would need to add each to path?
           'git+{}#egg={}'.format(url,name)]
    
    call_successful = _call_a_logged_process(cmd, logger_for_external)
    if not call_successful:
        # raise a more general exception for the CLI?
        pass

def get_other_packages(packages, log_folder):
    '''Python package dependencies for python packages'''
    if source == 'default_source':
        info = dependencies[name][dependencies[name][source]]
    else:
        info = dependencies[name][source]
    
    for package in packages:
        if info[package]['source'] == 'pypi':
            get_pypi(log_folder = log_folder, **dependencies[package])
        else:
            get_download(log_folder = log_folder, **dependencies[package])

def get_GATK():
    print("Visit https://www.broadinstitute.org/gatk/download/auth?package=GATK "\
    "and download after reviewing and agreeing license (and signing up to the "\
    "forum). Put jar file to external_programs/GenomeAnalysisTK/GenomeAnalysisTK.jar")


### general get function that calls specific get functions above

def get(name, log_folder, source = 'default_source'):
    '''get and install a dependency
    
    
    The information stored in the dependencies dictionary is used.
    
    Parameters
    ----------
    name : str
        Name of dependency corresponding to key in dependencies dictionary
    log_folder : str
        A usually time-stamped folder path to write log files into
    source : str
        One of "download", "repository" or "DVCS" depending on what information
        is available in the dependencies dictionary, else "default_source"
    '''
    if source == 'default_source':
        source = dependencies[name][source]
        
    info = dependencies[name][source]
    
    if source == 'repository' and 'pypi_name' in info:
        get_pypi(info['pypi_name'], name, log_folder = log_folder)
    elif source == 'DVCS' and 'pip_git' in info:
        get_git_pip(name, description = info['description'], source = info['git'], 
                url = info['git'], destination = destination_packages, 
                preparation = info['preparation'], checker = info['checker'], 
                log_folder = log_folder, commit = info['commit'])
    elif source == 'DVCS' and 'git' in info:
        get_git(name, info['git'], destination_programs, info['preparation'], 
                log_folder = log_folder, commit = info['commit'])
    elif source == 'download':
        get_download(info['name'], info['url'], info['destination'], 
                info['preparation'], info['checksum'], log_folder = log_folder)

def pysam_install(logger):
    _call_a_logged_process([_sys.executable, 'setup.py', 'build'], logger)
    for f in _os.listdir('build/lib.linux-x86_64-2.7/pysam'):
        if f[-3:] == '.so':
            _shutil.copy('build/lib.linux-x86_64-2.7/pysam/{}'.format(f), 'pysam/{}'.format(f))
    
    try:
        _shutil.rmtree('../pysam')
    except OSError:
        pass
    
    _os.rename('pysam', '../pysam')

def discosnp_install(logger):
    _call_a_logged_process(['bash','./compile_discoSnp++.sh'], logger)
    # ensure VCF generation works if system default python is version 3
    # VCF code is python 2 only
    fixed = open('run_VCF_creator.sh').read().replace('python ','python2 ')
    open('run_VCF_creator.sh','w').write(fixed)

def fix_python2(exepy, x = 1, **kwargs):
    '''replace first x instances of "python" to "python2"'''
    fixed = open(exepy).read().replace('python ' ,'python2 ', x)
    open(exepy,'w').write(fixed)

def fix_python2_shebang(exepy, x = 1, **kwargs):
    '''replace first x instances of "python" to "python2"'''
    fixme = open(exepy).readlines()
    fixme[0] = fixme[0].replace('python' ,'python2', x)
    open(exepy,'w').write(''.join(fixme))

def prep_python_install(extras, logger):
    _call_a_logged_process([_sys.executable, 'setup.py', 'install'] + extras, 
            logger)

def prep_simple_make(logger, path = False, configure = False, alt_command = False):
    if path:
        _os.chdir(_os.path.sep.join(path))
    
    if configure:
        _call_a_logged_process(['./configure'], logger)
    
    if alt_command:
        _call_a_logged_process([alt_command], logger)
    else:
        _call_a_logged_process(['make'], logger)
    
    if path:
        _os.chdir(_os.path.sep.join([_os.path.pardir]*len(path)))

def chmod_xr(logger_name, pattern):
    files = _glob(pattern)
    for f in files:
        st = _os.stat(f)
        _os.chmod(f, st.st_mode | _stat.S_IXUSR | _stat.S_IRUSR)



### specific check functions (called by general check functions below)

def check_no_error(name, path = '', 
                   log_folder = '.',
                   system = False, 
                   java_commands = False, 
                   extra = False, 
                   extra_pre = False,
                   success_returncode = 0):
    '''
    check a binary executable can be called.
    
    A failure is considered as OSError (usually FileNotFound error but could 
    also be an access problem), not a non-zero exit status. Some programs 
    return non-zero exit status if no commands given.
    
    Can check in default baga location or a specified path or the system path 
    or java. Additional commands and be appended with 'extra' list or pre-
    pended with 'extra-pre' list e.g., specifying 'python2' for when python3 
    defaults but python >3.3 is installed but not supported.
    
    Parameters
    ----------
    path : str
        Path to executable to try to run
    log_folder : str
        A usually time-stamped folder path to write log files into
    system : bool
        Check system install, not local install
    extra : list
        Extra arguments to append to system call
    extra_pre : list
        Extra arguments to prepend to system call
    success_returncode: int
        Return code that means the program was found
    '''
    if not _os.path.exists(log_folder):
        _os.makedirs(log_folder)
    
    
    main_logger, logger_for_external = _configure_logging(log_folder, name, 
            filename_suffix = '{}_check'.format(name))
    
    if java_commands:
        # have to extract name of the actual .jar for log file name
        path_to_jar = java_commands[java_commands.index('-jar')+1]
        external_program = path_to_jar.split(_os.path.sep)[-1][:-4]
    elif len(path) > 0:
        external_program = path[-1]
    
    if java_commands:
        # could do with a check for java install here
        p = _subprocess.Popen(['java'] + java_commands, stdout=_subprocess.PIPE, 
                stderr=_subprocess.STDOUT)
        output = p.stdout.read()
        output = output.decode('utf-8')
        logger_for_external.log(PROGRESS, 'STDOUT and STDERR for {} ({}) was:'.format(
                external_program, p.pid))
        for line in output.split('\n'):
            logger_for_external.log(PROGRESS, '| '+line)
        if 'Error:' in output:
            return(False)
        else:
            return(True)
    else:
        if system:
            # check the version that is in the system path:
            # just take name of executable, omit the path to it (if present)
            cmd = path[-1:]
        elif path[0][0] == _os.path.sep or path[0] == '':
            # absolute path to elsewhere provided
            cmd = [_os.path.sep.join(path)]
        else:
            # path is to local install path
            cmd = [_os.path.sep.join([destination_programs] + path)]
        
        if extra:
            cmd += extra
        
        if extra_pre:
            cmd = extra_pre + cmd
        
        # use _call_a_logged_process() here?
        returncode = _call_a_logged_process(cmd, logger_for_external, 
                return_returncode = True)
        if returncode == success_returncode:
            return(True)
        else:
            return(False)


def check_python_package(package, major_version = -1, minor_version = -1, 
        system = False, package_path = destination_packages, **kwargs):
    '''Check if a package can be imported from the locally installed collection
    
    Print the results to stderr so they can be collected separately if calling 
    baga_cli from the outside. Calling baga_cli from baga is required when a 
    package has just been installed because a running python instance does not 
    seem to be able to detect newly installed packages.
    
    Parameters
    ----------
    package : str
        Name of package to attempt to import
    major_version : int
        To compare with major version of imported module
    minor_version : int
        To compare with minor version of imported module
    package_path : string
        Path where package expected to be installed
    system : bool
        Enforce checking in system python installation not local (not 
        implemented)
    '''
    # set up logging
    # unlike other check() for externals, don't need additional log file
    # main logger name (which could be Task-level) will always match module 
    # name here (Dependencies) == __name__.split('.')[-1]
    
    main_logger = _logging.getLogger(__name__.split('.')[-1])
    main_logger = _logging.LoggerAdapter(main_logger, {'task': __name__.split('.')[-1]})
    main_logger.info('Checking whether {} package available'.format(package))
    # add a handler to stderr ?
    if system:
        # try system installed <== currently not implemented in CLI nor tested here
        x = _sys.path.pop(0)
        try:
            globals()['_'+package] = __import__(package)
            # this doesn't seem like a good idea if testing other system packages after this one . . .
            _sys.path.insert(0,destination_packages)
            
        except ImportError:
            print('Failed to import {} from system install'.format(package))
            _sys.path.insert(0,destination_packages)
            return(False)
    else:
        try:
            # import as with initial underscore to avoid name collisions
            # maybe should use importlib here
            globals()['_'+package] = __import__(package)
            
        except ImportError as e:
            main_logger.error('Failed to import {} from {}'.format(
                    package, package_path))
            main_logger.error(e)
            return(False)
    
    if hasattr(globals()['_'+package], '__version__'):
        v = getattr(globals()['_'+package], '__version__')
    elif hasattr(globals()['_'+package], 'version'):
        v = getattr(globals()['_'+package], 'version')
    else:
        v = 'unknown'
    
    origin = getattr(globals()['_'+package], '__file__')
    imported = 'Successfully imported {} version {} from {}\n'.format(
            package, v, origin)
    
    if major_version >= 0:
        failresult = "Could not check major version of this package: {} (either "\
                "package does not report version in a conventional way (unknown) "\
                "or this is a baga bug! Please raise an issue at "\
                "github.com/daveuu/baga.). For now, we'll assume this is OK if "\
                "you just installed this package".format(v)
        if isinstance(v, str):
            try:
                this_major_version = int(v.split('.')[0].split('_')[0])
            except ValueError:
                main_logger.warning(failresult)
                # imported but couldn't check version. Not usually fatal.
                return(True)
            if minor_version >= 0:
                try:
                    #this_minor_version = int(v.split('.')[1].split('_')[0])
                    this_minor_version = int(v.split('.')[1])
                except ValueError:
                    main_logger.warning(failresult.replace('major','minor'))
                    # imported but couldn't check version. Not usually fatal.
                    return(True)
        elif isinstance(v, tuple):
            this_major_version = v[0]
            if minor_version:
                this_minor_version = v[1]
        else:
            main_logger.error(failresult)
        main_logger.info(imported)
        if minor_version >= 0:
            requested = '{}.{}'.format(major_version,minor_version)
        else:
            requested = '{}.x'.format(major_version)
        
        mediumresult = 'version {} installed but older version {} '\
                'specified in Dependencies (probably OK if backwards-compatibility '\
                'has been maintained)\n'.format(v, requested)
        badresult = 'Version {} installed but Dependencies module specifies '\
                'version {}: baga is not tested with this version. "Dependencies '\
                '--get {}" should install a newer and compatible version.\n'\
                ''.format(v, requested, package)
        if major_version == this_major_version:
            if minor_version >= 0:
                if minor_version == this_minor_version:
                    return(True)
                elif minor_version < this_minor_version:
                    main_logger.warning(mediumresult)
                    return(True)
                else:
                    main_logger.error(badresult)
                    return(False)
            else:
                return(True)
        elif major_version < this_major_version:
            main_logger.info(mediumresult)
            return(True)
        else:
            main_logger.warning(badresult)
            return(False)
    else:
        main_logger.info(imported)
        main_logger.info("baga's Dependencies module doesn't specify a version "\
                "for this package so we'll assume it's OK\n")
        return(True)


def check_git():
    # set up logging
    # unlike other check() for externals, don't need additional log file
    # main logger name (which could be Task-level) will always match module 
    # name here (Dependencies)
    this_module = __name__.split('.')[-1]
    main_logger = _logging.getLogger(this_module)
    main_logger = _logging.LoggerAdapter(main_logger, {'task': this_module})
    main_logger.info('Writing output for {} to {}'.format(external_program, 
            filename))
    try:
        version = _subprocess.check_output(['git','--version'])
        ## p23
        version = version.decode('utf-8')
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

### general check functions that call specific check functions above

def check(name, log_folder, source = 'default_source'):
    '''check whether a dependency is available
    
    
    The information stored in the dependencies dictionary is used.
    
    Parameters
    ----------
    name : str
        Name of dependency corresponding to key in dependencies dictionary
    log_folder : str
        A usually time-stamped folder path to write log files into
    source : str
        One of "download", "repository" or "DVCS" depending on what information
        is available in the dependencies dictionary, else "default_source"
    '''
    if source == 'default_source':
        info = dependencies[name][dependencies[name][source]]
    else:
        info = dependencies[name][source]
    checker = info['checker']['function']
    checker_args = info['checker']['arguments']
    if checker == check_python_package:
        # its a python package
        result = checker(info['import_name'], log_folder = log_folder, **checker_args)
    else:
        result = checker(name, log_folder = log_folder, **checker_args)
    return(result)

def checkpackage(name, log_folder, **kwargs):
    '''call a conventional --check for checking new python packages'''
    this_module = __name__.split('.')[-1]
    task_logger = _logging.getLogger(this_module)
    task_logger = _logging.LoggerAdapter(task_logger, {'task': this_module})
    import subprocess
    task_logger.log(PROGRESS, 'Will check whether a Python package "{}" is available '\
            '(by calling baga again)'.format(name))
    if _PY3:
        cmd = ['python3', _sys.argv[0], '--nosplash', this_module, '--check', name]
    else:
        cmd = ['python2', _sys.argv[0], '--nosplash', this_module, '--check', name]
    task_logger.log(PROGRESS, 'Will call: {}'.format(' '.join(cmd)))
    # add try..except here?
    proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE )
    o,e = proc.communicate()
    if proc.returncode == 1:
        return(False)
    elif proc.returncode == 0:
        return(True)

## still getting false negative successful installs . . .
def checkget(checkgetthese, log_folder):
    '''check whether a list of dependencies are available and get if absent
    
    
    The information stored in the dependencies dictionary is used.
    
    Parameters
    ----------
    checkgetthese : list
        Of (name, source) tuples where name corresponding to key in 
        dependencies dictionary and source is a level one key, possibly
        "default_source" else "download", "repository" or "DVCS".
    log_folder : str
        A usually time-stamped folder path to write log files into
    '''
    this_module = __name__.split('.')[-1]
    task_logger = _logging.getLogger(this_module)
    task_logger = _logging.LoggerAdapter(task_logger, {'task': this_module})
    # store summary infor organised by package in a dict
    check_summary = {}
    # store last install status after finding or attempting install
    check_results = {}
    task_logger.log(PROGRESS, 'Will check, and if necessary get: {}'.format(
            ', '.join(['{} ({})'.format(name,source) for name,source in \
            checkgetthese])))
    for (name,source) in checkgetthese:
        task_logger.log(PROGRESS, 'Checking for {}'.format(name))
        check_summary[name] = []
        alreadygot = check(name, log_folder, source = source)
        if alreadygot:
            task_logger.info('Found: {}'.format(name))
            check_summary[name] += ['\t{}: found!'.format(name)]
            check_results[name] = alreadygot
        else:
            task_logger.log(PROGRESS, 'Not found: {}'.format(name))
            task_logger.log(PROGRESS, 'Attempting to install: {}'.format(name))
            check_summary[name] += ['\t{0}: not found . . . '.format(name),
                    'Attempting to install {}.'.format(name)]
            if source == 'default_source':
                source = dependencies[name][source]
            
            info = dependencies[name][source]
            
            get(name, log_folder, source = source)
            
            if info['checker']['function'] == check_python_package:
                # necessary to call baga_cli again if just installed a package
                gotnow = checkpackage(name, log_folder, source = source)
            else:                
                gotnow = check(name, log_folder, source = source)
            
            if gotnow:
                task_logger.info('Installed successfully! Found: {}'.format(name))
                check_summary[name] += ["Installed successfully: found!"]
            else:
                task_logger.error('Failed to install: {}. (missing "\
                        "dependencies? Bad download URL?)'\
                        ''.format(name))
                check_summary[name] += ["Failed to install (missing "\
                        "dependencies? Bad download URL?)"]
            
            # could still be False for a failed install
            check_results[name] = gotnow
        
    return(check_summary,check_results)



# these can be installed by BAGA

# where and how to get dependencies are described in the
# baga.Dependencies.dependencies dictionary:
# level one key is name of package, module or program
# level two keys are sources and the default to use
# there are three possible source types:
    # distributed version control system (DVCS):
        # git via github.com or bitbucket.org,
        # mercurial via bitbucket.org
    # repository:
        # pip via pypi
        # (could add apt-get, yum etc if they support local install?)
    # download:
        # binaries in a tarball
        # source+compile in a tarball
# at least one must be present to be useful!
# there is a default to use if not explicitly specified by the function call

dependencies = {}

dependencies['sickle'] = {
    'default_source':'DVCS',
    'DVCS': {
        'name': 'sickle',
        'description': 'per position quality score trimming of short reads',
        'destination': destination_programs,
        'git': 'https://github.com/najoshi/sickle',
        'commit': 'f3d6ae3460e94ed8162d314b74fb470c27624432',
        'preparation': [{'function': prep_simple_make,'arguments': 
                {'path': False, 'configure': False}}],
        'checker': {'function': check_no_error,'arguments': 
                {'path': ['sickle', 'sickle']}}
        }
    }

dependencies['cutadapt'] = {
    'default_source':'DVCS',
    ## pypi can install from a git repo: if this went into a virtualenv, where would the executable go? <== how to know when to pypi from a DVCS?
    'DVCS': {
        'name': 'cutadapt',
        'description': 'remove library preparation artifacts from short reads',
        # clone with git, use specific preparation scripts
        'git': 'https://github.com/marcelm/cutadapt',
        # install with pip (not yet implemented) can commit be specified?
        # 'pip_git': 'https://github.com/marcelm/cutadapt',
        'commit': '3bee1330cab4a50fb37b02cd0c82fa5203e7fbf5',
        'destination': destination_programs,
        # this one is awkward - could only install it to user's folder
        'preparation': [{'function': prep_python_install, 'arguments': 
                {'extras': ['--user']}}],
        'checker': {'function': check_no_error, 'arguments': {'path': 
                [_os.path.expanduser('~'), '.local','bin','cutadapt'],
                # usually this is not required: 0 return code
                # else is 1, the default success for non-zero
                # but occassionally is 2 as in this case
                'success_returncode':2}}
        },
    'repository': {
        'import_name': 'cutadapt',
        'pypi_name': 'cutadapt',
        'description': 'remove library preparation artifacts from short reads',
        'source': 'pypi',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':1,'minor_version':9,'micro_version':1}}
        }
    }

dependencies['bwa'] = {
    'default_source':'DVCS',
    'DVCS': {
        'name': 'bwa',
        'description': 'Burrows-Wheeler Aligner',
        'git': 'https://github.com/lh3/bwa.git',
        'commit': 'bb7533e32ca55d53e842c045721b4ff4b407a857',
        'destination': destination_programs,
        'preparation': [{'function': prep_simple_make,
                    'arguments': {'path': False, 'configure': False}}],
        'checker': {'function': check_no_error,
                    'arguments': {'path': ['bwa', 'bwa']}}
        }
    }

dependencies['htslib'] = {
    'default_source':'DVCS',
    'DVCS': {
        'name': 'htslib',
        'description': 'htslib for samtools etc.',
        'git': 'https://github.com/samtools/htslib',
        'commit': '1dee8d5a6fb6b10debd7562bf0abd6f5869693a3',
        'destination': destination_programs,
        'preparation': [{'function': prep_simple_make,
                    'arguments': {'path': False, 'configure': False}}],
        'checker': {'function': None,
                    'arguments': {'path': None}}
        }
    }

dependencies['samtools'] = {
    'default_source':'DVCS',
    'DVCS': {
        'name': 'samtools',
        'description': 'samtools for processing SAM and BAM files',
        'git': 'https://github.com/samtools/samtools',
        'commit': 'd49c73b9ce5f0b169f82d8c79e2e411a581dc5c7',
        'destination': destination_programs,
        'preparation': [
                    {'function': get_git,
                    # default_source is only way to specify source for dependencies
                     'arguments': dependencies['htslib'][dependencies['htslib']['default_source']]},
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
    }


dependencies['picard'] = {
    'default_source':'download',
    'download': {
        'name': 'picard',
        'description': 'manipulation of SAM and BAM files',
        'url': 'https://github.com/broadinstitute/picard/releases/download/'\
                '1.135/picard-tools-1.135.zip',
        'checksum': None,
        'destination': destination_programs,
        # typical for pre-compiled binaries
        'preparation': None,
        'checker': {'function': check_no_error,
                    'arguments': {'java_commands': ['-Xmx5g', '-jar', 
                            _os.path.sep.join([destination_programs,
                            'picard-tools-1.135','picard.jar']), 
                            'MarkDuplicates', '--version']}}
        }
    }

dependencies['seq-align'] = {
    'default_source':'DVCS',
    'DVCS': {
        'name': 'seq-align',
        'description': 'global pairwise sequence alignment',
        'git': 'https://github.com/noporpoise/seq-align',
        'commit': 'd0aa4804a93e46fdf267ef4b89f7cc66f2f0eeb8',
        'destination': destination_programs,
        'preparation': [{'function': prep_simple_make,
                    'arguments': {'path': False, 'configure': False}}],
        'checker': {'function': check_no_error,
                    'arguments':{'path': ['seq-align', 'bin', 
                            'needleman_wunsch']}}
        }
    }

dependencies['spades'] = {
    'default_source':'download',
    'download': {
        'name': 'spades',
        'description': 'short read de novo assembler',
        'url': 'http://spades.bioinf.spbau.ru/release3.6.1/'\
                'SPAdes-3.6.1-Linux.tar.gz',
        'checksum': None,
        'destination': destination_programs,
        'preparation': None,
        'checker': {'function': check_no_error, 'arguments':{
                'path': ['SPAdes-3.6.1-Linux', 'bin', 'spades.py'], 
                # for python2, returncode 2 is no file found,
                'extra_pre': ['python2']}}
        }
    }

# originally from i4m
dependencies['bowtie2'] = {
    'default_source':'download',
    'download': {
        'name': 'bowtie2',
        'description': 'short read aligner',
        'url': 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/'\
                '2.2.6/bowtie2-2.2.6-linux-x86_64.zip',
        'checksum': None,
        'destination': destination_programs,
        'preparation': [
                    {'function': chmod_xr,
                     'arguments': {'pattern':'bowtie2*'}}
        ],
        'checker': {'function': check_no_error,
                    'arguments': {'path': ['bowtie2-2.2.6', 'bowtie2']}},
        'exe': {
                'main':_os.path.sep.join([destination_programs,
                        'bowtie2-2.2.6', 'bowtie2'])},
        }
    }

dependencies['seqtk'] = {
    'default_source':'DVCS',
    'DVCS': {
        'name': 'seqtk',
        'description': 'manipulation of fastq files',
        'git': 'https://github.com/lh3/seqtk',
        'commit': '4feb6e81444ab6bc44139dd3a125068f81ae4ad8',
        'destination': destination_programs,
        'preparation': [{'function': prep_simple_make,
                    'arguments': {'path': False, 'configure': False}}],
        'checker': {'function': check_no_error,
                    'arguments':{'path': ['seqtk', 'seqtk']}},
        'exe': {
                'main':_os.path.sep.join([destination_programs, 
                        'seqtk', 'seqtk'])}
        }
    }


dependencies['biopython'] = {
    'default_source':'repository',
    'repository': {
        'pypi_name': 'biopython',
        'import_name': 'Bio',
        'description': 'BioPython for reading, writing and downloading '\
                'sequence data',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':1,'minor_version':6}}
        }
    }

dependencies['pysam'] = {
    'default_source':'repository',
    'repository': {
        'pypi_name': 'pysam',
        'import_name': 'pysam',
        'description': 'for processing SAM and BAM files',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':0, 'minor_version':8}}
        }
    }

dependencies['svgwrite'] = {
    'default_source':'repository',
    'repository': {
        'pypi_name': 'svgwrite',
        'import_name': 'svgwrite',
        'description': 'generate plots in SVG files',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':1}}
        }
    }

dependencies['dendropy'] = {
    'default_source':'repository',
    'repository': {
        'pypi_name': 'dendropy',
        'import_name': 'dendropy',
        'description': 'phylogenetic computation',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':4,'minor_version':1}}
        }
    }

dependencies['pytest'] = {
    'default_source':'repository',
    'repository': {
        'pypi_name': 'pytest',
        'import_name': 'pytest',
        'description': 'pytest',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':2,'minor_version':8}}
        }
    }

# may eventually want system numpy if better optimised?
# may need to set venv to look at system packages?
dependencies['numpy'] = {
    'default_source':'repository',
    'repository': {
        'pypi_name': 'numpy',
        'import_name': 'numpy',
        'description': 'numerical python',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':1,'minor_version':10}}
        }
    }

# biom format needs numpy already installed
# but knows it needs to install scipy via pip . . .
# but scipy requires a Fortran compiler for library dfftpack . . . send user warning?
dependencies['biom'] = {
    'default_source':'DVCS',
    'DVCS': {
        'import_name': 'biom',
        'description': 'The Biological Observation Matrix (BIOM) format',
        'git': 'https://github.com/biocore/biom-format',
        'pip_git': 'https://github.com/biocore/biom-format',
        'preparation': [{'function': get_other_packages,
                    'arguments': {'package_list':[['numpy']], # really an argument list of length one which is another list
                                  # this is the only thing that differentiates this prepare() from others that need some chdir in get_download()
                                  # this should be improved
                                  'just_packages':True}}],
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':2}}
        },
    'repository': {
        'import_name': 'biom',
        'pypi_name': 'biom-format',
        'description': 'The Biological Observation Matrix (BIOM) format',
        'preparation': [{'function': get_other_packages,
                    'arguments': {'package_list':[['numpy']], # really an argument list of length one which is another list
                                  # this is the only thing that differentiates this prepare() from others that need some chdir in get_download()
                                  # this should be improved
                                  'just_packages':True}}],
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':2}}
        }
    }

#h5py

dependencies['mash'] = {
    'default_source':'download',
    'download': {
        'name': 'mash',
        'description': 'fast genome distances',
        'url': 'https://github.com/marbl/Mash/releases/download/v1.0.2/'\
                'mash-Linux64-v1.0.2.tar.gz',
        'checksum': None,
        'preparation': None,
        'destination': destination_programs,
        'checker': {'function': check_no_error,
                    'arguments': {'path': ['mash']}},
        'exe': {
                'main':_os.path.sep.join([destination_programs,
                        'mash'])},
        }
    }
#wget --output-document=mash-Linux64-v1.0.2.tar.gz https://github.com/marbl/Mash/releases/download/v1.0.2/mash-Linux64-v1.0.2.tar.gz

dependencies['pda'] = {
    'default_source':'download',
    'download': {
        'name': 'pda',
        'description': 'Phylogenetic Diversity Analyzer',
        'url': 'http://www.cibiv.at/software/pda/download/pda-1.0.3/'\
                'pda-1.0.3-Linux.tar.gz',
        'checksum': None,
        'preparation': None,
        'destination': destination_programs,
        'checker': {'function': check_no_error,
                    'arguments': {'path': ['pda-1.0.3-Linux','bin','pda']}},
        'exe': {
                'main':_os.path.sep.join([destination_programs,
                        'pda-1.0.3-Linux','bin','pda'])},
        }
    }

dependencies['cd-hit'] = {
    'default_source':'download',
    'download': {
        'name': 'cd-hit',
        'description': 'short read aligner',
        'url': 'https://github.com/weizhongli/cdhit/releases/download/V4.6.5/'\
                'cd-hit-v4.6.5-2016-0304.tar.gz',
        'checksum': None,
        'destination': destination_programs,
        'preparation': [
                    {'function': prep_simple_make,
                    'arguments': {'path': False, 'configure': False}}
        ],
        'checker': {'function': check_no_error,
                    'arguments': {'path': ['cd-hit-v4.6.5-2016-0304', 'cd-hit'],
                            'success_returncode':1}},
        'exe': {
                'main':_os.path.sep.join([destination_programs,
                        'cd-hit-v4.6.5-2016-0304', 'cd-hit'])},
        }
    }

dependencies['roaringbitmap'] = {
    'default_source':'repository',
    'repository': {
        'pypi_name': 'roaringbitmap',
        'import_name': 'roaringbitmap',
        'description': 'roaringbitmap fast Set implementation',
        'checker': {'function': check_python_package,
                    'arguments': {'major_version':0,'minor_version':3}}
        }
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
'seq-align',
]

dependencies_by_task['CallVariants'] = [
# 'GATK' checked separately when path specified to GATK
'pysam'
]

dependencies_by_task['FilterVariants'] = [

]

dependencies_by_task['ComparativeAnalysis'] = [
'phyml',
'clonalframeml',
'dendropy',
'biopython',
'svgwrite',
]

dependencies_by_task['SimulateReads'] = [
'gemsim',
]

dependencies_by_task['Homology'] = [
'cd-hit',
'mash',
'pda',
'roaringbitmap'
]




# these are easier installed by a GNU/Linux package manager

sys_dependencies = {}

sys_dependencies['time'] = {
    'required_by': ['lmat'],
    'linux_flavour':'all'
}

sys_dependencies['build-essentials'] = {
    'required_by': ['Dependencies building from source'],
    'linux_flavour':'Debian'
}


def main():
    # get
    if dependencies[name]['source'] == 'git':
        if dependencies[name]['pip']:
            get_git_pip(**dependencies[name])
        else:
            get_git(**dependencies[name])
    elif dependencies[name]['source'] == 'download':
        get_download(**dependencies[name])
    elif dependencies[name]['source'] == 'pypi':
        get_pypi(**dependencies[name])

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


if __name__ == '__main__':
    main()
