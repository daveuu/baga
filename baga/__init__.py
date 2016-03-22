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

import gzip as _gzip
import json as _json
import os as _os
import subprocess as _subprocess
import multiprocessing as _multiprocessing
import sys as _sys
import re as _re
import tarfile as _tarfile
import json as _json
import time as _time
from glob import glob as _glob
from datetime import datetime as _datetime
import shlex as _shlex
from hashlib import md5 as _md5
from array import array as _array
import logging as _logging
import textwrap as _textwrap

PY3 = _sys.version_info > (3,)

if PY3:
    from urllib.error import URLError as _URLError
    from urllib.request import urlopen as _urlopen
    from io import BytesIO as _BytesIO
    import pickle as _pickle
else:
    from urllib import urlopen as _urlopen
    from urllib2 import URLError as _URLError
    from cStringIO import StringIO as _StringIO
    import cPickle as _pickle
  



## define formatters and filters etc for logging

try:
    # get width for wrapping
    info = _subprocess.check_output(['stty','size'])
    rows, text_width = map(int,info.rstrip().split())
except FileNotFoundError:
    # else pick a narrow one
    rows, text_width = 30, 70

# for main CLI logger, make a formatter for writing to files
# include date, time and level - omit year for space and because in file name
formatter_tofile = _logging.Formatter('%(asctime)s %(levelname)s: %(message)s', 
        '%m/%d %H:%M:%S')

# for INFO+ from Task to main
# as above but include the task and module name
formatter_tofile_task_to_main = _logging.Formatter(
        '%(asctime)s %(task)s %(name)s %(levelname)s: %(message)s', 
        '%m/%d %H:%M:%S')

# clever line wrapping formatter Class for stdout
class WrappedFixedIndentingLog(_logging.Formatter):
    def __init__(self, fmt=None, datefmt=None, style='%', width=70, indent=4):
        if PY3:
            super().__init__(
                    fmt=fmt, datefmt=datefmt, style=style)
        else:
            # no style key word in Python 2
            super(WrappedFixedIndentingLog, self).__init__(
                    fmt=fmt, datefmt=datefmt)
        self.wrapper = _textwrap.TextWrapper(width=width, 
                subsequent_indent=' '*indent)
    
    def format(self, record):
        return self.wrapper.fill(super(WrappedFixedIndentingLog, self).format(record))

# use a simplified formatter for printing to console at CLI level
wrapping_formatter_stdout = WrappedFixedIndentingLog(
        '%(asctime)s %(message)s', '%H:%M:%S', 
        width = text_width, indent = 4)

# include "level" in format for errors, warnings etc
wrapping_formatter_stderr = WrappedFixedIndentingLog(
        '%(asctime)s %(levelname)s: %(message)s', '%H:%M:%S', 
        width = text_width, indent = 4)

# for Task-level loggers to console, include task name and module
wrapping_formatter_stdout_task = WrappedFixedIndentingLog(
        '%(asctime)s %(task)s %(name)s: %(message)s', '%H:%M:%S', 
        width = text_width, indent = 4)

wrapping_formatter_stderr_task = WrappedFixedIndentingLog(
        '%(asctime)s %(task)s %(name)s %(levelname)s: %(message)s', '%H:%M:%S', 
        width = text_width, indent = 4)

# make a logging level between debug and info for detailed progress without
# usually useless detail
PROGRESS = 15
_logging.addLevelName(PROGRESS, 'PROGRESS')

# make a filter for sending lower level messages to stdout via baga_cli_to_stdout handler
class UpToFilter(_logging.Filter):
    '''
    provide message level integer up to and including which messages are passed
    
    retains "name" functionality of standard Filter
    '''
    def __init__(self, up_to_maximum, name=""):
        super(UpToFilter, self).__init__(name)
        self.max_level = up_to_maximum
    
    def filter(self, record):
        # non-zero return means we log this message
        if 0 < record.levelno <= self.max_level:
            return True
        else:
            return False

def configureLogger(sample_name, main_log_filename, console_verbosity_lvl, 
        logger_name = '-CLI-'):
    '''
    get and configure a logger for writing to console and a file
    
    Two types of loggers are used:
        (i) the main CLI logger to console and file ('-CLI-')
        (ii) a logger per task potentially used by several modules
    
    If logger_name is not '-CLI-' it should be that of a task e.g. 
    'CollectData'. Some modules are task-specific so in those cases 
    the Task and Module names may be the same. Loggers pass their task
    to Formatters and are included in some of the logs.
    
    Loggers log to stdout and stderr (if level is warning or above).
    
    Loggers also log to files with each Task and each externally called tool 
    having its own timestamped file. Task-level log files (where logger_name is 
    also the name of the Task) are written to the same folder as main_log_filename
    and each path is logged in the main log file.
    
    console_verbosity_lvl acts on console output and corresponds to 
    conventional Python logging level integers. Verbosity to other destinations
    can be tweaked directly on the resulting logger, if desired.
    
    Returns:
        logger: a logging.Logger instance
        log_folder: folder in which log files are written to
    '''
    logger = _logging.getLogger(logger_name)
    # always lowest: filter in handlers etc
    logger.setLevel(_logging.DEBUG)
    # decide filename for log
    if logger_name == '-CLI-':
        this_log_filename = main_log_filename
        this_log_folder = _os.path.sep.join(main_log_filename.split(_os.path.sep)[:-1])
    else:
        # task logger: determine filename to write to
        time_for_task_log_folder = _datetime.now().strftime("%y-%m-%d_%H-%M-%S")
        this_log_folder = _os.path.sep.join(['baga-'+sample_name+'_logs', 
                time_for_task_log_folder+'_'+logger_name])
        try:
            _os.makedirs(this_log_folder)
        except FileExistsError:
            pass
        # with a leading '-' the -main-.log will sometimes sort before the
        # timestamped logs which start with a number. However, this depends on
        # LANG and LC_COLLATE environment settings which cannot be controlled
        # here! LC_COLLATE='en_GB.utf8' ls ; LC_COLLATE='POSIX' ls
        # so just using 00 instead
        this_log_filename = _os.path.sep.join([this_log_folder, '00_main.log'])
        # will need this returned for passing to this_log_filename for each MetaSample.__init__()
        # self.log_folder: per task, formatter to include module name
    
    ## handler to log file
    to_log_file = _logging.FileHandler(this_log_filename, mode = 'a', 
            encoding = 'utf-8')
    to_log_file.setLevel(_logging.DEBUG)
    to_log_file.setFormatter(formatter_tofile)
    logger.addHandler(to_log_file)
    
    if logger_name != '-CLI-':
        # this is a task logger: write filename in main CLI log
        main_logger = _logging.getLogger('-CLI-')
        main_logger.info('Logger for {} will write to {}'.format(logger_name, 
                this_log_filename))
        # also, log INFO to main log file but not progress
        to_main_log_file = _logging.FileHandler(main_log_filename, mode = 'a', 
                encoding = 'utf-8')
        to_main_log_file.setLevel(_logging.INFO)
        # this formatter includes 'task' required in extra kwarg when logging
        to_main_log_file.setFormatter(formatter_tofile_task_to_main)
        logger.addHandler(to_main_log_file)
        
    ## handlers to console: stout and stderr
    to_stdout = _logging.StreamHandler(stream = _sys.stdout)
    to_stderr = _logging.StreamHandler(stream = _sys.stderr)
    if logger_name != '-CLI-':
        # this is a task logger: include Task in output
        to_stdout.setFormatter(wrapping_formatter_stdout_task)
        to_stderr.setFormatter(wrapping_formatter_stderr_task)
    else:
        to_stdout.setFormatter(wrapping_formatter_stdout)
        to_stderr.setFormatter(wrapping_formatter_stderr)
    
    # initial levels to stdout: logging.DEBUG to logging.INFO
    to_stdout.setLevel(_logging.DEBUG)
    to_stdout.addFilter(UpToFilter(_logging.INFO))
    # initial levels to stderr: logging.WARNINGS and logging.ERRORS
    to_stderr.setLevel(_logging.WARNING)
    # update console verbosity to console (stdout+stderr) from CL options
    # verbosity lower than INFO will silence stdout because of a Filter above
    to_stdout.setLevel(console_verbosity_lvl)
    # only update stderr handler if level higher than WARNINGS
    if console_verbosity_lvl > _logging.WARNING:
        to_stderr.setLevel(console_verbosity_lvl)
    
    logger.addHandler(to_stdout)
    logger.addHandler(to_stderr)
    
    if logger_name != '-CLI-':
        # use an adapter to include Task name if necessary
        logger = _logging.LoggerAdapter(logger, {'task': logger_name})
    
    return(logger, this_log_folder)



def bagasave(object_tuple,file_name):
    _pickle.dump(object_tuple, _gzip.open('%s.baga' % file_name,'wb'))

def bagaload(file_name):
    if file_name[-5:] == '.baga':
        return _pickle.load(_gzip.open('%s' % file_name,'rb'))
    else:
        return _pickle.load(_gzip.open('%s.baga' % file_name,'rb'))


def load(file_name):
    '''
    load metadata from previous baga objects into dict
    
    Will reconstruct dicts of arrays if each array was saved as:
        '__<dictname>__<key_as_arrayname>'
    Data type is inferred from start of dictname:
        'sequence*' for unicode
        'ratio*' for float
        all others for integer
    '''
    def loadArray(member_name,contents):
        if 'sequence' in member_name:
            if PY3:
                array_data = _array('b', contents.getvalue())
                # saved as bytes for space saving
                # converted to unicode for test convenience
                # but 'u' is deprecated since 3.3 . . .
                # https://docs.python.org/3/library/array.html
                array_data = _array('u', array_data.tobytes().decode('utf-8'))
            else:
                array_data = _array('c', contents.getvalue())
        elif member_name.startswith('ratio'):
            array_data = _array('f', contents.getvalue())
        else:
            array_data = _array('i', contents.getvalue())
        return(array_data)
    
    meta_data = {}
    arraydicts = {}
    with _tarfile.open(file_name, "r:gz") as tar:
        for member in tar:
            if PY3:
                contents = _BytesIO(tar.extractfile(member).read())
            else:
                contents = _StringIO(tar.extractfile(member).read())
            if member.name[:2] == '__' and sum(['__' == a+member.name[n+1] \
                    for n,a in enumerate(member.name[:-1])]) == 2:
                # if name is '__<dictname>__<arrayname>' build an array dict
                dictname, arrayname = member.name.split('__')[1:]
                array_data = loadArray(dictname, contents)
                try:
                    arraydicts[dictname][arrayname] = array_data
                except KeyError:
                    arraydicts[dictname] = {arrayname:array_data}
            else:
                # else it's an array or JSONable itself
                try:
                    # either json serialised conventional objects
                    contents = _json.loads(contents.getvalue().decode('utf-8'))
                except ValueError:
                    # or longer python array.array objects ### how did this work for integers? Didn't!
                    array_data = loadArray(member.name, contents)
                    if PY3:
                        contents = _array('u', contents.getvalue())
                    else:
                        contents = _array('c', contents.getvalue())
                meta_data[member.name] = contents
        # add dicts of arrays to the rest of meta_data
        for member_name,d in arraydicts.items():
            meta_data[member_name] = d
    return(meta_data)

def save(metadata, file_name):
    '''save metadata dict from previous baga object into an baga file'''
    omissions = []
    with _tarfile.open(file_name, "w:gz") as tar:
        for att_name, att in metadata.items():
            if isinstance(att, _array):
                if PY3:
                    ioob = _BytesIO(att.tostring())
                else:
                    ioob = _StringIO(att.tostring())
                ioob.seek(0, _os.SEEK_END)
                length = ioob.tell()
                ioob.seek(0)
                thisone = _tarfile.TarInfo(name = att_name)
                thisone.size = length
                tar.addfile(tarinfo = thisone, fileobj = ioob)
                #log.debug
                print('Stored in {}: "{}" (array)'.format(file_name, 
                        att_name))
            else:
                # try saving everything else here by jsoning
                try:
                    #log.debug
                    print('Attempting store in {}: "{}", {}'.format(file_name, 
                            att_name, type(att)))
                    if PY3:
                        ioob = _BytesIO()
                    else:
                        ioob = _StringIO()
                    att_json = _json.dumps(att)
                    ioob.write(att_json.encode('utf-8'))
                    ioob.seek(0, _os.SEEK_END)
                    length = ioob.tell()
                    ioob.seek(0)
                    thisone = _tarfile.TarInfo(name = att_name)
                    thisone.size = length
                    tar.addfile(tarinfo = thisone, fileobj = ioob)
                    #log.debug
                    print('Stored in {}: "{}", {}'.format(file_name, 
                            att_name, att))
                except TypeError:
                    #log.debug
                    print('Not JSONable: "{}", {}'.format(att_name, type(att)))
                    # ignore non-jsonable things like functions
                    # include unicodes, strings, lists etc etc
                    omissions += [att_name]
    return(omissions)


def cpu_count():
    try:
        return _multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass
    
    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')
        
        if res > 0:
            return res
    except IOError:
        pass
    
    # Windows
    try:
        res = int(_os.environ['NUMBER_OF_PROCESSORS'])
        
        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

def get_available_memory(info_file = '/proc/meminfo'):
    '''Get free memory (GB) in systems with a /proc/meminfo'''
    try:
        meminfo = open(info_file).readlines()
    except IOError:
        print('Could not open file for collecting memory status: {}'.format(info_file))
        return(None)
    
    total = -1
    free = -1
    buffers = -1
    cache = -1
    for line in meminfo:
        if 'MemTotal' in line:
            total = float(line.split(' ')[-2])/1048576
        elif 'MemFree' in line:
            free = float(line.split(' ')[-2])/1048576
        elif 'Buffers' in line:
            buffers = float(line.split(' ')[-2])/1048576
        elif 'Cached' in line:
            cache = float(line.split(' ')[-2])/1048576
    
    e = 'Problems parsing {}:\n{}\n'.format(info_file, ''.join(meminfo))
    assert -1 not in (total, free, buffers, cache), e
    
    available = free + buffers + cache
    print('Memory available: {:.2f} GB, ({:.2f} + {:.2f} + {:.2f}) of {:.2f} GB'.format(available, free, buffers, cache, total))
    return(available)


def decide_max_processes(max_cpus):
    '''
    decide how many CPUs to use
    basically one less than maximum if not specified
    maximum of more than available specified
    else the amount specified
    '''
    total_cpus = cpu_count()
    if max_cpus < 0:
        # negative number is number less than total on system
        # default is one less
        max_processes = total_cpus + max_cpus
        
        if max_processes <= 0:
            print('WARNING: max_cpus is set too low ({}), total CPUs is {}, so \
                   will us a single CPU'.format(max_cpus, total_cpus))
            
            max_processes = 1
        
    elif max_cpus > total_cpus:
        print('WARNING: max_cpus is set too high ({}), total CPUs is {}, so \
               will us all available CPUs'.format(max_cpus, total_cpus))
        
        max_processes = total_cpus
        
    else:
        max_processes = max_cpus
    
    return( max_processes )

def report_time(start_time, this_item_index, total_items):
    '''
    print time elapsed to nearest minute or second
    print an estimate for time remaining given item number and total items
    '''
    interval = _time.time() - start_time
    remaining =   (interval / float(this_item_index + 1)) \
                  * (total_items - (this_item_index + 1))
    
    def get_unit(interval_seconds):
        if interval_seconds < 60:
            n = '{:.0f}'.format(interval_seconds)
            u = 'second'
        else:
            n = '{:.0f}'.format(interval_seconds / 60)
            u = 'minute'
        
        if n == '1':
            p = ''
        else:
            p = 's'
        
        return('{} {}{}'.format(n,u,p))
    
    print('\n{} elapsed after {} of {}'.format(get_unit(interval), this_item_index+1, total_items))
    print('Estimated time left: {}\n'.format(get_unit(remaining)))


def get_exe_path(name):
    '''
    Use the Dependencies module to get the path to an executable.
    Allows keeping that information in one place.
    If returned value is a string, it is the executable path;
    if returned value is a list, it is the executable path with path
    to script.
    Therefore, output needs to be checked and handled appropriately 
    with subprocess methods.
    '''
    from baga.Dependencies import dependencies as _dependencies
    path_to_exe = []
    if _dependencies[name]['checker']['arguments']['path'][0][0] != '/' or \
       len(_dependencies[name]['checker']['arguments']['path'][0]) == 0 :
        # not absolute path e.g. for user install of cutadapt
        # so include absolute path to baga local installs for exe
        path_to_exe += [_dependencies[name]['destination']]
    
    path_to_exe += _dependencies[name]['checker']['arguments']['path']
    
    path_to_exe = _os.path.sep.join(path_to_exe)
    
    try:
        # will return list because executable + script
        # would need to be checked and handled appropriately with subprocess methods
        pre = _os.path.sep.join(_dependencies[name]['checker']['arguments']['extra_pre'])
        path_to_exe = [pre,path_to_exe]
    except KeyError:
        # nothing to prepend e.g., python2 for spades
        # i.e. python2 when python3.3 is last python3 supported but >3.3 installed
        pass
    
    return(path_to_exe)


def get_jar_path(name):
    '''
    Use the Dependencies module to get the path to a Java jar file.
    Allows keeping that information in one place.
    '''
    from baga.Dependencies import dependencies as _dependencies
    return(_dependencies[name]['checker']['arguments']['java_commands'][2])

def parse_new_arguments(arguments, cmd):
    '''
    Update cmd dict from string of arguments to pass to wrapped program
    
    Handles -/--, flags, quotes and multiple arguments per option
    '''
    lexer = _shlex.shlex(arguments)
    # prevent splitting floating point arguments
    lexer.wordchars += '.'
    new_arguments = {}
    collect = False
    pre = []
    for token in lexer:
        if token == '-':
            collect = False
            pre += [token]
        else:
            if not collect:
                k = ''.join(pre + [token])
                new_arguments[k] = None
                collect = True
                pre = []
            else:
                try:
                    new_arguments[k] += [token]
                except TypeError:
                    new_arguments[k] = [token]
    
    for opt,arg in new_arguments.items():
        if arg is not None:
            new = opt + ' ' + ' '.join(arg)
            # remove ""
            use_arg = [a.strip('"') for a in arg]
        else:
            new = opt
            use_arg = arg
        try:
            if cmd[opt] is not None:
                old = opt + ' ' + ' '.join(cmd[opt])
            else:
                old = opt
            if old != new:
                print('Using option: {} (instead of: {})'.format(new, old))
                cmd[opt] = [a.strip() for a in arg]
        except KeyError:
            print('Using option: {}'.format(new))
            cmd[opt] = use_arg
    
    return(cmd)


if PY3:
    use_class = FileNotFoundError
else:
    use_class = IOError

class InheritanceFileNotFoundError(use_class):
    '''
    Raised when a required metadata baga file is not found.
    
    These are generated by prior steps in the baga pipeline.
    '''

class NCBItaxonomyFileNotFoundError(use_class):
    '''
    Raised when the pre-parsed NCBItaxonomy file is not found.
    '''


class Base(object):
    '''
    Strange Base class to "absorb" arguments
    
    see:
    https://stackoverflow.com/questions/8972866/correct-way-to-use-super-argument-passing
    http://python-history.blogspot.co.uk/2010/06/inside-story-on-new-style-classes.html
    http://fuhm.net/super-harmful/
    http://rhettinger.wordpress.com/2011/05/26/super-considered-super/
    '''
    def __init__(self, *args, **kwargs): pass


class MetaSample(Base):
    '''
    MetaSample class is inherited in each baga module where different methods are added

    Common properties are __init__ where loading from prior analysis is dealt with 
    and logging for wrappers is set up and saveLocal() to save (meta)data to a file. Child
    classes include Genome and Reads in  the CollectData module.
    '''


    def __init__(self,  sample_name,
                        module_name,
                        task_name = False,
                        console_verbosity_lvl = False,
                        log_folder = False,
                        analysis_path = '.',
                        skip_analysis_path = False,
                        inherit_from = False,
                        omit_from_inheritance = [],
                        *args, **kwargs):
        '''
        Providing only sample name will instantiate a new object.
        
        Optionally providing "inherit_from = <object class>" will load metadata 
        from a saved object. If <object class> is same as this one (or for 
        convenience, 'self'), an attempt will be made to reload a previously 
        saved version of this object for this sample_name to resume analysis at 
        this stage of a pipeline.
        
        Usually, <object class> is from a previous, upstream, analysis step in a 
        pipeline.
        
        skip_analysis_path is for subclasses that don't need an analysis path 
        created e.g. CollectData.Taxonomy
        '''
        super(MetaSample, self).__init__(*args, **kwargs)
        # log what happens during .__init__()
        if task_name:
            # if this has been called via the CLI it is part of a task
            logger = _logging.getLogger(task_name)
            # need to add adaptor after getting (handlers etc are retained from configureLogger)
            # in this namespace logger is the adaptor but not beyond
            logger = _logging.LoggerAdapter(logger, {'task': task_name})
        else:
            # else conform to conventional Python logging and use module name
            logger = _logging.getLogger(module_name)
        
        self.logger = logger
        self.sample_name = sample_name
        self.console_verbosity_lvl = console_verbosity_lvl
        
        ## inheritance of metadata ##
        # will make own .i attribute of inherited metadata
        # don't inherit these universally
        inheritance_filter = set([
                'externals_for_logging',
                'log_for',
                'file_name',
                'log_folder',
                'i',
                ])
        
        # don't inherit these in a particular subclass
        inheritance_filter.update(omit_from_inheritance)
        
        if inherit_from == 'self':
            # restore from same kind of object
            # not tested task name i.e., class using type(self).__name__
            file_name = '{}.{}-{}.baga'.format(module_name, type(self).__name__, self.sample_name)
            path = _os.path.sep.join([analysis_path, file_name])
            logger.debug('__init__() loading from: {}'.format(path))
            metadata = load(path)
            for name,content in metadata.items():
                if name == 'sample_name':
                    if content != self.sample_name:
                        raise RuntimeError("mismatch between requested sample "\
                                "name ({}) and '{}' in local saved object: {}. "\
                                "Delete or rename it to continue".format(
                                self.sample_name, content, path))
                    else:
                        continue
                if name in ('externals_for_logging','log_for'):
                    # do not restore: set up below
                    # do restore: log_folder
                    continue
                # set as attributes
                setattr(self, name, content)
                logger.debug('__init__() loaded: {} ({})'.format(name, 
                        type(content)))
        elif isinstance(inherit_from, list):
            i = {}
            for Module_Datatype,Class_SubtaskTool in inherit_from:
                # module name describes data type
                # class name describes subtask and software tool used
                # CLI-based pipelines consist of tasks, tasks consist of subtasks
                # retrieve (meta)data from objects produced by upstream subtask(s)
                file_name = 'baga.{}.{}-{}.baga'.format(Module_Datatype, 
                        Class_SubtaskTool, self.sample_name)
                logger.debug('__init__() inheriting from: {}'.format(file_name))
                try:
                    metadata = load(file_name)
                except FileNotFoundError as e:
                    raise InheritanceFileNotFoundError(e)
                if len(metadata):
                    if Module_Datatype not in i:
                        i[Module_Datatype] = {}
                    i[Module_Datatype][Class_SubtaskTool] = {}
                    # how to log this from a specific module if defined in superclass?
                    # print('Inheriting from: {}'.format(file_name))
                    # store inheritance in dictionaries
                    for name,content in metadata.items():
                        if name == 'sample_name':
                            if content != self.sample_name:
                                raise RuntimeError("mismatch between requested "\
                                        "sample name ({}) and '{}' in local "\
                                        "saved object: {}. Delete or rename it "\
                                        "to continue".format(self.sample_name, 
                                        content, file_name))
                            else:
                                continue
                        if name in inheritance_filter:
                            # don't inherit module specific operational stuff
                            # and the previous inheritance i.e., 'i'
                            # (could change this later?)
                            continue
                        i[Module_Datatype][Class_SubtaskTool][name] = content
                        logger.debug('__init__() inherited: {}'.format(name))
            self.i = i
            
            inherited_paths = set()
            if inherit_from:
                # establish path for output: check inheritance
                # applies for 'self' or 'upstream' metadata
                for module, tasks in self.i.items():
                    for task, inheritance in tasks.items():
                        if 'analysis_path' in inheritance:
                            inherited_paths.add(inheritance['analysis_path'])
            
            ## determine analysis path ##
            if skip_analysis_path:
                ### I think this is redundant
                # e.g. CollectData --taxonomy doesn't need an analysis path created
                logger.debug('__init__() not creating analysis path (e.g. for CollectData --taxonomy)')
            else:
                logger.debug('__init__() determining analysis path for object in '\
                        'module: {}'.format(module_name))
                if len(inherited_paths) > 1:
                    # could have a --new option?
                    raise NotImplemented("inherited (meta)data objects from more than one "\
                            "previous pipline. Analyses were in: {}".format(
                            ' and '.join(sorted(inherited_paths))))
                elif len(inherited_paths) == 1:
                    # use the previously found path
                    self.analysis_path = inherited_paths.pop()
                    logger.debug('__init__() inherited data folder: {}'.format(self.analysis_path))
                    logger.info('Not using provided path for analysis output ({}) '\
                            'because inherited object provided one: {}'.format(analysis_path, 
                            self.analysis_path))
        else:
            # no inheritance/loading of metadata no previous output:
            # start a new one
            self.analysis_path = _os.path.sep.join([analysis_path,
                    'baga-{}_{}'.format(self.sample_name, 
                    _datetime.now().strftime("%y-%m-%d_%H-%M-%S"))])
            # if being run from baga_cli.py, this logger came via configureLogger() which
            # would have been passed console_verbosity_lvl, which would have been applied to
            # stdout and stderr handlers . . . .
            logger.debug('__init__() new data folder: {}'.format(self.analysis_path))
            
        try:
            if PY3:
                _os.makedirs(self.analysis_path, exist_ok = True)
            else:
                if not _os.path.exists(self.analysis_path):
                    _os.makedirs(self.analysis_path)
        except OSError:
            raise OSError('Could not make output path at: {}'.format(
                    self.analysis_path))
        
        ## determine file name to eventually save metadata to ##
        # (and for monitoring which files were inherited from if info later
        # restored from file)
        # if inheriting from self, will be the same as before
        self.file_name = '{}.{}-{}.baga'.format(module_name, type(self).__name__, 
                self.sample_name)
        logger.debug('__init__() will save (meta)data to: {}'.format(self.file_name))
        
        ## determine logger name ##
        # for logging this module
        # (name only known at inherited instantiation so __name__ will no work)
        self.module_name = module_name
        if task_name:
            # if called from the CLI in the context of a Task (subparser)
            self.logger_name = task_name
        else:
            self.logger_name = module_name
        
        if log_folder:
            # task level logging determines folder in baga_cli.py
            self.log_folder = log_folder
          

    def saveLocal(self, exclude = []):
        '''
        Save an baga object to a file
        
        Objects usually contain sample metadata information. The data is 
        saved as a compressed archive with attributes as JSON or raw Python 
        array.array objects as bytes or strings. Dictionaries of arrays are
        deconstructed and saved as array.array bytes and can be reconstrcuted
        by baga metadata inheritance in MetaSample's __init__() method.
        '''
        def saveArray(array_data, member_name, tar):
            if array_data.typecode == 'u':
                # 1 byte per character instead of 2 (4?)
                # 'u' Deprecated since version 3.3, will be removed in version 4.0
                # https://docs.python.org/3/library/array.html
                # so how to conveniently and efficiently put a string into an array?
                # implement own fetch methods for extracting region or writing to text file?
                # better way of doing this?
                array_data = _array('b', bytes(''.join(array_data.tolist()), encoding = 'utf-8'))
                #array_data = _array('b', array_data.encode('utf-8'))
            if PY3:
                ioob = _BytesIO(array_data.tobytes())
            else:
                ioob = _StringIO(array_data.tostring())
            ioob.seek(0, _os.SEEK_END)
            length = ioob.tell()
            ioob.seek(0)
            thisone = _tarfile.TarInfo(name = member_name)
            thisone.size = length
            tar.addfile(tarinfo = thisone, fileobj = ioob)
        
        omissions = []
        with _tarfile.open(self.file_name, "w:gz", compresslevel = 3) as tar:
            for att_name, att in self.__dict__.items():
                if att_name in exclude:
                    omissions += [att_name]
                    self.logger.debug('Excluding from {}: "{}" ({})'.format(
                            self.file_name, att_name, type(att)))
                if isinstance(att, _array):
                    saveArray(att, att_name, tar)
                    self.logger.debug('Stored in {}: "{}" (array)'.format(
                            self.file_name, att_name))
                elif isinstance(att, dict) and \
                        all([isinstance(v, _array) for v in att.values()]):
                    # can save dicts of arrays efficiently as array.tobytes()
                    # .tostring() in python2 which is what .tofile() does
                    for array_name,array_data in att.items():
                        member_name = '__{}__{}'.format(att_name,array_name)
                        saveArray(array_data, member_name, tar)
                        self.logger.debug('Stored in {}: "{}" (array) in {} '\
                                '(dict) as {}'.format(self.file_name, 
                                array_name, att_name, member_name))
                else:
                    # try saving everything else here by jsoning
                    try:
                        if PY3:
                            ioob = _BytesIO()
                        else:
                            ioob = _StringIO()
                        att_json = _json.dumps(att)
                        ioob.write(att_json.encode('utf-8'))
                        ioob.seek(0, _os.SEEK_END)
                        length = ioob.tell()
                        ioob.seek(0)
                        thisone = _tarfile.TarInfo(name = att_name)
                        thisone.size = length
                        tar.addfile(tarinfo = thisone, fileobj = ioob)
                        self.logger.debug('Stored in {}: "{}", {}'.format(
                                self.file_name, att_name, type(att)))
                    except TypeError:
                        # could be warning but expect functions to fail here
                        self.logger.debug('Not JSONable: "{}", {}'.format(
                                att_name, type(att)))
                        # ignore non-jsonable things like functions
                        # include unicodes, strings, lists etc etc
                        omissions += [att_name]
        return(omissions)


    def launch_external(self, cmd, 
                              external_program_name, 
                              log_folder, 
                              main_logger = False, 
                              stdout_filename = False,
                              shell = False):
        '''
        Using a subprocess.Popen object, launch a program, collect output 
        and log useful information like pid, issuing command, duration etc.
        If main_logger provided, which is usually a Task-level logger (not CLI), 
        key info but not details are also logged there.
        '''
        if not _os.path.exists(log_folder):
            _os.makedirs(log_folder)
        
        # set up handler for this subprocess
        # send to a unique time-stamped log file
        # also add a handler to the main log (at debug level?)
        time_start = _datetime.now()
        time_stamp = time_start.strftime("%d_%H-%M-%S-%f")
        filename = _os.path.join(log_folder, time_stamp+'__'+external_program_name+'.log')
        this_ext_file_handler = _logging.FileHandler(filename, mode='w', encoding='utf-8')
        this_ext_file_handler.setLevel(_logging.DEBUG)
        logger_name = self.module_name+'_{}'.format(external_program_name)
        logger_for_external = _logging.getLogger(logger_name)
        logger_for_external.setLevel(_logging.DEBUG)
        logger_for_external.addHandler(this_ext_file_handler)
        
        # log issuing command
        command_str = 'Command: "{}"'.format(' '.join(cmd))
        if main_logger:
            main_logger.log(PROGRESS, command_str)
        logger_for_external.log(PROGRESS, command_str)
        time_for_log = time_start.strftime("%H:%M:%S on %d/%m/%y (%a)")
        
        if stdout_filename:
            # for data generated and sent to stdout
            logger_for_external.info('Writing stdout to file: {}'.format(stdout_filename))
            with open(stdout_filename, 'wb') as for_stdout:
                proc = _subprocess.Popen(cmd, stdout = for_stdout, 
                        stderr = _subprocess.PIPE, shell = shell)
                # wait for completion
                starting = 'Started {} ({}) at {}'.format(external_program_name, proc.pid, time_for_log)
                logger_for_external.info(starting)
                if main_logger:
                    # info level also goes to CLI
                    main_logger.info(starting)
                    main_logger.info('Writing log for {} (STDERR, etc) to: {}'.format(
                            external_program_name, filename))
                logger_for_external.info('Writing STDOUT for {} ({}) to {}'.format(
                        external_program_name, proc.pid, stdout_filename))
                stdout_not_used, stderr = proc.communicate()
                #stderr = proc.stderr.read()
            
        else:
            proc = _subprocess.Popen(cmd, stdout = _subprocess.PIPE, 
                    stderr = _subprocess.PIPE, shell = shell)
            logger_for_external.log(PROGRESS, 'Will include stdout in this log')
            starting = 'Started {} ({}) at {}'.format(external_program_name, proc.pid, time_for_log)
            logger_for_external.info(starting)
            if main_logger:
                main_logger.info(starting)
                main_logger.info('Writing log for {} (STDERR, etc) to: {}'.format(
                        external_program_name, filename))
            # wait for completion
            stdout, stderr = proc.communicate()
        
        # work out how long it took
        time_end = _datetime.now()
        time_for_log = time_end.strftime("%H:%M:%S on %d/%m/%y (%a)")
        finished_at = 'Finished {} ({}) at {}'.format(external_program_name, proc.pid, 
                time_for_log)
        logger_for_external.info(finished_at)
        return_code = 'Return code for {} ({}) was: {}'.format(external_program_name, 
                proc.pid, proc.returncode)
        
        if proc.returncode == 0:
            logger_for_external.info(return_code)
            # some programs (e.g. LMAT) do not produce informative return codes
            if not stdout_filename:
                if 'error' in stdout.decode('utf-8').lower():
                    main_logger.warning('Despite return code 0, '\
                            '"Error" string detected in STDOUT: check logs, '\
                            'including thoses generated directly by {}'\
                            ''.format(external_program_name))
            if 'error' in stderr.decode('utf-8').lower():
                main_logger.warning('Despite return code 0, '\
                        '"Error" string detected in STDERR: check logs, '\
                        'including thoses generated directly by {}'\
                        ''.format(external_program_name))
        else:
            logger_for_external.warning(return_code)
        
        if main_logger:
            main_logger.log(PROGRESS, finished_at)
            if proc.returncode == 0:
                main_logger.log(PROGRESS, return_code)
            else:
                main_logger.warning('STDERR for {} ({}) was:'.format(external_program_name, 
                        proc.pid))
                for line in stderr.decode('utf-8').split('\n'):
                    main_logger.warning('| '+line)
                if not stdout_filename:
                    main_logger.warning('STDOUT for {} ({}) was:'.format(external_program_name, 
                            proc.pid))
                    for line in stdout.decode('utf-8').split('\n'):
                        main_logger.warning('| '+line)
        
        run_time = time_end - time_start
        seconds = int(run_time.total_seconds())
        hours, remainder = divmod(seconds,60*60)
        minutes, seconds = divmod(remainder,60)
        ran_for = 'Run time for {} ({}): {} hours, {} minutes, {} seconds'.format(
                external_program_name, proc.pid, hours, minutes, seconds)
        logger_for_external.log(PROGRESS, ran_for)
        if main_logger:
            main_logger.log(PROGRESS, ran_for)
        
        if stdout_filename:
            logger_for_external.info('STDOUT for {} ({}) was written '\
                    'to {} ({:,} bytes)'.format(external_program_name, proc.pid, 
                    stdout_filename, _os.path.getsize(stdout_filename)))
        else:
            # only write to log if not writing to file
            logger_for_external.log(PROGRESS, 'STDOUT for {} ({}) was:'.format(
                    external_program_name, proc.pid))
            for line in stdout.decode('utf-8').split('\n'):
                logger_for_external.log(PROGRESS, '| '+line)
        
        logger_for_external.log(PROGRESS, 'STDERR for {} ({}) was:'.format(
                external_program_name, proc.pid))
        for line in stderr.decode('utf-8').split('\n'):
            logger_for_external.log(PROGRESS, '| '+line)
        
        # remove handler (logger is per application, not per subprocess, but log 
        # files are per process)
        logger_for_external.removeHandler(this_ext_file_handler)
        
        return(proc.returncode)
