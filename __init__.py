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

import cPickle as _cPickle
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
import shlex as _shlex
from cStringIO import StringIO as _StringIO
from array import array as _array


# update this for objects that are mostly long lists of numbers:
    # open a gzipped tar file and add firstly json'd list of info like length, resolution
    # then array.tofile()
    # load using array.fromfile() after getting info out of first file

def bagasave(object_tuple,file_name):
    _cPickle.dump(object_tuple, _gzip.open('%s.baga' % file_name,'wb'))

def bagaload(file_name):
    if file_name[-5:] == '.baga':
        return _cPickle.load(_gzip.open('%s' % file_name,'rb'))
    else:
        return _cPickle.load(_gzip.open('%s.baga' % file_name,'rb'))

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
