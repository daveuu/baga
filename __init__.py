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
# Dr Michael A Brockhurst (University of York, UK)
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
    return(path_to_exe)

def get_jar_path(name):
    '''
    Use the Dependencies module to get the path to a Java jar file.
    Allows keeping that information in one place.
    '''
    from baga.Dependencies import dependencies as _dependencies
    return(_dependencies[name]['checker']['arguments']['java_commands'][2])
