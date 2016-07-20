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
Commandline interface to the Bacterial and Archaeal Genome Analyser 
(baga) Python package.
'''
import argparse
import textwrap
import sys
import re
import os
import subprocess
import logging
from datetime import datetime
from collections import defaultdict

# needed by logging bits in baga
try:
    # get width for wrapping
    info = subprocess.check_output(['stty','size'])
    rows, text_width = map(int,info.rstrip().split())
except FileNotFoundError:
    # else pick a narrow one
    rows, text_width = 30, 70

from baga import Dependencies
from baga import configureLogger
from baga import InheritanceFileNotFoundError
from baga import NCBItaxonomyFileNotFoundError

def check_files(files):
    have_these = set()
    if isinstance(files, dict):
        # dictionary of paired fastqs
        for run_acc, thesefiles in files.items():
            have = 0
            for n,f in thesefiles.items():
                try:
                    with open(f, 'rb') as filein:
                        if os.path.getsize(f):
                            # confirm it isn't an empty file
                            have += 1
                except IOError:
                    pass
            
            if have == len(thesefiles):
                have_these.add(run_acc)
    else:
        # list of SAMs or BAMs
        for alnfile in files:
            try:
                with open(alnfile, 'rb') as filein:
                    if os.path.getsize(alnfile):
                        # confirm it isn't an empty file
                        have_these.add(alnfile)
            except IOError:
                pass
            
    
    if have_these == set(files):
        return(True)
    else:
        return(False)


def delete_files(files, extra = False):
    total_size = 0
    if isinstance(files, dict):
        for run_acc, thesefiles in files.items():
            for n,f in thesefiles.items():
                total_size += os.path.getsize(f)
                print('Deleting {}'.format(f))
                os.unlink(f)
                if extra:
                    if extra in f:
                        presubsampled = f.replace(extra,'')
                        total_size += os.path.getsize(presubsampled)
                        print('Deleting {}'.format(presubsampled))
                        os.unlink(presubsampled)
    else:
        for alnfile in files:
            total_size += os.path.getsize(alnfile)
            print('Deleting {}'.format(alnfile))
            os.unlink(alnfile)
    
    return(total_size)

from baga import get_exe_path
def exe_fail(program_name):
    print('Could not find the {} executable at executable at {}.'.format(
                        program_name, get_exe_path(program_name)))
    
    print('You can check if it is installed using:')
    print('{} Dependencies --check {}'.format(sys.argv[0], program_name))
    print('You can install it locally using:')
    print('{} Dependencies --get {}'.format(sys.argv[0], program_name))
    sys.exit(1)

def check_baga_path(baga_prefix, supplied_name):
    '''
    to allow for either path and filename or just name of baga objects, 
    do some checks to see what was supplied at commandline. Return both
    the name and the path to the baga file
    '''
    if baga_prefix[-1] == '-':
        baga_prefix = baga_prefix[:-1]
    # try as path
    try:
        use_path = os.path.expanduser(supplied_name)
        use_name = use_path.split(os.path.sep)[-1].replace(baga_prefix+'-','').replace('.baga','')
        with open(use_path,'rb') as f:
            return(use_path,use_name)
    except IOError:
        pass
    
    # try as baga name
    try:
        use_path = '{}-{}.baga'.format(baga_prefix, supplied_name)
        print(use_path)
        with open(use_path,'rb') as f:
            return(use_path,supplied_name)
    except IOError:
        pass
    
    return(False,False)

def check_direct_arguments(arguments, wrapped_tools):
    '''
    Parse commandline arguments to be passed to wrapped tools.
    
    Checks whether specified target programs are recognised among task-specific 
    "wrapped_tools".
    '''
    if arguments:
        # to avoid conflicts on the command line, these direct arguments to 
        # wrapped tools need to use underscores instead of dashes, so replace
        # former with latter here
        def underscores(txt):
            '''convert __ to -- and _ to -'''
            a = re.sub('(_)(_)', r'-\1', txt)
            a = re.sub('([- ]|^)(_)([\w])', r'\1-\3', a)
            return(a)
        direct_arguments = {a:underscores(b) for a,b in zip(arguments[::2],
                                               arguments[1::2])}
        assert len(direct_arguments) != 0, 'Supplied direct arguments ({}) '\
                'not recognised. Should be: "--arguments <program name> '\
                '<\'_o 5 __option 2\'>" (use underscore, _o and __option, '\
                'instead of hyphen, -o and --option)'
        recognised = set(direct_arguments) & set(wrapped_tools)
        if len(recognised) == 0:
            recognised = ['none']
        not_recognised = set(direct_arguments) - set(wrapped_tools)
        print(direct_arguments, not_recognised, len(not_recognised))
        assert len(not_recognised) == 0, 'Some of the wrapped '\
                'program names specified are not recognised:\n'\
                'Recognised: {}\n'\
                'Not recognised: {}\n'.format(', '.join(recognised),
                                              ', '.join(not_recognised))
        return(direct_arguments)
    else:
        # CLI defaults to False for empty
        return({})
text_width = 70

title = 'Bacterial and Archaeal Genome Analyser'
subtitle = textwrap.fill('Novel analyses and wrapped tools pipelined for convenient processing of genome sequences', text_width)
version_date = 'December 20 2015'
version_num = 0.2
authors = 'David Williams'
email = 'david.williams.at.liv.d-dub.org.uk'
blurb = '''Work on this software was started at The University of Liverpool, UK 
with funding from The Wellcome Trust (093306/Z/10) awarded to:

Dr Steve Paterson (The University of Liverpool, UK)
Dr Craig Winstanley (The University of Liverpool, UK)
Dr Michael A Brockhurst (The University of York, UK)


Copyright (C) 2015 David Williams
License GPLv3+: GNU GPL version 3 or later
This is free software: you are free to change and redistribute it
There is NO WARRANTY, to the extent permitted by law

'''
splash = '\n{title}:\n\n{subtitle}\n\nVersion {version_num} ({version_date})\n\n{authors}\n{email}\n\n\n{blurb}\n'.format(
        title = title,
        subtitle = subtitle,
        version_date = version_date,
        version_num = version_num,
        authors = authors,
        email = email,
        blurb = blurb)


## get options from command line and decide what to do

parser = argparse.ArgumentParser(
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('{} {}:\n{}'.format(
                                                title, 
                                                version_num, 
                                                subtitle),
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('The {} (baga) will perform all necessary downloads and analysis in the current directory. Please ensure there is sufficient disk space for downloads of short read data (many large fastq files) and subsequent analyses. Third party programs can either be downloaded and compiled by baga (except git, make and GATK) or installed system wide before-hand by the user or installed locally with paths to executables supplied in the dependencies.conf file.\n\nExample usage: %(prog)s CollectData ...\n'.format(title), text_width, replace_whitespace = False))


parser.add_argument('--version', action='version', 
                    version='%(prog)s {} ({})'.format(version_num, version_date))

group = parser.add_mutually_exclusive_group(required=False)

#parser.add_argument('--nosplash',
group.add_argument('--nosplash',
    help = "Supress printing of start-up splash info.",
    action = 'store_true',
    default = False)

group.add_argument('--splash',
    help = "Print start-up splash info.",
    action = 'store_true',
    default = True)

parser.add_argument('-v', '--verbosity', 
    type = int, 
    help = "set verbosity output to console: 0 for silent except errors, 1 to "\
            "include warnings, 2 to include updates on progress, 3 to report "\
            "general progress, 4 to include lots of detail (maximum verbosity "\
            "for debugging). More detail is always provided in log files.",
    default = 2,
    choices = [0,1,2,3,4])

subparser_adder = parser.add_subparsers(title = 'Analyses', dest="subparser")

# arguments common to all modules (except Dependencies)
common_module_arguments = subparser_adder.add_parser('common', add_help=False)

common_module_arguments.add_argument('-F', "--force", 
    help = "even if output files present before an analysis, repeat analysis "\
    "and overwrite them. Default behaviour is to continue to next step with "\
    "repeating. If e.g., repeating with different parameters, you may want "\
    "'--force' to overwrite with new output.",
    action = 'store_true')


common_module_arguments.add_argument('-c', "--max_cpus", 
    help = "maximum number of cpus (or threads) to use where parallelisation "\
            "possible. Negative values are less than total CPUs. Default is -1 "\
            "i.e., one less than total.",
    type = int,
    default = -1)

### not sure which are universal in baga yet

# arguments required for all data processing modules (--sample_name)
# (optional for CollectData so specified separately)
sample_name_argument = subparser_adder.add_parser('sample_name_argument', 
        add_help=False)
sample_name_argument.add_argument('-n', "--reads_name", 
    help = "name of sample for which data will be loaded (originally used in a "\
            "CollectData command). This is required.",
    required = True)

####### single letters inconsistant with CallVariants <==

# consider this as CollectData only?
# or retain option to alter --analysis_path? <== not yet implemented
# usually inherited from upstream.
# common_module_arguments.add_argument('-p', "--analysis_path", 
    # help = "name of path in which to write input data. Ideally fast storage "\
    # "such as a hard disk installed in the device doing the computations. Avoid "\
    # "network storage (because it is slower), solid state disks are typically "\
    # "fastest. Defaults to current working directory",
    # type = str,
    # default = '.',
    # metavar = 'PATH')


parser_Dependencies = subparser_adder.add_parser('Dependencies',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Check required external programs \
are available and optionally get them.',
                                            text_width,
                                            replace_whitespace = False))
                # epilog = textwrap.fill('Groups of read data sets from the \
# CollectData option are loaded by providing the file name.\n\n\
# Example usage: "%(prog)s PrepareReads --reads_name ERR953490plus5others --adaptors --trim --align FM209186"\n'\
# .format(title), text_width, replace_whitespace = False))

group_check_or_get = parser_Dependencies.add_mutually_exclusive_group()

group_check_or_get.add_argument('-c', "--check", 
    help = "check BAG Analyser has access to a dependency either in the system path or locally",
    type = str,
    choices = sorted(Dependencies.dependencies),
    nargs = '+')

group_check_or_get.add_argument('-p', "--checkpackage", 
    help = "check BAG Analyser has access to a python package dependency. Sometimes required if baga just installed a package and a fresh Python instance is required to check it",
    type = str,
    #choices = sorted(Dependencies.dependencies),
    choices = sorted([name for name,info in Dependencies.dependencies.items() if \
                info['checker']['function'] == Dependencies.check_python_package]))

group_check_or_get.add_argument('-g', "--get", 
    help = "get (or explain how to get) a dependency for BAG Analyser",
    type = str,
    choices = sorted(Dependencies.dependencies),
    nargs = '+')

group_check_or_get.add_argument('-C', "--checkget", 
    help = "check a dependency for BAG Analyser and attempt to get if not available",
    type = str,
    choices = sorted(Dependencies.dependencies),
    nargs = '+')

group_check_or_get.add_argument('-f', "--checkgetfor", 
    help = "check a set of dependencies for a BAG Analyser task and attempt to get those that are not available",
    type = str,
    choices = sorted(Dependencies.dependencies_by_task),
    nargs = '+')

parser_Dependencies.add_argument("-V", "--versions_file", 
    help="specify file containing versions of software dependencies to use. Defaults to versions.yaml in current folder, falls back to versions.yaml in same folder as baga_cli.py (which might be the same one). If no yaml file specified or found, will use a set of versions built into Dependencies.py that were current in late 2015.", 
    type = str,
    default = 'versions.yaml')


parser_CollectData = subparser_adder.add_parser('CollectData',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        parents = [common_module_arguments],
        description = textwrap.fill('Download and parse genomes or download '\
        'short reads for analysis by {}'.format(title), text_width, 
        replace_whitespace = False),
        epilog = textwrap.fill('Genomes can be loaded from genbank files with '\
        'user provided paths or downloaded from the National Center for '\
        'Biotechnology Information (NCBI) with Assembly database search terms, '\
        'usually accession numbers or strain, species or genus names. '\
        'Short reads can be downloaded from the European Nucleotide Archive '\
        'with user provided Run Accession numbers\n\nExample usage: '\
        '"%(prog)s -r ERR953490,ERR953494,ERR953501"\n'.format(title), 
        text_width, replace_whitespace = False))

mutually_exclusive_group = parser_CollectData.add_mutually_exclusive_group()

# could use nargs='+' or similar here to get a list directly,
# but not sure about spaces in paths
mutually_exclusive_group.add_argument("-g", "--genomes", 
    help="(down)load and parse genomes for analysis. Requires one or more "\
    "filename paths ending .gbk or .gbff or search terms to query the NCBI "\
    "Assembly database. Alternatively the path to a single text file with an "\
    "Assembly accession or other single-result search term per line (that will "\
    "also be used as sample_name and used in a filename). A typical file is "\
    "generated after downloading from Entrez to help with reproducibility (more "\
    "general search terms can return inconsistent results over time). All results "\
    "will be downloaded unless exclusions specified e.g., --complete_only. Genomes "\
    "are saved to an fast internal format unless --keep_gbk is specified in which "\
    "case the full .gbk files are also saved to disk. \n E.g., '--genomes "\
    "Campylobacter --complete_only' will download all complete bacterial genomes "\
    "of the genus 'Campylobacter'",
    type = str,
    nargs = '+')

parser_CollectData.add_argument('-G', "--keep_gbk", 
    help = "retain the original genbank file from NCBI and save to disk in "\
    "addition to the fast internal format",
    action = 'store_true')

parser_CollectData.add_argument('-S', "--refseq_only", 
    help = "only download if sequence is part of the RefSeq collection, "\
    "ignore if in GenBank only",
    action = 'store_true')

assembly_status_filters = parser_CollectData.add_mutually_exclusive_group()

assembly_status_filters.add_argument('-C', "--complete_only", 
    help = "only download complete genome sequences",
    action = 'store_true')

assembly_status_filters.add_argument('-D', "--draft_only", 
    help = "only download partial genome sequences in contigs or scaffolds",
    action = 'store_true')



mutually_exclusive_group.add_argument("-t", "--taxonomy", 
    help="download and update local version of the NCBI taxonomy. "\
    "An alternative URL can be provided but should not be necessary.", 
    type = str,
    nargs = '?',
    default = None,
    const = 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz',
    metavar = 'URL')



mutually_exclusive_group.add_argument("-r", "--reads_download", 
    help="download short reads for analysis", 
    type = str,
    nargs = '+',
    metavar = 'ACCESSION_NUMBER')



mutually_exclusive_group.add_argument("-R", "--reads_path", 
    help="path to local short reads in fastq files for analysis. All read files will be collcted. Alternatively, a pattern containing '*' '?' and other shell expansion characters or an explicit list of read files can also be supplied.", 
    type = str,
    nargs = '+')

parser_CollectData.add_argument('-n', "--reads_group_name", 
    help = "optional for downloading from NCBI, required for loading from local path")



parser_CollectData.add_argument('-e', "--email_address", 
    type = str,
    help = "required for downloading from NCBI")


parser_PrepareReads = subparser_adder.add_parser('PrepareReads',
                parents = [common_module_arguments],
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Prepare reads for alignment to \
genome sequence by removing adaptor sequences and trimming by position specific \
quality scores. Align reads to a reference genome sequence.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('Groups of read data sets from the \
CollectData option are loaded by providing the file name.\n\n\
Example usage: "%(prog)s --reads_name ERR953490plus5others --adaptors --trim --align FM209186"\n', 
text_width, replace_whitespace = False))



parser_PrepareReads.add_argument('-n', "--reads_name", 
    help = "name of read datasets group generated by the CollectData option",
    type = str,
    required = True)

parser_PrepareReads.add_argument('-s', "--subsample_to_cov", 
    help = "subsample reads to a requested average coverage depth of a given \
genome size. This provides smaller files of a consistent size saving storage \
space and processing time later and benefitting some analyses like de novo \
assembly. E.g.\n\
'--subsample_to_cov 80 5000000' for 80x coverage of a 5 Mb genome.",
    type = int,
    nargs = 2,
    metavar = ('COVERAGE_DEPTH','GENOME_LENGTH'))

parser_PrepareReads.add_argument('-a', "--adaptors", 
    help = "cut residual adaptor sequences from reads using CutAdapt. Defaults to \
using read sets previously subsampled to adequate estimated reference genome \
coverage. Optionally use original reads, even if subsampling has been performed, \
with '--adaptors fullsample'.",
    type = str,
    nargs='?',
    const = 'subsample',
    choices = ['fullsample','subsample'])

parser_PrepareReads.add_argument('-t', "--trim", 
    help = "trim read ends based on quality scores using Sickle",
    action = 'store_true')

parser_PrepareReads.add_argument('-D', "--delete_intermediates", 
    help = "delete intermediate fastq files to save space. Files are only deleted if those for next stage are found",
    action = 'store_true')


parser_AlignReads = subparser_adder.add_parser('AlignReads',
                parents = [common_module_arguments],
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Align reads to a reference genome sequence.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('Reads prepared by the PrepareReads module are loaded. \
A genome prepared by the CollectData option are loaded by providing the file name.\n\n\
Example usage: "%(prog)s --reads_name ERR953490plus5others --genome_name FM209186 --align"\n', 
text_width, replace_whitespace = False))

parser_AlignReads.add_argument('-n', "--reads_name", 
    help = "name of read datasets group prepared by the PrepareReads option",
    type = str,
    required = True)

parser_AlignReads.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option",
    type = str,
    required = True)

parser_AlignReads.add_argument('-a', "--align", 
    help = "align reads to a genome using Burrows Wheeler Aligner (BWA). \
Requires name of genome prepared by CollectData option. Convert to BAM, sort \
and index",
    action = 'store_true')

parser_AlignReads.add_argument('-d', "--deduplicate", 
    help = "Remove duplicates using Picard",
    action = 'store_true')

parser_AlignReads.add_argument('-r', "--indelrealign", 
    help = "realign read alignments near potential indels using GATK",
    action = 'store_true')

parser_AlignReads.add_argument('-p', "--prepared", 
    help = "if your reads were trimmed and cleaned etc. already, without using BAGA's PrepareReads options, this option allows you to continue directly after the CollectData options.",
    action = 'store_true')

parser_AlignReads.add_argument('-P', "--GATK_jar_path", 
    help = "path to Genome Analysis Toolkit (GATK) jar file. See also --JRE_path if system Java version is not compatible with your GATK version",
    type = str)

parser_AlignReads.add_argument('-J', "--JRE_path", 
    help = "path to JAVA runtime environment binary file for use with GATK, Picard and other software written in Java",
    type = str)

parser_AlignReads.add_argument('-D', "--delete_intermediates", 
    help = "delete intermediate SAM and BAM files to save space. Files are only deleted if those for next stage are found",
    action = 'store_true')


parser_SimulateReads = subparser_adder.add_parser('SimulateReads',
                parents = [common_module_arguments],
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Simulate reads with optional '\
                        'variants from a reference genome sequence.',
                        text_width,
                        replace_whitespace = False),
                epilog = textwrap.fill('Example usage: "%(prog)s --genome_name FM209186 --num_SNPs 100"\n', 
text_width, replace_whitespace = False))

parser_SimulateReads.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option from which to generate reads from",
    type = str,
    required = True)

parser_SimulateReads.add_argument('-G', "--gemsim", 
    help = "generate reads using GemSIM",
    action = 'store_true')

parser_SimulateReads.add_argument('-n', "--num_individuals", 
    help = "genome population size",
    type = int,
    default = 1)

parser_SimulateReads.add_argument('-D', "--large_deletions", 
    help = "ranges of large deletions e.g., prophage in reference missing in "\
            "samples. If specified, a second set of individuals are generated "\
            "with these deletions. Currently, total variants are shared among "\
            "this double sized population. In the future, the second set will "\
            "probably share same variants as first set, the only difference "\
            "being the large deletions. e.g., '--large_deletions 10000 20000 "\
            "80000 90000' will omit regions 10000-20000 bp and 80000-90000 bp.",
    type = int,
    nargs = '+',
    metavar = 'CHROMOSOME_POSITION')

parser_SimulateReads.add_argument('-r', "--random_seed", 
    help = "set the random seed for the pseudo-random number generator for "\
            "reproducible variants.",
    type = int,
    default = 684651)

parser_SimulateReads.add_argument('-s', "--num_SNPs", 
    help = "total single nucleotide polymorphisms in population",
    type = int)

parser_SimulateReads.add_argument('-d', "--num_deletions", 
    help = "total small deletion polymorphisms in population",
    type = int)

parser_SimulateReads.add_argument('-i', "--num_insertions", 
    help = "total small deletion polymorphisms in population",
    type = int)


parser_Repeats = subparser_adder.add_parser('Repeats',
                parents = [common_module_arguments],
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Detect repetitive regions in a \
genome sequence. Plot these regions as chromosome map with percent \
identity and some features. The intended use is to mark the repetitive regions \
for exclusion from a short read mapping experiment. Repetitive regions are \
more likely to contain ambiguous variants that may be caused by divergence \
between duplicated regions in the reference genome or sample and not by \
variation at orthologous regions.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('The Repeats finder is expected to \
operate on a genome obtained via the CollectData option.\n\n\
Example usage: %(prog)s --genome_name FM209186.1 --find\n', 
text_width, replace_whitespace = False))

parser_Repeats.add_argument('-g', "--genome_name", 
    help = "name of genome obtained and prepared by the CollectData option",
    type = str,
    required = True)

parser_Repeats.add_argument('-f', "--find", 
    help = "find repeats using Burrows-Wheeler Aligner (BWA) alignments",
    action = 'store_true')

parser_Repeats.add_argument('-m', "--method", 
    help = "The baga method (default) can be used as part of a variant "\
    "calling pipeline by providing a variant filter at regions that are "\
    "so similar, they might cause ambiguous read alignment. The nucmer "\
    "method (nucmer_check) is to provide a means to compare a previous BAGA "\
    "repeats analysis with the 'classic' nucmer method from the MUMmer "\
    "package. The BAGA method is much slower but performs optimal global "\
    "alignments between repeats and performs codon alignments where "\
    "appropriate, so is more accurate than nucmer but much slower.",
    type = str,
    default = "baga",
    choices = ["baga","nucmer_check"])

parser_Repeats.add_argument('-i', "--minimum_percent_identity", 
    help = "the lowest nucleotide percent identity permitted repeat over regions",
    type = int,
    default = 98)

parser_Repeats.add_argument('-l', "--minimum_repeat_length", 
    help = "shortest repeat length: should be similar to insert size of "\
            "paired end reads as smaller (non-tandem) repeats should be "\
            "resolvable and have unambiguous mappings",
    type = int,
    default = 400)

parser_Repeats.add_argument('-p', "--plot", 
    help = "plot repeats using svgwrite library found using '--find'",
    action = 'store_true')

parser_Repeats.add_argument('-s', "--summarise", 
    help = "summarise repeats found using '--find' as a .csv file and printed to screen",
    action = 'store_true')


parser_Structure = subparser_adder.add_parser('Structure',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Detect regions with rearrangements \
between sample genome and reference genome sequence by examining distribution of \
inferred insert sizes from paired end short reads of sample aligned to reference. \
Proportion of paired reads not designated as "proper_pair" by Burrows-Wheeler \
Aligner is expected to increase over regions of sequence rearrangements. Alignments \
in these regions violate assumptions of read mapping method and and any variants \
called at these regions should be excluded from down-stream analyses. \n\nPlot \
these regions as chromosome map with ratio of non-proper pair to proper pair \
classifications and threshold at which deemed to deviate from expected.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('The Structure checker is expected to \
operate on a genome obtained via the CollectData option; and either a set of BAM \
files in a user provided path or a reads group processed with the AlignReads option. \
\n\n\
Example usages:\n\
%(prog)s --genome_name FM209186.1 --reads_name Liverpool --check\n\
%(prog)s --genome_name FM209186.1 --alignments_paths path/to/my/bams --check\n', 
text_width, replace_whitespace = False))

mutually_exclusive_group = parser_Structure.add_mutually_exclusive_group(required=True)

mutually_exclusive_group.add_argument('-n', "--reads_name", 
    help = "name of read datasets group if processed by PrepareReads and AlignReads options. Should match --reads_name option used previously",
    type = str)

mutually_exclusive_group.add_argument('-a', "--alignments_paths", 
    help = "path to paired-end short read alignments to a reference genome. If a directory path(s) is provided, all *.BAM and *.bam files will be included. A list of filenames with full path or with *, ? or [] characters can also be provided (with unix-style pathname pattern expansion for the latter)",
    type = str,
    nargs = '+')

mutually_exclusive_group.add_argument('-b', "--checkinfos_path", 
    help = "path to baga.Structure.CheckerInfo-<sample_name>__<genome_name>.baga files for summarising regions indicated by rearrangements filter with --summarise. Alternatively use --summarise with --reads_name and --genome_name if available (typically as part of a baga short read analysis pipeline). If a directory path(s) is provided, all baga.Structure.CheckerInfo-*__*.baga files will be included. A list of filenames with full path or with *, ? or [] characters can also be provided (with unix-style pathname pattern expansion for the latter)",
    type = str,
    nargs = '+')

parser_Structure.add_argument('-g', "--genome_name", 
    help = "required with --reads_name option. Name of genome obtained by the CollectData option and used with AlignReads option",
    type = str)

parser_Structure.add_argument('-c', "--check", 
    help = "check sequence rearrangements using paired end insert sizes inferred by Burrows-Wheeler Aligner (BWA) alignments",
    action = 'store_true')

parser_Structure.add_argument('-p', "--plot", 
    help = "plot regions affected by structural rearrangements found using '--check' Uses the svgwrite library",
    action = 'store_true')

parser_Structure.add_argument('-r', "--plot_range", 
    help = "plot a specific region to see how it was affected by structural rearrangements found using '--check'. Uses the svgwrite library",
    type = int,
    nargs = 2,
    metavar = 'CHROMOSOME_POSITION')

parser_Structure.add_argument('-t', "--ratio_threshold", 
    help = "When checking for rearrangements, the ratio of non-proper pair to proper pair assigned reads above this value cause a region of read alignments to be considered rearranged. This ratio tends to zero rhen the distance between aligned paired reads is close to the expectation according to the estimated mean fragment size. It increases to around one adjacent to rearrangements e.g., within a fragment's length of a large deletion in the sample/insertion in the reference. Lower values are more sensitive to rearrangements but might include false positive rearrangements, but these can be examined by local de novo assembly of reads and pairwise alignment of contig with reference. If used to filter regions affected by unreliable short read alignments for variant calling, lower values are more conservative (will exclude more false positive variants) but might cause omission of true positive variants. Default = 0.15",
    type = float,
    default = 0.15,
    metavar = 'FLOAT')

## these two should probably be in a parent parser inherited for most tasks
parser_Structure.add_argument('-S', "--include_samples", 
    help = "With --plot or --plot_range, restrict plotting to these samples else if omitted, plots for all samples are produced.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_Structure.add_argument('-e', "--exclude_samples", 
    help = "exclude these samples from the analysis. If omitted, no samples are excluded.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_Structure.add_argument('-s', "--summarise", 
    help = "summarise regions of putative rearrangements found using '--check' as a .csv file and printed to screen. Requires --genome and --reads_name, optionally with --include_samples. Or just --genome_name and --include_samples to specify samples with corresponding 'baga.Structure.CheckerInfo-<sample_name>__<genome_name>.baga' available in current folder.",
    action = 'store_true')

parser_Structure.add_argument('-C', "--collect", 
    help = "extract short reads aligned to regions of putative rearrangements found using '--find', write to a .fastq file, assemble each of them using SPAdes and align contigs back to reference chromosome. Requires --reads_name, optionally with --include_samples to limit to specific samples and --collect_range to specify a range to collect reads for if not those already found by '--find'. Alternatively use --checkinfos_path to specify samples with corresponding 'baga.Structure.CheckerInfo-<sample_name>__<genome_name>.baga'.",
    action = 'store_true')

parser_Structure.add_argument('-R', "--collect_ranges", 
    help = "extract short reads aligned at specified region or regions, write to a .fastq file, assemble de novo using SPAdes and align contigs back to reference chromosome. --num_padding_positions will be set to zero. If more than one region set e.g., --collect_ranges 10000 20000 80000 90000, the reads in range 10000-20000 bp and 80000-90000 bp along with unmapped and poorly mapped reads will be assembled together.",
    type = int,
    nargs = '+',
    metavar = 'CHROMOSOME_POSITION')

parser_Structure.add_argument('-F', "--force", 
    help = "overwrite existing assemblies when using --collect/-C.",
    action = 'store_true',
    default = False)

parser_Structure.add_argument('-P', "--num_padding_positions", 
    help = "for the de novo assembly of collected reads that were aligned/mapped to the reference genome at regions of putative rearrangements, optionally specify additional padding (bp) around each region to collect more reads and assmeble longer contigs for pairwise alignment back to the reference",
    type = int,
    default = 5000,
    metavar = 'NUM_BASEPAIRS')

parser_Structure.add_argument('-m', "--max_memory", 
    help = "maximum memory to use in gigabytes for each assembly. If not specified, total available at launch time will be used.",
    type = int)

parser_Structure.add_argument('-l', "--min_align_region", 
    help = "when using --collect, set minimum region to align among those reported as potentially rearranged (by --check).",
    type = int,
    default = 200)



parser_CallVariants = subparser_adder.add_parser(
        'CallVariants',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = textwrap.fill(
                'Call variants with Genome Analysis Tool Kit (GATK) from each '\
                'of a group of read sets previously aligned to a genome via the '\
                'PrepareReads option.',
                 text_width,
                 replace_whitespace = False),
        epilog = textwrap.fill('Groups of read sets from the AlignReads option '\
                'are loaded by providing the name supplied previously.\n\n'\
                'Example usage: %(prog)s --reads_name ERR953490plus5others '\
                '--calleach --calljoint --hardfilter\n', 
text_width, replace_whitespace = False))

parser_CallVariants.add_argument('-r', "--reads_name", 
    help = "name of read datasets group processed by PrepareReads and "\
    "AlignReads options. For GATK a single group can be processed. For "\
    "DiscoSNP++, one or more groups can be processed",
    type = str,
    required = False,
    nargs = '+')

parser_CallVariants.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option",
    type = str,
    required = False)

parser_CallVariants.add_argument('-F', "--force", 
    help = "overwrite existing per isolate data: required when repeating an "\
    "analysis else previous versions retained. Retention of previous versions "\
    "is convenient for resuming an interrupted analysis in which only some "\
    "read sets were processed.",
    action = 'store_true',
    default = False)

parser_CallVariants.add_argument('-C', "--check", 
    help = "check previously called variants in a VCF by de novo assembly of "\
    "relevent shorts reads and aligning the contigs to each variant part of "\
    "the reference genome.",
    action = 'store_true',
    default = False)

parser_CallVariants.add_argument('-p', "--vcfs_paths", 
    help = "path to vcf files. If a directory path(s) is provided, all *.VCF "\
    "and *.vcf files will be included. A list of filenames with full path or "\
    "with *, ? or [] characters can also be provided (with unix-style pathname "\
    "pattern expansion for the latter)",
    type = str,
    nargs = '+')

parser_CallVariants.add_argument('-a', "--alignments_paths", 
    help = "path to paired-end short read alignments to a reference genome. If a "\
    "directory path(s) is provided, all *.BAM and *.bam files will be included. A "\
    "list of filenames with full path or with *, ? or [] characters can also be "\
    "provided (with unix-style pathname pattern expansion for the latter)",
    type = str,
    nargs = '+',
    metavar = 'PATH_TO_BAMs')

parser_CallVariants.add_argument('-N', "--new", 
    help = "start new variant calling: required when starting whole \
CallVariants process again from prepared read alignments.",
    action = 'store_true',
    default = False)

parser_CallVariants.add_argument('-n', "--max_cpus", 
    help = "maximum number of cpu cores used when parallel processing",
    type = int,
    default = -1)

parser_CallVariants.add_argument('-m', "--max_memory", 
    help = "maximum memory to use in gigabytes. If not "\
    "specified, GATK will use 8GB, SPAdes will use total "\
    "available at launch time.",
    type = int)

parser_CallVariants.add_argument('-s', "--callsingles", 
    help = "call variants in each alignment on a per sample basis (not for a joint "\
"analysis, see --calleach and --calljoint for that). "\
"Called 1st time on uncalibrated alignments, called 2nd after base quality score "\
"recalibration.",
    action = 'store_true')

parser_CallVariants.add_argument('-c', "--calleach", 
    help = "call variants in each alignment in preparation for a joint analysis. "\
"Called 1st time on uncalibrated alignments, called 2nd after base quality score "\
"recalibration.",
    action = 'store_true')

parser_CallVariants.add_argument('-j', "--calljoint", 
    help = "call variants in all alignments in a joint analysis. Called 1st time "\
"on uncalibrated alignments, called 2nd after base quality score recalibration.",
    action = 'store_true')

parser_CallVariants.add_argument('-f', "--hardfilter", 
    help = "apply 'hard filtering' thresholds on called variants to decrease "\
"false positives using GATK",
    action = 'store_true')

parser_CallVariants.add_argument('-R', "--recalibrate", 
    help = "apply read base quality score recalibration using GATK",
    action = 'store_true')

parser_CallVariants.add_argument('-P', "--GATK_jar_path", 
    help = "path to Genome Analysis Toolkit (GATK) jar file. See also "\
"--JRE_path if system JAVA is version 1.8",
    type = str)

parser_CallVariants.add_argument('-J', "--JRE_path", 
    help = "path to JAVA runtime environment version 1.7 binary file for use "\
"with GATK versions 3.3 or 3.4 (not compatible with JRE 1.8!)",
    type = str)

parser_CallVariants.add_argument('-d', "--calldisco", 
    help = "call variants de novo from short reads using DiscoSNP++.",
    action = 'store_true')

parser_CallVariants.add_argument('-e', "--use_existing_graph", 
    help = "Use previously generated DiscoSNP++ graph. Be sure the last graph generated matches the specified reads!",
    action = 'store_true')

parser_CallVariants.add_argument('-A', "--arguments", 
    help = "Send direct arguments to a wrapped tool. E.g. --arguments DiscoSNP++ '_k 41'. The leading dashes (-) must be replaced with underscores (_) and the arguments must be in quotations",
    nargs = '*')


parser_FilterVariants = subparser_adder.add_parser('FilterVariants',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Apply filters determined by \
the Repeats and Structure options on variant calls and report tables of \
effects of filters on different classes of variants. VCF files will be \
copied with updated information marking certain variants to be excluded. \
Vaiants can be inferred using the CallVariants option.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('Filter regions and VCF files are \
loaded by providing the read set names or VCF filenames.\n\n\
Example usage: %(prog)s --reads_name ERR953490plus5others \
--genome FM209186 --filters genome_repeats rearrangements\n', 
text_width, replace_whitespace = False))

mutually_exclusive_group = parser_FilterVariants.add_mutually_exclusive_group(required=True)

mutually_exclusive_group.add_argument('-n', "--reads_name", 
    help = "name of read datasets group if processed by PrepareReads and AlignReads options. Should match --reads_name option used previously",
    type = str,
    nargs = '+')

mutually_exclusive_group.add_argument('-p', "--vcfs_paths", 
    help = "path to vcf files. If a directory path(s) is provided, all *.VCF "\
    "and *.vcf files will be included. A list of filenames with full path or "\
    "with *, ? or [] characters can also be provided (with unix-style pathname "\
    "pattern expansion for the latter)",
    type = str,
    nargs = '+')

parser_FilterVariants.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option",
    type = str,
    required = True)

parser_FilterVariants.add_argument('-f', "--filters", 
    help = "names of filters to apply. One or more of: genome_repeats provided by the Repeats option; rearrangements provided by the Structure option.",
    type = str,
    nargs = '+',
    # these choices must correspond to known_filters in CallVariants.Filter.__init__()
    # or to other_filters which are listed by variant calling program, currently just GATK
    choices = ['genome_repeats','rearrangements','GATK'],
    metavar = 'FILTER_NAME',
    required = True)

parser_FilterVariants.add_argument('-s', "--include_samples", 
    help = "restrict filtering to these samples. If omitted, plots for all samples are produced.",
    type = str,
    nargs = '+')

parser_FilterVariants.add_argument('-i', "--path_to_rearrangements_info", 
    help = "optionally supply path to where rearrangement filter information is for all samples (if not in current directory; they look like baga.Structure.CheckerInfo-<samplename>__<genomename>.baga). These must be generated using the 'Structure --check' option and are each generated from the same .bam alignment file as the corresponding, supplied VCF files.",
    type = str)


parser_SummariseVariants = subparser_adder.add_parser('SummariseVariants',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = textwrap.fill('Generate various .csv text files for '\
                'convenient viewing and downstream analysis.', text_width,
                replace_whitespace = False),
        epilog = textwrap.fill('Example usage: %(prog)s --simple --vcfs_paths '\
                'path/to/*.vcf\n', text_width, replace_whitespace = False))

mutually_exclusive_group = parser_SummariseVariants.add_mutually_exclusive_group(required=True)

# either: reads_name and genome_name to infer CallVariants.Caller object path
mutually_exclusive_group.add_argument('-n', "--reads_name", 
        help = "name of read datasets group if processed by CallVariants options. "\
        "Should match --reads_name option used previously",
        type = str,
        nargs = '+')

mutually_exclusive_group.add_argument('-p', "--vcfs_paths", 
        help = "path to vcf files. If a directory path(s) is provided, all *.VCF "\
        "and *.vcf files will be included. A list of filenames with full path or "\
        "with *, ? or [] characters can also be provided (with unix-style pathname "\
        "pattern expansion for the latter)",
        type = str,
        nargs = '+')

# even without reads_name, genome is still useful to provide additional annotations
parser_SummariseVariants.add_argument('-g', "--genome_names", 
        help = "name of genome obtained by the CollectData option",
        type = str,
        nargs = '+')

parser_SummariseVariants.add_argument('-f', "--filters", 
        help = "names of filters to apply. One or more of: genome_repeats provided "\
        "by the Repeats option; rearrangements provided by the Structure option.",
        type = str,
        nargs = '+',
        # these choices must correspond to known_filters in CallVariants.Filter.__init__()
        # or to other_filters which are listed by variant calling program, currently just GATK
        choices = ['genome_repeats','rearrangements','GATK'],
        metavar = 'FILTER_NAME')

parser_SummariseVariants.add_argument('-S', "--simple", 
        help = "generate a simple .csv table corresponding to the VCF file rows.",
        action = 'store_true')

parser_SummariseVariants.add_argument('-L', "--lists", 
        help = "generate a .csv table listing allele frequencies in whole dataset, "\
        "between samples and reference and among samples excluding reference.",
        action = 'store_true')

parser_SummariseVariants.add_argument('-C', "--cumulative", 
        help = "summarise the cumulative effect of filters in .csv table.",
        action = 'store_true')

parser_SummariseVariants.add_argument('-s', "--include_samples", 
        help = "restrict filtering to these samples. If omitted, plots for all "\
        "samples are produced.",
        type = str,
        nargs = '+')


parser_CheckLinkage = subparser_adder.add_parser('CheckLinkage',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Check for alleles on the same \
read or read pair (fragment; and therefore the same chromosomes) called at \
polymorphisms in a sample of pooled genomic DNA. Infrequent co-incidence of \
variants on the same read in nearby polymorphisms implies variants occuring in \
different genomes in the sample (separate lineages) and has been described as a \
"multidiverse" signature in:\n\n\
Lieberman, T. D., Flett, K. B., Yelin, I., Martin, T. R., McAdam, A. J., Priebe, G. P. & Kishony, R. \n\
Genetic variation of a bacterial pathogen within individuals with cystic fibrosis provides a record of selective pressures.\n\
Nature Genetics, 2013, 46, 82-87.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('The linkage checker is expected to \
operate on a short reads aligned to a genome and the corresponding variants. \
These can be provided by a previous BAGA analysis or by providing paths to BAM \
and VCF files. \
\n\n\
Example usages:\n\
%(prog)s --genome_name FM209186.1 --reads_name Liverpool --check\n\
%(prog)s --alignments_paths path/to/my/bams --vcfs_paths path/to/my/vcfs --check --genome_name FM209186.1\n', 
text_width, replace_whitespace = False))

mutually_exclusive_group = parser_CheckLinkage.add_mutually_exclusive_group(required=True)

mutually_exclusive_group.add_argument('-n', "--reads_name", 
    help = "name of read datasets group if processed by PrepareReads and AlignReads options. Should match --reads_name option used previously",
    type = str)

mutually_exclusive_group.add_argument('-a', "--alignment_paths", 
    help = "path to paired-end short read alignments to a reference genome. If a directory path(s) is provided, all *.BAM and *.bam files will be included. A list of filenames with full path or with *, ? or [] characters can also be provided (with unix-style pathname pattern expansion for the latter)",
    type = str,
    nargs = '+',
    metavar = 'PATH_TO_BAMs')

parser_CheckLinkage.add_argument('-g', "--genome_name", 
    help = "Name of genome obtained by the CollectData option and if the BAGA variant calling pipeline was used, the same genome used with AlignReads and CallVariants options",
    type = str,
    required=True)

parser_CheckLinkage.add_argument('-p', "--vcfs_paths", 
    help = "path to vcf files. If a directory path(s) is provided, all *.VCF and *.vcf files will be included. A list of filenames with full path or with *, ? or [] characters can also be provided (with unix-style pathname pattern expansion for the latter)",
    type = str,
    nargs = '+',
    metavar = 'PATH_TO_VCFs')

parser_CheckLinkage.add_argument('-c', "--check", 
    help = "check variant linkage on paired reads (fragments) in a pooled dataset",
    action = 'store_true')

parser_CheckLinkage.add_argument('-S', "--include_samples", 
    help = "Restrict checking to these samples else if omitted, all samples are checked.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_CheckLinkage.add_argument('-F', "--force", 
    help = "overwrite existing files when using --collect/-C.",
    action = 'store_true',
    default = False)


parser_ComparativeAnalysis = subparser_adder.add_parser('ComparativeAnalysis',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Build multiple sequence alignments \
from SNPs in VCFs with a reference genome, infer phylogenies and homologous \
recombination events.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('Example usage: %(prog)s --buildMSA \
--reads_name Liverpool --genome FM209186.1  --include_invariants\n', 
text_width, replace_whitespace = False))

mutually_exclusive_group1 = parser_ComparativeAnalysis.add_mutually_exclusive_group(required=True)

mutually_exclusive_group1.add_argument('-m', "--build_MSA", 
    help = "build a multiple sequence alignment from a reference genome and SNPs listed in VCF files",
    action = 'store_true')

mutually_exclusive_group1.add_argument('-i', "--infer_phylogeny", 
    help = "infer a phylogeny from a multiple sequence alignment",
    action = 'store_true')

mutually_exclusive_group1.add_argument('-r', "--infer_recombination", 
    help = "infer recombination from a phylogeny and multiple sequence alignment",
    action = 'store_true')

mutually_exclusive_group1.add_argument('-p', "--plot_phylogeny", 
    help = "plot a phylogeny, possibly including recombination inferences",
    action = 'store_true')

# for build_MSA
parser_ComparativeAnalysis.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option on which to base a MSA along with SNPs",
    type = str)

mutually_exclusive_group2 = parser_ComparativeAnalysis.add_mutually_exclusive_group(required=False)

mutually_exclusive_group2.add_argument('-n', "--reads_name", 
    help = "name of read datasets groups to include in an MSA, if processed by PrepareReads and AlignReads options. Should match --reads_name option used previously",
    type = str,
    nargs = '+')

mutually_exclusive_group2.add_argument('-v', "--vcfs_paths", 
    help = "path to vcf files. If a directory path(s) is provided, all *.VCF and *.vcf files will be included. A list of filenames with full path or with *, ? or [] characters can also be provided (with unix-style pathname pattern expansion for the latter)",
    type = str,
    nargs = '+')

parser_ComparativeAnalysis.add_argument('-s', "--include_samples", 
    help = "restrict MSA to these samples. If omitted, all samples are included.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_ComparativeAnalysis.add_argument('-e', "--exclude_samples", 
    help = "exclude these samples from MSA. If omitted, no samples are excluded.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_ComparativeAnalysis.add_argument('-l', "--include_invariants", 
    help = "include invariant sites to make a Long alignment. Invariant sites are required to accurately estimate some parameters",
    action = 'store_true')

parser_ComparativeAnalysis.add_argument('-c', "--core_only", 
    help = "only include sites present in all samples",
    action = 'store_true')

parser_ComparativeAnalysis.add_argument('-B', "--sample_bams", 
    help = "path to file linking sample names to BAM files to collect reference genome coverage information. Format: prer line <sample name><tab><path to BAM>",
    type = str)

# for infer_phylo
parser_ComparativeAnalysis.add_argument('-M', "--path_to_MSA", 
    help = "path to multiple nucleotide sequence alignment (e.g. as generated by --build_MSA). Amino acid sequences will be implemented in the future.",
    type = str)

parser_ComparativeAnalysis.add_argument('-P', "--program", 
    help = "software to use for inference: currently only phyML!",
    type = str,
    choices = ['phyML'],
    default = 'phyML')

parser_ComparativeAnalysis.add_argument('-o', "--out_group", 
    help = "outgroup for rooting. Can be single sample name or several but must be monophyletic: not checked!",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

# other parameters for phyML to follow

# infer phylogeny and recombination
parser_ComparativeAnalysis.add_argument('-b', "--num_bootstraps", 
    help = "number of bootstrap replicates for --infer_phylo or --infer_recombination.",
    type = int,
    default = 0)

# infer recombination
parser_ComparativeAnalysis.add_argument('-T', "--plot_transfers", 
    help = "plot a phylogeny, also plot transfers inferred by ClonalFrame (in <treename>.importation_status.txt)",
    action = 'store_true')

parser_ComparativeAnalysis.add_argument('-t', "--path_to_tree", 
    help = "path to starting phylogeny for inferring recombination over and correct.",
    type = str)


# plot phylogeny
# parser_ComparativeAnalysis.add_argument('-t', "--path_to_CFML_tree", 
    # help = "For plotting phylogeny with transfers, path to phylogeny updated by ClonalFrameML when inferring recombination.",
    # type = str)

parser_ComparativeAnalysis.add_argument('-N', "--use_names", 
    help = "For plotting phylogeny, path to file with tab delimited list of actual tip label and desired label, for translating labels.",
    type = str)

parser_ComparativeAnalysis.add_argument('-G', "--genome_length", 
    help = "For plotting transfers on phylogeny, reference genome length (alternatively supply --genome_name of CollectData saved object).",
    type = int)

parser_AssembleReads = subparser_adder.add_parser('AssembleReads',
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = textwrap.fill('Assemble reads into contiguous chromosome sequences.',
                                text_width,
                                replace_whitespace = False),
    epilog = textwrap.fill('Example usage: %(prog)s --denovo spades '\
    '--reads_name Liverpool\n', 
    text_width, replace_whitespace = False))

parser_AssembleReads.add_argument('-n', "--reads_name", 
    help = "name of read datasets groups on which to perform assembly. Reads "\
    "should have collected and prepared with BAGA's CollectData and "\
    "PrepareReads options. Should match --reads_name option used previously",
    required = True,
    type = str,
    nargs = '+')

parser_AssembleReads.add_argument('-s', "--include_samples", 
    help = "Implement Me! restrict assembly to these samples within the read "\
    "datasets groups provided by --reads_name/-n. If omitted, all samples are "\
    "included.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_AssembleReads.add_argument('-e', "--exclude_samples", 
    help = "Implement Me! exclude these samples from assembly from among samples "\
    "within the read datasets groups provided by --reads_name/-n. If omitted, no "\
    "samples are excluded.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_AssembleReads.add_argument('-P', "--program", 
    help = "software to use for assembly: currently only spades!",
    type = str,
    choices = ['spades'],
    default = 'spades')

parser_AssembleReads.add_argument('-m', "--max_memory", 
    help = "maximum memory to use in gigabytes for each assembly. If not "\
    "specified, total available at launch time will be used.",
    type = int)

if '--nosplash' not in sys.argv:
    print(splash)

args = parser.parse_args()

task_name = args.subparser
    

# generate folder for logs and path to main log file
# includes CLI-provided sample_name in log filename path:
# ensure no awkward characters for filename path
if hasattr(args, 'sample_name') and args.sample_name is not None:
    use_sample_name, sanitised = sanitize_filename(args.sample_name)
elif task_name == 'Dependencies':
    # allow for CLI use without a sample name e.g. checking Dependencies
    # logging file paths will include 'dependencycheck'
    use_sample_name, sanitised = 'dependencycheck', ''
elif task_name == 'CollectData' and args.taxonomy is not None:
    use_sample_name, sanitised = 'taxonomy', ''
elif task_name == 'CollectData' and args.genomes is not None:
    use_sample_name, sanitised = 'getgenomes', ''
else:
    # fall back
    # logging file paths will include 'nosample'
    use_sample_name, sanitised = 'nosample', ''

# timestamped folder so sorting by time places it ahead of log for module usage
time_for_main_log = datetime.now().strftime("%y-%m-%d_%H-%M-%S")
main_log_folder = os.path.sep.join([
        'baga-{}_logs'.format(use_sample_name),
        # 00 so it sorts above other log folders made immediately after
        '{}_0_CLI'.format(time_for_main_log)])
try:
    os.makedirs(main_log_folder)
except FileExistsError:
    pass

main_log_filename = os.path.sep.join([main_log_folder,'main.log'])

# NB CL options increase with verbosity, logging level decreases!
# make a dict to convert CLI verbosity amount 0-4 to logging level 
# threshold 0-50
from baga import PROGRESS
verbosities = {}
verbosities[0] = logging.ERROR
verbosities[1] = logging.WARNING
verbosities[2] = logging.INFO
verbosities[3] = PROGRESS
verbosities[4] = logging.DEBUG

baga_cli, main_log_folder = configureLogger(use_sample_name, main_log_filename, 
        verbosities[args.verbosity])


# start logging
time_start = datetime.now()
time_stamp = time_start.strftime("%H:%M:%S on %a %d %b, %Y")
baga_cli.info('')
baga_cli.info('|=== Starting baga analysis at {} ===|'.format(time_stamp))
baga_cli.info('')

# now main logging set up, inform if sample name contained any funky characters
# that had to be omitted
if len(sanitised):
    baga_cli.warning('Sample name changed for use in filenames: {} changed '\
            'to {} by replacing {} with _'.format(args.sample_name, 
            use_sample_name, ','.join(sanitised)))



if hasattr(args, 'GATK_jar_path') and args.GATK_jar_path:
    # check Java is installed and version
    use_java = 'java'
    if hasattr(args, 'JRE_path') and args.JRE_path:
        #'-j', "--JRE_path"
        use_java = args.JRE_path
        print('Using provided JAVA JRE at {}'.format(use_java))
    else:
        print('Using system JAVA JRE')
    
    try:
        p = subprocess.Popen([use_java, '-version'], stdout = subprocess.PIPE, 
                stderr = subprocess.STDOUT)
    except OSError:
        sys.exit('Could not find the Java runtime environment, needed for GATK '\
                'and Picard. Please ensure you have Java installed and that the '\
                '"java" executable is in the system path or that you provided '\
                'the correct path for --JRE_path')
    
    try:
        java_version_string = p.stdout.read().split('\n')[0]
        java_ver_maj, java_ver_min1, java_ver_min2 = re.findall(
                '([0-9])\.([0-9])\.([0-9_]+)', java_version_string)[0]
        java_ver_maj, java_ver_min1 = map(int,[java_ver_maj, java_ver_min1])
        if java_ver_maj == 1 and java_ver_min1 == 7:
            print('Java 1.7: found!')
        elif java_ver_maj == 1 and java_ver_min1 == 8:
            print('Java 1.8: found!')
        else:
            sys.exit('GATK v3.3 to v3.5 requires Java v1.7, v3.6 can run with '\
                    'Java v1.8 but your version is {}.{}. Please install Java '\
                    'v1.7 to continue or use --JRE_path to specify Java v1.7 '\
                    'binary to use.'.format(java_ver_maj,java_ver_min1))
    except IndexError:
        print(output)
        sys.exit('There was a problem checking your Java version. Please '\
                'report this as a baga bug at '\
                'https://github.com/daveuu/baga/issues.')
    
    # check GATK jar file
    p = subprocess.Popen([use_java, '-Xmx8g', '-jar', args.GATK_jar_path, 
            '--help'], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    output = p.stdout.read()
    if 'Error: Unable to access jarfile' in output:
        sys.exit('Java could not find jar file for GATK at: {}\nPlease provide '\
                'path to a GATK v3.3 to v3.6 jar file with --GATK_jar_path.'\
                ''.format(args.GATK_jar_path))
    else:
        try:
            GATK_ver_maj, GATK_ver_min1, GATK_ver_min2, GATK_ver_min3 = \
                    re.findall('\(GATK\) v([0-9])\.([0-9]+)-([0-9]+)-([a-z0-9]+),', 
                    output.split('\n')[1])[0]
            GATK_ver_maj, GATK_ver_min1, GATK_ver_min2 = map(int, 
                    [GATK_ver_maj, GATK_ver_min1, GATK_ver_min2])
            if GATK_ver_maj == 3 and GATK_ver_min1 >= 3:
                print('Found GATK v{}.{}. Should be OK with this version of '\
                        'baga . . .'.format(GATK_ver_maj, GATK_ver_min1))
            elif GATK_ver_maj == 3 and GATK_ver_min1 < 3:
                print('WARNING: Found GATK v{}.{}. Not tested in this baga '\
                        'version (v3.3 and later were) . . . use at your own '\
                        'risk!')
            else:
                sys.exit('Expecting GATK version v3, minor version 3 or '\
                        'greater but found {}.{}\n'\
                        'Please provide path to GATK v3.3 jar file with '\
                        '--GATK_jar_path.'.format(GATK_ver_maj, GATK_ver_min1))
        except IndexError:
            if 'UnsupportedClassVersionError' in output:
                sys.exit('There seems to be a conflict between your version of '\
                        'Java and this version of GATK. GATK v3.6 requires '\
                        'Java JRE 1.8. GATK v3.3 to v3.5 work with either '\
                        'Java JRE 1.7 or v1.8).\n'\
                        'The full error message was:\n{}'.format(output))
            elif 'Invalid or corrupt jarfile' in output:
                sys.exit('Please check the path provided with "--GATK_jar_path". '\
                        'An error was returned:\n{}'.format(output))
            else:
                sys.exit('There was a problem checking your --GATK_jar_path argument. '\
                        'Please report this as a baga bug at '\
                        'https://github.com/daveuu/baga/issues. The full error message '\
                        'from the attempt to run GATK was:\n{}'.format(output))
    
    ## if we got this far, java version and GATK versions are compatible


### Check Dependencies ###

if args.subparser == 'Dependencies':
    print('\n-- Dependencies check/get module --')
    
    if args.versions_file:
        try:
            use_versions_file = args.versions_file
            versions = open(use_versions_file).read()
        except IOError:
            use_versions_file = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1] + ['versions.yaml'])
            try:
                versions = open(use_versions_file).read()
            except IOError:
                use_versions_file = False
        
        updated = False
        if use_versions_file:
            # parse it
            use_versions = {}
            main_key = None
            for line in versions.split('\n'):
                if not line.startswith('#') and len(line.rstrip('\n')) > 0:
                    if not line.startswith(' '):
                        main_key,value = line.rstrip('\n').split(':',1)
                    else:
                        this_key,value = line.rstrip('\n').lstrip(' ').split(':',1)
                        try:
                            use_versions[main_key][this_key] = value.strip(' ')
                        except KeyError:
                            use_versions[main_key] = {this_key : value.strip(' ')}
            for external_name, info in use_versions.items():
                for key, new in info.items():
                    if new == 'None':
                        new = None
                    try:
                        current = Dependencies.dependencies[external_name][key]
                        updated = True
                        if current != new:
                            Dependencies.dependencies[external_name][key] = new
                            print('{} for {} updated from {} to {}'.format(key, 
                                    external_name, current, new))
                    except KeyError:
                        print('Warning: ignoring unrecognised entry in {}: '\
                                '{} => {}'.format(use_versions_file,external_name,key))
        if updated:
            print('Used {} for version information'.format(os.path.abspath(use_versions_file)))
            # update checker exe path from download url where needed
            # this is clunky
            new_exe = Dependencies.dependencies['picard']['url'].split('/')[-1][:-4]
            Dependencies.dependencies['picard']['checker']['arguments']['java_commands'][2] = \
                    os.path.sep.join(['external_programs', new_exe, 'picard.jar'])
            new_exe = Dependencies.dependencies['spades']['url'].split('/')[-1][:-7]
            Dependencies.dependencies['spades']['checker']['arguments']['path'][0] = new_exe
            new_exe = Dependencies.dependencies['mummer']['url'].split('/')[-1][:-7]
            Dependencies.dependencies['mummer']['checker']['arguments']['path'][0] = new_exe
    
    def get(name):
        if Dependencies.dependencies[name]['source'] == 'git':
            Dependencies.get_git(**Dependencies.dependencies[name])
        elif Dependencies.dependencies[name]['source'] == 'download':
            Dependencies.get_download(**Dependencies.dependencies[name])
    
    def check(name):
        # this would need changing with the dependencies dict in Dependencies
        if 'local_packages' in Dependencies.dependencies[name]['destination']:
            checker = Dependencies.dependencies[name]['checker']['function']
            checker_args = Dependencies.dependencies[name]['checker']['arguments']
            result = checker(Dependencies.dependencies[name]['name'], **checker_args)
            #checker(dependencies[name]['name'], system = True)
        elif 'external_programs' in Dependencies.dependencies[name]['destination']:
            checker = Dependencies.dependencies[name]['checker']['function']
            checker_args = Dependencies.dependencies[name]['checker']['arguments']
            result = checker(**checker_args)
        else:
            sys.exit('The destination for this package or program is unknown: {}\n'.format(
                                            Dependencies.dependencies[name]['destination'])
                                            )
        return(result)
    
    def checkpackage(name):
        '''this calls a conventional --check for checking new python packages'''
        import subprocess
        cmd = [sys.argv[0], '--nosplash', 'Dependencies', '--check', name]
        proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE )
        o,e = proc.communicate()
        #print(e)
        if proc.returncode == 1:
            return(False)
        elif proc.returncode == 0:
            return(True)
    
    def checkget(checkgetthese):
        check_summary = []
        check_results = []
        for name in checkgetthese:
            print('\nChecking for {}:\n'.format(name))
            alreadygot = check(name)
            
            if alreadygot:
                check_summary += ['\n{}: found!'.format(name)]
                check_results += [alreadygot]
            else:
                check_summary += ['\n{0}: not found . . . \nAttempting to install.\n'.format(name, sys.argv[0])]
                get(name)
                if 'local_packages' in Dependencies.dependencies[name]['destination']:
                    # it's a python package so need to forget the import if it was an older version and re-import
                    # imported as _name to avoid collisions . . .
                    try:
                        del sys.modules[name]
                        #del globals()['_'+name]
                    except KeyError:
                        pass
                
                if Dependencies.dependencies[name]['checker']['function'] == Dependencies.check_python_package:
                    gotnow = checkpackage(name)
                else:                
                    gotnow = check(name)
                
                if gotnow:
                    check_summary[-1] += "Installed successfully: found!"
                else:
                    check_summary[-1] += "Failed to install . . . there may be dependencies missing or some other problem . . ."
                
                check_results += [gotnow]
            
        return(check_summary,check_results)
    
    if args.get:
        for name in args.get:
            get(name)
    
    if args.check:
        check_summary = []
        check_results = []
        for name in args.check:
            print('\nChecking for {}:\n'.format(name))
            check_results += [check(name)]
            if check_results[-1]:
                check_summary += ['\n{}: found!'.format(name)]
            else:
                check_summary += ['\n{0}: not found . . . \nTry "{1} Dependencies --get {0}"'.format(name, sys.argv[0])]
            
            print(check_summary[-1])
        
        if len(check_summary) > 1:
            print(''.join(['\n\nSummary:\n'] + sorted(check_summary))+'\n')
        
        if not all(check_results):
            sys.exit('\n\nOne or more dependencies are unavailable . . . \n')
    
    if args.checkpackage:
        # call baga to do a conventional check on a python package
        result = checkpackage(args.checkpackage)
        if result:
            print('\n{}: found!'.format(args.checkpackage))
        else:
            print('\n{0}: not found . . . \nTry "{1} Dependencies --get {0}"'.format(args.checkpackage, sys.argv[0]))
            
    
    if args.checkget:
        check_summary,check_results = checkget(args.checkget)
        
        if len(check_summary):
            print(''.join(['\n\nSummary:\n'] + sorted(check_summary))+'\n')
        
        if not all(check_results):
            sys.exit('\n\nOne or more baga dependencies are unavailable and could not be installed . . . \n')
    
    if args.checkgetfor:
        for task in args.checkgetfor: 
            checkthese = sorted(Dependencies.dependencies_by_task[task])
            print('Checking on dependencies for {} ({})'.format(task, ', '.join(checkthese)))
            check_summary,check_results = checkget(checkthese)
            
            if len(check_summary):
                print(''.join(['\n\nSummary for {}:\n'.format(task)] + sorted(check_summary))+'\n')
            
            if not all(check_results):
                sys.exit('\n\nOne or more baga dependencies are unavailable and could not be installed . . . \n')

### Download Genomes ###

if task_name == 'CollectData':
    baga_cli.info('\t-- Data collection task --')
    from baga import CollectData
    # configure logger for this Task
    task_logger, task_log_folder = configureLogger(use_sample_name, 
            main_log_filename, verbosities[args.verbosity], 
            logger_name = task_name)
    ### Download Genomes ###
    if args.genomes:
        # first check if we're downloading from a list in a file
        if len(args.genomes) == 1 and os.path.exists(args.genomes[0]) and \
                not re.search('\.gbk|\.GBK|\.gbff|\.GBFF', args.genomes[0]):
            baga_cli.info('Found {} - will attempt to parse for search terms and '\
                    'URLs'.format(args.genomes[0]))
            raw = open(args.genomes[0]).readlines()
            accs_urls_dict = {}
            for n,line in enumerate(raw):
                k_v = re.split('[\t ]+', line.rstrip())
                if len(k_v) == 2:
                    accs_urls_dict[k_v[0]] = k_v[1]
                elif len(k_v) == 1:
                    accs_urls_dict[k_v[0]] = ''
                else:
                    baga_cli.error('Failed to parse {} at line {}: {}'.format(
                            args.genomes[0], n+1, line))
                    sys.exit(1)
            baga_cli.log(PROGRESS, 'Found {} genomes to download'.format(
                    len(accs_urls_dict)))
            try:
                genome = CollectData.Genome(task_name = task_name, 
                        console_verbosity_lvl = verbosities[args.verbosity],
                        log_folder = task_log_folder)
            except ValueError as e:
                baga_cli.error(e)
                sys.exit(1)
            downloaded = genome.downloadFromList(accs_urls_dict, force = args.force, 
                    retain_gbk = args.keep_gbk, user_email = args.email_address, 
                    path = args.analysis_path)
            if set(downloaded) != set(accs_urls_dict):
                failed = sorted(set(accs_urls_dict) - set(downloaded))
                baga_cli.error('Failed to download these: {}. See logs for '\
                        'details in {}'.format(', '.join(failed), 
                        task_log_folder))
                sys.exit(1)
        else:
            # work out what accessions and/or paths to genbank files were provided
            if ',' in args.genomes:
                # allow for single string with ',' delimiter
                use_genome = args.genomes.split(',')
                ### what if its a single single and not a list?
            else:
                use_genome = args.genomes
            load_gbks = []
            search_terms = []
            for g in use_genome:
                genome_no_quotes = g.strip('"').strip("'")
                if re.search('\.gbk|\.GBK|\.gbff|\.GBFF', genome_no_quotes):
                    load_gbks += [genome_no_quotes]
                    task_logger.debug("Will attempt to load genbank file: {}"\
                            "".format(genome_no_quotes))
                else:
                    # anything that isn't *.gbk or *.gbff is assumed to be an
                    # accession number
                    search_terms += [genome_no_quotes]
                    task_logger.debug("Will query NCBI Assembly database with: {}"\
                            "".format(genome_no_quotes))
            
            # check email address provided for Entrez, if any accession found
            if args.email_address is None and len(search_terms) > 0:
                task_logger.error('User email address is required for '\
                'downloading from NCBI with accession numbers. Detected {} items '\
                'without ".gbk" or ".GBK" and assumed to be search terms'\
                ''.format(len(search_terms)))
                sys.exit(1)
            
            loaded_genomes = {}
            genome = CollectData.Genome(task_name = task_name, 
                    console_verbosity_lvl = verbosities[args.verbosity],
                    log_folder = task_log_folder)
            for gbk in load_gbks:
                task_logger.info('Loading {}'.format(gbk))
                genome.loadFromGBK(gbk)
                # attempt to find a sensible sample_name
                patt = '(GC[FA]_[0-9]{9}\.[0-9]{1,2})_'
                m = re.match(patt,gbk)
                if m is not None and len(m.groups()) == 1:
                    task_logger.debug('Matched {} in {}'.format(patt,gbk))
                    # Assembly accession in filename
                    sample_name = m.groups()[0]
                else:
                    task_logger.debug('Did not match {} in {}'.format(patt,gbk))
                    # else just filename without a conventional extensions
                    # if present
                    sample_name = gbk.replace('.gbk','').replace('.gbff','')
                
                task_logger.debug('will use sample name: {}'.format(sample_name))
                genome.sample_name = sample_name
                use_filename = 'baga.{}.Genome-{}.baga'.format(task_name, sample_name)
                task_logger.debug('will use filename: {}'.format(use_filename))
                genome.file_name = use_filename
                genome.source = gbk
                genome.saveLocal()
                task_logger.info('Saved from {} to {}'.format(gbk, use_filename))
            
            download_urls = {}
            for search_term in search_terms:
                task_logger.info('Searching with: {}'.format(search_term))
                genome.queryEntrezAssembly(search_term, args.email_address)
                include_complete = True
                include_scaffolds = True
                include_contigs = True
                if args.complete_only:
                    include_scaffolds = False
                    include_contigs = False
                elif args.draft_only:
                    include_complete = False
                
                total = len(genome.assemblies_info)
                # download results of query
                filenames = []
                genome_identifiers = []
                not_downloaded = []
                c = 0
                try:
                    while c < total:
                        c += 1
                        task_logger.info('Fetching {} of {}'.format(c, total))
                        res = genome.downloadNextAssembly(
                                include_complete = include_complete, 
                                include_scaffolds = include_scaffolds, 
                                include_contigs = include_contigs, 
                                refseq_only = args.refseq_only, 
                                force = args.force, 
                                retain_gbk = args.keep_gbk,
                                path = '.') # args.analysis_path
                        if genome.file_name:
                            # False if not DLd
                            genome_identifiers += [res]
                            filenames += [genome.file_name]
                            download_urls[genome.sample_name] = genome.source
                            genome.saveLocal(exclude = ["assemblies_info", 
                                    "assemblies_problems"])
                            task_logger.info('saved genome to: {}'.format(
                                    genome.file_name))
                            task_logger.info('IMPORTANT: use "--genome_name {}" to '\
                                    'use this genome for baga analyses'.format(
                                    genome.sample_name))
                            # no option to save to other folder yet
                            # if args.analysis_path != '.':
                                # new_path = os.path.sep.join([args.analysis_path,
                                        # genome.file_name])
                                # genome.logger.log(PROGRESS, 
                                        # 'Moving to new path: {}'.format(new_path))
                                # os.rename(genome.file_name, new_path)
                        else:
                            not_downloaded += [res]
                except StopIteration:
                    task_logger.info('Completed downloads for search term: {}. {} '\
                            'downloaded, {} not downloaded'.format(search_term, 
                            len(filenames), len(not_downloaded)))
                    pass
            
            if len(search_terms):
                # save a record of downloaded genomes
                filename = '{}{}genome_accessions_downloaded.txt'.format(
                        task_log_folder, os.path.sep)
                with open(filename,'w') as fout:
                    for accession,url in sorted(download_urls.items()):
                        fout.write('{}\t{}\n'.format(accession,url))
                    task_logger.info('Put a list of accessions with download '\
                            'URLs in {}'.format(filename))
    ### old genomes ###
    # if args.genomes is not None:
        # from baga import CollectData
        # # work out what accessions and/or paths to genbank files were provided
        # if ',' in args.genomes:
            # # allow for single string with ',' delimiter
            # use_genome = args.genomes.split(',')
        # else:
            # use_genome = args.genomes
        
        # load_gbks = []
        # load_bagas = []
        # collect_accessions = []
        # for g in use_genome:
            # genome_no_quotes = g.strip('"').strip("'")
            # if re.search('\.gbk|\.GBK', genome_no_quotes):
                # load_gbks += [genome_no_quotes]
            # elif re.search('\.baga|\.baga', genome_no_quotes):
                # load_bagas += [genome_no_quotes]
            # else:
                # # anything that isn't .gbk or .baga is assumed to be an accession number
                # collect_accessions += [genome_no_quotes]
        
        # # check email address provided for Entrez, if any accession found
        # if args.email_address is None and len(collect_accessions) > 0:
            # print(textwrap.fill('User email address is required for downloading from NCBI \
# with accession numbers. Detected %s items without ".gbk" or ".GBK" and assumed to be \
# accession numbers' % (len(collect_accessions)), text_width))
            # sys.exit(1)
        
        # loaded_genomes = {}
        
        # for gbk in load_gbks:
            # print('Loading {}'.format(gbk))
            # genome = CollectData.Genome(local_path = gbk, format = 'genbank')
            # print('Storing for future use . . .')
            # genome.saveLocal()
            # #loaded_genomes[genome.id] = genome
        
        # for baga in load_bagas:
            # print('Loading {}'.format(baga))
            # genome = CollectData.Genome(local_path = gbk, format = 'baga')
            # print('Storing for future use . . .')
            # genome.saveLocal()
            # #loaded_genomes[genome.id] = genome
        
        # for genome_accession in collect_accessions:
            # print('Fetching {}'.format(genome_accession))
            # genome = CollectData.Genome(accession = genome_accession, user_email = args.email_address)
            # print('Storing for future use . . .')
            # genome.saveLocal()
            #loaded_genomes[genome.id] = genome
        
        # # download genomes from NCBI
        # if len(collect_accessions) > 0:
            # print('Found accessions for collection from NCBI:\n%s' % '\n'.join(collect_accessions))
            # for genome_accession in collect_accessions:
                # genome = CollectData.Genome()
                # print('Downloading: %s from NCBI' % genome_accession)
                # genome.getFromNCBI(genome_accession, args.email_address)
                # loaded_genomes[genome.genome_genbank_record.name.replace(' ','_')] = genome
        
        # # extract ORF IDs and chromosome ranges
        # for name, loaded_genome in loaded_genomes.items():
            # print('Extracting loci from %s' % name)
            # loaded_genome.extractLoci()
            # loaded_genome.saveLocal(name)

# no option to not extractLoci for now

def sanitize_filename(proposed_name):
    invalidChars = '!"#$%&\'()*+,/:;<=>?@[\\]^`{|}~ '
    sanitised = set()
    use_name = []
    for l in proposed_name:
        if l in invalidChars:
            use_name += ['_']
            sanitised.add(l)
        else:
            use_name += [l]
    
    use_name = ''.join(use_name)
    return(use_name, sanitised)


### (down)load Reads ###

if args.subparser == 'CollectData':
    if args.reads_download is not None:
        from baga import CollectData
        # only download so only accession numbers?
        if ',' in args.reads_download:
            # allow for single string with ',' delimiter
            use_reads = args.reads_download.split(',')
        else:
            use_reads = args.reads_download
        
        read_accessions = [r.strip('"').strip("'") for r in use_reads]
        read_accessions.sort()
        reads = CollectData.Reads()
        reads.getFromENA(read_accessions)
        # make a name from run accessions . . .
        if args.reads_group_name is None:
            use_reads_group_name = read_accessions[0]+'plus%sothers' % (len(read_accessions) - 1)
            print('No reads group name provided. Using: {}'.format(use_reads_group_name))
        else:
            use_reads_group_name, sanitised = sanitize_filename(args.reads_group_name)
            
            if len(sanitised) > 0:
                print('\nWarning: replaced {} with _ in provided reads group name.'.format(', '.join(sanitised)))
                print('Provided: {}'.format(args.reads_group_name))
                print('Using: {}\n\n'.format(use_reads_group_name))
            
        reads.saveLocal(use_reads_group_name)

### Load Reads from path ###

if args.subparser == 'CollectData':
    if args.reads_path is not None:
        # check read group name provided
        if args.reads_group_name is None:
            print(textwrap.fill('--reads_group_name is required for loading reads from a local path', text_width))
            sys.exit(1)
        
        use_reads_path = []
        for path in args.reads_path:
            use_reads_path += [path.strip('"').strip("'")]
        
        use_reads_group_name, sanitised = sanitize_filename(args.reads_group_name)
        
        if len(sanitised) > 0:
            print('\nWarning: replaced {} with _ in provided reads group name.'.format(', '.join(sanitised)))
            print('Provided: {}'.format(args.reads_group_name))
            print('Using: {}\n\n'.format(use_reads_group_name))
        
        from baga import CollectData
        reads = CollectData.Reads()
        reads.getFromPath(use_reads_path)
        reads.saveLocal(use_reads_group_name)

### Prepare Reads ###

if args.subparser == 'PrepareReads':
    print('\n-- Short Reads Preparation module --')
    import baga
    if args.reads_name is not None:
        
        use_path_reads,use_name_reads = check_baga_path('baga.CollectData.Reads', args.reads_name)
        e = 'Could not locate a saved baga.CollectData.Reads-<reads_name>.baga for reads group given: {}'.format(args.reads_name)
        assert all([use_path_reads,use_name_reads]), e
        
        from baga import PrepareReads
        
        if args.subsample_to_cov is not None:
            # make a new CallVariants.Reads object
            print('Loading reads group %s' % use_name_reads)
            downloaded_reads = baga.bagaload('baga.CollectData.Reads-{}'.format(use_name_reads))
            reads = PrepareReads.Reads(downloaded_reads)
            read_cov_depth, genome_size = args.subsample_to_cov
            print('Subsampling reads group {} to {}x coverage for a {:,} bp genome'.format(use_name_reads, read_cov_depth, genome_size))
            reads.subsample(genome_size, read_cov_depth, force = args.force)
            reads.saveLocal(use_name_reads)
        
        if args.adaptors is not None:
            # need to test parallelism with lots of stdout reports
            if args.subsample_to_cov is None:
                # didn't subsample in current analysis
                if args.adaptors == 'fullsample':
                    print('Not using previously subsampled reads because "--adaptors fullsample" supplied.')
                    print('Loading collected reads group %s' % use_name_reads)
                    downloaded_reads = baga.bagaload('baga.CollectData.Reads-%s' % use_name_reads)
                    reads = PrepareReads.Reads(downloaded_reads)
                else:
                    print('Loading subsampled reads group %s' % use_name_reads)
                    reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name_reads)
            
            try:
                reads.cutAdaptors(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('cutadapt')
            
            reads.saveLocal(use_name_reads)
        
        if args.trim:
            # could check whether adaptor cut read files exist
            # need to test parallelism with lots of stdout reports
            if not args.adaptors:
                # load a previously adaptor cut reads set
                print('Loading processed reads group %s' % use_name_reads)
                reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name_reads)
            
            print('\nTrimming reads . . .')
            try:
                reads.trim(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('sickle')
            
            reads.saveLocal(use_name_reads)
        
        if args.delete_intermediates:
            print('Checking on intermediate fastq files to delete . . .')
            if not args.adaptors and not args.trim:
                # load a previously adaptor cut reads set
                print('Loading processed reads group %s' % use_name_reads)
                reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name_reads)
            
            total_size = 0
            # check stage 1 files and stage 2 files
            if hasattr(reads,'read_files') and hasattr(reads,'adaptorcut_read_files'):
                stage1s = check_files(reads.read_files)
                stage2s = check_files(reads.adaptorcut_read_files)
                if stage2s:
                    if stage1s:
                        # delete stage 1 files if have all <== not for now . . .
                        # print('Deleting original or subsampled fastq files . . .')
                        # total_size += delete_files(reads.read_files, extra = '_subsmp')
                        print('Retaining original fastq files even though processed versions exist because re-downloading is time consuming!')
                    else:
                        print('Some or all of original or subsampled fastq files seem to have been deleted')
                else:
                    print('Missing some cutadapt-processed files: not deleting originals or subsampled')
            
            if hasattr(reads,'adaptorcut_read_files') and hasattr(reads,'trimmed_read_files'):
                stage2s = check_files(reads.adaptorcut_read_files)
                stage3s = check_files(reads.trimmed_read_files)
                
                if stage3s:
                    if stage2s:
                        # delete stage 2 files if have all
                        print('Deleting cutadapt-processed fastq files . . .')
                        total_size += delete_files(reads.adaptorcut_read_files)
                    else:
                        print('Some or all of cutadapt-processed fastq files seem to have been deleted')
                else:
                    print('Missing some sickle-processed files: not deleting cutadapt-processed')
            
            if total_size:
                print('Saved {:.2f} Gb by deleting intermediate files'.format(total_size/1000000000.0))
            else:
                print('Nothing deleted.')

### Align Reads ###

if task_name == 'AlignReads':
    print('\n-- Read Aligning module --')
    # configure logger for this Task
    task_logger, task_log_folder = configureLogger(use_sample_name, main_log_filename, 
            verbosities[args.verbosity], logger_name = task_name)
    if args.reads_name is not None and args.genome_name is not None:
        # first check whether GATK path is needed
        if args.indelrealign and not args.GATK_jar_path:
            
            print('''Please supply:

--GATK_jar_path

if using:

--indelrealign

''')
            sys.exit(1)
        
        # ensure upstream files are available: genome
        use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
        e = 'Could not locate a saved baga.CollectData.Genome-<genome_name>.baga for name given: {}'.format(args.genome_name)
        assert all([use_path_genome,use_name_genome]), e
        
        import baga
        from baga import AlignReads
        from baga import CollectData
        # ensure upstream files are available: reads
        if args.prepared:
            # not prepared by BAGA, load direct from a CollectData.Reads file
            use_path_reads,use_name_reads = check_baga_path(
                    'baga.CollectData.Reads', args.reads_name)
            e = 'Could not locate a saved baga.CollectData.Reads-<reads_name>.baga '\
                    'for reads group given: {}'.format(args.reads_name)
            assert all([use_path_reads,use_name_reads]), e
            from baga import PrepareReads
            print('Loading reads group %s' % use_name_reads)
            downloaded_reads = baga.bagaload('baga.CollectData.Reads-{}'\
                    ''.format(use_name_reads))
            reads = PrepareReads.Reads(downloaded_reads)
            # generate a PreparedReads.Reads file for use below
            # reads already trimmed etc
            reads.trimmed_read_files = reads.read_files
            reads.saveLocal(use_name_reads)
        
        # BAGA prepared or just generated
        use_path_reads,use_name_reads = check_baga_path(
                'baga.PrepareReads.Reads', args.reads_name)
        e = 'Could not locate a saved baga.PrepareReads.Reads-<reads_name>.baga '\
                'for reads group given: {}'.format(args.reads_name)
        assert all([use_path_reads,use_name_reads]), e
        alns_name = '__'.join([use_name_reads, use_name_genome])
        
        if args.align:
            print('Loading processed reads group %s' % use_name_reads)
            prepared_reads = baga.bagaload(use_path_reads)
            print('Loading genome %s' % use_name_genome)
            genome = CollectData.Genome(sample_name = use_name_genome, 
                    inherit_from = 'self')
            # create alignment object
            alignments = AlignReads.SAMs(sample_name = alns_name, 
                    reads = prepared_reads, genome = genome,
                    task_name = task_name, 
                    console_verbosity_lvl = verbosities[args.verbosity],
                    log_folder = task_log_folder)  # use inherit from here?, inherit_from = upstreams)
            print('\nAligning reads . . .')
            # let BWA estimate insert size and assign proper_pairs
            try:
                alignments.align(max_cpus = args.max_cpus, force = args.force)
            except OSError:
                exe_fail('bwa')
            
            try:
                alignments.toBAMs(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('samtools')
            
            # need to include the genome name for aligning a group of reads sets to more than one genome
            alignments.saveLocal()
        
        
        if args.deduplicate:
            if not args.align:
                # add an exception here and inform to use --align first
                print('Loading previously processed read alignments: %s' % alns_name)
                alignments = AlignReads.SAMs(sample_name = alns_name, 
                        task_name = task_name, 
                        console_verbosity_lvl = verbosities[args.verbosity],
                        log_folder = task_log_folder, inherit_from = 'self')
                
            try:
                alignments.removeDuplicates(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('picard')
            
            try:
                alignments.sortIndexBAMs(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('samtools')
            
            alignments.saveLocal()
        
        if args.indelrealign:
            if not args.deduplicate:
                # add an exception here and inform to use --align first
                print('Loading previously processed read alignments: %s' % alns_name)
                alignments = AlignReads.SAMs(sample_name = alns_name, 
                        task_name = task_name, 
                        console_verbosity_lvl = verbosities[args.verbosity],
                        log_folder = task_log_folder, inherit_from = 'self')
            
            try:
                os.makedirs('genome_sequences')
            except OSError:
                pass
            
            # seems not to save .genome in alignments?
            genome_fna = 'genome_sequences/%s.fna' % alignments.genome_name
            #genome_fna = 'genome_sequences/{}.fna'.format(use_name_genome)
            
            if not os.path.exists(genome_fna):
                from Bio import SeqIO
                from Bio.Seq import Seq
                from Bio.SeqRecord import SeqRecord
                records_to_write = []
                for seq_id,seq_array in alignments.genome_sequence.items():
                    records_to_write += [SeqRecord(Seq(seq_array.tostring()), 
                            id = seq_id, description = self.genome_names[seq_id])]
                
                SeqIO.write(records_to_write, genome_fna, 'fasta')
            
            alignments.IndelRealignGATK(
                        jar = args.GATK_jar_path.split(os.path.sep), 
                        use_java = use_java,
                        force = args.force, 
                        max_cpus = args.max_cpus)
            
            alignments.saveLocal()
        
        if args.delete_intermediates:
            print('Checking on intermediate alignment files to delete . . .')
            if not any([args.indelrealign, args.deduplicate, args.align]):
                # load a previously adaptor cut reads set
                print('Loading previously processed read alignments: %s' % alns_name)
                alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}.baga'.format(alns_name))
            
            attrs = ('aligned_read_files', 'paths_to_BAMs', 'paths_to_BAMs_dd', 'paths_to_BAMs_dd_si', 'ready_BAMs')
            descr = {   'aligned_read_files':'initial SAMs', 
                        'paths_to_BAMs':'initial BAMs', 
                        'paths_to_BAMs_dd':'deduplicated BAMs', 
                        'paths_to_BAMs_dd_si':'sorted, deduplicated BAMs', 
                        'ready_BAMs':'indel-realigned BAMs'}
            
            total_size = 0
            for attr1,attr2 in zip(attrs[:-1],attrs[1:]):
                if hasattr(alignments,attr1) and hasattr(alignments,attr2):
                    have1, have2 = False, False
                    if attr1 == 'aligned_read_files':
                        files_to_check1 = alignments.__dict__[attr1].values()
                    else:
                        files_to_check1 = alignments.__dict__[attr1]
                    
                    if attr2 == 'ready_BAMs':
                        files_to_check2 = alignments.__dict__[attr2][0]
                    else:
                        files_to_check2 = alignments.__dict__[attr2]
                    
                    have1 = check_files(files_to_check1)
                    have2 = check_files(files_to_check2)
                    
                    if have2:
                        if have1:
                            # delete stage 1 files if have all
                            print('Deleting {} . . .'.format(descr[attr1]))
                            total_size += delete_files(files_to_check1)
                        else:
                            print('Some or all of {} seem to have been deleted'.format(descr[attr1]))
                    else:
                        print('Missing some {}: not deleting {}'.format(descr[attr2],descr[attr1]))
            
            if total_size:
                print('Saved {:.2f} Gb by deleting intermediate files'.format(total_size/1000000000.0))
            else:
                print('Nothing deleted.')

if args.subparser == 'SimulateReads':
    print('\n-- Read Simulation module --')
    use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
    e = 'Could not locate a saved baga.CollectData.Genome-<genome_name>.baga for name given: {}'.format(args.genome_name)
    import baga
    from baga import CollectData
    from baga import SimulateReads
    print('Loading genome %s' % use_name_genome)
    genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
    
    large_deletions = {}
    x = 1
    for i in range(len(args.large_deletions))[::2]:
        large_deletions['Deletion_{}'.format(x)] = tuple(args.large_deletions[i:i+2])
        x += 1
    
    simulator = SimulateReads.Simulator(genome = genome, 
                   num_individuals = args.num_individuals,
                   large_deletions = large_deletions,
                   random_seed = args.random_seed)
    
    simulator.do(num_SNPs = args.num_SNPs, num_deletions = args.num_deletions,
            num_insertions = args.num_insertions)
    
    if args.gemsim:
        simulator.generateReads(max_cpus = args.max_cpus)


### Repeats ###

if task_name == 'Repeats':
    print('\n-- Chromosome repeats detection module --')
    # configure logger for this Task
    task_logger, task_log_folder = configureLogger(use_sample_name, 
            main_log_filename, verbosities[args.verbosity], 
            logger_name = task_name)
    
    e = '-i/--minimum_percent_identity must be between 0 and 100 percent '\
            '(low values not recommended!)'
    assert 0 < args.minimum_percent_identity <= 100, e
    
    import baga
    from baga import Repeats
    from baga import CollectData
    
    use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', 
            args.genome_name)
    e = 'Could not locate a saved baga.CollectData.Genome-<genome_name>.baga '\
            'for name given: {}'.format(args.genome_name)
    assert all([use_path_genome,use_name_genome]), e
    
    if args.find:
        print('Loading genome %s' % use_name_genome)
        genome = CollectData.Genome(sample_name = use_name_genome, 
                inherit_from = 'self')
        finder = Repeats.Finder(genome = genome,
                task_name = task_name, 
                console_verbosity_lvl = verbosities[args.verbosity],
                log_folder = task_log_folder)
        
        if args.method == 'baga':
            # minimum_percent_identity defaults to 98%, argument takes 0.98 so *0.01
            # minimum_repeat_length defaults to 400
            finder.findRepeats(
                    minimum_percent_identity = args.minimum_percent_identity * 0.01, 
                    minimum_repeat_length = args.minimum_repeat_length,
                    max_extensions = 25)
            # save to file
            finder.saveLocal(serialiser = 'pickle')
            # also save just the ranges for filtering <== check if this is still used by FilterVariants etc
            # might just load full Finder object
            baga.bagasave(finder.ambiguous_ranges, 'baga.Repeats.filter_regions-{}'.format(use_name_genome))
        elif args.method == 'nucmer_check':
            finder.findRepeatsNucmer(minimum_percent_identity = args.minimum_percent_identity * 0.01, 
                          minimum_repeat_length = args.minimum_repeat_length)
            
            finder.compareRepeatRegions()
    
    if args.plot:
        # if not args.find:
            # print('Loading repeats and genome: %s' % use_name_genome)
            # finder = baga.bagaload('baga.Repeats.Finder-%s' % use_name_genome)
        
        ## this bit needs reworking to work properly with the pattern:
        # object = _MetaSample(genome = genome,
                # task_name = task_name, 
                # console_verbosity_lvl = verbosities[args.verbosity],
                # log_folder = task_log_folder)
        # ==> how is verbosity inherited for restored objects?
        # need debug levels and same logger (?) in those cases . . .
        # specifically here --find always instantiates a new one, --plot always loads what --find made
        # [but a Repeat module-level function is called, not a Plotter instantiated <== this is target]
        # while --find always loads a Genome created by CollectData --genome
        # (is logging etc transferred there or is there separation between restoring data (state)
        # and restoring a complete object with its own logger etc . . . ?)
        # e.g., should a restored object's methods ever be called?
        # Perhaps pipeline restore points should only be after a new instantiation plus method-based analysis
        # . . . restoring later is for data only?
        
        Repeats.plotRepeats(use_name_genome, outdir = ['plots_repeats'], force = args.force)
    
    if args.summarise:
        Repeats.summariseRepeats(use_name_genome)

### Structure ###
if task_name == 'Structure':
    print('\n-- Chromosome sequence rearrangement detection module --')
    # configure logger for this Task
    task_logger, task_log_folder = configureLogger(use_sample_name, 
            main_log_filename, verbosities[args.verbosity], 
            logger_name = task_name)
    import baga
    from baga import Structure
    from baga import CollectData
    
    if not(args.check or args.plot or args.plot_range or args.summarise or args.collect):
        parser_Structure.error('Need at least one of --check/-c, --plot/-p or --plot_range/-r or --summarise/-s or --collect/-C')
    
    if args.check or args.plot or args.plot_range or args.collect or args.collect_ranges:
        # collect BAMs
        if args.reads_name:
            # baga pipeline information provided
            if not args.genome_name:
                parser.error('--genome_name/-g is required with --reads_name/-n. (The baga CollectData-processed genome used with the AlignReads option)')
            
            use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
            e = 'Could not locate a saved baga.CollectData.Genome-<genome_name>.baga for name given: {}'.format(args.genome_name)
            assert all([use_path_genome,use_name_genome]), e
            
            # in case full filename provided
            use_name_group = args.reads_name.replace('baga.AlignReads.SAMs-', '' , 1).replace('.p.gz', '').replace('.baga', '')
            
            print('Loading alignments information for: {}__{} from AlignReads output'.format(use_name_group, use_name_genome))
            
            from baga import AlignReads
            alns_name = '{}__{}'.format(use_name_group, use_name_genome)
            alignments = AlignReads.SAMs(sample_name = alns_name, 
                    inherit_from = 'self')
            
            e = 'the reads for "--reads_name/-n {}" seem to not have been '\
                    'fully processed by the AlignReads module: they are '\
                    'missing the "ready_BAMs" attribute. Please ensure the '\
                    'AlignReads commands "--align --deduplicate --indelrealign" '\
                    'have been performed.'.format(args.reads_name)
            assert hasattr(alignments, 'ready_BAMs'), e
            ### shouldn't all BAMs have headers parsed and stored in dict with (sample,genome) from the start?
            BAMs = alignments.ready_BAMs[-1]
            sample_names = sorted(alignments.read_files)
            
        elif args.alignments_paths:
            # list of folders or files provided
            BAMs = []
            for path in args.alignments_paths:
                if os.path.isdir(path):
                    path_contents = os.listdir(path)
                    theseBAMs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('BAM', 'bam')]
                    e = 'No BAM files (*.bam or *.BAM) found in:\n{}'.format(args.alignments_paths)
                    assert len(theseBAMs), e
                    BAMs += theseBAMs
                else:
                    # add file
                    BAMs += [path]
        
        # check on requested samples: crop BAMs accordingly
        if args.include_samples or args.exclude_samples:
            if not args.reads_name:
                # need sample names so will parse BAMs
                import pysam
                sample_names = {}
                for BAM in BAMs:
                    sample_names[pysam.Samfile(BAM, 'rb').header['RG'][0]['ID']] = BAM
            
            found_labels = sorted(sample_names)
            if args.include_samples:
                missing_labels = sorted(set(args.include_samples) - set(sample_names))
                found_labels = sorted(set(args.include_samples) & set(sample_names))
                e = ['None of the requested sample labels were found among the '\
                        'previously checked reads.']
                e += ['Requested: {}'.format(', '.join(args.include_samples))]
                e += ['Available: {}'.format(', '.join(sorted(sample_names)))]
                assert len(found_labels), '\n'.join(e)
                print('Found {} samples to use after --include_samples: {}'\
                        ''.format(len(found_labels), ', '.join(found_labels)))
                
                if len(missing_labels):
                    print('WARNING: could not find the following requested '\
                            'samples among previously checked reads.')
                    print(', '.join(missing_labels))
                
            if args.exclude_samples:
                print('Found {} samples to use, --exclude_samples provides {} '\
                        'to remove.'.format(len(found_labels), 
                        len(args.exclude_samples)))
                found_labels = sorted(set(found_labels) - set(args.exclude_samples))
                print('{} samples remain to analyse: {}'.format(len(found_labels), 
                        ', '.join(found_labels)))
                
            # update BAMs for args.check
            #print(BAMs)
            if not args.reads_name:
                BAMs = [sample_names[sample] for sample in found_labels]
            else:
                BAMs = [BAM for BAM in BAMs if \
                        BAM.split(os.path.sep)[-1].split('__')[0] in found_labels]
            
            #print(BAMs)
            # update sample_names for arg.plot
            sample_names = sorted(found_labels)
        
        if args.check:
            # check these genome-aligned read sets
            checkers = Structure.checkStructure(BAMs, min_mapping_quality = 5, 
                    smoothed_resolution = 10, 
                    ratio_threshold = args.ratio_threshold, 
                    genome_name = use_name_genome, force = args.force, 
                    task_name = task_name, 
                    console_verbosity_lvl = verbosities[args.verbosity], 
                    log_folder = task_log_folder
                    )
            
            # regarding pipelined inheritance of saved states:
            # this would be a new style object . . BUT should be able to use inherit_from = 'AlignReads.Aligner'
            # that would make CLI much simpler . . . or should 'AlignReads.Aligner' be a 
            # default to be over-ridden by 'self'?
            
            # would also need to be more OO with a do function callable from here . .
            # each do() should have a good doc string describing the order and why
        
        if args.plot or args.plot_range:
            if args.plot_range:
                if args.plot_range[0] > args.plot_range[1]:
                    sys.exit('--plot_range values must be ascending!')
            # need genome for plotting
            if not args.genome_name:
                parser.error('--genome_name/-g is required for plotting. (A baga CollectData-processed genome)')
            else:
                use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
                print('Loading genome %s' % use_name_genome)
                genome = CollectData.Genome(local_path = 'baga.CollectData.Genome-{}.baga'.format(use_name_genome), format = 'baga')
            
            if not args.include_samples and not args.reads_name:
                # need to get sample names from BAMs because not supplied
                # and didn't collect yet because args.include_samples not provided
                try:
                    import pysam
                except ImportError:
                    sys.exit('Need pysam to get sample names if not provided for plotting. Use Dependencies --get pysam to install locally')
                
                sample_names = []
                for BAM in BAMs:
                    sample_names += [pysam.Samfile(BAM, 'rb').header['RG'][0]['ID']]
                
                sample_names.sort()
            
            plot_folder = 'plots_structure'
            filter_name = 'high non-proper pairs'
            
            for sample in sample_names:
                # get information for plotting
                filein = 'baga.Structure.CheckerInfo-{}__{}.baga'.format(sample, use_name_genome)
                
                try:
                    checker_info = Structure.loadCheckerInfo(filein)
                except IOError:
                    print('Could not find: {}'.format(filein))
                
                e = 'Genome name for checker {} ({}) does not match name of supplied genome ({})'.format(
                                                                            filein, 
                                                                            checker_info['genome_name'], 
                                                                            genome.id)
                assert checker_info['genome_name'] == genome.id, e
                
                outdir = os.path.sep.join([plot_folder, checker_info['genome_name']])
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                
                print('Plotting filter regions for {} reads aligned to {}'.format(sample, genome.id))
                
                if args.plot_range:
                    # just the requested range
                    do_ranges = [args.plot_range]
                    
                elif args.plot:
                    # all ranges for filtering for this sample
                    
                    # currently range selection is rearrangements filter, not including the extensions
                    # do_ranges = checker_info['suspect_region']['rearrangements_extended']
                    
                    do_ranges = []
                    for s,e in checker_info['suspect_regions']['rearrangements']:
                        # select consistant plotting region for comparison between samples
                        plot_chrom_start = int(round(s - 500 - 100, -3))
                        if s + 2500 > e:
                            plot_chrom_end = plot_chrom_start + 2500
                        else:
                            plot_chrom_end = int(round(e + 500 + 100, -3))
                        
                        do_ranges += [(plot_chrom_start,plot_chrom_end)]
                
                for plot_chrom_start,plot_chrom_end in do_ranges:
                    
                    if args.reads_name:
                        plot_filename = '{:07d}_{:07d}_{}__{}__{}.svg'.format(  plot_chrom_start, 
                                                                                plot_chrom_end, 
                                                                                checker_info['genome_name'], 
                                                                                use_name_group, 
                                                                                sample)
                    else:
                        plot_filename = '{:07d}_{:07d}_{}__{}.svg'.format(      plot_chrom_start, 
                                                                                plot_chrom_end, 
                                                                                checker_info['genome_name'],
                                                                                sample)
                    
                    plot_output_path = [outdir, plot_filename]
                    plot_output_path = os.path.sep.join(plot_output_path)
                    print(plot_output_path)
                    plotter = Structure.Plotter(checker_info, genome, plot_output_path)
                    plotter.doPlot(plot_chrom_start, plot_chrom_end, panel = ((1,1),(1,1)), label = sample)
    
    ## check for allowed combinations for --summarise
    if args.checkinfos_path:
        if not args.summarise:
            parser.error('--summarise/-s is required with --checkinfos_path/-b.')
        if args.genome_name or args.reads_name:
            parser.error('--genome_name/-g and --reads_name/-n cannot be used with --checkinfos_path/-b.')
    elif args.reads_name and not args.genome_name:
        parser.error('--genome_name/-g is required with --reads_name/-n. (The baga CollectData-processed genome used with the AlignReads option)')
    
    if args.summarise or args.collect or args.collect_ranges:
        # both tasks share some requirements: deal with these first
        
        if args.reads_name:
            # baga pipeline information provided
            if not args.genome_name:
                parser.error('--genome_name/-g is required with --reads_name/-n. (The baga CollectData-processed genome used with the AlignReads option)')
            
            use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
            assert all([use_path_genome,use_name_genome]), 'Could not locate genome given: {}'.format(args.genome_name)
            
            # in case full filename provided
            use_name_group = args.reads_name.replace('baga.AlignReads.SAMs-', '' , 1).replace('.p.gz', '').replace('.baga', '')
            
            ## this was already done above <== need to refactor this whole section
            baga_file = 'baga.AlignReads.SAMs-{}__{}.baga'.format(use_name_group, use_name_genome)
            print('Loading alignments information for: {} aligned to {} from {} output'.format(use_name_group, use_name_genome, baga_file))
            from baga import AlignReads
            alns_name = '{}__{}'.format(use_name_group, use_name_genome)
            alignments = AlignReads.SAMs(sample_name = alns_name, 
                    inherit_from = 'self')
            sample_names = sorted(alignments.read_files)
        
        if args.include_samples and args.reads_name:
            missing_labels = sorted(set(args.include_samples) - set(sample_names))
            found_labels = sorted(set(args.include_samples) & set(sample_names))
            e = ['None of the requested sample labels were found among the rearrangements filter-checked reads.']
            e += ['Requested: {}'.format(', '.join(args.include_samples))]
            e += ['Available: {}'.format(', '.join(sorted(sample_names)))]
            assert len(found_labels), '\n'.join(e)
            
            if len(missing_labels):
                print('WARNING: could not find the following requested samples among previously checked reads.')
                print(', '.join(missing_labels))
            
            # update sample_names
            sample_names = sorted(found_labels)
        
        if args.checkinfos_path:
            # treat include_samples as list of files to deal with
            # list of folders or files provided
            print(args.checkinfos_path)
            baga_filenames = {}
            for path in args.checkinfos_path:
                if os.path.isdir(path):
                    path_contents = os.listdir(path)
                    for f in path_contents:
                        if f[-5:] == '.baga' and f[:27] == 'baga.Structure.CheckerInfo-' and '__' in f[27:-5]:
                            baga_filenames[tuple(f[27:-5].split('__'))] = f
                    
                    e = 'No baga.Structure.CheckerInfo-*__*.baga files found in:\n{}'.format(args.checkinfos_path)
                    assert len(baga_filenames), e
                else:
                    # add file
                    f = path.split(os.path.sep)[-1]
                    if f[-5:] == '.baga' and f[:27] == 'baga.Structure.CheckerInfo-' and '__' in f[27:-5]:
                        baga_filenames[tuple(f[27:-5].split('__'))] = f
            
            e = 'Could not find valid baga.Structure.CheckerInfo files at {}'.format(', '.join(args.checkinfos_path))
            assert len(baga_filenames) > 0, e
            
        else:
            try:
                baga_filenames = dict([(tuple([sample, use_name_genome]),
                        'baga.Structure.CheckerInfo-{}__{}.baga'.format(sample, use_name_genome)) \
                        for sample in sample_names])
            except NameError:
                print('need --checkinfos_path or --reads_name')
                sys.exit(1)
        
        checker_info = {}
        for (sample, genome_name), filein in sorted(baga_filenames.items()):
            try:
                print('Loading from {}'.format(filein))
                filein.replace('baga.Structure.CheckerInfo-','')
                this_checker_info = Structure.loadCheckerInfo(filein)
            except IOError as e:
                print('Cannot access provided baga.Structure.CheckerInfo file: {}'.format(filein))
                print(e)
            
            e = 'Genome name for checker {} ({}) does not match supplied genome name ({})'.format(
                                                                        filein, 
                                                                        this_checker_info['genome_name'], 
                                                                        genome_name)
            assert this_checker_info['genome_name'] == genome_name, e
            # else proceed
            checker_info[sample, genome_name] = this_checker_info
        
        if args.summarise:
            #### summarise only
            if args.reads_name:
                if args.genome_name:
                    foutname = 'rearrangements_regions_{}__{}.csv'.format(use_name_group,use_name_genome)
                else:
                    foutname = 'rearrangements_regions_{}.csv'.format(use_name_group)
            else:
                foutname = 'rearrangements_regions.csv'
            
            print('Writing to {}'.format(foutname))
            with open(foutname, 'w') as fout:
                fout.write('"genome","chromosome","sample","filter","start","end"\n')
                for (sample, genome_name), info in sorted(checker_info.items()):
                    print('Writing out {} regions'.format(sample))
                    for replicon_id, ranges in info['suspect_regions']['rearrangements'].items():
                        for start, end in ranges:
                            fout.write('"{}","{}","rearrangements1",{},{}\n'\
                                    ''.format(genome_name, replicon_id, sample, start, end))
                        t = len(ranges)
                        s = sum([(e - s) for s,e in ranges])
                        print('  in replicon (or contig) {}, {} regions spanning {:,} basepairs '\
                                'are affected by rearrangements thus having ambiguous 1:1 '\
                                'orthology.'.format(replicon_id, t, s))
                    for replicon_id, ranges in info['suspect_regions']['rearrangements_extended'].items():
                        for start, end in ranges:
                            fout.write('"{}","{}","rearrangements2",{},{}\n'\
                                    ''.format(genome_name, replicon_id, sample, start, end))
                        t = len(ranges)
                        s = sum([(e - s) for s,e in ranges])
                        print('  in replicon (or contig) {}, {} additional regions spanning {:,} '\
                                'basepairs are adjacent to the above regions but have a >50% zero '\
                                'read depth over a moving window (i.e., any aligned reads have a '\
                                'patchy distribution, are usually rare). These are typically large '\
                                'deletions including missing prophage and genomic islands\n'\
                                ''.format(replicon_id, t, s))
        
        # probably better in a .do_collect() method in the module
        if args.collect or args.collect_ranges:
            use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
            assert all([use_path_genome,use_name_genome]), 'Could not locate genome given: {}'.format(args.genome_name)
            genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
            ## which BAMs?
            ## this is a bit clunky . . probably best to parse BAM headers once, above
            ## move into baga function from cli?
            import pysam
            BAMs_by_ids = {}
            for BAM in BAMs:
                header = pysam.Samfile(BAM, 'rb').header
                BAMs_by_ids[(header['RG'][0]['ID'],header['SQ'][0]['SN'])] = BAM
            
            from baga import AssembleReads
            if args.max_memory:
                use_mem_gigs = args.max_memory
            else:
                # round down available GBs
                use_mem_gigs = int(baga.get_available_memory())
            
            for (sample, genome_name), info in sorted(checker_info.items()):
                path_to_fastq_folder = os.path.sep.join(['read_collections', genome_name])
                if not os.path.exists(path_to_fastq_folder):
                    os.makedirs(path_to_fastq_folder)
                
                print('Extracting reads aligned near rearrangements between {} and genome {} . . .'.format(sample, genome_name))
                collector = Structure.Collector(BAMs_by_ids[(sample, genome_name)])
                
                e = 'mismatch between BAM genome ({}) and genome used by BAGA ({})'.format(
                                        collector.reads.references[0],
                                        genome_name)
                assert genome_name == collector.reads.references[0], e
                
                if args.collect_ranges:
                    # single assembly of reads aligned to one or more ranges in reference
                    # versus separate assemblies of multiple ranges i.e., those with putative rearrangements
                    single_assembly = True
                    # ensure pairs of start-end ranges given
                    e = 'Odd number of ranges provided. Required: start-end, '\
                    'start-end integers as --collect_ranges start end start end'
                    assert len(args.collect_ranges) % 2 == 0, e
                    e = 'Ranges must be non-overlapping and in ascending order'
                    assert sorted(args.collect_ranges) == args.collect_ranges, e
                    use_regions = zip(args.collect_ranges[::2],args.collect_ranges[1::2])
                    use_num_padding_positions = 0
                else:
                    single_assembly = False
                    # join main rearrangement zones with extended regions if found
                    # to get contiguous blocks for investigation
                    all_regions = sorted(
                                  info['suspect_regions']['rearrangements'] + \
                                  info['suspect_regions']['rearrangements_extended'])
                    
                    from collections import Counter
                    use_regions = [a for b in all_regions for a in b]
                    c = Counter(use_regions)
                    use_regions = [a for a in use_regions if c[a] == 1]
                    use_regions = zip(use_regions[::2],use_regions[1::2])
                    use_num_padding_positions = args.num_padding_positions
                
                collector.getUnmapped()
                r1_out_path_um, r2_out_path_um, rS_out_path_um = collector.writeUnmapped(path_to_fastq_folder)
                
                # assemble poorly/unmapped alone first
                reads_path_unmapped = {}
                output_folder_um = '_'.join(r1_out_path_um.split('_')[:-1]).split(os.path.sep)[-1]
                reads_path_unmapped[output_folder_um] = r1_out_path_um, r2_out_path_um, rS_out_path_um
                path_to_bad_unmapped_contigs = os.path.sep.join(['read_collections',
                                                    genome_name, 
                                                    output_folder_um,
                                                    'contigs.fasta'])
                if os.path.exists(path_to_bad_unmapped_contigs) and \
                   os.path.getsize(path_to_bad_unmapped_contigs) > 0 and \
                   not args.force:
                    print('Found assembly at {}\nUse --force/-F to overwrite. Skipping . . .'.format(path_to_bad_unmapped_contigs))
                else:
                    if not args.force:
                        print('Nothing found at {}. Doing assembly.'.format(path_to_bad_unmapped_contigs))
                    
                    reads = AssembleReads.DeNovo(paths_to_reads = reads_path_unmapped)
                    reads.SPAdes(output_folder = ['read_collections', genome_name], mem_num_gigs = use_mem_gigs,
                            only_assembler = True, careful = False)
                
                # assemble read from each region with poorly/unmapped
                reads_paths = {}
                # make a second dict of reads for assembly, all values for unmapped reads
                # that need to be included in each assembly
                reads_path_unmapped = {}
                assemblies_by_region = {}
                for (s,e) in use_regions:
                    collector.makeCollection(s, e, use_num_padding_positions)
                    r1_out_path, r2_out_path, rS_out_path = collector.writeCollection(path_to_fastq_folder)
                    if not r1_out_path:
                        # if no reads found, False returned
                        # do not add to reads_paths dict for assembly
                        print('debug: no reads found by collector')
                        continue
                    # put assembly in folder with same name as read files
                    output_folder = '_'.join(r1_out_path.split('_')[:-1]).split(os.path.sep)[-1]
                    path_to_contigs = os.path.sep.join(['read_collections', 
                                                        genome_name, 
                                                        output_folder, 
                                                        'contigs.fasta'])
                    assemblies_by_region[s,e] = path_to_contigs
                    if os.path.exists(path_to_contigs) and \
                            os.path.getsize(path_to_contigs) > 0 and \
                            not args.force:
                        print('Found assembly at {}\nUse --force/-F to overwrite. Skipping . . .'.format(path_to_contigs))
                    else:
                        reads_paths[output_folder] = (r1_out_path, r2_out_path, rS_out_path)
                        reads_path_unmapped[output_folder] = r1_out_path_um, r2_out_path_um, rS_out_path_um
                
                print('debug: len(assemblies_by_region) == {}'.format(len(assemblies_by_region)))
                reads = AssembleReads.DeNovo(paths_to_reads = reads_paths, 
                        paths_to_reads2 = reads_path_unmapped)
                reads.SPAdes(output_folder = ['read_collections', genome_name], 
                        mem_num_gigs = use_mem_gigs, single_assembly = single_assembly, 
                        only_assembler = True, careful = False)
                # a dict of paths to contigs per region
                aligner = Structure.Aligner(genome)
                unmappedfasta = os.path.sep.join(['read_collections', 
                                                  genome_name, 
                                                  output_folder_um, 
                                                  'contigs.fasta'])
                if os.path.exists(unmappedfasta) and os.path.getsize(unmappedfasta) > 0:
                    if len(assemblies_by_region) > 0:
                        # provide dict of range tuples
                        aligner.alignRegions(assemblies_by_region, 
                                use_num_padding_positions, 
                                path_to_omit_sequences = unmappedfasta, 
                                single_assembly = single_assembly,
                                min_region_length = args.min_align_region)
                        aligner.reportAlignments()
                    else:
                        print('WARNING: no assembled regions found. Either there are '\
                                'none which is fine, or SPAdes assemblies failed '\
                                'to finish. You could check SPAdes log files in '\
                                'folders in {}'.format(os.path.sep.join([
                                'read_collections', genome_name])))
                else:
                    print('WARNING: no assembled unmapped and poorly mapped reads found at:\n{}'.format(unmappedfasta))
                    try:
                        r1_size = os.path.getsize(r1_out_path_um)
                        r2_size = os.path.getsize(r2_out_path_um)
                        print('but reads, {} ({:,} bytes) and {} ({:,} bytes), exist . . check SPAdes assembly log in {}'.format(
                                                                r1_out_path_um,
                                                                r1_size,
                                                                r2_out_path_um,
                                                                r2_size,
                                                                unmappedfasta.replace('contigs.fasta','')))
                    except IOError:
                        print('WARNING: could not find unmapped and poorly '\
                        'aligned reads at:\n{}\n{}\nthis is unexpected but '\
                        'conceivable (if ALL reads really did map to reference!).'.format(
                                r1_out_path_um,r2_out_path_um))
                    print('proceeding with alignment of assembled putatively rearranged regions to reference nonetheless')
                    if len(assemblies_by_region) > 0:
                        aligner.alignRegions(assemblies_by_region, use_num_padding_positions, single_assembly = single_assembly)
                        aligner.reportAlignments()
                    else:
                        print('WARNING: no assembled regions found. Either there are '\
                                'none which is fine, or SPAdes assemblies failed '\
                                'to finish. You could check SPAdes log files in '\
                                'folders in {}'.format(os.path.sep.join([
                                'read_collections', genome_name])))

### Call Variants ###

if args.subparser == 'CallVariants':
    print('\n-- Variant Calling module --\n')
    # check whether GATK path is needed
    if any([args.callsingles,
            args.calleach, 
            args.calljoint, 
            args.hardfilter, 
            args.recalibrate
            ]):
        assert args.reads_name, '--reads_name is required for calling with GATK'
        if len(args.reads_name) != 1:
            sys.exit('Only one reads group can be processed by GATK per '\
            'analysis (supplied: {})'.format(', '.join(args.reads_name)))
        if args.calldisco:
            sys.exit('--calldisco cannot be used with any GATK options!')
        elif not args.GATK_jar_path:
            sys.exit('''Please supply:

--GATK_jar_path

if using any of:

--callsingles
--calleach
--calljoint
--hardfilter
--recalibrate

''')
        elif not args.genome_name:
            sys.exit('--genome_name is required for GATK')
    if args.genome_name:
        use_path_genome,use_name_genome = check_baga_path(
                'baga.CollectData.Genome', args.genome_name)
        assert all([use_path_genome,use_name_genome]), 'Could not locate a '\
                'saved baga.CollectData.Genome-<genome_name>.baga for name '\
                'given: {}'.format(args.genome_name)
        from baga import CollectData
        genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
    
    import baga
    from baga import CallVariants
    if args.check:
        assert (args.vcfs_paths and args.genome_name and args.alignments_paths), ''\
        '--check needs --vcfs_paths to know which variants to check for and '\
        '--genome_name to align contigs and --alignments_paths for the BAM files '\
        'against.'
        def collectfiles(paths, file_extensions):
            files = []
            for path in paths:
                if os.path.isdir(path):
                    path_contents = os.listdir(path)
                    thesefiles = [os.path.sep.join([path,f]) for f in path_contents if \
                                        f.split(_os.path.extsep)[-1] in file_extensions]
                    glob_extensions = ['*.{}'.format(e) for e in file_extensions]
                    assert len(thesefiles), 'No files ({}) found in:\n{}'\
                    ''.format(' or '.join(glob_extensions), args.alignments_paths)
                    files += thesefiles
                else:
                    # add file
                    files += [path]
            return(files)
        BAMs = collectfiles(args.alignments_paths, ('BAM', 'bam'))
        VCFs = collectfiles(args.vcfs_paths, ('VCF', 'vcf'))
        checker = CallVariants.Checker(VCFs, BAMs, genome)
        checker.doCheck(num_padding = 1000, max_memory = args.max_memory, 
                force = args.force)
    else:
        if args.calldisco:
            assert args.reads_name, '--reads_name is required for calling with DiscoSNP++'
            # check if direct commands make sense
            direct_arguments = check_direct_arguments(args.arguments, 
                    wrapped_tools = ['DiscoSNP++'])
            # allows for multiple tools to have direct arguments processed together
            try:
                use_arguments = direct_arguments['DiscoSNP++']
            except KeyError:
                use_arguments = False
            # load baga reads
            from baga import PrepareReads
            use_these = []
            use_these_names = []
            for reads_group in args.reads_name:
                use_path_reads,use_name_reads = check_baga_path(
                        'baga.PrepareReads.Reads', reads_group)
                use_these_names += [use_name_reads]
                assert all([use_path_reads,use_name_reads]), 'Could not locate a saved '\
                        'baga.PrepareReads.Reads-<reads_name>.baga for reads group '\
                        'given: {}'.format(reads_group)
                use_these += [PrepareReads.Reads(path_to_baga = use_path_reads)]
                print('Loaded: {}'.format(use_path_reads))
            if args.genome_name:
                caller = CallVariants.CallerDiscoSNP(reads = use_these, genome = genome)
            else:
                caller = CallVariants.CallerDiscoSNP(reads = use_these)
            if args.genome_name:
                add_prefix = use_name_genome + '_' + '+'.join(use_these_names)
            else:
                add_prefix = 'noref_' + '+'.join(use_these_names)
            caller.call(use_existing_graph = args.use_existing_graph, add_prefix = add_prefix, arguments = use_arguments)
        elif any([args.callsingles,
                args.calleach, 
                args.calljoint, 
                args.hardfilter, 
                args.recalibrate
                ]):
            # check if direct commands make sense
            direct_arguments = check_direct_arguments(args.arguments, 
                    wrapped_tools = ['HaplotypeCaller', 'GenotypeGVCFs'])
            # GATK only other option implemented.
            # check_baga_path() can handle full object names
            # e.g. baga.CollectData.Genome-mygenome.baga as well as actual names
            # i.e. mygenome . . . but not for this reads__genome compound name
            # so args.reads must be given correctly. Feedback given is a file can't
            # be found so should still be fine for ease-of-use.
            reads_genome_name = '__'.join([args.reads_name[0], use_name_genome])
            alns_name = reads_genome_name
            # see what is already available (existing analysis of data from
            # previous baga stage or data from this stage)
            use_path_alns,use_name_alns = check_baga_path(
                    'baga.AlignReads.SAMs', reads_genome_name)
            use_path_caller,use_name_caller = check_baga_path(
                    'baga.CallVariants.CallerGATK', reads_genome_name)
            
            if args.max_memory:
                max_memory = args.max_memory
            else:
                max_memory = 8
            if args.new:
                print('Starting new variant calling analysis because --new '\
                        'requested. Will overwrite any previous analyses.')
                from baga import AlignReads
                assert all([use_path_alns,use_name_alns]), 'Could not locate a saved '\
                        'baga.AlignReads.SAMs-<reads_name>__<genome_name>.baga for '\
                        'reads group and genome combination given: {}'.format(
                        use_name_alns)
                print('Loading alignments information for: {} from AlignReads output'.format(use_name_alns))
                alignments = AlignReads.SAMs(sample_name = alns_name, 
                        inherit_from = 'self')
                caller = CallVariants.CallerGATK(alignments = alignments)
            elif use_path_caller:
                # attempt to resume (if use_path_caller not False)
                print('Loading existing variants call analysis for: {}'.format(use_name_alns))
                print('(use --new to start variant calling again)')
                caller = CallVariants.CallerGATK(
                        baga = use_path_caller)
            else:
                # start new as previous not found
                assert all([use_path_alns,use_name_alns]), 'Could not locate a saved '\
                        'baga.AlignReads.SAMs-<reads_name>__<genome_name>.baga for '\
                        'reads group and genome combination given: {}'.format(
                        use_name_alns)
                print('Starting new variant calling analysis (could not find '\
                        'previous baga.CallVariants.CallerGATK-{}.baga)'.format(
                        reads_genome_name))
                from baga import AlignReads
                print('Loading alignments information for: {} from AlignReads '\
                        'output'.format(use_name_alns))
                alignments = AlignReads.SAMs(sample_name = alns_name, 
                        inherit_from = 'self')
                caller = CallVariants.CallerGATK(alignments = alignments)
            
            # because --new, setting --force to true to overwrite each output file
            if args.new:
                print('Because --new, also setting --force to overwrite each output file')
                args.force = True
            
            if args.callsingles and (args.calleach or args.calljoint):
                print('--callsingles for calling variants in individual samples '\
                        'cannot be used with the --calleach plus --calljoint '\
                        'combination that is used for joint variant calling for a cohort')
                sys.exit(1)
            
            if args.callsingles:
                try:
                    use_arguments = direct_arguments['HaplotypeCaller']
                except KeyError:
                    use_arguments = False
                caller.CallVCFsGATK(
                            mem_num_gigs = max_memory,
                            jar = args.GATK_jar_path.split(os.path.sep),
                            use_java = use_java,
                            force = args.force, 
                            max_cpus = args.max_cpus,
                            arguments = use_arguments)
                
                caller.saveLocal(use_name_alns)
                
            if args.calleach:
                try:
                    use_arguments = direct_arguments['HaplotypeCaller']
                except KeyError:
                    use_arguments = False
                caller.CallgVCFsGATK(
                            mem_num_gigs = max_memory,
                            jar = args.GATK_jar_path.split(os.path.sep),
                            use_java = use_java,
                            force = args.force, 
                            max_cpus = args.max_cpus,
                            arguments = use_arguments)
                
                caller.saveLocal(use_name_alns)
            
            if args.calljoint:
                try:
                    use_arguments = direct_arguments['GenotypeGVCFs']
                except KeyError:
                    use_arguments = False
                caller.GenotypeGVCFsGATK(
                            reads_genome_name,
                            jar = args.GATK_jar_path.split(os.path.sep), 
                            use_java = use_java,
                            # ultimately, scale this by the number of samples involved
                            # needed >8 for a set of 40
                            mem_num_gigs = max_memory, 
                            force = args.force,
                            arguments = use_arguments)
                
                caller.saveLocal(use_name_alns)
            
            if args.hardfilter:
                caller.hardfilterSNPsGATK(
                            jar = args.GATK_jar_path.split(os.path.sep),
                            use_java = use_java,
                            force = args.force)
                
                caller.hardfilterINDELsGATK(
                            jar = args.GATK_jar_path.split(os.path.sep),
                            use_java = use_java,
                            force = args.force)
                
                caller.saveLocal(use_name_alns)
            
            if args.recalibrate:
                # this is slow!
                caller.recalibBaseScoresGATK(
                            jar = args.GATK_jar_path.split(os.path.sep),
                            use_java = use_java,
                            force = args.force, 
                            mem_num_gigs = max_memory, 
                            max_cpus = args.max_cpus)
                
                caller.saveLocal(use_name_alns)
    

### Filter Variants ###

if args.subparser == 'FilterVariants':
    print('\n-- Filter Variants (part of the Variant Calling module) --\n')
    
    ## to apply variants, provide one reads group name
    if args.reads_name:
        if len(args.reads_name) > 1:
            sys.exit('Filters can only be applied to one group of reads at a time. Multiple sets can be handled with the --report option, though.')
    
    ## to report effects of filters, provide one or more read group names along with --report 
    
    from baga import CallVariants
    
    use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
    e = 'Could not locate a saved baga.CollectData.Genome-<genome_name>.baga for name given: {}'.format(args.genome_name)
    assert all([use_path_genome,use_name_genome]), e
    
    VCFs_for_report = {}
    
    # collect VCFs
    if args.reads_name:
        for these_reads in args.reads_name:
            # baga pipeline information provided
            # allow for multiple rounds of recalibration at end of CallVariants i.e., 1 or 2 and select 2 if available
            import baga
            # sometimes the baga from the previous step in the pipeline is not actually needed
            # so this name and file check could be relaxed
            use_path_reads,use_name_reads = check_baga_path('baga.PrepareReads.Reads', these_reads)
            e = 'Could not locate a saved baga.PrepareReads.Reads-<reads_name>.baga for reads group given: {}'.format(these_reads)
            assert all([use_path_reads,use_name_reads]), e
            alns_name = '__'.join([use_name_reads, use_name_genome])
            
            filein = 'baga.CallVariants.CallerGATK-{}.baga'.format(alns_name)
            caller = CallVariants.CallerGATK(baga = filein)
            
            if hasattr(caller, 'path_to_hardfiltered_SNPs') and hasattr(caller, 'path_to_hardfiltered_INDELs'):
                if isinstance(caller.path_to_hardfiltered_SNPs[-1], list):
                    # --callsingles
                    VCFs = caller.path_to_hardfiltered_SNPs[-1] + caller.path_to_hardfiltered_INDELs[-1]
                else:
                    # --calleach --calljoint
                    VCFs = [caller.path_to_hardfiltered_SNPs[-1], caller.path_to_hardfiltered_INDELs[-1]]
                # more than one can be handled with --report though
                # only implemented for separate VCFs currently
                VCFs_for_report[these_reads] = {
                        'SNPs':caller.path_to_hardfiltered_SNPs[-1], 
                        'InDels':caller.path_to_hardfiltered_INDELs[-1]}
            elif hasattr(caller, 'path_to_unfiltered_VCF'):
                print('WARNING: path to GATK hardfiltered variants not found in {}'.format(filein))
                print('It is recommended to complete the GATK variant calling with the CallVariants module')
                if isinstance(caller.path_to_unfiltered_VCF[-1],list):
                    VCFs = caller.path_to_unfiltered_VCF[-1]
                else:
                    VCFs = [caller.path_to_unfiltered_VCF[-1]]
            elif hasattr(caller, 'paths_to_raw_gVCFs'):
                print('WARNING: path to GATK joint called variants not found in {}'.format(filein))
                print('It is recommended to complete the GATK variant calling with the CallVariants module')
                if isinstance(caller.paths_to_raw_gVCFs[-1],list):
                    VCFs = caller.paths_to_raw_gVCFs[-1]
                else:
                    VCFs = [caller.paths_to_raw_gVCFs[-1]]
            else:
                print('WARNING: path to GATK called variants not found in {}'.format(filein))
                sys.exit('It seems the analysis described in {} is incomplete. Try completing or rerunning using the CallVariants module'.format(filein))
        
    else:
        # list of folders or files provided in args.vcfs_paths
        VCFs = []
        for path in args.vcfs_paths:
            if os.path.isdir(path):
                path_contents = os.listdir(path)
                theseVCFs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('VCF', 'vcf')]
                e = 'No VCF files (*.vcf or *.VCF) found in:\n{}'.format(args.alignments_paths)
                assert len(theseVCFs), e
                VCFs += theseVCFs
            else:
                # add file
                VCFs += [path]
    
    print('Loaded VCF locations:\n{}'.format('\n'.join(VCFs)))
    
    # check accessible and collect contained chromsome/sequence accessions
    # and sample name from columns
    # ... could at this point merge all samples called into one VCF?
    sample_names = {}
    for VCF in VCFs:
        try:
            with open(VCF, 'r') as filein:
                header, header_section_order, colnames, variants = CallVariants.parseVCF(VCF)
                contigs = {}
                for contiginfo in header['contig']:
                    bits = contiginfo.split('<')[-1].split('>')[0].split(',')
                    contiginfodict = dict([bit.split('=') for bit in bits])
                    #contig_id,contig_len = 
                    contigs[contiginfodict['ID']] = contiginfodict['length']
                for sample_name in colnames[9:]:
                    sample_names[sample_name] = contiginfo
            #sample_names[VCF] = {} ## delete?
        except IOError as e:
            print(e)
            sys.exit('Failed to open: {}'.format(VCF))
            
    # genome is only needed when generating csv summary with ORFs affected
    from baga import CollectData
    print('Loading genome %s' % use_name_genome)
    #genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
    genome = CollectData.Genome(use_name_genome, inherit_from = 'self')
    
    print('Loading filter information . . .')
    
    filters = {}
    if 'genome_repeats' in args.filters:
        from baga import Repeats
        finder = Repeats.Finder(use_name_genome, inherit_from = 'self')
        filters['genome_repeats'] = finder.ambiguous_ranges
        # print a summary
        for replicon_id,ambiguous_ranges in sorted(filters['genome_repeats'].items()):
            t = len(ambiguous_ranges)
            s = sum([(e - s) for s,e in ambiguous_ranges])
            print('Reference genome {}, chromosome sequence {}, contains {} '\
                    'repeated regions spanning {:,} basepairs'.format(use_name_genome, 
                    replicon_id, t, s))
    
    if 'rearrangements' in args.filters:
        
        from baga import Structure
        
        filters['rearrangements'] = {}
        #for sample,checker in checkers.items():
        for sample,genomeinfo in sorted(sample_names.items()):
            
            if args.include_samples:
                if sample not in args.include_samples:
                    continue
            
            filein = 'baga.Structure.CheckerInfo-{}__{}.baga'.format(sample, 
                    use_name_genome)
            if args.path_to_rearrangements_info:
                filein = os.path.sep.join([args.path_to_rearrangements_info,filein])
            
            checker_info = Structure.loadCheckerInfo(filein)
            
            e = 'Genome name for checker {} ({}) does not match supplied '\
                    'genome name ({})'.format(filein, 
                    checker_info['genome_name'], genome.sample_name)
            name_match = checker_info['genome_name'] == genome.sample_name
            genome_lengths = sorted([len(seq) for \
                    seq in genome.sequence.values()])
            checker_genome_lengths = sorted([length for \
                    length in checker_info['genome_lengths'].values()])
            # checker_info['genome_lengths'] = len(genome.sequence)
            length_match = genome_lengths == checker_genome_lengths
            
            if not name_match:
                print('WARNING: genome name used for rearrangements filter '\
                        'for {} ({}) does not match requested genome to use '\
                        '({}).'.format(sample, checker_info['genome_name'], 
                        genome.sample_name))
                if length_match:
                    print('Genome lengths match so proceeding and assuming '\
                            'alternative names or accession numbers for same '\
                            'genome ( bp).'.format('+'.join(str(l) for l in \
                            checker_genome_lengths)))
                else:
                    sys.exit('ERROR: genome length mismatch. Different genome '\
                            'used for rearrangments filter analysis to one '\
                            'requested to use to apply filter? ({} bp vs '\
                            '{:,} bp).'.format('+'.join(str(l) for l in genome_lengths),
                            '+'.join(str(l) for l in checker_genome_lengths)))
            
            # report brief summary for this sample
            print('For sample {}, relative to {}:'.format(sample, genome.sample_name))
            for replicon_id, ranges in checker_info['suspect_regions']['rearrangements'].items():
                t = len(ranges)
                s = sum([(e - s) for s,e in ranges])
                print('  in replicon (or contig) {}, {} regions spanning {:,} basepairs '\
                        'are affected by rearrangements thus having ambiguous 1:1 '\
                        'orthology.'.format(replicon_id, t, s))
            
            for replicon_id, ranges in checker_info['suspect_regions']['rearrangements_extended'].items():
                t = len(ranges)
                s = sum([(e - s) for s,e in ranges])
                print('  in replicon (or contig) {}, {} additional regions spanning {:,} '\
                        'basepairs are adjacent to the above regions but have a >50% zero '\
                        'read depth over a moving window (i.e., any aligned reads have a '\
                        'patchy distribution, are usually rare). These are typically large '\
                        'deletions including missing prophage and genomic islands\n'\
                        ''.format(replicon_id, t, s))
    
    filter_applier = CallVariants.Filter(VCFs, genome)  #, use_name_reads)
    filter_applier.doFiltering(filters)

    

### Summarise Variants ###

if args.subparser == 'SummariseVariants':
    print('\n-- Summarise Variants (part of the Variant Calling module) --\n')
    
    from baga import CallVariants
    ## to apply variants, provide one reads group name
    if args.reads_name:
        assert args.genome_names, "--genome_names is required with --reads_name"
    
    VCFs_for_report = {}
    
    # collect genomes
    # plural genomes is a bit of a hack at the moment and only works here with named VCFs . . .
    if args.genome_names:
        from baga import CollectData
        use_genomes = []
        for genome_name in args.genome_names:
            use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', genome_name)
            assert all([use_path_genome,use_name_genome]), 'Could not locate genome: {}'.format(genome_name)
            use_genomes += [CollectData.Genome(use_name_genome, inherit_from = 'self')]
    else:
        use_genomes = False
    
    # collect VCFs
    if args.reads_name:
        for these_reads in args.reads_name:
            # baga pipeline information provided
            # allow for multiple rounds of recalibration at end of CallVariants i.e., 1 or 2 and select 2 if available
            import baga
            # sometimes the baga from the previous step in the pipeline is not actually needed
            # so this name and file check could be relaxed
            use_path_reads,use_name_reads = check_baga_path('baga.PrepareReads.Reads', these_reads)
            e = 'Could not locate a saved baga.PrepareReads.Reads-<reads_name>.baga for reads group given: {}'.format(these_reads)
            assert all([use_path_reads,use_name_reads]), e
            alns_name = '__'.join([use_name_reads, use_name_genome])
            
            filein = 'baga.CallVariants.CallerGATK-{}.baga'.format(alns_name)
            caller = CallVariants.CallerGATK(baga = filein)
            
            if hasattr(caller, 'path_to_hardfiltered_SNPs') and hasattr(caller, 'path_to_hardfiltered_INDELs'):
                if isinstance(caller.path_to_hardfiltered_SNPs[-1], list):
                    # --callsingles
                    VCFs = caller.path_to_hardfiltered_SNPs[-1] + caller.path_to_hardfiltered_INDELs[-1]
                else:
                    # --calleach --calljoint
                    VCFs = [caller.path_to_hardfiltered_SNPs[-1], caller.path_to_hardfiltered_INDELs[-1]]
                # more than one can be handled with --report though
                # only implemented for separate VCFs currently
                VCFs_for_report[these_reads] = {'SNPs':caller.path_to_hardfiltered_SNPs[-1], 'InDels':caller.path_to_hardfiltered_INDELs[-1]}
            elif hasattr(caller, 'path_to_unfiltered_VCF'):
                print('WARNING: path to GATK hardfiltered variants not found in {}'.format(filein))
                print('It is recommended to complete the GATK variant calling with the CallVariants module')
                if isinstance(caller.path_to_unfiltered_VCF[-1],list):
                    # called with --callsingle: one sample per VCF
                    VCFs = caller.path_to_unfiltered_VCF[-1]
                else:
                    # called with --calleach --calljoint: multi-sample VCFs
                    VCFs = [caller.path_to_unfiltered_VCF[-1]]
            elif hasattr(caller, 'paths_to_raw_gVCFs'):
                print('WARNING: path to GATK joint called variants not found in {}'.format(filein))
                print('It is recommended to complete the GATK variant calling with the CallVariants module')
                if isinstance(caller.paths_to_raw_gVCFs[-1],list):
                    # called with --callsingle: one sample per VCF
                    VCFs = caller.paths_to_raw_gVCFs[-1]
                else:
                    # called with --calleach --calljoint: multi-sample VCFs
                    VCFs = [caller.paths_to_raw_gVCFs[-1]]
            else:
                print('WARNING: path to GATK called variants not found in {}'.format(filein))
                sys.exit('It seems the analysis described in {} is incomplete. Try completing or rerunning using the CallVariants module'.format(filein))
        
    else:
        # list of folders or files provided in args.vcfs_paths
        VCFs = []
        for path in args.vcfs_paths:
            if os.path.isdir(path):
                path_contents = os.listdir(path)
                theseVCFs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('VCF', 'vcf')]
                e = 'No VCF files (*.vcf or *.VCF) found in:\n{}'.format(args.alignments_paths)
                assert len(theseVCFs), e
                VCFs += theseVCFs
            else:
                # add file
                VCFs += [path]

    print('Loaded VCF locations:\n{}'.format('\n'.join(VCFs)))
    
    message = 'Need --filters FILTER_NAME [FILTER_NAME] for report type "{}"'
    if args.cumulative:
        assert args.filters, message.format('cumulative')
        print('Reporting cumulative variant totals by class and group as '\
                'filters applied : {}'.format('+'.join(args.filters)))
        for genome in use_genomes:
            CallVariants.reportCumulative(args.filters, genome.id, VCFs_for_report)
    
    if args.lists:
        assert args.filters, message.format('lists')
        print('Reporting lists of variants by class and group with filters'\
                ': {}'.format('+'.join(args.filters)))
        for genome in use_genomes:
            CallVariants.reportLists(args.filters, genome.id, VCFs_for_report)
    
    if args.simple:
        summariser = CallVariants.Summariser(VCFs, genomes = use_genomes)
        summariser.simple()

### Check Linkage ###

if args.subparser == 'CheckLinkage':
    print('\n-- Check Linkage (part of the Variant Calling module) --\n')
    # required input: paths to corresponding VCFs and BAMs
    # which baga objects contain that information?
    from baga import CollectData
    use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
    from baga import CallVariants
    if args.reads_name:
        # allegedly works in 2 and 3
        raise NotImplementedError('Waiting for pooled samples to be implemented in baga.CallVariants. Use --vcfs_paths --alignments_paths instead')
        # part of baga pipeline
        if not args.genome_name:
            parser.error('--genome_name/-g is required with --reads_name/-n. (The baga CollectData-processed genome used with the AlignReads option)')
        
        # in this case, don't need to load the genome, just have its sanitised name
        # e = 'Could not locate a saved baga.CollectData.Genome-<genome_name>.baga for name given: {}'.format(args.genome_name)
        # assert all([use_path_genome,use_name_genome]), e
        # allow for multiple rounds of recalibration at end of CallVariants i.e., 1 or 2 and select 2 if available
        # sometimes the baga from the previous step in the pipeline is not actually needed
        # so this name and file check could be relaxed
        use_path_reads,use_name_reads = check_baga_path('baga.PrepareReads.Reads', args.reads_name)
        e = 'Could not locate a saved baga.PrepareReads.Reads-<reads_name>.baga for reads group given: {}'.format(these_reads)
        assert all([use_path_reads,use_name_reads]), e
        alns_name = '__'.join([use_name_reads, use_name_genome])
        
        filein = 'baga.CallVariants.CallerGATK-{}.baga'.format(alns_name)
        caller = CallVariants.CallerGATK(baga = filein)
        
        if hasattr(caller, 'path_to_hardfiltered_SNPs') and hasattr(caller, 'path_to_hardfiltered_INDELs'):
            # only one for VCFs so no overwriting
            VCFs = [caller.path_to_hardfiltered_SNPs[-1], caller.path_to_hardfiltered_INDELs[-1]]
        elif hasattr(caller, 'path_to_unfiltered_VCF'):
            print('WARNING: path to GATK hardfiltered variants not found in {}'.format(filein))
            print('It is recommended to complete the GATK variant calling with the CallVariants module')
            VCFs = [caller.path_to_unfiltered_VCF[-1]]
        elif hasattr(caller, 'paths_to_raw_gVCFs'):
            print('WARNING: path to GATK joint called variants not found in {}'.format(filein))
            print('It is recommended to complete the GATK variant calling with the CallVariants module')
            VCFs = caller.paths_to_raw_gVCFs[-1]
        else:
            print('WARNING: path to GATK called variants not found in {}'.format(filein))
            sys.exit('It seems the analysis described in {} is incomplete. Try completing or rerunning using the CallVariants module'.format(filein))
        
        print('Loading alignments information for: {}__{} from AlignReads output'.format(use_name_group, use_name_genome))
        from baga import AlignReads
        alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}__{}.baga'.format(use_name_group, use_name_genome))
        
        e = 'the reads for "--reads_name/-n {}" seem to not have been fully processed by the AlignReads module: they are missing the "ready_BAMs" attribute. Please ensure the AlignReads commands "--align --deduplicate --indelrealign" have been performed.'.format(args.reads_name)
        assert hasattr(alignments, 'ready_BAMs'), e
        ### shouldn't all BAMs have headers parsed and stored in dict with (sample,genome) from the start?
        BAMs = alignments.ready_BAMs[-1]
        sample_names = sorted(alignments.read_files)
        # check on requested samples: crop BAMs accordingly <== and VCFs, and what about sample matching . . .
        if args.include_samples:
            if not args.reads_name:
                # need sample names so will parse BAMs
                import pysam
                sample_names = {}
                for BAM in BAMs:
                    sample_names[pysam.Samfile(BAM, 'rb').header['RG'][0]['ID']] = BAM
            
            missing_labels = sorted(set(args.include_samples) - set(sample_names))
            found_labels = sorted(set(args.include_samples) & set(sample_names))
            e = ['None of the requested sample labels were found among the previously checked reads.']
            e += ['Requested: {}'.format(', '.join(args.include_samples))]
            e += ['Available: {}'.format(', '.join(sorted(sample_names)))]
            assert len(found_labels), '\n'.join(e)
            
            if len(missing_labels):
                print('WARNING: could not find the following requested samples among previously checked reads.')
                print(', '.join(missing_labels))
            
            # update BAMs for args.check
            print(BAMs)
            if not args.reads_name:
                BAMs = [sample_names[sample] for sample in found_labels]
            else:
                BAMs = [BAM for BAM in BAMs if BAM.split(os.path.sep)[-1].split('__')[0] in found_labels]
            
            print(BAMs)
            # update sample_names for arg.plot
            sample_names = sorted(found_labels)
            
    else:
        # list of folders or files provided in args.vcfs_paths
        VCFs = []
        for path in args.vcfs_paths:
            if os.path.isdir(path):
                path_contents = os.listdir(path)
                theseVCFs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('VCF', 'vcf')]
                e = 'No VCF files (*.vcf or *.VCF) found in:\n{}'.format(args.alignments_paths)
                assert len(theseVCFs), e
                VCFs += theseVCFs
            else:
                # add file
                VCFs += [path]
        print('Loaded VCF locations:\n{}'.format('\n'.join(VCFs)))
        BAMs = []
        for path in args.alignments_paths:
            if os.path.isdir(path):
                path_contents = os.listdir(path)
                theseBAMs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('BAM', 'bam')]
                e = 'No BAM files (*.bam or *.BAM) found in:\n{}'.format(args.alignments_paths)
                assert len(theseBAMs), e
                BAMs += theseBAMs
            else:
                # add file
                BAMs += [path]
        print('Loaded BAM locations:\n{}'.format('\n'.join(BAMs)))
        print('Loading genome %s' % use_name_genome)
        genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
    
    if args.check:
        linkage_checker = CallVariants.Linkage(genome = genome, vcf_paths = VCFs, alignment_paths = BAMs)
        linkage_checker.doLinkageCheck()
    else:
        print('use --check to actually fo the checking . . .')

if args.subparser == 'ComparativeAnalysis':
    print('\n-- Comparative Analyses --\n')
    
    from baga import AlignReads
    from baga import CallVariants
    from baga import ComparativeAnalysis
    from baga import CollectData
    
    if args.build_MSA:
        # required: reference genome
        # required: reads name OR path to VCFs (. . optional path to BAMs to properly include gaps)
        # ... need to link VCFs to BAMs . . . . 
        # check appropriate combination of options provided
        if args.reads_name:
            if args.sample_bams:
                print('If --reads_name/-n provided, --sample_bams/-B is not necessary: ignoring latter')
        elif args.vcfs_paths:
            # these two added as "add_mutually_exclusive_group"
            if args.include_invariants:
                if not args.sample_bams:
                    print('WARNING: making a full-length multiple-sequence alignment without checking read alignments for missing pieces of chromosome is a risky assumption! Only proceed if you know there are no missing pieces of chromosome among your samples relative to the reference chromosome and/or BAMs are unavailable.')
        
        # VCFs contain sample names, require file sample_name\tpath_to_bam\n
        use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
        e = 'Could not locate a saved baga.CollectData.Genome-<genome_name>.baga for name given: {}'.format(args.genome_name)
        assert all([use_path_genome,use_name_genome]), e
        if args.reads_name:
            # baga pipeline information provided <== only way currently implemented
            # allow for multiple rounds of recalibration at end of CallVariants i.e., 1 or 2 and select 2 if available
            
            import baga
            path_to_SNPs_VCFs = []
            path_to_InDels_VCFs = []
            paths_to_BAMs = []
            for these_reads in args.reads_name:
                print('Collecting VCFs for {}'.format(these_reads))
                # sometimes the baga from the previous step in the pipeline is not actually needed
                # so this name and file check could be relaxed
                use_path_reads,use_name_reads = check_baga_path('baga.PrepareReads.Reads', these_reads)
                e = 'Could not locate a saved baga.PrepareReads.Reads-<reads_name>.baga for reads group given: {}'.format(these_reads)
                assert all([use_path_reads,use_name_reads]), e
                alns_name = '__'.join([use_name_reads, use_name_genome])
                
                filein = 'baga.CallVariants.CallerGATK-{}.baga'.format(alns_name)
                caller = CallVariants.CallerGATK(baga = filein)
                
                if hasattr(caller, 'path_to_hardfiltered_SNPs') and \
                        hasattr(caller, 'path_to_hardfiltered_INDELs'):
                    if isinstance(caller.path_to_hardfiltered_SNPs[-1],list):
                        # called with --callsingle: one sample per VCF
                        VCFs = caller.path_to_hardfiltered_SNPs[-1]
                    else:
                        # called with --calleach --calljoint: multi-sample VCFs
                        VCFs = [caller.path_to_hardfiltered_SNPs[-1]]
                    for checkthis in VCFs:
                        try:
                            with open(checkthis) as fin:
                                #### list
                                path_to_SNPs_VCFs += [checkthis]
                                print('Found: {}'.format(checkthis))
                        except IOError:
                            print('Could not find: {}'.format(checkthis))
                            sys.exit('You may need to rerun the analysis that '\
                                    'should have generated that file')
                    
                    if isinstance(caller.path_to_hardfiltered_SNPs[-1],list):
                        # called with --callsingle: one sample per VCF
                        VCFs = caller.path_to_hardfiltered_INDELs[-1]
                    else:
                        # called with --calleach --calljoint: multi-sample VCFs
                        VCFs = [caller.path_to_hardfiltered_INDELs[-1]]
                    for checkthis in VCFs:
                        try:
                            with open(checkthis) as fin:
                                path_to_InDels_VCFs += [checkthis]
                                print('Found: {}'.format(checkthis))
                        except IOError:
                            print('Could not find: {}'.format(checkthis))
                            sys.exit('You may need to rerun the analysis that should have generated that file')
                    
                else:
                    print('ERROR: path to GATK hardfiltered variants not found in {}'.format(filein))
                    print('It seems the analysis described in {} is incomplete. Try completing or rerunning using the CallVariants module'.format(filein))
                    NotImplementedError('Building multiple alignments from one VCF per sample is not yet implemented: Coming soon!')
                    # see code in ApplyFilters section to collect other per sample VCFs
                
                # greedily collect version of each vcf with the most filters applied
                import re
                from glob import glob
                patt = re.compile('(__F_)')
                path_to_SNPs_VCFs_use = []
                for path in path_to_SNPs_VCFs:
                    numfilters = {}
                    for path2 in glob(path.replace('.vcf','*.vcf').replace('.VCF','*.VCF')):
                        numfilters[path2] = len(re.findall(patt,path2))
                    
                    path_to_SNPs_VCFs_use += [sorted(numfilters, key = numfilters.get)[-1]]
                
                print('Using:\n{}'.format('\n'.join(path_to_SNPs_VCFs_use)))
                path_to_SNPs_VCFs = path_to_SNPs_VCFs_use
                
                print('Collecting BAMs for {}'.format(these_reads))
                alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}.baga'.format(alns_name))
                for BAM in alignments.ready_BAMs[-1]:
                    checkthis = BAM
                    try:
                        with open(checkthis) as fin:
                            paths_to_BAMs += [checkthis]
                            print('Found: {}'.format(checkthis))
                    except IOError:
                        print('Could not find: {}'.format(checkthis))
                        sys.exit('You may need to rerun the analysis that should have generated that file')
            
            MSA_filename = '{}__{}_SNPs'.format(use_name_genome,'_'.join(args.reads_name))
            paths_to_VCFs = path_to_SNPs_VCFs + path_to_InDels_VCFs
            
        else:
            # list of folders or files provided in args.vcfs_paths
            # not part of a baga pipeline, so need BAMs linked to samples separately in --sample_bams
            paths_to_VCFs = []
            paths_to_BAMs = []
            for path in args.vcfs_paths:
                if os.path.isdir(path):
                    path_contents = os.listdir(path)
                    theseVCFs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('VCF', 'vcf')]
                    e = 'No VCF files (*.vcf or *.VCF) found in:\n{}'.format(args.alignments_paths)
                    assert len(theseVCFs), e
                    paths_to_VCFs += theseVCFs
                else:
                    # add file
                    paths_to_VCFs += [path]
            
            # could forbid full size alignments without BAMs to inform missing regions
            # if any([args.include_invariants, args.sample_bams]):
                # e = '--include_invariants and --sample_bams must be used together.'
                # assert all([args.include_invariants, args.sample_bams]), e
            
            # but will allow it with a warning for flexibility in case BAMs unavailable etc
            if args.sample_bams:
                # path_to_VCFs and file linking samples to supplied instead of full baga pipeline info
                try:
                    BAMs = dict([line.rstrip().split('\t') for line in open(args.sample_bams).readlines()])
                except IOError:
                    print('there was a problem reading file: {}'.format(args.sample_bams))
                except ValueError:
                    print('there was a problem parsing file: {}'.format(args.sample_bams))
                
                # make a list for .getCoverageRanges()
                paths_to_BAMs = sorted(BAMs.values())
            
            path_to_InDels_VCFs = False
            # could generate better name here? 
            if len(paths_to_VCFs) == 1:
                vcf_name = '1_VCF'
            else:
                vcf_name = '{}_VCFs'.format(len(paths_to_VCFs))
            
            MSA_filename = '{}__{}_SNPs'.format(use_name_genome,vcf_name)
        
        print('Loaded VCF locations:\n{}'.format('\n'.join(paths_to_VCFs)))
        
        ### now collected required info: build MSA
        print('Loading genome %s' % use_name_genome)
        genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
        
        MSA_builder = ComparativeAnalysis.MultipleSequenceAlignment(paths_to_VCFs)
        MSA_builder.collectVariants(samples_to_include = args.include_samples,
                                    samples_to_exclude = args.exclude_samples)
        
        
        if len(paths_to_BAMs):
            print('Loaded BAM locations:\n{}'.format('\n'.join(paths_to_BAMs)))
            MSA_builder.getCoverageRanges(paths_to_BAMs)
        
        MSA_builder.writeMSA(   MSA_filename, 
                                strict_core = args.core_only, 
                                include_invariants = args.include_invariants, 
                                genome = genome)
    
    if args.infer_phylogeny:
        # args
        assert args.path_to_MSA is not None, '--path_to_MSA is required for --infer_phylogeny'
        phylo_analyser = ComparativeAnalysis.Phylogenetics(args.path_to_MSA)
        if args.program == 'phyML':
            #currently only and default choice!
            phylo_analyser.estimate_phylogeny_PhyML(num_bootstraps = args.num_bootstraps, collect_previous = False)
            phylo_analyser.load_tree()
            phylo_analyser.reroot_to_outgroup(args.out_group)
            phylo_analyser.write_tree()
    
    if args.infer_recombination:
        assert args.path_to_MSA is not None, '--path_to_MSA is required for --infer_recombination (should be that used to estimate tree at --path_to_tree)'
        assert args.path_to_tree is not None, '--path_to_tree is required for --infer_recombination (should be that estimated from alignment at --path_to_MSA)'
        phylo_analyser = ComparativeAnalysis.Phylogenetics(args.path_to_MSA, path_to_tree = args.path_to_tree)
        # bit of a fudge dealing with rooted tree . . . to be improved
        phylo_analyser.collectPhyMLstats(path_to_stats_file = args.path_to_tree.replace('_rooted','').replace('_phyml_tree','') + '_phyml_stats')
        phylo_analyser.infer_recombination(bootstraps = args.num_bootstraps) #, output_suffix = '_rooted')
    
    if args.plot_phylogeny:
        # should check either or etc here
        if args.genome_name:
            use_path_genome,use_name_genome = check_baga_path('baga.CollectData.Genome', args.genome_name)
            print('Loading genome %s' % use_name_genome)
            genome = CollectData.Genome(local_path = use_path_genome, format = 'baga')
            genome_length = len(genome.sequence)
        elif args.genome_length:
            genome_length = args.genome_length
        else:
            print('Provide --genome_name or --genome_length for a scale bar unit of actual substitutions')
            genome_length = False
        
        if args.plot_transfers:
            bits = args.path_to_tree.split(os.path.extsep)
            bits[-2] = bits[-2]+'_transfers'
            plot_output_path = os.path.extsep.join(bits[:-1] + ['svg'])
        else:
            plot_output_path = os.path.extsep.join(args.path_to_tree.split(os.path.extsep)[:-1] + ['svg'])
        
        print('Plotting to {}'.format(plot_output_path))
        
        phylo_plotter = ComparativeAnalysis.Plotter(
                            plot_output_path,
                            args.path_to_tree.replace('_rooted',''),   # deal with rootedness at some point
                            # smaller values keeps edges in smaller central
                            # zone allowing for longer tip labels
                            plot_width_prop = 0.75, 
                            plot_height_prop = 0.75)
        
        # for translating labels: supply tab delimited list of tip label and desired label
        # args.use_names
        ## either provide string that is path to raw ClonalFrameML output
        #plot_transfers = 'Both_all_core_positions_rooted.importation_status.txt'
        ## or provide a SNPs_by_homoplasies processed dictionary from a .summarise_recombined_variants() analysis
        #plot_transfers = phylo_analyser.SNPs_by_homoplasies
        if args.plot_transfers:
            plot_transfers = args.path_to_tree.replace('labelled_tree.newick','importation_status.txt').replace('_rooted','')   # deal with rootedness at some point
        else:
            plot_transfers = False
        
        if args.out_group:
            outgroup_label_list = args.out_group
        else:
            outgroup_label_list = []
        
        phylo_plotter.doPlot(outgroup_label_list = outgroup_label_list, 
                             stroke_width = 3, 
                             label_size = 15, 
                             plotinnerlabels = False,
                             plottiplabels = True,
                             plot_extra_lines = False,
                             direction = 'right',
                             plot_transfers = plot_transfers,
                             use_names = args.use_names,
                             scale_bar = True,
                             genome_length = genome_length)


### Assemble Reads ###

if args.subparser == 'AssembleReads':
    print('\n-- Reads Assembly module --')
    import baga
    from baga import AssembleReads
    for this_reads_name in args.reads_name:
        use_path_reads,use_name_reads = check_baga_path('baga.PrepareReads.Reads', this_reads_name)
        e = 'Could not locate a saved baga.PrepareReads.Reads-<reads_name>.baga '\
                'for reads group given: {}'.format(this_reads_name)
        assert all([use_path_reads,use_name_reads]), e
        print('Loading processed reads group %s' % use_name_reads)
        prepared_reads = baga.bagaload(use_path_reads)
        
        if args.program == 'spades':
            reads = AssembleReads.DeNovo(baga = prepared_reads)
            if args.max_memory:
                use_mem_gigs = args.max_memory
            else:
                # round down available GBs
                use_mem_gigs = int(baga.get_available_memory())
            
            # for more reliable: only_assembler = True, careful = False
            reads.SPAdes(mem_num_gigs = use_mem_gigs, only_assembler = False, 
                    careful = True)
    
    # if args.delete_intermediates:
        # print('Checking on intermediate fastq files to delete . . .')
        # if not args.adaptors and not args.trim:
            # # load a previously adaptor cut reads set
            # print('Loading processed reads group %s' % use_name_reads)
            # reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name_reads)
        
        # total_size = 0
        # # check stage 1 files and stage 2 files
        # if hasattr(reads,'read_files') and hasattr(reads,'adaptorcut_read_files'):
            # stage1s = check_files(reads.read_files)
            # stage2s = check_files(reads.adaptorcut_read_files)
            # if stage2s:
                # if stage1s:
                    # # delete stage 1 files if have all <== not for now . . .
                    # # print('Deleting original or subsampled fastq files . . .')
                    # # total_size += delete_files(reads.read_files, extra = '_subsmp')
                    # print('Retaining original fastq files even though processed versions exist because re-downloading is time consuming!')
                # else:
                    # print('Some or all of original or subsampled fastq files seem to have been deleted')
            # else:
                # print('Missing some cutadapt-processed files: not deleting originals or subsampled')
        
        # if hasattr(reads,'adaptorcut_read_files') and hasattr(reads,'trimmed_read_files'):
            # stage2s = check_files(reads.adaptorcut_read_files)
            # stage3s = check_files(reads.trimmed_read_files)
            
            # if stage3s:
                # if stage2s:
                    # # delete stage 2 files if have all
                    # print('Deleting cutadapt-processed fastq files . . .')
                    # total_size += delete_files(reads.adaptorcut_read_files)
                # else:
                    # print('Some or all of cutadapt-processed fastq files seem to have been deleted')
            # else:
                # print('Missing some sickle-processed files: not deleting cutadapt-processed')
        
        # if total_size:
            # print('Saved {:.2f} Gb by deleting intermediate files'.format(total_size/1000000000.0))
        # else:
            # print('Nothing deleted.')


# work out how long it took
time_end = datetime.now()
time_stamp = time_end.strftime("%H:%M:%S on %a %d %b, %Y")

baga_cli.info('')
baga_cli.info('|====== Finished tasks at {} =======|'.format(time_stamp))
baga_cli.info('')

run_time = time_end - time_start
seconds = int(run_time.total_seconds())
hours, remainder = divmod(seconds,60*60)
minutes, seconds = divmod(remainder,60)
baga_cli.info('Total run time: {} hours, {} minutes, {} seconds  (PID:{})'.format(
        hours, minutes, seconds, os.getpid()))

logging.shutdown()
