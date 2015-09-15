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
from collections import defaultdict

# a basic installation keeps the commandline interface in the package folder
path_to_baga_package = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-2])
# so importing the module needs to be done from there
sys.path.append(path_to_baga_package)

from baga import Dependencies

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


    
text_width = 70

title = 'Bacterial and Archaeal Genome Analyser'
subtitle = textwrap.fill('Novel analyses and wrapped tools pipelined for convenient processing of genome sequences', text_width)
version_date = 'July 10 2015'
version_num = 0.1
authors = 'David Williams'
email = 'david.williams.at.liv.d-dub.org.uk'
blurb = '''Work on this software was started at The University of Liverpool, UK 
with funding from The Wellcome Trust (093306/Z/10) awarded to:

Dr Steve Paterson (The University of Liverpool, UK)
Dr Craig Winstanley (The University of Liverpool, UK)
Dr Michael A Brockhurst (University of York, UK)


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

print(splash)

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

# group = parser.add_mutually_exclusive_group()
# group.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
# group.add_argument("-q", "--quiet", help="suppress output", action="store_true")

subparser_adder = parser.add_subparsers(title = 'Analyses', dest="subparser")


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

parser_CollectData = subparser_adder.add_parser('CollectData',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Download and parse genomes or \
download short reads for analysis by {}'.format(title),
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('Genomes can be loaded from genbank \
files with user provided paths or downloaded from the National Center for \
Biotechnology Information (NCBI) with user provided accession numbers. Short \
reads can be downloaded from the European Nucleotide Archive with user provided \
Run Accession numbers\n\nExample usage: "%(prog)s -r \
ERR953490,ERR953494,ERR953501,ERR953509,ERR953491,ERR953513"\n'\
.format(title), text_width, replace_whitespace = False))


mutually_exclusive_group = parser_CollectData.add_mutually_exclusive_group()

# could use nargs='+' or similar here to get a list directly,
# but not sure about spaces in paths
mutually_exclusive_group.add_argument("-g", "--genomes", 
    help="(down)load and parse genomes for analysis", 
    type = str,
    nargs = '+')

mutually_exclusive_group.add_argument("-r", "--reads_download", 
    help="download short reads for analysis", 
    type = str,
    nargs = '+')

mutually_exclusive_group.add_argument("-R", "--reads_path", 
    help="path to local short reads in fastq files for analysis. All read files will be collcted. Alternatively, a pattern containing '*' '?' and other shell expansion characters or an explicit list of read files can also be supplied.", 
    type = str,
    nargs = '+')


parser_CollectData.add_argument('-e', "--email_address", 
    help = "required for downloading from NCBI")

parser_CollectData.add_argument('-n', "--reads_group_name", 
    help = "optional for downloading from NCBI, required for loading from local path")

parser_PrepareReads = subparser_adder.add_parser('PrepareReads',
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

parser_PrepareReads.add_argument('-N', "--max_cpus", 
    help = "maximum number of cpu cores used when parallel processing",
    type = int,
    default = -1)

parser_PrepareReads.add_argument('-F', "--force", 
    help = "overwrite existing data: required when repeating an analysis else \
previous versions retained. Retension of previous versions is convenient for \
resuming an interrupted analysis in which only some read sets were processed.",
    action = 'store_true',
    default = False)

parser_PrepareReads.add_argument('-s', "--subsample_to_cov", 
    help = "subsample reads to a requested coverage of a given genome size. This \
provides smaller files of a consistent size saving storage space and processing \
time later and benefitting some analyses like de novo assembly. E.g.\n\
'--subsample_to_cov 80 5000000' for 80x coverage of a 5 Mb genome.",
    type = int,
    nargs = 2)

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

parser_AlignReads.add_argument('-N', "--max_cpus", 
    help = "maximum number of cpu cores used when parallel processing",
    type = int,
    default = -1)

parser_AlignReads.add_argument('-F', "--force", 
    help = "overwrite existing data: required when repeating an analysis else \
previous versions retained. Retension of previous versions is convenient for \
resuming an interrupted analysis in which only some read sets were processed.",
    action = 'store_true',
    default = False)

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

parser_AlignReads.add_argument('-P', "--GATK_jar_path", 
    help = "path to Genome Analysis Toolkit (GATK) jar file",
    type = str)

parser_AlignReads.add_argument('-D', "--delete_intermediates", 
    help = "delete intermediate SAM and BAM files to save space. Files are only deleted if those for next stage are found",
    action = 'store_true')


parser_Repeats = subparser_adder.add_parser('Repeats',
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

parser_Repeats.add_argument('-i', "--minimum_percent_identity", 
    help = "the lowest nucleotide percent identity permitted repeat over regions",
    type = int,
    default = 98)

parser_Repeats.add_argument('-l', "--minimum_repeat_length", 
    help = "shortest repeat length: should be similar to insert size of paired end reads as smaller (non-tandem) repeats should be resolvable and have unambiguous mappings",
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
%(prog)s --genome_name FM209186.1 --alignments_path path/to/my/bams --check\n', 
text_width, replace_whitespace = False))

mutually_exclusive_group = parser_Structure.add_mutually_exclusive_group(required=True)

mutually_exclusive_group.add_argument('-n', "--reads_name", 
    help = "name of read datasets group if processed by PrepareReads and AlignReads options. Should match --reads_name option used previously",
    type = str)

mutually_exclusive_group.add_argument('-a', "--alignments_path", 
    help = "path to paired-end short read alignments to a reference genome. If a directory path(s) is provided, all *.BAM and *.bam files will be included. A list of filenames with full path or with *, ? or [] characters can also be provided (with unix-style pathname pattern expansion for the latter)",
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

parser_Structure.add_argument('-s', "--include_samples", 
    help = "restrict plotting to these samples. Applies to --plot and --plot_range. If omitted, plots for all samples are produced.",
    type = str,
    nargs = '+')


parser_CallVariants = subparser_adder.add_parser('CallVariants',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Call variants with Genome Analysis \
Tool Kit (GATK) from each of a group of read sets previously aligned to a \
genome via the PrepareReads option.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('Groups of read sets from the \
AlignReads option are loaded by providing the name supplied previously.\n\n\
Example usage: %(prog)s --reads_name ERR953490plus5others \
--calleach --calljoint --hardfilter\n', 
text_width, replace_whitespace = False))

parser_CallVariants.add_argument('-r', "--reads_name", 
    help = "name of read datasets group processed by PrepareReads and AlignReads options",
    type = str,
    required = True)

parser_CallVariants.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option",
    type = str,
    required = True)

parser_CallVariants.add_argument('-F', "--force", 
    help = "overwrite existing per isolate data: required when repeating an \
analysis else previous versions retained. Retension of previous versions is \
convenient for resuming an interrupted analysis in which only some read sets \
were processed.",
    action = 'store_true',
    default = False)

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
    help = "maximum number of memory to use in gigabytes",
    type = int,
    default = 8)

parser_CallVariants.add_argument('-c', "--calleach", 
    help = "call variants in each alignment in preparation for a joint analysis.\
Called 1st time on uncalibrated alignments, called 2nd after base quality score \
recalibration.",
    action = 'store_true')

parser_CallVariants.add_argument('-j', "--calljoint", 
    help = "call variants in all alignments in a joint analysis. Called 1st time \
on uncalibrated alignments, called 2nd after base quality score recalibration.",
    action = 'store_true')

parser_CallVariants.add_argument('-f', "--hardfilter", 
    help = "apply 'hard filtering' thresholds on called variants to decrease \
false positives using GATK",
    action = 'store_true')

parser_CallVariants.add_argument('-R', "--recalibrate", 
    help = "apply read base quality score recalibration using GATK",
    action = 'store_true')

parser_CallVariants.add_argument('-P', "--GATK_jar_path", 
    help = "path to Genome Analysis Toolkit (GATK) jar file",
    type = str)



parser_ApplyFilters = subparser_adder.add_parser('ApplyFilters',
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = textwrap.fill('Apply filters determined by \
the Repeats and Structure options on variant calls. VCF files will be \
copied with updated information marking certain variants to be excluded. \
Vaiants can be inferred using the CallVariants option.',
                                            text_width,
                                            replace_whitespace = False),
                epilog = textwrap.fill('Filter regions and VCF files are \
loaded by providing the read set names (or VCF filenames, to be implemented).\n\n\
Example usage: %(prog)s --reads_name ERR953490plus5others \
--genome FM209186 --filters genome_repeats rearrangements\n', 
text_width, replace_whitespace = False))

mutually_exclusive_group = parser_ApplyFilters.add_mutually_exclusive_group(required=True)

mutually_exclusive_group.add_argument('-n', "--reads_name", 
    help = "name of read datasets group if processed by PrepareReads and AlignReads options. Should match --reads_name option used previously",
    type = str,
    nargs = '+')

mutually_exclusive_group.add_argument('-p', "--vcfs_path", 
    help = "path to vcf files. If a directory path(s) is provided, all *.VCF and *.vcf files will be included. A list of filenames with full path or with *, ? or [] characters can also be provided (with unix-style pathname pattern expansion for the latter)",
    type = str,
    nargs = '+')

parser_ApplyFilters.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option",
    type = str,
    required = True)

parser_ApplyFilters.add_argument('-f', "--filters", 
    help = "names of filters to apply. One or more of: genome_repeats provided by the Repeats option; rearrangements provided by the Structure option.",
    type = str,
    nargs = '+',
    # these choices must correspond to known_filters in CallVariants.Filter.__init__()
    # or to other_filters which are listed by variant calling program, currently just GATK
    choices = ['genome_repeats','rearrangements','GATK'],
    required = True)

parser_ApplyFilters.add_argument('-s', "--include_samples", 
    help = "restrict filtering to these samples. If omitted, plots for all samples are produced.",
    type = str,
    nargs = '+')

parser_ApplyFilters.add_argument('-r', "--report", 
    help = "summarise the effect of filters. Only choice currently implemented is cumulative effect on total variants.",
    type = str,
    nargs = '+',
    choices = ['cumulative'])

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

mutually_exclusive_group = parser_ComparativeAnalysis.add_mutually_exclusive_group(required=True)

mutually_exclusive_group.add_argument('-m', "--build_MSA", 
    help = "build a multiple sequence alignment from a reference genome and SNPs listed in VCF files",
    action = 'store_true')

mutually_exclusive_group.add_argument('-i', "--infer_phylogeny", 
    help = "infer a phylogeny from a multiple sequence alignment",
    action = 'store_true')

mutually_exclusive_group.add_argument('-r', "--infer_recombination", 
    help = "infer recombination from a phylogeny and multiple sequence alignment",
    action = 'store_true')

mutually_exclusive_group.add_argument('-p', "--plot_phylogeny", 
    help = "plot a phylogeny, possibly including recombination inferences",
    action = 'store_true')

# for build_MSA
parser_ComparativeAnalysis.add_argument('-g', "--genome_name", 
    help = "name of genome obtained by the CollectData option on which to base a MSA along with SNPs",
    type = str)

parser_ComparativeAnalysis.add_argument('-n', "--reads_name", 
    help = "name of read datasets groups to include in an MSA. Should match --reads_name option used previously in baga GATK-based read calling pipline",
    type = str,
    nargs = '+')

parser_ComparativeAnalysis.add_argument('-s', "--include_samples", 
    help = "restrict MSA to these samples. If omitted, plots for all samples are included.",
    type = str,
    nargs = '+',
    metavar = 'SAMPLE_NAME')

parser_ComparativeAnalysis.add_argument('-l', "--include_invariants", 
    help = "include invariant sites to make a Long alignment. Invariant sites are required to accurately estimate some parameters",
    action = 'store_true')

parser_ComparativeAnalysis.add_argument('-c', "--core_only", 
    help = "only include sites present in all samples",
    action = 'store_true')

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

args = parser.parse_args()



# if args.verbose:
    # print("verbosity turned on")

# if args.quiet:
    # print("output suppressed")


if hasattr(args, 'GATK_jar_path') and args.GATK_jar_path:
    # check Java is installed and version
    try:
        p = subprocess.Popen(['java', '-version'], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    except OSError:
        sys.exit('Could not find the Java runtime environment, needed for GATK and Picard. Please ensure you have Java installed and that the "java" executable is in the system path')
    
    try:
        java_version_string = p.stdout.read().split('\n')[0]
        maj, min1, min2 = re.findall('([0-9])\.([0-9])\.([0-9_]+)', java_version_string)[0]
        maj, min1 = map(int,[maj, min1])
        if maj == 1 and min1 == 7:
            print('Java 1.7: found!')
        else:
            sys.exit('GATK v3.3 requires Java v1.7 but your version is {}.{}. Please install Java v1.7 to continue'.format(maj,min1))
    except IndexError:
        print(output)
        sys.exit('There was a problem checking your Java version. Please report this as a baga bug.')
    
    # check GATK jar file
    p = subprocess.Popen(['java', '-Xmx8g', '-jar', args.GATK_jar_path, '--help'], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    output = p.stdout.read()
    if 'Error: Unable to access jarfile' in output:
        sys.exit('Java could not find jar file for GATK at: {}\nPlease provide path to GATK v3.3 jar file with --GATK_jar_path.'.format(args.GATK_jar_path))
    else:
        try:
            maj, min1, min2, min3 = re.findall('\(GATK\) v([0-9])\.([0-9]+)-([0-9]+)-([a-z0-9]+),', output.split('\n')[1])[0]
            maj, min1, min2 = map(int, [maj, min1, min2])
            if maj == 3 and min1 == 3:
                print('GATK v3.3: found!')
            elif maj == 3 and min1 > 3:
                print('Found GATK v{}.{}. Should be OK but this version of baga was tested with v3.3 . . .'.format(maj,min1))
            elif maj == 3 and min1 < 3:
                print('WARNING: Found GATK v{}.{}. Not tested in this baga version (v3.3 was) . . . use at your own risk!')
            else:
                sys.exit('Expecting GATK version 3.3 but found {}.{}\nPlease provide path to GATK v3.3 jar file with --GATK_jar_path.'.format(maj,min1))
        except IndexError:
            print(output)
            sys.exit('There was a problem checking your --GATK_jar_path argument. Please report this as a baga bug.')

### Check Dependencies ###


if args.subparser == 'Dependencies':
    print('\n-- Dependencies check/get module --')
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
    
    if args.checkget:
        check_summary = []
        check_results = []
        for name in args.checkget:
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
                
                gotnow = check(name)
                if gotnow:
                    check_summary[-1] += "Installed successfully: found!"
                else:
                    check_summary[-1] += "Failed to install . . . there may be dependencies missing or some other problem . . ."
                    check_summary[-1] += ". . . currently python modules like dendropy cannot be checked straight after getting . . . so rerunning this command may succeed for those: try that."
                
                check_results += [gotnow]
        
        if len(check_summary): # > 1:
            print(''.join(['\n\nSummary:\n'] + sorted(check_summary))+'\n')
        
        if not all(check_results):
            sys.exit('\n\nOne or more baga dependencies are unavailable and could not be installed . . . \n')
        

### Download Genomes ###

if args.subparser == 'CollectData':
    print('\n-- Data collection module --')
    if args.genomes is not None:
        from baga import CollectData
        # work out what accessions and/or paths to genbank files were provided
        if ',' in args.genomes:
            # allow for single string with ',' delimiter
            use_genome = args.genomes.split(',')
        else:
            use_genome = args.genomes
        
        load_gbks = []
        load_bagas = []
        collect_accessions = []
        for g in use_genome:
            genome_no_quotes = g.strip('"').strip("'")
            if re.search('\.gbk|\.GBK', genome_no_quotes):
                load_gbks += [genome_no_quotes]
            elif re.search('\.baga|\.baga', genome_no_quotes):
                load_bagas += [genome_no_quotes]
            else:
                # anything that isn't .gbk or .baga is assumed to be an accession number
                collect_accessions += [genome_no_quotes]
        
        # check email address provided for Entrez, if any accession found
        if args.email_address is None and len(collect_accessions) > 0:
            print(textwrap.fill('User email address is required for downloading from NCBI \
with accession numbers. Detected %s items without ".gbk" or ".GBK" and assumed to be \
accession numbers' % (len(collect_accessions)), text_width))
            sys.exit(1)
        
        loaded_genomes = {}
        
        for gbk in load_gbks:
            print('Loading {}'.format(gbk))
            genome = CollectData.Genome(local_path = gbk, format = 'genbank')
            print('Storing for future use . . .')
            genome.saveLocal()
            #loaded_genomes[genome.id] = genome
        
        for baga in load_bagas:
            print('Loading {}'.format(baga))
            genome = CollectData.Genome(local_path = gbk, format = 'baga')
            print('Storing for future use . . .')
            genome.saveLocal()
            #loaded_genomes[genome.id] = genome
        
        for genome_accession in collect_accessions:
            print('Fetching {}'.format(genome_accession))
            genome = CollectData.Genome(accession = genome_accession, user_email = args.email_address)
            print('Storing for future use . . .')
            genome.saveLocal()
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
        
        use_name = args.reads_name.replace('baga.CollectData.Reads-', '' , 1).replace('.baga', '')
        
        from baga import PrepareReads
        
        if args.subsample_to_cov is not None:
            # make a new CallVariants.Reads object
            print('Loading reads group %s' % use_name)
            downloaded_reads = baga.bagaload('baga.CollectData.Reads-%s' % use_name)
            reads = PrepareReads.Reads(downloaded_reads)
            read_cov_depth, genome_size = args.subsample_to_cov
            print('Subsampling reads group {} to {}x coverage for a {:,} bp genome'.format(use_name, read_cov_depth, genome_size))
            reads.subsample(genome_size, read_cov_depth, force = args.force)
            reads.saveLocal(use_name)
        
        if args.adaptors is not None:
            # need to test parallelism with lots of stdout reports
            if args.subsample_to_cov is None:
                # didn't subsample in current analysis
                if args.adaptors == 'fullsample':
                    print('Not using previously subsampled reads because "--adaptors fullsample" supplied.')
                    print('Loading collected reads group %s' % use_name)
                    downloaded_reads = baga.bagaload('baga.CollectData.Reads-%s' % use_name)
                    reads = PrepareReads.Reads(downloaded_reads)
                else:
                    print('Loading subsampled reads group %s' % use_name)
                    reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name)
            
            try:
                reads.cutAdaptors(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('cutadapt')
            
            reads.saveLocal(use_name)
        
        if args.trim:
            # could check whether adaptor cut read files exist
            # need to test parallelism with lots of stdout reports
            if not args.adaptors:
                # load a previously adaptor cut reads set
                print('Loading processed reads group %s' % use_name)
                reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name)
            
            print('\nTrimming reads . . .')
            try:
                reads.trim(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('sickle')
            
            reads.saveLocal(use_name)
        
        if args.delete_intermediates:
            print('Checking on intermediate fastq files to delete . . .')
            if not args.adaptors and not args.trim:
                # load a previously adaptor cut reads set
                print('Loading processed reads group %s' % use_name)
                reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name)
            
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

if args.subparser == 'AlignReads':
    print('\n-- Read Aligning module --')
    if args.reads_name is not None and args.genome_name is not None:
        # first check whether GATK path is needed
        if args.indelrealign and not args.GATK_jar_path:
            
            print('''Please supply:

--GATK_jar_path

if using:

--indelrealign

''')
            sys.exit(1)
        
        use_name_reads = args.reads_name.replace('baga.PrepareReads.Reads-', '' , 1).replace('.baga', '')
        use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
        alns_name = '__'.join([use_name_reads, use_name_genome])
        
        import baga
        
        from baga import AlignReads
        from baga import CollectData
        
        if args.align:
            print('Loading processed reads group %s' % use_name_reads)
            prepared_reads = baga.bagaload('baga.PrepareReads.Reads-%s' % use_name_reads)
            print('Loading genome %s' % use_name_genome)
            genome = CollectData.Genome(local_path = 'baga.CollectData.Genome-{}.baga'.format(use_name_genome), format = 'baga')
            # create alignment object
            alignments = AlignReads.SAMs(reads = prepared_reads, genome = genome)
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
            alignments.saveLocal(alns_name)
        
        
        if args.deduplicate:
            if not args.align:
                # add an exception here and inform to use --align first
                print('Loading previously processed read alignments: %s' % alns_name)
                alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}.baga'.format(alns_name))
                
            try:
                alignments.removeDuplicates(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('picard')
            
            try:
                alignments.sortIndexBAMs(force = args.force, max_cpus = args.max_cpus)
            except OSError:
                exe_fail('samtools')
            
            alignments.saveLocal(alns_name)
        
        if args.indelrealign:
            if not args.deduplicate:
                # add an exception here and inform to use --align first
                print('Loading previously processed read alignments: %s' % alns_name)
                alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}.baga'.format(alns_name))
            
            try:
                os.makedirs('genome_sequences')
            except OSError:
                pass
            
            genome_fna = 'genome_sequences/%s.fna' % alignments.genome_id
            
            if not os.path.exists(genome_fna):
                from Bio import SeqIO
                from Bio.Seq import Seq
                from Bio.SeqRecord import SeqRecord
                SeqIO.write(SeqRecord(Seq(alignments.genome_sequence.tostring()), id = alignments.genome_id), 
                            genome_fna, 
                            'fasta')
            
            alignments.IndelRealignGATK(
                        jar = args.GATK_jar_path.split(os.path.sep), 
                        force = args.force, 
                        max_cpus = args.max_cpus)
            
            alignments.saveLocal(alns_name)
        
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

### Repeats ###

if args.subparser == 'Repeats':
    print('\n-- Chromosome repeats detection module --')
    
    e = '-i/--minimum_percent_identity must be between 0 and 100 percent (low values not recommended!)'
    assert 0 < args.minimum_percent_identity <= 100, e
    
    import baga
    from baga import Repeats
    from baga import CollectData
    
    use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
    
    if args.find:
        print('Loading genome %s' % use_name_genome)
        genome = CollectData.Genome(local_path = 'baga.CollectData.Genome-{}.baga'.format(use_name_genome), format = 'baga')
        finder = Repeats.Finder(genome)
        # minimum_percent_identity defaults to 98%, argument takes 0.98 so *0.01
        # minimum_repeat_length defaults to 400
        finder.findRepeats(minimum_percent_identity = args.minimum_percent_identity * 0.01, minimum_repeat_length = args.minimum_repeat_length)
        finder.saveLocal(use_name_genome)
        # also save just the ranges for filtering
        baga.bagasave(finder.ambiguous_ranges, 'baga.Repeats.filter_regions-{}'.format(use_name_genome))
    
    if args.plot:
        # if not args.find:
            # print('Loading repeats and genome: %s' % use_name_genome)
            # finder = baga.bagaload('baga.Repeats.Finder-%s' % use_name_genome)
        
        Repeats.plotRepeats(use_name_genome, outdir = ['plots_repeats'], force = True)
    
    if args.summarise:
        loaded = False
        try:
            ambiguous_ranges = baga.bagaload('baga.Repeats.filter_regions-{}'.format(use_name_genome))
            loaded = True
            
        except IOError:
            print('Could not load baga.Repeats.filter_regions-{}'.format(use_name_genome))
            print('You probably need to find the repeats first with the --find option')
        
        if loaded:
            total = 0
            with open('baga.Repeats.filter_regions-{}.csv'.format(use_name_genome), 'w') as fout:
                fout.write('"length (bp)","start","end"\n')
                for s,e in ambiguous_ranges:
                    print('{:,} bp from {:,} to {:,}'.format(e-s+1, s, e))
                    total += e-s+1
                    fout.write('"{}","{}","{}"\n'.format(e-s+1, s, e))
            
            print('\n{} repetitive, ambiguous regions spanning {:,} bp of chromosome\n'.format(len(ambiguous_ranges), total))

### Structure ###

if args.subparser == 'Structure':
    print('\n-- Chromosome sequence rearrangement detection module --\n')
    import baga
    from baga import Structure
    from baga import CollectData
    
    if not(args.check or args.plot or args.plot_range):
        parser_Structure.error('Need at least one of --check/-c, --plot/-p or --plot_range/-r')
    
    if args.check or args.plot or args.plot_range:
        # collect BAMs
        if args.reads_name:
            # baga pipeline information provided
            if not args.genome_name:
                parser.error('--genome_name/-g is required with --reads_name/-n. (The baga CollectData-processed genome used with the AlignReads option)')
            
            use_name_group = args.reads_name.replace('baga.AlignReads.SAMs-', '' , 1).replace('.p.gz', '')
            use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
            print('Loading alignments information for: {}__{} from AlignReads output'.format(use_name_group, use_name_genome))
            from baga import AlignReads
            alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}__{}.baga'.format(use_name_group, use_name_genome))
            sample_names = sorted(alignments.read_files)
            BAMs = alignments.ready_BAMs[-1]
            
        elif args.alignments_path:
            # list of folders or files provided
            BAMs = []
            for path in args.alignments_path:
                if os.path.isdir(path):
                    path_contents = os.listdir(path)
                    theseBAMs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('BAM', 'bam')]
                    e = 'No BAM files (*.bam or *.BAM) found in:\n{}'.format(args.alignments_path)
                    assert len(theseBAMs), e
                    BAMs += theseBAMs
                else:
                    # add file
                    BAMs += [path]
    
    # check on requested samples: crop BAMs accordingly
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
    
    if args.check:
        # check these genome-aligned read sets
        checkers = Structure.checkStructure(BAMs, mean_param = 5, min_mapping_quality = 5, resolution = 10, step = 1)
    
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
            
            if args.plot:
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
                
            else:
                # just the requested range
                do_ranges = [args.plot_range]
            
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

### Call Variants ###

if args.subparser == 'CallVariants':
    print('\n-- Variant Calling module --\n')
    if args.reads_name is not None and args.genome_name:
        # first check whether GATK path is needed
        if any([args.calleach, 
                args.calljoint, 
                args.hardfilter, 
                args.recalibrate
                ]) and not args.GATK_jar_path:
            
            print('''Please supply:

--GATK_jar_path

if using any of:

--calleach
--calljoint
--hardfilter
--recalibrate

''')
            sys.exit(1)
        
        use_name_reads = args.reads_name.replace('baga.AlignReads.Reads-', '' , 1).replace('.baga', '')
        use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
        alns_name = '__'.join([use_name_reads, use_name_genome])
        
        import baga
        from baga import CallVariants
        
        if os.path.exists('baga.CallVariants.Caller-{}.baga'.format(alns_name)) and not args.new:
            print('Loading existing variants call analysis for: {}'.format(alns_name))
            print('(use --new to start variant calling again)')
            caller = CallVariants.Caller(baga = 'baga.CallVariants.Caller-{}.baga'.format(alns_name))
        elif os.path.exists('baga.CallVariants.Caller-{}.baga'.format(alns_name)) and args.new:
            print('Starting new variant calling analysis bacause --new requested. Will overwrite any previous analyses.')
            from baga import AlignReads
            print('Loading alignments information for: {} from AlignReads output'.format(alns_name))
            alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}.baga'.format(alns_name))
            caller = CallVariants.Caller(alignments = alignments)
        elif not os.path.exists('baga.CallVariants.Caller-{}.baga'.format(alns_name)):
            print('Starting new variant calling analysis. Will overwrite any previous analyses.')
            print('(could not find baga.CallVariants.Caller-{}.baga)'.format(alns_name))
            from baga import AlignReads
            print('Loading alignments information for: {} from AlignReads output'.format(alns_name))
            alignments = AlignReads.SAMs(baga = 'baga.AlignReads.SAMs-{}.baga'.format(alns_name))
            caller = CallVariants.Caller(alignments = alignments)
        
        # because --new, setting --force to true to overwrite each output file
        if args.new:
            print('Because --new, also setting --force to overwrite each output file')
            args.force = True
        
        if args.calleach:
            caller.CallgVCFsGATK(
                        mem_num_gigs = args.max_memory,
                        jar = args.GATK_jar_path.split(os.path.sep),
                        force = args.force, 
                        max_cpus = args.max_cpus)
            
            caller.saveLocal(alns_name)
        
        if args.calljoint:
            caller.GenotypeGVCFsGATK(
                        alns_name,
                        jar = args.GATK_jar_path.split(os.path.sep), 
                        # ultimately, scale this by the number of samples involved
                        # needed >8 for a set of 40
                        mem_num_gigs = args.max_memory, 
                        force = args.force)
            
            caller.saveLocal(alns_name)
        
        if args.hardfilter:
            caller.hardfilterSNPsGATK(
                        jar = args.GATK_jar_path.split(os.path.sep),
                        force = args.force)
            
            caller.hardfilterINDELsGATK(
                        jar = args.GATK_jar_path.split(os.path.sep),
                        force = args.force)
            
            caller.saveLocal(alns_name)
        
        if args.recalibrate:
            # this is slow!
            caller.recalibBaseScoresGATK(
                        jar = args.GATK_jar_path.split(os.path.sep),
                        force = args.force, 
                        mem_num_gigs = args.max_memory, 
                        max_cpus = args.max_cpus)
            
            caller.saveLocal(alns_name)

### Apply Filters ###

if args.subparser == 'ApplyFilters':
    print('\n-- Apply Filters (part of the Variant Calling module) --\n')
    
    ## to apply variants, provide one reads group name
    if len(args.reads_name) > 1 and not args.report:
        sys.exit('Filters can only be applied to one group of reads at a time. Multiple sets can be handled with the --report option, though.')
    
    ## to report effects of filters, provide one or more read group names along with --report 
    
    from baga import CallVariants
    
    use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
    
    VCFs_for_report = {}
    
    # collect VCFs
    if args.reads_name:
        for these_reads in args.reads_name:
            # baga pipeline information provided
            # allow for multiple rounds of recalibration at end of CalVariants i.e., 1 or 2 and select 2 if available
            import baga
            use_name_reads = these_reads.replace('baga.AlignReads.Reads-', '' , 1).replace('.baga', '')
            alns_name = '__'.join([use_name_reads, use_name_genome])
            
            filein = 'baga.CallVariants.Caller-{}.baga'.format(alns_name)
            caller = CallVariants.Caller(baga = filein)
            
            if hasattr(caller, 'path_to_hardfiltered_SNPs') and hasattr(caller, 'path_to_hardfiltered_INDELs'):
                # only one for VCFs so no overwriting
                VCFs = [caller.path_to_hardfiltered_SNPs[-1], caller.path_to_hardfiltered_INDELs[-1]]
                # more than one can be handled with --report though
                # only implemented for separate VCFs currently
                VCFs_for_report[these_reads] = {'SNPs':caller.path_to_hardfiltered_SNPs[-1], 'InDels':caller.path_to_hardfiltered_INDELs[-1]}
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
            
            # old code for manually collecting baga pipelined VCFs
            # collected = defaultdict(list)
            # pattern = re.compile('{}__{}_([0-9])_samples_hardfiltered_[A-Za-z]+.vcf$'.format(use_name_reads, use_name_genome))
            # files = os.listdir(os.path.sep.join(['variants', use_name_genome]))
            # for f in files:
                # m = re.match(pattern, f)
                # if m is not None:
                    # collected[m.groups()[0]] += [os.path.sep.join(['variants', use_name_genome, f])]
            
            # eg_path = "{}{}__{}_samples_hardfiltered_SNPs.vcf".format(os.path.sep.join(['variants', use_name_genome,'']), use_name_reads, use_name_genome)
            # e = "Could not find any VCFs in the expected location e.g., {}".format(eg_path)
            # assert len(collected) > 0, e
            # VCF_paths = [files for n,files in sorted(collected.items())][-1]
        
    else:
        # list of folders or files provided in args.vcf_paths
        VCFs = []
        for path in args.vcfs_path:
            if os.path.isdir(path):
                path_contents = os.listdir(path)
                theseVCFs = [os.path.sep.join([path,f]) for f in path_contents if f[-3:] in ('VCF', 'vcf')]
                e = 'No VCF files (*.vcf or *.VCF) found in:\n{}'.format(args.alignments_path)
                assert len(theseVCFs), e
                VCFs += theseVCFs
            else:
                # add file
                VCFs += [path]

    print('Loaded VCF locations:\n{}'.format('\n'.join(VCFs)))
    
    # either do a report or go on and actually do filtering
    if args.report:
        for reporttype in args.report:
            print('{} == cumulative'.format(reporttype))
            if reporttype == 'cumulative':
                filter_order = args.filters
                print('Reporting . . .: {}, {}'.format(args.report,'+'.join(filter_order)))
                CallVariants.reportCumulative(filter_order, use_name_genome, VCFs_for_report)
    else:
        # check accessible and collect sample names
        sample_names = set()
        for VCF in VCFs:
            try:
                with open(VCF, 'r') as filein:
                    header, header_section_order, colnames, variants = CallVariants.parseVCF(VCF)
                    sample_names.update(colnames[9:])
                    pass
            except IOError as e:
                print(e)
                sys.exit('Failed to open: {}'.format(VCF))
        
        sample_names = sorted(sample_names)
        
        # genome is only needed when generating csv summary with ORFs affected
        from baga import CollectData
        print('Loading genome %s' % use_name_genome)
        genome = CollectData.Genome(local_path = 'baga.CollectData.Genome-{}.baga'.format(use_name_genome), format = 'baga')
        
        print('Loading filter information . . .')
        
        filters = {}
        if 'genome_repeats' in args.filters:
            from baga import Repeats
            filein = 'baga.Repeats.FinderInfo-{}.baga'.format(use_name_genome)
            finder_info = Repeats.loadFinderInfo(filein)
            filters['genome_repeats'] = finder_info['ambiguous_ranges']
            # print a summary
            t = len(filters['genome_repeats'])
            s = sum([(e - s) for s,e in filters['genome_repeats']])
            print('Reference genome sequence {} contains {} repeated regions spanning {:,} basepairs'.format(use_name_genome, t, s))
        
        if 'rearrangements' in args.filters:
            
            from baga import Structure
            
            filters['rearrangements'] = {}
            #for sample,checker in checkers.items():
            for sample in sample_names:
                
                if args.include_samples:
                    if sample not in args.include_samples:
                        continue
                
                filein = 'baga.Structure.CheckerInfo-{}__{}.baga'.format(sample, use_name_genome)
                checker_info = Structure.loadCheckerInfo(filein)
                
                e = 'Genome name for checker {} ({}) does not match supplied genome name ({})'.format(filein, checker_info['genome_name'], genome.id)
                assert checker_info['genome_name'] == genome.id, e
                
                # report brief summary for this sample
                print('For sample {}, relative to {}:'.format(sample, genome.id))
                #t = len(checker_info['suspect_regions']['high non-proper pairs'])
                t = len(checker_info['suspect_regions']['rearrangements'])
                s = sum([(e - s) for s,e in checker_info['suspect_regions']['rearrangements']])
                print('  {} regions spanning {:,} basepairs are affected by rearrangements thus having ambiguous 1:1 orthology.'.format(t, s))
                
                #t = len(checker_info['suspect_regions_extensions']['high non-proper pairs'])
                t = len(checker_info['suspect_regions']['rearrangements_extended'])
                s = sum([(e - s) for s,e in checker_info['suspect_regions']['rearrangements_extended']])
                print('  {} additional regions spanning {:,} basepairs are adjacent to the above regions but have a >50% zero read depth over a moving window (i.e., any aligned reads have a patchy distribution, are usually rare). These are typically large deletions including missing prophage and genomic islands\n'.format(t, s))
                
                filters['rearrangements'][sample] = checker_info['suspect_regions']
        
        filter_applier = CallVariants.Filter(VCFs, genome)  #, use_name_reads)
        filter_applier.doFiltering(filters)


                
        

if args.subparser == 'ComparativeAnalysis':
    print('\n-- Comparative Analyses --\n')
    
    from baga import AlignReads
    from baga import CallVariants
    from baga import ComparativeAnalysis
    from baga import CollectData
    
    if args.build_MSA:
        # baga pipeline information provided <== only way currently implemented
        # allow for multiple rounds of recalibration at end of CalVariants i.e., 1 or 2 and select 2 if available
        use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
        import baga
        path_to_SNPs_VCFs = []
        path_to_InDels_VCFs = []
        paths_to_BAMs = []
        for these_reads in args.reads_name:
            print('Collecting VCFs for {}'.format(these_reads))
            use_name_reads = these_reads.replace('baga.AlignReads.Reads-', '' , 1).replace('.baga', '')
            alns_name = '__'.join([use_name_reads, use_name_genome])
            
            filein = 'baga.CallVariants.Caller-{}.baga'.format(alns_name)
            caller = CallVariants.Caller(baga = filein)
            
            if hasattr(caller, 'path_to_hardfiltered_SNPs') and hasattr(caller, 'path_to_hardfiltered_INDELs'):
                checkthis = caller.path_to_hardfiltered_SNPs[-1]
                try:
                    with open(checkthis) as fin:
                        path_to_SNPs_VCFs += [checkthis]
                        print('Found: {}'.format(checkthis))
                except IOError:
                    print('Could not find: {}'.format(checkthis))
                    sys.exit('You may need to rerun the analysis that should have generated that file')
                
                checkthis = caller.path_to_hardfiltered_INDELs[-1]
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
        
        print('Loading genome %s' % use_name_genome)
        genome = CollectData.Genome(local_path = 'baga.CollectData.Genome-{}.baga'.format(use_name_genome), format = 'baga')
        
        MSA_builder = ComparativeAnalysis.MultipleSequenceAlignment(path_to_SNPs_VCFs, path_to_InDels_VCFs)
        MSA_builder.parseVCFs()
        #MSA_builder.parseVCFs(filters_include = [])
        #MSA_builder.parseVCFs(filters_exclude = ['genome_repeats'], filters_include = [])
        #MSA_builder.getCoverageRanges(paths_to_BAMs, BAM_suffix = 'realn.bam')
        MSA_builder.getCoverageRanges(paths_to_BAMs)
        MSA_builder.writeMSA(   '{}__{}_SNPs'.format(use_name_genome,'_'.join(args.reads_name)), 
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
        phylo_analyser = ComparativeAnalysis.Phylogenetics(args.path_to_MSA, args.path_to_tree)
        # bit of a fudge dealing with rooted tree . . . to be improved
        phylo_analyser.collectPhyMLstats(path_to_stats_file = args.path_to_tree.replace('_rooted','').replace('_phyml_tree','') + '_phyml_stats')
        phylo_analyser.infer_recombination(bootstraps = args.num_bootstraps) #, output_suffix = '_rooted')
    
    if args.plot_phylogeny:
        # should check either or etc here
        if args.genome_name:
            use_name_genome = args.genome_name.replace('baga.CollectData.Genome-', '' , 1).replace('.baga', '')
            print('Loading genome %s' % use_name_genome)
            genome = CollectData.Genome(local_path = 'baga.CollectData.Genome-{}.baga'.format(use_name_genome), format = 'baga')
            genome_length = len(genome.sequence)
        elif args.genome_length:
            genome_length = genome_length
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
                            args.path_to_tree.replace('_rooted',''),   # deal with rootedness at some point
                            plot_output_path,
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
        
        phylo_plotter.doPlot(outgroup_label_list = ['FM209186.1'], 
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
        
