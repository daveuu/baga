# Bacterial and Archaeal Genome Analyser

## Novel analyses and wrapped tools pipelined for convenient processing of genome sequences

#### David Williams

## Introduction

**The Bacterial and Archaeal Genome Analyser (BAGA, pronounced "baga") is a commandline application and Python 2 package (3 coming soon) for diverse analyses of genome sequence data designed to facilitate reproducible research.**

**Input data** can be complete genome sequences and/or paired end short reads from Illumina sequencers, typically whole genome shotgun libraries. **Tasks** might include variant calling and resolving population structure of resequenced pathogen isolates, analysis of evolution experiments and comparative genomics including phylogenomics. Click [here](https://baga.readthedocs.org/en/latest/) to go straight to the documentation.

BAGA is a wrapper for proven third party tools<sup>1</sup>, but also includes novel algorithms for identifying chromosomal rearrangements and sequence repeats known to increase the likelihood of false positive variant calls<sup>2</sup>, the means to filter those probable false positive variants out of a dataset, the means to create custom pipelines for reproducible analyses<sup>3</sup>, and can generate various informative plots<sup>4</sup>. It is under active development: new features and much more documentation will be appearing shortly.

1. e.g. [BWA](http://bio-bwa.sourceforge.net/) for short read alignment to longer sequences, [GATK](https://www.broadinstitute.org/partnerships/education/broade/best-practices-variant-calling-gatk-1) for variant calling and [ClonalFrameML](https://github.com/xavierdidelot/clonalframeml) for homologous recombination inference
2. Variant calls in such regions are unreliable and should be filtered because conventional variant calling algorithms would be unaware of potential misalignments caused by the loss of homology and might therefore report false positive variant calls e.g., near chromosomal rearrangements caused by mobile genetic elements. Detailed characterisation of those regions can be made by local de novo assemblies of reads and alignment of resulting contigs to the reference sequence
3. researchers can make use of version-control and digital object identifiers to generate citable and reproducible analyses for peer review publication
4. BAGA can plot all automatically indicated regions such as those prone to misalignment of short reads because of structural differences between a reference sequence and a sampled genome, e.g. a missing prophage (see point 2 above) ![](docs/images/2689000_2691500_NC_011770.1__Liverpool__ERR953478.png)

Please see [the documentation](https://baga.readthedocs.org/en/latest/) for more details and step by step guides for performing various analysis pipelines and making your research more easily reproducible.

## Funding

Work on this software was started at [**The University of Liverpool**](https://www.liv.ac.uk), UK with funding from [**The Wellcome Trust**](http://www.wellcome.ac.uk/) (093306/Z/10) awarded to:

* Dr **Steve Paterson** (The University of Liverpool, UK)
* Dr **Craig Winstanley** (The University of Liverpool, UK)
* Dr **Michael A Brockhurst** (The University of York, UK)

License GPLv3+: GNU GPL version 3 or later. This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law.
