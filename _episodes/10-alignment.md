---
Title: "Read Alignment and Small Variant Calling"
---
## Getting Started
Log in to your instances through Jupyter notebooks and launch a terminal.

Next, find the places in the directories where you are going to work.

All of the data you will need for these exercises are in
```bash
/data/alignment
```

There are multiple subdirectories inside that directory that you will access for various parts of this session.
Part of the aim of these exercises is for you to become more familiar with working with files and commands in addition to learning the processing steps, 
so for the most part there will not be full paths to the files you need, but all of them will be under the alignment directory.

You CAN put your output anywhere you want to. That's also part of the workshop.
I would recommend you use the following space:

```bash
/workshop/output
```

I would recommend making a directory under that for everything to keep it separate from what you did yesterday or will do in the following sessions.
You could do something like

```bash
mkdir /workshop/output/alignment
```

However, if you think you will find it confusing to have "alignment" in both in the input and output directory paths, 
you could try something different like align_out or alnout or whatever you want. It's only important that you remember it.

Part of this work will require accessing all the various data and feeding it into the programs you will run and directing your output.
Certain steps will require that you be able to find the output of previous steps to use it for the input to the next steps.
Keeping your work organized this way is a key part of running multi-step analysis pipelines.
Yesterday we discussed absolute and relative paths. During this session you will need to keep track of where you are in the filesystem (pwd),
how to access your input data from there, and how to access your output files and directories from there.

It is up to you if you want to do everything with abslute paths, with relative paths, or a mix. It is up to you if you want to cd to the output directory,
the input directory, or a third place. Do whichever feels most comfortable to you. Note that if you choose to do everything with absolute paths,
your working directory does not matter, but you will type more.

Tab completion is your friend. If there are only three things I want you to remember from this workshop, they are tab completion, tab completion, and tab completion.

## Copy and Paste
Copy and paste is not your friend. It may be possible to copy and paste commands from this markdown into the Jupyter terminal.
Please note that this more often than not fails in general because text in html or pdf or Word documents has non-printing characters or non-standard 
punctuation (like em dashes and diretional quotation marks) which are not processed by the unix shell and will generate errors.
The same is true with copying text that spans multiple lines. Because of this, most of the code blocks in this section are intentionally constructed to
prevent copying and pasting. That is, the commands will be constructed, but full paths will not be given and real names of files will be substituted with placeholders.
As a convention, if a name is in *UPPERCASE ITALIC* or as \_UPPERCASE\_, it is the placeholder for a filename or directory name you need to provide. If you copy and paste this,
it will fail to work.

## The exercises we will do
What we will do today is take some data from the Genome in a Bottle sample, also part of 1000 Genomes, NA12878, and this person's parents, NA12891 and NA12892,
and align them to the human genome. Aligning a whole human genome would take hours, and analyzing it even longer, so I have extracted a subset of reads that
align to only chromosome 20, which is about 2% of the genome and will align in a reasonable amount of time. We will still align to the entire genome,
because you should always align to the whole genome even if you think you have only targeted a piece of it because there are always off target reads,
and using the whole genome for alignment improves alignment quality. We will walk through the following steps:

1. Build an index for the bwa aligner (on a small part of the genome)
2. Align reads to the genome with bwa mem
3. Merge multiple alignment files and sort them in coordinate order (with gatk or samtools)
4. Mark duplicate reads with gatk
5. Call variants using gatk HaplotypeCaller (on a small section of chromosome 20)
6. Filter variants using gatk VariantFiltration

We will also skip a couple of steps, indel realignment and base quality score recalibration (BQSR), but we will discuss them.

Along the way we will discuss the formatting and uses of several file formats:
- fasta
- fastq
- sam
- bam
- cram
- gatk dup metrics
- vcf
- gvcf

## Running bwa
bwa is the Burrow-Wheeler Aligner. We're going to use it to align short reads.

First we're going to see an example of building the bwt index.
Indexing an entire human genome would again take an hour or more,
so we will just index chromosome 20. Look for it in the references subdirectory of /data/alignment
and make a copy of it to your workspace using cp then index it with bwa. This should take about 1 minute.

```bash
cp /data/alignment/references/chr20.fa _MYCHR20COPY_
bwa index -a bwtsw _MYCHR20COPY_
```

The -a bwtsw tells bwa which indexing method to use. It also supports other indexing schemes that are faster for small genomes
(up to a few megabases, like bacteria or most fungi) but do not scale to vertebrate sized genomes. You can see more detail about this
by running bwa index with no other arguments, which produces the usage.

If you now look in the directory where you put _MYCHR20COPY_, you will now see a bunch of other files
with the same base or root name as your chromosome fa, but with different extensions.
When you use the aligner, you tell it the base file name of your fasta and it will find all the indices
automatically.

Normally we do not want to make copies of files to do work, but in this case you had to because indexing
requires creating new files and you do not have write permission to /data/alignment/references, so you cannot make the index files (try it if you want).

We will not actually use this index, so feel from to delete these index files or move them out of the way. This was just so you can see how that step works.

Now we will actually run bwa on the chr20 data, but against the whole reference genome, which is in

```bash
/data/alignment/references/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

The reads from chromosome 20 are in

```bash
/data/alignment/chr20/NA12878.fq.dir/
```

We will stop here and take a look at each of these file formats.

```bash
less /data/alignment/references/GRCh38_full_analysis_set_plus_decoy_hla.fa
less /data/alignment/chr20/NA12878.fq.dir/NA12878_TTGCCTAG-ACCACTTA_HCLHLDSXX_L001_1.fastq
```
Now we will align the reads from one pair of fastqs to the reference. The command is below, but please read the rest of this section before you start running this.

```bash
bwa mem -K 100000000 -Y -t 8 -o _SAMFILE1_ _REFERENCE_ _READ1FASTQ_ _READ2FASTQ_
```

To break down that command line:
- bwa mem launches the bwa minimal exact match algorithm
- -K 100000000 tells bwa to process reads in chunks of 100M, regardless of the number of threads, which mainly forces it to produce deterministic output
- -Y says to use soft instead of hard clipping for supplementary alignments
- -t 8 tells it to use 8 threads, which maximally utilizes your 8 vCPU cloud instances
- -o _SAMFILE1_ tells it to put its output (in SAM format) in _SAMFILE1_, which by convention ends in ".sam"
- _READ1FASTQ_ and _READ2FASTQ_ are the paired fastq files

By default, bwa sends output to stdout (which is your screen, unless you redirected it), which is not very useful by itself,
but allows you to pipe it into another program directly so that you do not have write an intermediate alignment file. In our case,
we want to be able to look at and manipulate that intermediate file, but it can be more efficient not to write all that data to
disk just to read it back in and throw the sam file away later.

bwa can take either unpaired reads, in which there is only one read file argument, or paired reads, in which case there are two as here,
or if you have interleaved fastq (with read 1s and read 2s alternating in the same file), it will take one fastq, with the -p option.
Note that the fastqs come after the reference because there is always only one reference but could be one or two fastqs.

When you pick the fastqs, pick a pair that have the same name _except_ for _1.fastq and _2.fastq at the end.

There is another option you can add. We do not strictly need it here, but you would want to do this in a production environment,
and files you get from a sequencing center already aligned ought to have this done. This is to add read group information to the sam. _Read groups_
are collection of reads that saw exactly the sam lab process (extaction, library construction, sequencing flow cell and lane)
and thus should share the same probabilities of sequencing errors or artifacts. Labeling these during alignment allows us to keep track of them downstream
if we want to do read group specific processing. The argument to bwa mem is -R _RG_INFO_ and looks like this:

```bash
-R @RG\tID:NA12878_TTGCCTAG-ACCACTTA_HCLHLDSXX_L001\tPL:illumina\tPM:Unknown\tLB:NA12878\tDS:GRCh38\tSM:NA12878\tCN:NYGenome\tPU:HCLHLDSXX.1.TTGCCTAG
```

If you want to format this for your reads, you can reconstruct the information from the filename.
- ID is the name of the read group as is the fastq filename up to but _not_ including the _1/_2.fastq
- PL is illumina
- PM is Unknown
- LB is NA12878
- DS is GRCh38
- SM is NA12878
- CN is NYGenome
- PU is the ID with NA12878_ chopped off (it's formatted slightly differently above, but for our purposes, it's the same information)

Adding the read groups is obviously a lot of effort and also error prone (I wrote a perl script to launch all the jobs and format this argument),
so if you do not want to do it, you can skip it.

The bwa command should run for about 3 minutes. Once you finish one set, please pick a _different_ set of fastqs and align those to the _same reference_
with a _different output file_. You will need at least 2 aligned sam files for the next steps. Your output sam file should be about 700Mb (for each run you do,
but will vary slightly by read group). If it is much smaller than this or the program ran for much less time, you have done something wrong.

We will pause here and take a look at the sam file format.

## Merging and Sorting
As you noticed, you now have multiple sam files. This is typical that we do one alignment for each lane of sequencing we do,
and we often do many lanes of sequencing because we pool many samples and run them across many flow cells.
