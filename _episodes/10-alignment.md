---
title: "Read Alignment and Small Variant Calling"
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
punctuation (like em dashes and directional quotation marks) which are not processed by the unix shell and will generate errors.
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
/data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

The reads from chromosome 20 are in

```bash
/data/alignment/chr20/NA12878.fq.dir/
```

We will stop here and take a look at each of these file formats.

```bash
less /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
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
They are also in sam file format, which is bulky, and we want bam instead. Lastly, they are not all sorted in coordinate order,
which is needed in order to index the files for easy use.

Thus, we want to merge and sort these files. There are two programs that will allow us to do this.
One is samtools. Samtools is written in C using htslib, a C library for interacting with sam/bam/cram format data.
It has the advantage of being fairly easy to use on the command line, running efficiently in C, and allowing data to
be piped in and out. It has the disadvantage of generally having less functionality than gatk, requiring more (but simpler)
steps to do something, and occasionally generating sam/bam files that are mildly non-compliant with the spec.
The other is GATK (the Genome Analysis ToolKit), written in Java using htsjdk.
GATK has the advantage of allowing you do multiple steps in one go
(as you may see here), offering many different functions, and being more generally compliant with specifications and
interoperable with other tools. It has the drawbacks of long and convoluted command lines that do not easily support
piping data (although with version 4 they have abstracted the interface to the JVM, which helps a lot with the complex usage).

First, we will look at how to do this with GATK:

```bash
gatk MergeSamFiles [-I _INPUTSAM_]... -O _MERGEDBAM_ -SO coordinate --CREATE_INDEX true
```

Note the [-I _INPUTSAM_]... portion. Because GATK does not know how many sam/bam files you want to merge in advance,
it makes you specify each one with the -I argument. This means you either need to tediously type this out by hand
(a good time to perhaps consider cd to the directory with the sam files) or programmatically assemble the command
line using a script outside the unix command line.

If you are behind on the previous step or want to make a bigger bam file, all the aligned sams can be found in

```bash
/data/alignment/chr20/NA12878.aln
```

They differ slightly in the flow cell names and the lane numbers. Tab completion will be your friend typing these out.

At the end of this, we will get one bam with all the sams merged, sorted in genomic coordinate order, and then indexed
(the .bai file next to the merged bam. The -SO coordinate argument to gatk tells it that in addition to whatever
else it is doing, it should sort the reads in coordinate order before spitting them out. The --CREATE_INDEX true
tells it to make an index for the new bam.

If you run this with all 12 sams, it should take 3-4 minutes (not counting the time to type the command).

We can also do this using samtools, but it takes multiple steps. We have to merge first, then sort, the index.

```bash
samtools merge _MERGEDBAM_ _INPUTSAM_...
```

Note that the syntax is a lot simpler. Also, because all the input files go at the end without any option/argument
flags, you can use a glob (\*) to send all the sam files in a directory without typing them out (e.g.,
/data/alignment/chr20/NA12878.aln/NA12878*sam), which is a lot easier to type. However, you will note there is 
no index, and if you looked inside, the reads not sorted. So we still need to do:

```bash
samtools sort -o _SORTEDBAM_ _INPUTBAM_
samtools index _SORTEDBAM_
```

Note that the -o is back. samtools sort again spits the output to stdout if you don't say otherwise. The index command
does not take an output because it just creates the index using the base name of the input. It does, however,
use the extension .bam.bai instead of just .bai. Both are valid.

The whole samtools thing takes about the same time as GATK's one action. In my test it also made a sorted bam that was about 80%
the size of GATK's, although both programs offer different levels of compression and the defaults may be different.

## Marking Duplicates
Sometimes you get identical read pairs (or single end reads with identical alignments, but those are more
likely to occur by random chance than identical pairs). These are usually assumed to be the result of artifical
duplication during the sequencing process, and we want to mark all but one copy as duplicates so they
can be ignored for downstream analyses. Again, both GATK and samtools do this, but even the author of samtools
says that GATK's marker is superior and you should not use samtools. The GATK argument looks like this:

```bash
gatk MarkDuplicates -I _INPUTBAM_ -O _DEDUPBAM_ -M _METRICS_ --CREATE_INDEX true
```

We give it the name of input and the name of our output. We also have to give a metrics file. This is a text file
and typically get the extension .txt. It contains informtion about the duplicate marking that can be useful.
We use --CREATE_INDEX because we will get a new bam with the dups marked, but we do not need -SO because we
are not reordering the bam.

If you are stuck up to this point, you can find a merged and sorted bam in /data/alignment/chr20/NA12878.aln.

We can take a look at the duplicate marked bam and the metrics file.

## NOT Realigning Indels
Some of you may have heard of a process called indel realignment (although maybe not). This used to be a part
of the GATK pipeline that would attempt to predict the presence of indels and adjust the alignments to prevent
false SNVs from being called due to misalignment of reads with indels. The HaplotypeCaller now does this internally
(but by default does not alter the bam), so we do not do it separately now. (It is also very time consuming.)

## Not Running Base Quality Score Recalibration (BQSR)
There is much discussion about whether this still has value at all. Mainly if you are running a large
project over several years, you may be forced to change sequencing platforms, and recalibration can help
with compatibility. Otherwise, it is probably not worthwhile.

We will not run it for several reasons:
- These reads went through a pipeline which binned quality, which means the recalibration will be very coarse
- These reads have already gone through BQSR once. Generally, running a second time is not terrible, but in the limit, infinitely recalibrating will eventually collapse all bases to a single quality score
- It is unclear whether quality recalibration will remain relevant going forward
- It takes a long time
- The command lines are really tedious to type (yes, even compared to MergeSamFiles)

## Calling Variants
We will call SNPs and indels using GATK's HaplotypeCaller. HaplotypeCaller takes a long time to run, so we can't even
do all of chromosome 20 in a reasonable amount of time. Instead, we will run a small piece of chromosome 20, from 61,723,571-62,723,570
(yes there is a reason I picked that). There are two ways you can do this. You can extract that region from the bam using samtools

```bash
samtools view -b -o _SMALLBAM_ _DEDUPBAM_ chr20:61723571-62723570
```

You can then run HaplotypeCaller on the _SMALLBAM_.

Alternatively, you can run HaplotypeCaller on the entire chr 20 bam, but tell it to only look in that region.

```bash
gatk HaplotypeCaller -O _VCF_ -L chr20:61723571-62723570 -I _BAM_ -R /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

Calling variants is only half the process. Now we need to filter them. We also do that with GATK. For a whole genome,
we typically use Variant Quality Score Recalibration (VQSR), but there are not enough variants in the tiny
region we called to build a model for that, so we will use the best practices filtering. It's a long tedious command again:

```bash
gatk VariantFiltration -R /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-O _FILTEREDVCF_ -V _INPUTVCF_ \
--filter-expression 'QD < 2.0' --filter-name QDfilter \
--filter-expression 'MQ < 40.0' --filter-name MQfilter \
--filter-expression 'FS > 60.0' --filter-name FSfilter \
--filter-expression 'SOR > 3.0' --filter-name SORfilter \
--filter-expression 'MQRankSum < -12.5' --filter-name MQRSfilter \
--filter-expression 'ReadPosRankSum < -8.0' --filter-name RPRSfilter
```

You can type the \ characters to spread this across several lines so it is easier to read.

Now we can look at the filtered VCF. If you did not get a filtered vcf, you can use one from /data/alignment/vcfs/NA12878.chr20.61.final.vcf

There is one last thing we will do if we have time, which is make a gvcf.
When we do large projects, we typically do not call the samples one at a time. We preprocess each sample with 
haplotype caller than jointly call them using GenotypeGVCFs. We will not do that since we are only looking at one
sample, but we will make the GVCF.

```bash
gatk HaplotypeCaller -O _GVCF_ -L chr20:61723571-62723570 -I _BAM_ -R /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa -ERC GVCF
```

We can take a look at the GVCF and then end here.

Feel free to go back and run any of these steps again with different programs or different inputs.
