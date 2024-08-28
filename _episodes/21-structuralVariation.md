---
title: "Structural Variation in short reads"
teaching: 90
exercises: 0
questions:
- What is a structural variant?
- Why is structural variantion important?
objectives:
- Explain the difference between SNVs, INDELs, and SVs.
- Explain the different types of SVs.
- Explain the evidence we use to discover SVs.
keypoints:
- Structural variants are more difficult to identify.
- Discovery of SVs usually requires multiple types of read evidence.
- There is often significant disagreement between SV callers.
---

## What are structural variants

Structural variation is most typically defined as variation affecting larger fragments of the genome
than SNVs and InDels; for our purposes those 50 base pairs or greater. This is an admittedly 
arbitrary definition, but it provides us a useful cutoff between InDels and SVs. 

> ## Importance of SVs
>
> SVs affect an order of magnitude more bases in the human genome in comparison to SNVs (Pang et al, 2010) 
> and are more likely to associate with disease.
{: .keypoints}

Structural variation encompases several classes of variants including deletions, insertions, 
duplications, inversions, translocations, and copy number variations (CNVs). CNVs are a subset of 
structural variations, specifically deletions and duplications, that affect large (>10kb) segments 
of the genome.

> ## Breakpoints
>
> The term breakpoint is used to denote a boundry between a structural variation and the reference.
{: .keypoints}

> ## Examples
> ##### Deletion
> ![Deletion]({{ page.root }}/fig/SV.deletion.png)
> 
> ##### Insertion
> ![Insertion]({{ page.root }}/fig/SV.insertion.png)
> 
> ##### Duplication
> ![Duplication]({{ page.root }}/fig/SV.duplication.png)
> 
> ##### Inversion
> ![Inversion]({{ page.root }}/fig/SV.inversion.png)
> 
> ##### Translocation
> ![Translocation]({{ page.root }}/fig/SV.translocation.png)
{: .solution}

## Simple Indels
![Simple Alignment]({{ page.root }}/fig/SV.simpleAln.png)
![Simple SVs]({{ page.root }}/fig/SV.simpleSV.png)


## Detecting structural variants in short-read data

Because structural variants are most often larger than the individual reads we must use different 
types of read evidence than those used for SNVs and InDels which can be called by simple read alignment.
We use three types of read evidence to discover structural variations: discordant read pairs, 
split-reads, and read depth. 

Discordant read pairs have insert sizes that fall significantly outside the normal distribution of 
insert sizes.

##### Insert size distribution
![Insert size distribution]({{ page.root }}/fig/SV.insertSize.png)

Split reads are those where part of the read aligns to the reference on one side of the breakpoint 
and the other part of the read aligns to the other side of the deletion breakpoint or to the 
inserted sequence. Read depth is where increases or decreases in read coverage occur versus the 
average read coverage of the genome.

##### Reads aligned to sample genome
![Reads aligned to sample]({{ page.root }}/fig/SV.readsVsample.png)

##### Reads aligned to reference genome
![Reads aligned to reference]({{ page.root }}/fig/SV.readsVref.png)

Coverage comes in two variants, sequence coverage and physical coverage. Sequence coverage is the 
number of times a base was read while physical coverage is the number of times a base was read or 
spanned by paired reads.

##### Sequence coverage
![Sequence coverage]({{ page.root }}/fig/SV.sequenceCov.png)

When there are no paired reads, sequence coverage equals the physical coverage. However, when
paired reads are introduced the two coverage metrics can vary widely. 

##### Physcial coverage
![Sequence coverage vs physical coverage]({{ page.root }}/fig/SV.physicalCov.png)

##### Read depth
![Read depth]({{ page.root }}/fig/SV.readDepth.png)

> ## Read signatures
> 
> ##### Deletion read signature
> ![Deletion read signature]({{ page.root }}/fig/SV.deletionSig.png)
> 
> ##### Inversion read signature
> ![Inversion read signature]({{ page.root }}/fig/SV.inversionSig.png)
> 
> ##### Tandem duplication read signature
> ![Tandem duplication read signature]({{ page.root }}/fig/SV.tandemDupSig.png)
> 
> ##### Translocation read signature
> ![Translocation read signature]({{ page.root }}/fig/SV.translocationSig.png)
{: .solution}

> ## Challenge
>
> What do you think the read signature of an insertion might look like?
>
> > ## Solution
> > Insert image
> >
> {: .solution}
{: .challenge}

## Copy number analysis

Calling of copy number variation from WGS data is done using read depth, where reads are counted
in bins or windows across the entire genome. The counts need to have some normalization applied to
them in order to account for sequencing irregularities such as mappability and GC content. These
normalized counts can then be converted into their copy number equivalents using a process called
segmentation. Read coverage is, however, inheirently noisy. It changes based on genomic regions, 
DNA quality, and other factors. This makes calling CNVs difficult and is why many CNV callers focus 
on large variants where it is easier to normalize away smaller confounding changes in read depth.

![CNV analysis]({{ page.root }}/fig/SV.cnvAnalysis.png)


## Caller concordance

Because SV callers 

![SV caller comparison]({{ page.root }}/fig/SV.algComparison.png)
