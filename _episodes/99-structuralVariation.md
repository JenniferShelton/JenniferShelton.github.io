---
title: "Structural Variation Discovery and Genotyping"
teaching: 1.5 hr
exercises: 3 hr
questions:
- What is a structural variant?
- Why is structural variantion important?
objectives:
- Explain the difference between SNVs, INDELs, and SVs.
- Explain the different types of SVs.
- Explain the evidence we use to discover SVs.
- Investigate how long-read data can improve SV calling.
keypoints:
- Structural variants are more difficult to identify.
- Discovery of SVs usually requires multiple types of read evidence.
- There is often significant disagreement between SV callers.
- Long-reads offer significant advantages over short-reads for SV calling.
- Genotyping of long-read discovered SVs in short-read data improves completeness but has limitations.
---

## What are structural variants

Structural variation is most typically defined as variation 50 base pairs or greater. This 
is an admittedly arbitrary definition, but it provides us a useful cutoff between INDELs
and SVs. Structural variation encompases several classes of variants including deletions,
insertions, duplications, inversions, translocations, and copy number variations (CNVs).
CNVs are a subset of structural variations, specifically deletions and duplications, that
affect large (>10kb) segments of the genome.

![Deletion]({{ page.root }}/fig/SV.deletion.png)
![Insertion]({{ page.root }}/fig/SV.insertion.png)
![Duplication]({{ page.root }}/fig/SV.duplication.png)
![Inversion]({{ page.root }}/fig/SV.inversion.png)
![Translocation]({{ page.root }}/fig/SV.translocation.png)

## Detecting structural variants in short-read data

We use three types of read evidence to determine structural variations: split-reads, discordant 
read pairs, and read depth. Split reads are those where part of the read aligns to the reference
on one side of the breakpoint and the other part of the read aligns to the other side of the 
deletion breakpoint or to the inserted sequence. Discordant read pairs have an insert size that
falls significantly outside the normal distribution of insert sizes. Read depth is where increases
or decreases in read coverage occur versus the average read coverage of the genome.

Coverage comes in two variants, sequence coverage and physical coverage. Sequence coverage is the 
number of times a base was read while physical coverage is the number of times a base was read or 
spanned by paired reads.

![Coverage]({{ page.root }}/fig/SV.sequenceCov.png)

When there are no paired reads, sequence coverage equals the physical coverage. However, when
paired reads are introduced the two coverage metrics can vary widely. 

![Coverage2]({{ page.root }}/fig/SV.physicalCov.png)

For detection of structural variation physical coverage is very important.

