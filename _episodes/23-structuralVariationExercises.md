---
title: "SV Exercises"
teaching: 0
exercises: 120
questions:
- How do the calls from short read data compare to those from the long read data?
- How do the 
objectives:
- Align long read data with minimap2
- Call SVs in short and long read data
- Regenotype LR SVs in a short read data set
---

## Exercises
1. Run minimap2 on chromosome 20
2. Run manta on chromosomes 1, 6, and 20
3. Run samtools merge for LR data
4. Run sniffles on merged LR bam


~~~
zcat /data/SV/long_read/inputs/NA12878_NRHG.chr20.fq.gz \
| minimap2 \
    -ayYL \
    --MD \
    --cs \
    -z 600,200 \
    -x map-ont \
    -t 8 \
    -R "@RG\tID:NA12878\tSM:NA12878\tPL:ONT\tPU:PromethION" \
    /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    /dev/stdin \
| samtools \
    sort \
    --threads \
    -M \
    -l 0 \
    -m 4G \
    -O bam \
> NA12878.minimap2.bam
~~~


> ## Super fun exercise
>
> ~~~
> wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
> chmod a+x bedtools.static
> mv bedtools.static miniconda3/envs/siw/bin/bedtools
> ~~~
> {: .source}
{: .challenge}

> ## Challenge
>
> How many variants do we call using manta?
>
> > ## Solution
> > Insert code block here
> {: .solution}
>
> What is the breakdown by type?
>
> > ## Solution
> > Insert code block here
> {: .solution}
{: .challenge}

> ## Challenge
>
> How many variants do we call using sniffles? 
>
> > ## Solution
> > ~~~
> > zgrep -v "^#" NA12878_NRHG.sniffles.vcf.gz | wc -l
> > ~~~
> > {: .source}
> > ~~~
> > 4209
> > ~~~
> > {: .output}
> {: .solution}
>
> What is the breakdown by type? 
>
> > ## Solution
> > ~~~
> > zgrep -v "^#" NA12878_NRHG.sniffles.vcf.gz | cut -f 8 | cut -f 2 -d ";" | sort | uniq -c 
> > ~~~
> > {: .source}
> > ~~~
> >     11 SVTYPE=BND
> >   1727 SVTYPE=DEL
> >      1 SVTYPE=DUP
> >   2464 SVTYPE=INS
> >      6 SVTYPE=INV
> > ~~~
> > {: .output}
> {: .solution}
>
> > ## Discussion
> >  
{: .challenge}
