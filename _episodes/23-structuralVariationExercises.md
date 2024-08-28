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
1. Run manta on our SR data for chromosomes 1, 6, and 20
2. Run minimap2 on chromosome 20
3. Run samtools merge for LR data
4. Run sniffles on merged LR bam


## Before we really start

> ## Super fun exercise
>
> ~~~
> cd
> wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
> chmod a+x bedtools.static
> mv bedtools.static miniconda3/envs/siw/bin/bedtools
> ~~~
> {: .source}
{: .challenge}

## First

We'll make our output directory and change our current directory to it.

~~~
mkdir /workshop/output/sv
cd /workshop/output/sv
~~~

## Manta

There are two steps to running manta on our data. First, we tell Manta what the 
inputs are by running `configManta.py`. 

Run the following code pieces one at a time.

~~~
/software/manta-1.6.0.centos6_x86_64/bin/configManta.py \
 --bam /data/alignment/combined/NA12878.dedup.bam \
 --referenceFasta /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa \
 --runDir manta_NA12878
~~~

~~~
./manta_NA12878/runWorkflow.py \
 --mode local \
 --jobs 8 \
 --memGb unlimited 
~~~

> ## Challenge
>
> Run manta on NA12878's parents (NA12891 and NA12892)
>
{: .challenge}

`ls -lh manta_*/results/variants/diploidSV.vcf.gz`
~~~
-rw-r--r-- 1 student student 137K Aug 27 22:54 manta_NA12878/results/variants/diploidSV.vcf.gz
-rw-r--r-- 1 student student 141K Aug 27 23:04 manta_NA12891/results/variants/diploidSV.vcf.gz
-rw-rw-r-- 1 student student 141K Aug 27 23:07 manta_NA12892/results/variants/diploidSV.vcf.gz
~~~
{: .output}

## Minimap2

> ## Challenge
>
> If we take the first 1000 lines of a fastq file how many reads do we get?
>
{: .challenge}

~~~
zcat /data/SV/long_read/inputs/NA12878_NRHG.chr20.fq.gz \
| head 20000 \
| minimap2 \
    -ayYL \
    --MD \
    --cs \
    -z 600,200 \
    -x map-ont \
    -t 8 \
    -R "@RG\tID:NA12878\tSM:NA12878\tPL:ONT\tPU:PromethION" \
    /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.mmi \
    /dev/stdin \
| samtools \
    sort \
    -M \
    -m 4G \
    -O bam \
> NA12878.minimap2.bam
~~~

#### Flags
* Minimap
    * -a   : output in the SAM format
    * -y   : Copy input FASTA/Q comments to output
    * -Y   : use soft clipping for supplementary alignments
    * -L   : write CIGAR with >65535 ops at the CG tag
    * --MD : output the MD tag
    * --cs : output the cs tag
    * -z   : Z-drop score and inversion Z-drop score
    * -x   : preset
    * -t   : number of threads
    * -R   : SAM read group line
* Samtools
    * -M        : Use minimiser for clustering unaligned/unplaced reads
    * -m        : Set maximum memory per thread
    * -O        : Specify output format

## Sniffles
~~~
sniffles \
 --threads 8 \
 --input /data/SV/long_read/minimap2/NA12878_NRHG.minimap2.chr1-6-20.bam \
 --vcf NA12878_NRHG.sniffles.vcf
~~~

## Regenotyping


> ## Another super fun exercise
>
> We need to make a copy of a file and change some of its contents.
>
> ~~~
> cd 
> cp /data/SV/inputs/NA12878.paragraph_manifest.txt .
> ~~~
> {: .source}
{: .challenge}

Open it in the editor on the side panel and change `/data/SV/bams/NA12878.chr1-6-20.bam` to
`/data/alignment/combined/NA12878.dedup.bam`. 

~~~ 15 min
cd /workshop/output/sv

~/paragraph-v2.4a/bin/multigrmpy.py \
 -i /data/SV/inputs/HGSVC_NA12878.chr1-6-20.vcf.gz \
 -m ./NA12878.paragraph_manifest.txt \
 -r /data/alignment/references/GRCh38_1000genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa \
 -o paragraph_NA12878 \
 --threads 8 \
 -M 400
~~~
{: .source}

## Secondary challenges

> ## Manta
>
> How many variants do we call using manta?
>
> > ## Solution
> > ~~~
> > zgrep -vc "^#" 
> > ~~~
> {: .solution}
>
> What is the breakdown by type?
>
{: .challenge}

> ## Sniffles
>
> How many variants do we call using sniffles? 
>
> > ## Solution
> > ~~~
> > zgrep -vc "^#" NA12878_NRHG.sniffles.vcf.gz 
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
