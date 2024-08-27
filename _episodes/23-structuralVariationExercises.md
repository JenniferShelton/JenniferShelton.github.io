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
> > Insert code block here
> {: .solution}
>
> > ## Discussion
> >  
{: .challenge}
