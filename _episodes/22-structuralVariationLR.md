---
title: "Structural Variation in long reads"
teaching: 60
exercises: 0
questions:
- What are the advantages/disadvantages of long reads?
- How might we leverage a combination of long and short read data?
objectives:
- Investigate how long-read data can improve SV calling.
keypoints:
- Long-reads offer significant advantages over short-reads for SV calling.
- Genotyping of long-read discovered SVs in short-read data allows for some scalability.
---

## Long read platforms

The two major platforms in long read sequencing are PacBio and Oxford Nanopore. 

PacBio's flagship is the Revio, which produces reads in the 5kb to 35kb range with very high accuracy.


![PacBio]({{ page.root }}/fig/SV.pbSeq.png)
![PacBio read length and quality]({{ page.root }}/fig/SV.pbRL.png)

Oxford Nanopore produces sequencers that range in size from the MinION, which is roughly smart
phone sized to the PromethION, the high throughput version that we have at NYGC. There are some
differences in the read outputs of the various platforms but the MinION has been shown to produce 
N50 read lengths over 100kb with maximum read lengths greater than 800kb using ONT's ultra-long
sequencing prep. The PromethION can produce even greater N50 values and can produce megabase long
reads. Typically these reads are lowe overall base quality than PacBio but ONT has steadily been
improving the base quality for their data.

![ONT MinION]({{ page.root }}/fig/SV.ontMin.png)
![ONT PromethION]({{ page.root }}/fig/SV.ontProm.png)
![ONT PromethION Read Length]({{ page.root }}/fig/SV.ontRL.png)
![ONT PromethION Quality]({{ page.root }}/fig/SV.ontQual.png)


## Advantages of long reads

The advantage of long reads is they map much more uniquely to the genome and can often span
repetitive elements in the genome that cause mapping quality issues with short reads. In long 
reads we are able to detect much larger events and in cases where the event is entirely inside
a read we are able to determine the breakpoints with much higher accuracy.

## SV calling in long reads

### Alignment

This works similarly to calling InDels in short reads where our primary information is the raw
gapped alignment to the reference. 

![Sniffles]({{ page.root }}/fig/SV.sniffles.png)

### Assembly

## Drawbacks

* We are still limited to some extent by the length of our reads and our ablility to span an entire 
  event with one or more reads and some regions of the genome are still very difficult to sequence 
  and align to.
* Long read sequencing is becoming more affordable but is still much more expensive than short read
  sequencing.
* Throughput is lower, reducing the turn around time for projects with large numbers of samples.
* Sequencing prep, especially for ultra-long protocols is tedious and difficult to perform 
  consitently.

## Genotyping LR SVs in SR data


> ## Question
> Does anyone know of any other technologies being used for structural variation?
{: .objectives}
