---
title: "Cancer Genomics"
teaching: 60
exercises: 0
---

## Slides

You can view Nico's talk [here](https://docs.google.com/presentation/d/1nJU9T-r8qJBmPblQ2TJouOvBDUCVRVBBxEMTs2RLzk0/edit?usp=sharing)

## Generate somatic variant calls with Mutect2

Mutect2 is a somatic mutation caller in the Genome Analysis Toolkit ([GATK])(https://gatk.broadinstitute.org/hc/en-us), designed for detecting single-nucleotide variants (SNVs) and INDELs (insertions and deletions) in cancer genomes. The input for Mutect2 includes preprocessed alignment files (CRAM files). The CRAMs include sequencing reads that have been aligned to a reference genome, sorted, duplicate-marked, and base quality score recalibrated. Mutect2 can call variants from paired tumor and normal sample CRAM files, and also has a tumor-only mode. 
Mutect2 uses the matched normal to additionally exclude rare germline variation not captured by the germline resource and individual-specific artifacts.

> ## Executing Mutect2
>
> ~~~
> gatk Mutect2 \
>   -R --- (reference.fasta) \
>   -I COLO-829_2B.variantRegions.cram \
>   -I COLO-829_829BL_1B.variantRegions.cram \
>   -germline-resource --- (gnomad file) \
>   -O --- (raw variants output vcf)
> ~~~
> {: .source}
{: .challenge}

## Filter somatic variant calls

After calling the variants, we filter out low-quality calls using the FilterMutectCalls tool which applies various filters like minimum variant-supporting read depth and mapping quality to distinguish true somatic mutations from artifacts. The full set of filters is described in the [Mutect2 GitHub repository](https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf).

## Filtering Mutect2 calls

~~~
gatk FilterMutectCalls \
   -V --- (raw variants output vcf) \
   -R --- (reference fasta) \
   -O --- (filtered variants  vcf)
~~~
{: .source}

## Anatomy of a VCF file & bcftools
[VCF file spec](https://samtools.github.io/hts-specs/VCFv4.5.pdf)
[MAF file spec](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format)
[bcftools](https://samtools.github.io/bcftools/bcftools.html)

# What is a VCF?
Variant Call Format (VCF) files are a widely used file format for representing genetic variation, with an official specification maintained by the Global Alliance for Genomics & Health ([GA4GH])(https://www.ga4gh.org). This organization also maintains the specifications for the SAM, BAM, and CRAM file formats. The VCF specification has evolved over time, and can now represent SNVs, INDELs, structural variants, and copy number variants, as well as any annotations. A large ecosystem of software exists for parsing and manipulating VCFs, the most useful of which is BCFtools. 

BCFtools is maintained by the same folks as SAMtools, and shares a lot of underlying code (namely the HTSlib library). For most routine tasks, such as sorting, filtering, annotating, and summarizing, BCFtools provides utilities so you don’t have to worry about your own implementation. For more complicated queries and analysis, you may need to write your own code, or use a command line tool with more sophisticated expressions like [vcfexpress](https://github.com/brentp/vcfexpress).

Another commonly used file format for representing variants is the Mutation Annotation Format (MAF). This format’s origins lie in The Cancer Genome Atlas (TCGA) project. It’s only meant to represent SNVs and INDELs, and the specification is quite rigid with respect to the allowed annotations. However, the data is represented in a convenient tabular format and the maftools R package exists for easy plotting, analysis, and comparison to TCGA. 

## File structures
Like SAM/BAM files, VCFs are split into two parts: the header, followed by the variant calls.
Exercise: Print the help menu for BCFtools
~~~
bcftools -h
~~~
{: .source}

Exercise: Print just the header with BCFtools, print just the variants with bcftools view
~~~
1. bcftools view -h ${vcf} OR bcftools head ${vcf}
2. bcftools view -H ${vcf}
~~~
{: .source}

#insert image here: header v variants

VCFs can contain information about multiple samples. In a cancer context, this is usually the paired tumor and normal, but you can imagine a situation in which we have multiple tumors from the same patient (e.g. a primary and metastasis). Each variant call can be thought of as having two parts, variant-level information (e.g., position, reference allele, impact on gene coding sequence), and more fine-grained sample-level information for a given variant (e.g. VAF, sequencing depth). 

#insert image here: variant-level vs sample-level info

Typically, we want to filter somatic variants to reduce false positives. The FILTER column indicates the relevant filters for a given variant call. Before filtering, this field is blank or contains just a period (“.”). In the earlier filtering step, Mutect2 “soft-filtered” the calls. That is, it filled in the FILTER column for us, but did not yet remove the calls. Sometimes it’s helpful to inspect these flagged calls to debug some downstream issue. To “hard filter”, we want to only keep calls with a PASS in the FILTER column. 

Exercise: Filter out the Mutect2 calls with bcftools filter
>
>> ## Solution
>>
>> ```
>> bcftools filter -i "FILTER='PASS'" -o ${out_vcf} ${in_vcf}
>> ```
> {: .solution}
{: .challenge}

Note: bcftools filter can also utilize other fields for filtering (e.g. VAF). Read more about filtering expressions [here](https://samtools.github.io/bcftools/howtos/filtering.html).

#insert screenshot with FILTER highlighted

The INFO field contains variant level annotations such as gene impact and population frequency. These are typically added by tools such as [Ensembl VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html), or bcftools annotate.

Each annotation should have a corresponding entry in the header indicating how the annotation should be parsed for bcftools and other tools (e.g. python’s Pysam or R’s VariantAnnotation), and a plaintext description for you, the user. 

#insert screenshot with INFO highlighed

The FORMAT column describes the order of sample-level annotations. Each sample is given a column after the FORMAT column, with the order described in FORMAT, and each annotation detailed in the header. 

#insert screenshot with FORMAT highlighted

# Exercise: 
What samples are contained in the Mutect2 VCF?
>
>> ## Hint
>>
>> Try using bcftools head
>> 
>> 
> {: .solution}
{: .challenge}

Finally, the reference genome contigs used in variant calling are also listed in the header. 

## Bgzip and tabix indexing

[Tabix: fast retrieval of sequence features from generic TAB-delimited files](https://pmc.ncbi.nlm.nih.gov/articles/PMC3042176/) 

Oftentimes, we want to access variants based on their position in the genome. Doing this in a reasonable amount of time requires an index. BCFtools comes with the utilities **bgzip** and **tabix** for accomplishing this. **Bgzip** implements a modified version of the gzip compression algorithm, compressing the data in blocks of a predetermined size. These bgzipped files are perfectly standards-compliant gzip files and will still work with zcat/gunzip/etc. Tabix is a general purpose algorithm for indexing **coordinate-sorted** genomic coordinate data in tabular format that relies on bgzip compression to work.  

Exercise: Sort the Mutect2 VCF, bgzip, and tabix-index VCF with bcftools sort, bgzip, and tabix
>
>> ## Solution
>>
>> ```
>> bcftools sort -o ${out_vcf} ${in_vcf}
>> bgzip ${out_vcf} 
>> tabix ${out_vcf}.gz
>> ```
>> Shorter answer (done in a single bcftools short command):
>> ```
>> bcftools sort -O z -Wtbi -o ${out_vcf_gz} ${in_vcf}
>> ```
> {: .solution}
{: .challenge}

Exercise: Print the variants overlapping the interval “chr1:1-1000000” in the uncompressed VCF with **bcftools view**. Try the same thing with the bgzipped+indexed VCF
>
>> ## Solution
>> 
>> ```
>> bcftools view ${vcf_gz} 'chr1:1-1000000'
>> ```

## Other BCFtools examples
* **isec**: perform set operations on multiple VCFs (e.g. intersection, union, set differences)
* **merge**: merge multiple VCFs
* **norm**: normalize indel representation, specify how multi-allelic sites should be handled 
   * **NOTE**: Normalizing INDELs is required before annotation, and when comparing output from multiple variant callers
* **reheader**: Rename samples and/or change contigs in the header 
* **query**: Extract information from a VCF in a user-defined format

## Extracting variants and plotting the VAF distribution

Exercise: Extract chromosome, position, and tumor VAF from the VCF. In order to extract the correct VAF values, you’ll need to provide BCFtools with the tumor sample name listed in the header. 
> 
>> ## Solution
>> 
>> ```
>> bcftools query -s ${tumor_sample_name} -f "%CHROM\t%POS\t[%AF]\n" inputFiles/COLO-829_2B--COLO-829BL_1B.snv.indel.high_confidence.v7.annotated.vcf > COLO-829_2B--COLO-829BL_1B.vafs.tsv
>> ```
> {: .solution}
{: .challenge}
[extracting information from VCFs](https://samtools.github.io/bcftools/howtos/query.html)

Exercise:  Plot VAF plot for chr21 and chr22. Identify something different on chr22. (for loop)
> 
>> ## Solution
>> 
>> ```
>> pdf(vaf_histogram.pdf')
>> f = 'COLO-829_2B--COLO-829BL_1B.vafs.tsv'
>> x = read.table(f, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr', 'pos', 'vaf')) 
>>
>> for (chr in unique(x$chr)) {
>>    hist(x$vaf[x$chr == chr], breaks=20, main=chr)
>> }
>> dev.off()
>> ```

