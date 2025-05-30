---
title: "Cancer Genomics"
teaching: 60
exercises: 0
---

## Slides

You can view Nico's talk [here](https://docs.google.com/presentation/d/1nJU9T-r8qJBmPblQ2TJouOvBDUCVRVBBxEMTs2RLzk0/edit?usp=sharing)

## Generate somatic variant calls with Mutect2

Mutect2 is a somatic mutation caller in the Genome Analysis Toolkit ([GATK](https://gatk.broadinstitute.org/hc/en-us)), designed for detecting single-nucleotide variants (SNVs) and INDELs (insertions and deletions) in cancer genomes. The input for Mutect2 includes preprocessed alignment files (CRAM files). The CRAMs include sequencing reads that have been aligned to a reference genome, sorted, duplicate-marked, and base quality score recalibrated. Mutect2 can call variants from paired tumor and normal sample CRAM files, and also has a tumor-only mode. 
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

> ## Filtering Mutect2 calls
>
> ~~~
> gatk FilterMutectCalls \
>   -V --- (raw variants output vcf) \
>   -R --- (reference fasta) \
>   -O --- (filtered variants  vcf)
> ~~~
> {: .source}
{: .challenge}

## Anatomy of a VCF file & bcftools
[VCF file spec](https://samtools.github.io/hts-specs/VCFv4.5.pdf)

[MAF file spec](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format)

[bcftools](https://samtools.github.io/bcftools/bcftools.html)

# What is a VCF?
Variant Call Format (VCF) files are a widely used file format for representing genetic variation, with an official specification maintained by the Global Alliance for Genomics & Health ([GA4GH](https://www.ga4gh.org)). This organization also maintains the specifications for the SAM, BAM, and CRAM file formats. The VCF specification has evolved over time, and can now represent SNVs, INDELs, structural variants, and copy number variants, as well as any annotations. A large ecosystem of software exists for parsing and manipulating VCFs, the most useful of which is BCFtools. 

BCFtools is maintained by the same folks as SAMtools, and shares a lot of underlying code (namely the HTSlib library). For most routine tasks, such as sorting, filtering, annotating, and summarizing, BCFtools provides utilities so you don’t have to worry about your own implementation. For more complicated queries and analysis, you may need to write your own code, or use a command line tool with more sophisticated expressions like [vcfexpress](https://github.com/brentp/vcfexpress).

Another commonly used file format for representing variants is the Mutation Annotation Format (MAF). This format’s origins lie in The Cancer Genome Atlas (TCGA) project. It’s only meant to represent SNVs and INDELs, and the specification is quite rigid with respect to the allowed annotations. However, the data is represented in a convenient tabular format and the maftools R package exists for easy plotting, analysis, and comparison to TCGA. 

## File structures
Like SAM/BAM files, VCFs are split into two parts: the header, followed by the variant calls.


**Exercise**: Print the help menu for BCFtools
~~~
bcftools -h
~~~
{: .source}

**Exercise**: Print just the header with BCFtools, print just the variants with bcftools view

~~~
1. bcftools view -h ${vcf} OR bcftools head ${vcf}
2. bcftools view -H ${vcf}
~~~

<img width="655" alt="Screenshot 2025-05-30 at 3 37 08 PM" src="https://github.com/user-attachments/assets/31c410c4-84d6-45e1-adbf-808fb5f2584a" />


VCFs can contain information about multiple samples. In a cancer context, this is usually the paired tumor and normal, but you can imagine a situation in which we have multiple tumors from the same patient (e.g. a primary and metastasis). Each variant call can be thought of as having two parts, variant-level information (e.g., position, reference allele, impact on gene coding sequence), and more fine-grained sample-level information for a given variant (e.g. VAF, sequencing depth). 

<img width="778" alt="Screenshot 2025-05-30 at 3 37 13 PM" src="https://github.com/user-attachments/assets/3868438d-622c-4758-b3c2-8585ca772259" />


Typically, we want to filter somatic variants to reduce false positives. The FILTER column indicates the relevant filters for a given variant call. Before filtering, this field is blank or contains just a period (“.”). In the earlier filtering step, Mutect2 “soft-filtered” the calls. That is, it filled in the FILTER column for us, but did not yet remove the calls. Sometimes it’s helpful to inspect these flagged calls to debug some downstream issue. To “hard filter”, we want to only keep calls with a PASS in the FILTER column. 

> **Exercise**: Filter out the Mutect2 calls with `bcftools filter`
>> ## Solution
>>
>> ```
>> bcftools filter -i "FILTER='PASS'" -o ${out_vcf} ${in_vcf}
>> ```
> {: .solution}
{: .challenge}

**Note**: bcftools filter can also utilize other fields for filtering (e.g. VAF). Read more about filtering expressions [here](https://samtools.github.io/bcftools/howtos/filtering.html).

<img width="782" alt="Screenshot 2025-05-30 at 4 02 09 PM" src="https://github.com/user-attachments/assets/c52bc50d-0e6e-43fe-91f7-9a8ac2f81138" />


The INFO field contains variant level annotations such as gene impact and population frequency. These are typically added by tools such as [Ensembl VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html), or `bcftools annotate`.

Each annotation should have a corresponding entry in the header indicating how the annotation should be parsed for bcftools and other tools (e.g. python’s Pysam or R’s VariantAnnotation), and a plaintext description for you, the user. 

<img width="723" alt="Screenshot 2025-05-30 at 4 02 39 PM" src="https://github.com/user-attachments/assets/aa5fc2fc-78c8-46da-8547-b19fa87775bc" />


The FORMAT column describes the order of sample-level annotations. Each sample is given a column after the FORMAT column, with the order described in FORMAT, and each annotation detailed in the header. 

<img width="741" alt="Screenshot 2025-05-30 at 4 02 58 PM" src="https://github.com/user-attachments/assets/3e7350f8-ff8a-4849-bd67-472462229ca1" />


> **Exercise**: What samples are contained in the Mutect2 VCF?
>
>> ## Hint
>>
>> Try using `bcftools head`
>> 
>> 
> {: .solution}
{: .challenge}

Finally, the reference genome contigs used in variant calling are also listed in the header. 

## Bgzip and tabix indexing

[Tabix: fast retrieval of sequence features from generic TAB-delimited files](https://pmc.ncbi.nlm.nih.gov/articles/PMC3042176/) 

Oftentimes, we want to access variants based on their position in the genome. Doing this in a reasonable amount of time requires an index. BCFtools comes with the utilities **bgzip** and **tabix** for accomplishing this. **Bgzip** implements a modified version of the gzip compression algorithm, compressing the data in blocks of a predetermined size. These bgzipped files are perfectly standards-compliant gzip files and will still work with zcat/gunzip/etc. Tabix is a general purpose algorithm for indexing **coordinate-sorted** genomic coordinate data in tabular format that relies on bgzip compression to work.  

> **Exercise**: Sort the Mutect2 VCF, bgzip, and tabix-index VCF with bcftools sort, bgzip, and tabix
>
>> ## Solution
>>
>> ```
>> bcftools sort -o ${out_vcf} ${in_vcf}
>> bgzip ${out_vcf} 
>> tabix ${out_vcf}.gz
>> ```
>> Shorter answer (done in a single short bcftools command):
>> ```
>> bcftools sort -O z -Wtbi -o ${out_vcf_gz} ${in_vcf}
>> ```
> {: .solution}
{: .challenge}

> **Exercise**: Print the variants overlapping the interval “chr1:1-1000000” in the uncompressed VCF with **bcftools view**. Try the same thing with the bgzipped+indexed VCF
>
>> ## Solution
>>
>> ```
>> bcftools view ${vcf_gz} 'chr1:1-1000000'
>> ```
> {: .solution}
{: .challenge}

### Other BCFtools examples
* **isec**: perform set operations on multiple VCFs (e.g. intersection, union, set differences)
* **merge**: merge multiple VCFs
* **norm**: normalize indel representation, specify how multi-allelic sites should be handled 
   * **NOTE**: Normalizing INDELs is required before annotation, and when comparing output from multiple variant callers
* **reheader**: Rename samples and/or change contigs in the header 
* **query**: Extract information from a VCF in a user-defined format

## Extracting variants and plotting the VAF distribution

> **Exercise**: Extract chromosome, position, and tumor VAF from the VCF. In order to extract the correct VAF values, you’ll need to provide BCFtools with the tumor sample name listed in the header. 
> [Extracting information from VCFs](https://samtools.github.io/bcftools/howtos/query.html)
>> ## Solution
>> 
>> ```
>> bcftools query -s ${tumor_sample_name} -f "%CHROM\t%POS\t[%AF]\n" inputFiles/COLO-829_2B--COLO-829BL_1B.snv.indel.high_confidence.v7.annotated.vcf > COLO-829_2B--COLO-829BL_1B.vafs.tsv
>> ```
> {: .solution}
{: .challenge}

> **Exercise**:  Plot VAF plot for chr21 and chr22. Identify something different on chr22. (for loop)
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
> {: .solution}
{: .challenge}

## Inspecting alignments with IGV

**Good candidates**:
1. chr22:22949473 - a straightforward clonal SNV
2. chr1:3243220 - a straightforward deletion
3. chr2:48028500 - a straightforward insertion

**Bad candidates**:
1. chr22:24,098,035 and chr22:24,098,086 - supporting reads have multiple mismatches, soft-clip on the left
2. chr22:29289914-29290094 - Five clustered variants all supported by the same set of reads, with soft-clip on the left. "Group by base at ..." is our friend here
3. chr22:19132585 - 7bp insertion with flanking mismatches and soft-clipping

> **Bonus exercise**: What organism is the contamination coming from?
> 
>> ## Solution
>> Right-click one of the artifact-supporting reads > Copy read sequence > paste into BLAST
>> <img width="586" alt="Screenshot 2025-05-30 at 4 04 32 PM" src="https://github.com/user-attachments/assets/05caabf1-db4e-40c4-ba0b-9b63592c4d9d" />
> {: .solution}
{: .challenge}

## Mutational Signatures (working in Python)

Resources:
* [**COSMIC SBS signature database**](https://cancer.sanger.ac.uk/signatures/sbs)
* [**SigProfilerAssignment documentation**](https://osf.io/mz79v/wiki/home/)
* [Paper](https://www.nature.com/articles/s41586-020-1943-3)
* Pre-run test inputs/outputs: `mutational_signatures/pre_run_files`

### Step-by-step version:
1. Create input vcf folder
> ~~~
> tar -xvzf input_vcfs.tar.gz
> ~~~

OR

> ~~~
> mkdir input_vcfs
>
> cp COLO-829_2B--COLO-829BL_1B.snv.indel.high_confidence.v7.annotated.vcf input_vcfs/
> ~~~

2. Install genome

> ~~~
> from SigProfilerMatrixGenerator import install as genInstall
> 
> genInstall.install('GRCh38', rsync=False, bash=True)
> ~~~

3. Generate mutational spectrum count matrix

> ~~~
> from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
> 
> matrices = matGen.SigProfilerMatrixGeneratorFunc("SIW_2025", "GRCh38", "input_vcfs", plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
> ~~~

Compare the SBS96 count matrix to the COSMIC database here. Looks closest to SBS7a (UV exposure) (input_vcfs/output/plots/SBS_96_plots_SIW_2025.pdf)
<img width="736" alt="Screenshot 2025-05-30 at 4 06 22 PM" src="https://github.com/user-attachments/assets/d1bf02c2-5877-434c-9f99-5de24fd92444" />



4. Deconvolve mutational signatures
> ~~~
> from SigProfilerAssignment import Analyzer as Analyze
>
> Analyze.cosmic_fit(samples="input_vcfs/output/SBS/SIW_2025.SBS96.all", output="output", input_type="matrix", context_type="96", genome_build="GRCh38")
> ~~~

5. Look at plots of mutational signatures
(output/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf)
>> ## Solution
>> SBS5: Unknown clock-like signature
>> SBS7a: UV-light exposure
>> SBS7b: UV-light expsure
>> SBS38: Unknown. Found only in ultraviolet light associated melanomas suggesting potential indirect damage from UV-light.
>> <img width="207" alt="Screenshot 2025-05-30 at 4 08 40 PM" src="https://github.com/user-attachments/assets/c68423c2-6a31-436d-9a28-231dc4a6182f" />
> {: .solution}
{: .challenge}

## Cohort-level analysis

Most cancer genomics studies involve the study of a whole cohort of samples, not just one or two. Furthermore, we often want to compare to previously-published data to uncover new biology, and as a way to sanity-check our data – e.g, did something go wrong during the many steps to get from biopsy to variant calls?

To dip our toe into cohort-level analysis, we’ll be using the R maftools package and curated TCGA data from its companion package, TCGAmutations. To start, let’s spin up an interactive R session. 

[Maftools user guide](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html)

[TCGAmutations user guide](https://github.com/PoisonAlien/TCGAmutations)

Start an interactive R session

Next, let's load up the R packages we'll need. Check what TCGA datasets are available and load up the breast cancer cohort.

>> ## Solution
>> ~~~
>> library(maftools) 
>> library(TCGAmutations)
>>
>> tcga_available()
>>
>> brca = tcga_load(study = "BRCA")
>> ~~~
> {: .solution}
{: .challenge}

In an interactive R session, just typing the variable will print its contents. Let’s check what’s inside our **brca** variable.

Next, let’s plot a basic summary of the cohort, including the variant type, coding consequences, reference and alternate alleles, tumor mutation burden, and the top 10 most mutated genes.

>> ## Solution
>> ~~~
>> pdf('maf_summary.pdf'
plotmafSummary(maf=brca, rmOutlier=TRUE, addStat='median', dashboard=TRUE, titvRaw=FALSE)
>> dev.off()
>> ~~~
>> <img width="474" alt="Screenshot 2025-05-30 at 4 10 02 PM" src="https://github.com/user-attachments/assets/cde6e498-7ba9-4584-9d80-0aaaf64bbf67" />
> {: .solution}
{: .challenge}



We can expand out the top 10 most mutated genes into an oncoprint (also called an oncoplot), showing distribution of mutations in each sample. Oncoprints can contain information from SNVs and INDELs, but also copy number variants. 

> **Exercise**: Generate an oncoprint of the top 10 most mutated genes in the breast cancer cohort
>> ## Solution
>> ~~~
>> pdf('oncoprint.pdf')
>> oncoplot(maf=brca, top=10) 
>> dev.off()
>> ~~~
>> <img width="486" alt="Screenshot 2025-05-30 at 4 10 29 PM" src="https://github.com/user-attachments/assets/01991cec-4a2b-450c-9aa2-df6a448cb04c" />
> {: .solution}
{: .challenge}



It can also be useful to zoom in on specific genes, to understand how the mutations are distributed within the genes – are they clustered in a few “hotspots”?
> **Exercise**: Generate lollipop plots of top 3 most frequently mutated genes in the oncoprint
>> ## Solution
>> ~~~
>> pdf('lollipop.pdf')
>> lollipopPlot(maf=brca, gene=’PIK3CA’, AACol=HGVSp_Short, showMutationRate=TRUE)
>> lollipopPlot(maf=brca, gene=’TTN’, AACol=’HGVSp_Short’, showMutationRate=TRUE)
>> lollipopPlot(maf=brca, gene=’TP53’, AACol=’HGVSp_Short’, showMutationRate=TRUE)
>> dev.off()
>> ~~~
>> <img width="395" alt="Screenshot 2025-05-30 at 4 11 21 PM" src="https://github.com/user-attachments/assets/7365f5af-e8f5-497d-87a9-9eec8b226d3a" />
>> <img width="416" alt="Screenshot 2025-05-30 at 4 11 09 PM" src="https://github.com/user-attachments/assets/f2d67cc5-9117-4609-b11b-ffb52908cdeb" />
>> <img width="459" alt="Screenshot 2025-05-30 at 4 11 00 PM" src="https://github.com/user-attachments/assets/2212fcd8-1c17-43be-b81c-e1e0227a4c37" />
> {: .solution}
{: .challenge}




Notice anything different when comparing the distribution of mutations in TTN versus PIK3CA and TP53? TTN is a very long gene and accumulates mutations by chance. Therefore, mutations in this gene are uniformly distributed, rather than clustered in regions of functional importance. The distribution of mutations in a gene’s coding sequence is a signal utilized in some driver gene discovery tools such as [OncodriveCLUST](https://www.google.com/url?q=https://academic.oup.com/bioinformatics/article/29/18/2238/240376&sa=D&source=docs&ust=1748636888210148&usg=AOvVaw2fjg456YNeSQgXfruUIEKs).

As mentioned earlier, it can be helpful to comapre a cohort against previously-published data. Let's compare this cohort's tumor mutation burden (TMB) to the rest of TCGA.
>> ## Solution
>> ~~~
>> pdf('tcga_comparison.pdf')
>> tcgaCompare(maf=brca, cohortName=‘Workshop', logscale=TRUE, capture_size=50)
>> dev.off()
>> ~~~
> {: .solution}
{: .challenge}

<img width="398" alt="Screenshot 2025-05-30 at 4 11 51 PM" src="https://github.com/user-attachments/assets/6ab3eba7-7855-44fd-b46e-b967708a9a2f" />


The oncoprint hinted at patterns of mutual exclusivity. We can check this in a more statistically rigorous manner with the `somaticInteractions` function.

> **Exercise**: Check patterns of co-occurence and mutual exclusivity
>> ## Solution
>> ~~~
>> pdf('interactions.pdf')
>> somaticInteractions(maf=brca, top=20, pvalue=c(0.05, 0.1))
>> dev.off()
>> ~~~
> {: .solution}
{: .challenge}

<img width="379" alt="Screenshot 2025-05-30 at 4 12 08 PM" src="https://github.com/user-attachments/assets/fab25087-9df0-422e-b233-08dbfdfa312a" />






