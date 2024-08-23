---
title: "Bulk RNA differential expression"
teaching: 60
exercises: 120
questions:
- How do I normalize a bulk RNA dataset?
- How can I summarize the relationship between samples?
- How can I compare expression based on biological or technical variables?
- How can I visualize the results of this comparison?
- How can I prepare differential expression results for input into downstream tools for biological insight?
- How can I interpret the results of functional enrichment tools?
objectives:
- Apply standard bulk workflow pre-processing methods such as DESeq2 normalization given a gene-sample count matrix.
- Execute dimensional reduction (PCA) and corresponding visualization.
- Compute differential expression between two groups, controlling for a secondary variable.
- Summarize differential expression (DE) results using a heatmap visualization.
- Interpret data visualizations from the above steps.
- Use base R commands to output text files from DE results for input into functional enrichment (GO/GSEA).
- Run GO/GSEA using a web-based tool, and describe the meaning of the resulting output.
keypoints:
- Will add this later.
---

# About this tutorial

This workflow is based on material from the rnaseqGene page on Bioconductor [(link)](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html).

This data is from a published dataset of an RNA-Seq experiment on airway smooth muscle (ASM) cell lines. From the abstract:

“Using RNA-Seq, a high-throughput sequencing method, we characterized transcriptomic changes in four primary human ASM cell lines that were treated with dexamethasone - a potent synthetic glucocorticoid (1 micromolar for 18 hours).”

Citation:

Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, Whitaker RM, Duan Q, Lasky-Su J, Nikolos C, Jester W, Johnson M, Panettieri R Jr, Tantisira KG, Weiss ST, Lu Q. “RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.” PLoS One. 2014 Jun 13;9(6):e99625. PMID: 24926665. GEO: GSE52778

DOI link [here.](https://doi.org/10.1371/journal.pone.0099625)

# Locating the data on the file system

The data has already been downloaded, and is in CSV format for easy reading into R.

The main files we will be working with today are a gene x sample count matrix, and a metadata file with design information for the samples.

Data is available in this directory:

```
# List path here.
```

If we list that path like so:

```
ls /path/to/dir
```
{: .language-r}

We find the following files:

```
airway_raw_counts.csv.gz
airway_sample_metadata.csv
```
{: .output}

# Start an R session

Start an R or RStudio session as you did in this previous lesson on R.

Then, we will run all of the following commands within that session.

# Bulk RNA differential expression workflow in R with DESeq2

## Load libraries.

Load all the libraries you will need for this tutorial using the `library` command. Today we will load `DESeq2`, `ggplot2`, and `pheatmap`.

```
library(DESeq2)
library(ggplot2)
library(pheatmap)
```
{: .language-r}

## Reading in files and basic data exploration.

First, before we read anything into R, what is the full path to the raw counts file?

>
> > ## Solution
> >
> > ```
> > /path/to/airway_raw_counts.csv.gz
> > ```
> {: .solution}
{: .challenge}

Use this to read into R using the `read.csv` command.

Let's read the help message for this command.

```
?read.csv
```
{: .language-r}

Main argument here is `file`. Let's try filling in the path from above.

>
> > ## Solution
> >
> > ```
> > read.csv(file=airway_raw_counts.csv.gz)
> > ```
> > {: .language-r}
> >
> > ```
> > Error: object 'airway_raw_counts.csv.gz' not found
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

Oops, let's fix that.

>
> > ## Solution
> >
> > ```
> > read.csv(file="airway_raw_counts.csv.gz")
> > ```
> >
> > ```
> >               gene_id SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
> > 1     ENSG00000000003        679        448        873        408       1138
> > 2     ENSG00000000005          0          0          0          0          0
> > 3     ENSG00000000419        467        515        621        365        587
> > 4     ENSG00000000457        260        211        263        164        245
> > 5     ENSG00000000460         60         55         40         35         78
> > 6     ENSG00000000938          0          0          2          0          1
> > 7     ENSG00000000971       3251       3679       6177       4252       6721
> > 8     ENSG00000001036       1433       1062       1733        881       1424
> > 9     ENSG00000001084        519        380        595        493        820
> > 10    ENSG00000001167        394        236        464        175        658
> > ...
> > [ reached 'max' / getOption("max.print") -- omitted 52566 rows ]
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

Still not quite right! Let's save the result of this command in an object called `raw.counts`.

>
> > ## Solution
> >
> > ```
> > raw.counts = read.csv(file="airway_raw_counts.csv.gz")
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

Use the `head` command to look at the first few lines.

```
head(raw.counts)
```
{: .language-r}

```
          gene_id SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
1 ENSG00000000003        679        448        873        408       1138
2 ENSG00000000005          0          0          0          0          0
3 ENSG00000000419        467        515        621        365        587
4 ENSG00000000457        260        211        263        164        245
5 ENSG00000000460         60         55         40         35         78
6 ENSG00000000938          0          0          2          0          1
  SRR1039517 SRR1039520 SRR1039521
1       1047        770        572
2          0          0          0
3        799        417        508
4        331        233        229
5         63         76         60
6          0          0          0
```
{: .output}

Looks like the first column is the gene ids.

It would be better if we could read in such that it would automatically make those the row names, so that all the values in the table could be numeric.

Let's look at the read.csv help message again to see if there is a way to do that.

```
?read.csv
```
{: .language-r}

If we scroll down we find the following:

```
row.names: a vector of row names.  This can be a vector giving the
          actual row names, or a single number giving the column of the
          table which contains the row names, or character string
          giving the name of the table column containing the row names.
```
{: .output}

Let's use this argument to read in the file again using read.csv, this time taking the first column as the row names.

>
> > ## Solution
> >
> > ```
> > raw.counts = read.csv(file="airway_raw_counts.csv.gz",row.names=1)
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

Look at the first few lines of raw.counts again.

```
head(raw.counts)
```
{: .language-r}

```
                SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
ENSG00000000003        679        448        873        408       1138
ENSG00000000005          0          0          0          0          0
ENSG00000000419        467        515        621        365        587
ENSG00000000457        260        211        263        164        245
ENSG00000000460         60         55         40         35         78
ENSG00000000938          0          0          2          0          1
                SRR1039517 SRR1039520 SRR1039521
ENSG00000000003       1047        770        572
ENSG00000000005          0          0          0
ENSG00000000419        799        417        508
ENSG00000000457        331        233        229
ENSG00000000460         63         76         60
ENSG00000000938          0          0          0
```
{: .output}

Looks good! Think we are ready to proceed with this object now.

Check how many rows there are.

```
nrow(raw.counts)
```

```
[1] 63677
```
{: .output}

So, we have 63,677 genes, and 8 samples (as we saw when we did head).

Next, read in the file with the study design.

File name is this, under the same path as the other file.

```
airway_sample_metadata.csv
```

We can once again use the csv command, and take the first column as row names.

Read into an object called `expdesign`.

>
> > ## Solution
> >
> > ```
> > expdesign = read.csv("airway_sample_metadata.csv",row.names=1)
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

Look at the first few rows of this object.

```
head(expdesign)
```
{: .language-r}

```
              cell   dex avgLength
SRR1039508  N61311 untrt       126
SRR1039509  N61311   trt       126
SRR1039512 N052611 untrt       126
SRR1039513 N052611   trt        87
SRR1039516 N080611 untrt       120
SRR1039517 N080611   trt       126
```
{: .output}

Actually, we only have a few samples here, so let's just print the whole thing.

```
expdesign
```
{: .language-r}

```
              cell   dex avgLength
SRR1039508  N61311 untrt       126
SRR1039509  N61311   trt       126
SRR1039512 N052611 untrt       126
SRR1039513 N052611   trt        87
SRR1039516 N080611 untrt       120
SRR1039517 N080611   trt       126
SRR1039520 N061011 untrt       101
SRR1039521 N061011   trt        98
```
{: .output}

We are mainly interested in columns `cell` (says which of the four cell lines the sample is) and `dex` (says whether or not the sample had drug treatment) here.

## Making a DESeq object using DESeq2

We are going to take the raw counts matrix, and the experimental design matrix, and make a special type of object called a `DESeqDataSet`.

We are going to use a function `DESeqDataSetFromMatrix`, from the `DESeq2` library that we already loaded, to do this.

Let's view the help message for that function.

```
?DESeqDataSetFromMatrix
```
{: .language-r}

Start of the help message looks like this:

```
DESeqDataSetFromMatrix(
       countData,
       colData,
       design,
       tidy = FALSE,
       ignoreRank = FALSE,
       ...
     )
```
{: .output}

We also find the following if we scroll down to the example:

```
countData <- matrix(1:100,ncol=4)
condition <- factor(c("A","A","B","B"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
```
{: .output}

If we explicitly stated the arguments in the example (and named the objects less confusingly), it would look like this:

```
mycounts = matrix(1:100,ncol=4)
condition <- factor(c("A","A","B","B"))
condition = data.frame(condition = condition)
#condition now looks like this:

dds = DESeqDataSetFromMatrix(countData = mycounts,
  colData = condition,
  design = ~condition)
```

Let's run this for our data set, but using `raw.counts` as the count matrix and `expdesign` as the design matrix.

For the design, let's look at the column names of expdesign again.

```
colnames(expdesign)
```
{: .language-r}

```
[1] "cell"      "dex"       "avgLength"
```
{: .output}

Here, we want to include both the cell line (`cell`) and treatment status (`dex`) in the design.

Save to an object called `myDESeqObject`.

>
> > ## Solution
> >
> > ```
> > myDESeqObject = DESeqDataSetFromMatrix(countData = raw.counts,
> >   colData = expdesign,
> >   design = cell,dex)
> > ```
> > {: .language-r}
> >
> > ```
> > Error: object 'dex' not found
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

Oops, looks like we did something wrong.

Should we be putting `cell` and `dex` in quotes maybe? That worked the last time we got this kind of error.

>
> > ## Solution
> >
> > ```
> > myDESeqObject = DESeqDataSetFromMatrix(countData = raw.counts,
> >   colData = expdesign,
> >   design = c("cell","dex"))
> > ```
> > {: .language-r}
> >
> > ```
> > Error in DESeqDataSet(se, design = design, ignoreRank) : 
> > 'design' should be a formula or a matrix
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

Hm, still not right.

I think what we are looking for is a design formula here.

If you search for how to make a design formula, you can find [this guide](https://genomicsclass.github.io/book/pages/expressing_design_formula.html), where we seem to see the following rules:

- Design formula should always start with `~`.
- If more than one variable, put a `+` in between.
- Nothing in quotes.

For example, their formula for the variables `diet` and `sex` was this:

```
~diet + sex
```

Let's see if we can finally fix our command to properly set the design to include `cell` and `dex` here.

>
> > ## Solution
> >
> > ```
> > myDESeqObject = DESeqDataSetFromMatrix(countData = raw.counts,
> >   colData = expdesign,
> >   design = ~cell + dex)
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

When we run the above, we get the following output:

```
Warning message:
In DESeqDataSet(se, design = design, ignoreRank) :
some variables in design formula are characters, converting to factors
```
{: .output}

This is just a warning message, doesn't mean we did anything wrong.

Looks like the command finally ran OK!

## Normalization on the DESeq object

We will next run the estimateSizeFactors command on this object, and save the output back to myDESeqObject.

This command prepares the object for the next step, where we normalize the data by library size (total counts).

```
myDESeqObject = estimateSizeFactors(myDESeqObject)
```

Next, we will run a normalization method called `rlogTransformation` (regularized-log transformation) on myDESeqObject, and save to a new object called `rlogObject`.

Basic syntax for how this works:

```
transformedObject = transformationCommand(object = oldObject)
```

Replace `transformedObject` here with the name we want for the new object, and `transformationCommand` with `rlogTransformation`. Here the `oldObject` is the DESeq object we already made.

>
> > ## Solution
> >
> > ```
> > rlogObject = rlogTransformation(object = myDESeqObject)
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

Let's also save just the normalized values from `rlogObject` in a separate object called `rlogMatrix`.

```
rlogMatrix = assay(rlogObject)
```

## Principal components analysis (PCA)

We will run DESeq2's plotPCA function on the object produced by the rlogTransformation command, to create a plot where we can see the relationship between samples by reducing into two dimensions.

Again, let's view the help message for this function.

```
?plotPCA
```
{: .language-r}

```
Help on topic ‘plotPCA’ was found in the following packages:

  Package               Library
  BiocGenerics          /nfs/sw/miniconda3/miniconda3-3.22.0/envs/R-4.2.3/lib/R/library
  DESeq2                /nfs/sw/miniconda3/miniconda3-3.22.0/envs/R-4.2.3/lib/R/library
```
{: .output}

Select that we want the message for the DESeq2 version of this command.

We find the following documentation.

```
Usage:

     ## S4 method for signature 'DESeqTransform'
     plotPCA(object, intgroup = "condition", ntop = 500, returnData = FALSE)
     
Arguments:

  object: a ‘DESeqTransform’ object, with data in ‘assay(x)’, produced
          for example by either ‘rlog’ or
          ‘varianceStabilizingTransformation’.

intgroup: interesting groups: a character vector of names in
          ‘colData(x)’ to use for grouping
```
{: .output}

Let's start by plotting based on `cell` variable.

>
> > ## Solution
> >
> > ```
> > plotPCA(object = rlogMatrix, intgroup = "cell")
> > ```
> > {: .language-r}
> >
> > ```
> > Error in (function (classes, fdef, mtable)  : 
> > unable to find an inherited method for function ‘plotPCA’ for signature ‘"matrix"’
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

Looks like we input the wrong object - should use the object that comes directly out of the rlogTransformation command.

Which object is that?

>
> > ## Solution
> >
> > ```
> > rlogObject
> > ```
> {: .solution}
{: .challenge}

Let's try the plotPCA command again using this object.

>
> > ## Solution
> >
> > ```
> > plotPCA(object = rlogObject, intgroup = "cell")
> > ```
> > {: .language-r}
> >
> > ![bulkRNA_plotPCA_intgroup_cell]({{ page.root }}/fig/bulkRNA_plotPCA_intgroup_cell.png)
> > 
> {: .solution}
{: .challenge}

Next, repeat but this time color by `dex`.

>
> > ## Solution
> >
> > ```
> > plotPCA(object = rlogObject, intgroup = "dex")
> > ```
> > {: .language-r}
> >
> > ![bulkRNA_plotPCA_intgroup_dex]({{ page.root }}/fig/bulkRNA_plotPCA_intgroup_dex.png)
> {: .solution}
{: .challenge}

> ## Interpreting the PCA plot
>
> Let's look at the two plots we made (can regenerate the first plot again if needed), and answer the following.
>
> - What principal component separates the four different cell lines? Which one separates treated vs. untreated?
> - How much variance is explained by each principal component (hint: see axis labels)?
{: .callout}

>
> > ## Solution
> >
> > PC2 separates the four different cell lines. PC1 separates treated vs. untreated. 
> >
> > PC1 explains 39% of the variance. PC2 explains 27% of the variance.
> {: .solution}
{: .challenge}
