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

## Differential expression testing

Now that we have created the DESeq object, and done some initial data exploration, it is time to run differential expression testing to see which genes are significantly different between conditions (here, drug-treated vs. not drug-treated).

One step we need to run to get these results is the command `DESeq`. Let's pull up the help message for this.

```
?DESeq
```
{: .language-r}

Start of this help message:

```
Usage:

     DESeq(
       object,
       test = c("Wald", "LRT"),
       fitType = c("parametric", "local", "mean", "glmGamPoi"),
       sfType = c("ratio", "poscounts", "iterate"),
       betaPrior,
       full = design(object),
       reduced,
       quiet = FALSE,
       minReplicatesForReplace = 7,
       modelMatrixType,
       useT = FALSE,
       minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
       parallel = FALSE,
       BPPARAM = bpparam()
     )
     
Arguments:

  object: a DESeqDataSet object, see the constructor functions
          ‘DESeqDataSet’, ‘DESeqDataSetFromMatrix’,
          ‘DESeqDataSetFromHTSeqCount’.
```
{: .output}

Seems like we can just input our DESeqDataSet object to the `object` argument?

Let's do that, and save as a new object called `dds`.

```
dds = DESeq(object = myDESeqObject)
```

If all goes well, you should get the following output as the command runs.

```
using pre-existing size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
```
{: .output}

Next, we will need to run the `results` function on `dds`.

Now, pull up the help message for that function.

```
?results
```
{: .language-r}

```
Usage:

     results(
       object,
       contrast,
       name,
       lfcThreshold = 0,
...)

Arguments:

  object: a DESeqDataSet, on which one of the following functions has
          already been called: ‘DESeq’, ‘nbinomWaldTest’, or
          ‘nbinomLRT’

contrast: this argument specifies what comparison to extract from the
          ‘object’ to build a results table. one of either:

            • a character vector with exactly three elements: the name
              of a factor in the design formula, the name of the
              numerator level for the fold change, and the name of the
              denominator level for the fold change (simplest case)

...

Examples:

     ## Example 1: two-group comparison
     
     dds <- makeExampleDESeqDataSet(m=4)
     
     dds <- DESeq(dds)
     res <- results(dds, contrast=c("condition","B","A"))
```
{: .output}

To highlight the example:

```
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","B","A"))
```
{: .output}

Replace `condition`, `B`, and `A` with the appropriate values for this data set.

We are looking to compare the effect of treatment with the drug.

Let's look at our experimental design matrix again.

```
expdesign
```

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

> ## Filling in the contrast argument to `results` from expdesign
>
> Some questions we need to answer here are:
>
> - What is the variable name for drug treatment here?
> - What are the two levels of this variable?
> - What order should we put these two levels in, if we want untreated to be the baseline?
{: .callout}

>
> > ## Solution
> >
> > - The variable name for drug treatment is `dex`.
> > - The two levels are `trt` and `untrt`.
> > - We should put `trt` first and `untrt` second so that the fold-change will be expressed as treated/untreated.
> {: .solution}
{: .challenge}

Now that we figured that out, we just need to plug everything into the function.

Output to an object called `res`.

>
> > ## Solution
> >
> > ```
> > res = results(object = dds, contrast = c("dex","trt","untrt"))
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

Let's look at the first few rows of the results using `head`.

```
head(res)
```
{: .language-r}

```
log2 fold change (MLE): dex trt vs untrt 
Wald test p-value: dex trt vs untrt 
DataFrame with 6 rows and 6 columns
                  baseMean log2FoldChange     lfcSE      stat      pvalue
                 <numeric>      <numeric> <numeric> <numeric>   <numeric>
ENSG00000000003 708.602170     -0.3812540  0.100654 -3.787752 0.000152016
ENSG00000000005   0.000000             NA        NA        NA          NA
ENSG00000000419 520.297901      0.2068126  0.112219  1.842943 0.065337292
ENSG00000000457 237.163037      0.0379204  0.143445  0.264356 0.791505742
ENSG00000000460  57.932633     -0.0881682  0.287142 -0.307054 0.758801924
ENSG00000000938   0.318098     -1.3782270  3.499873 -0.393793 0.693733530
                      padj
                 <numeric>
ENSG00000000003 0.00128292
ENSG00000000005         NA
ENSG00000000419 0.19646985
ENSG00000000457 0.91141962
ENSG00000000460 0.89500478
ENSG00000000938         NA
```
{: .output}

Run the following to get an explanation of what each column in this output means.

```
mcols(res)$description
```
{: .language-r}

```
[1] "mean of normalized counts for all samples"
[2] "log2 fold change (MLE): dex trt vs untrt" 
[3] "standard error: dex trt vs untrt"         
[4] "Wald statistic: dex trt vs untrt"         
[5] "Wald test p-value: dex trt vs untrt"      
[6] "BH adjusted p-values"
```
{: .output}

Ready to start to explore and interpret these results.

## Visualizing differential expression as a heatmap

Use the `which` function to get the indices (row numbers) of the genes with adjusted p-value < .01 (false discovery rate of 1%).

Then, save as an object called `degs`.

First, which column index or name contains the adjusted p-values?

>
> > ## Solution
> > 
> > Column "padj", or column number 6.
> {: .solution}
{: .challenge}

Then, we will need to extract this column. A few possible ways listed below.

>
> > ## Solution
> >
> > res$padj
> > res[,6]
> > res[,"padj"]
> {: .solution}
{: .challenge}

Next, help message for the `which` command (from base package) gives the following:

```
Usage:

     which(x, arr.ind = FALSE, useNames = TRUE)
     arrayInd(ind, .dim, .dimnames = NULL, useNames = FALSE)
     
Arguments:

       x: a ‘logical’ vector or array.  ‘NA’s are allowed and omitted
          (treated as if ‘FALSE’).

Examples:

     which(LETTERS == "R")
..
which(1:10 > 3, arr.ind = TRUE)
```

So, we just need to add a statement to get only values less than 0.01 from the column.

>
> > ## Solution
> >
> > ```
> > degs = which(res$padj < .01)
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

How many differentially expressed genes do we find?

```
length(degs)
```

```
[1] 2901
```
{: .output}

Next, let's subset the normalized expression matrix (`rlogMatrix`) to only these genes, and output to a new matrix rlogMatrix.degs.

Syntax of matrix subsetting:

```
mydat.subset = mydat[indices_to_subset,]
```

>
> > ## Solution
> >
> > ```
> > rlogMatrix.degs = rlogMatrix[degs,]
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

Let's run the pheatmap command on this matrix to output a heatmap to the console.

```
pheatmap(object = rlogMatrix.degs)
```
{: .language-r}

```
Error in pheatmap(object = rlogMatrix.degs) : 
  argument "mat" is missing, with no default
```
{: .output}

Well, that's why we should read the help message! Seems the main argument for the input object to pheatmap is `mat`, not `object`.

```
pheatmap(mat = rlogMatrix.degs)
```
{: .language-r}

![DE_heatmap_no_arguments]({{ page.root }}/fig/DE_heatmap_no_arguments.png)

Does not seem super informative! 

Let's look at the pheatmap help message again and see if we can add some arguments to make it better.

```
?pheatmap
```
{: .language-r}

Scrolling down to arguments, we see:

```
scale: character indicating if the values should be centered and
          scaled in either the row direction or the column direction,
          or none. Corresponding values are ‘"row"’, ‘"column"’ and
          ‘"none"’
```
{: .output}

Let's scale across the genes, so by row.

And also this:

```
show_rownames: boolean specifying if row names are be shown.
```
{: .output}

Boolean means either TRUE or FALSE. Let's turn this off.

>
> > ## Solution
> >
> > ```
> > pheatmap(mat = rlogMatrix.degs,
> >     scale="row",
> >     show_rownames=FALSE)
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

![DE_heatmap_with_arguments]({{ page.root }}/fig/DE_heatmap_with_arguments.png)

Think we mostly have what we want now. 

Except, the "SRR" sample names are not very informative.

It would be good to have an annotation bar at the top of the heatmap that said what cell line each sample is from, and whether it is treated or untreated.

Going back to the pheatmap help message:

```
annotation_row: data frame that specifies the annotations shown on left
          side of the heatmap. Each row defines the features for a
          specific row. The rows in the data and in the annotation are
          matched using corresponding row names. Note that color
          schemes takes into account if variable is continuous or
          discrete.

annotation_col: similar to annotation_row, but for columns.
```
{: .output}

We already have a data frame with the annotation we are looking for, so just need to give it as an argument here.

>
> > ## Solution
> >
> > ```
> > pheatmap(mat = rlogMatrix.degs,
> >     scale="row",
> >     show_rownames=FALSE,
> >     annotation_row=expdesign)
> > ```
> > {: .language-r}
> >
> > ```
> > Error in seq.int(rx[1L], rx[2L], length.out = nb) : 
> > 'from' must be a finite number
> > In addition: Warning messages:
> > 1: In min(x) : no non-missing arguments to min; returning Inf
> > 2: In max(x) : no non-missing arguments to max; returning -Inf
> > ```
> > {: .output}
> {: .solution}
{: .challenge}

Oops - looks like we input as an argument to annotation_row, when should have been an argument to annotation_col.

>
> > ## Solution
> >
> > ```
> > pheatmap(mat = rlogMatrix.degs,
> >     scale="row",
> >     show_rownames=FALSE,
> >     annotation_col=expdesign)
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

![DE_heatmap_with_annotation]({{ page.root }}/fig/DE_heatmap_with_annotation.png)

Looks like the samples separate primarily by treatment, as we would expect here since these are the differentially expressed genes.

We also see some variability by cell line, as we would have expected from the PCA plot.

Finally, we see a mix of genes upregulated vs. downregulated in the treated condition.
