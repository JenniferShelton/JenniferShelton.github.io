---
title: "Single-cell RNA"
teaching: 3h
exercises: 10
questions:
- How can I prepare a single-cell dataset for input into analysis?
- How can I evaluate the quality of a single-cell RNA-Seq dataset?
- How can I visualize a single-cell dataset?
- How can I annotate cell types in a single-cell dataset?
objectives:
- Move from a downloaded gene-cell count matrix to a filtered, normalized version of this dataset.
- Create visualizations of quality metrics and dimensional reduction to summarize the data.
- Annotate each cell with a cluster and a cell type label.
keypoints:
- Single-cell RNA analysis often starts with a gene-cell count matrix.
- Each dot in the scatter plots from a single-cell RNA analysis, or column in a heatmap, often represents a cell.
- Clustering is a key step in understanding biology.
- We can use known marker genes or run differential expression between clusters to annotate cell types after clustering.
---

# About this tutorial

This Single-cell RNA (scRNA) workflow is based on [this vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) from the Satija lab, the authors of the Seurat software package.

10X Genomics (which makes a popular technology for single-cell sequencing), regularly publishes datasets to showcase the quality of their kits. This is commonly done on peripheral blood mononuclear cells (PBMCs), as these are a readily available source of human samples (blood) compared to other tissues. These samples also contain a set of easily distinguishable, well-annotated cell types (immune cells).

The 10X website provides more information and downloadable files for this dataset [here](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0).

# Download and unpack the data

We will be using the filtered gene-cell matrix downloaded from the following link:

```
https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
```
{: .output}

Download the tarball with all the files for the expression matrix using the `wget` command.

```bash
wget https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
```

List the current working directory and see the download.

```bash
ls -F
```

```
pbmc3k_filtered_gene_bc_matrices.tar.gz
```
{: .output}

Unpack the archive with `tar`.

```bash
tar -xvf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

The command `tar` will print the path of the directories, subdirectories and files it unpacks on the screen.
```
filtered_gene_bc_matrices/
filtered_gene_bc_matrices/hg19/
filtered_gene_bc_matrices/hg19/matrix.mtx
filtered_gene_bc_matrices/hg19/genes.tsv
filtered_gene_bc_matrices/hg19/barcodes.tsv
```
{: .output}

Now list the current working directory.

```bash
ls -F
```

You should see a new directory (`filtered_gene_bc_matrices/`) that was not there before. List that directory (remember to try <kbd>Tab</kbd> to autocomplete after typing a few characters of the path).

```bash
ls -F filtered_gene_bc_matrices/
```

```
hg19/
```
{: .output}

List the sub-directory under that (remember to bring back the prior command with <kbd>↑</kbd>).

```bash
ls -F filtered_gene_bc_matrices/hg19/
```

```
barcodes.tsv  genes.tsv  matrix.mtx
```
{: .output}


We now have the relative path to the directory with the files we need. We will use this path to load the data into R.

# Start an R session

Start an R or RStudio session as you did in this previous lesson on R.

Then, we will run all of the following commands within that session.

# Single-cell workflow in R with Seurat

> ## Load libraries.
>
> Load all the libraries you will need for this tutorial using the `library` command. Today we will load `dplyr`, `Seurat`, `patchwork`. 
> ```
> library(dplyr)
> library(Seurat)
> library(patchwork)
> ```
> {: .language-r}
> {: .challenge}

> ## Read in counts and create a Seurat object.
> 
> ### Prepare the relative path to the directory with the files needed for the count matrix.
>
> Based on what we saw in the data download section,  what is the relative path to the directory with the files matrix.mtx, genes.tsv, and features.tsv files?
> > ## Solution
> > ```
> > filtered_gene_bc_matrices/hg19
> > ```
> {: .solution}
> 
> ### Read in based on the relative path.
> 
> Next, we will use the `Read10X` command to read in the downloaded counts, and `CreateSeuratObject` to create a Seurat object from this count matrix.
>
> Let’s look at the help message for Read10X.
> ```
> ?Read10X
> ```
> {: .language-r}
> If we scroll down to the examples section, we get the following:
> ```
> data_dir <- 'path/to/data/directory'
> expression_matrix <- Read10X(data.dir = data_dir)
> ```
> {: .output}
> Let's do something similar here, but replace 'path/to/data/directory' with the appropriate path.
> > ## Solution
> > ```
> > data_dir = 'filtered_gene_bc_matrices/hg19'
> > expression_matrix <- Read10X(data.dir = data_dir)
> > ```
> > {: .language-r}
> {: .solution}
>
> ### Creating the Seurat object
>
> From the Read10X help, the next line in the example was this:
>
> ```
> seurat_object = CreateSeuratObject(counts = expression_matrix)
> ```
> {: .language-r}
>
> If you ran the previous step as written (creating an object called expression_matrix), you should be able to run this line as-is.
>
> Then, we can proceed with the next steps in this workflow based on our Seurat object being called seurat_object.
> 
{: .challenge}

> ## Quality control metrics calculation
> 
> ### Calculating mitochondrial rates
> 
> One important metric correlated with cell viability (dead cells have higher rates) is mitochondrial rates, or the percent of reads going to mitochondrial genes.
>
> Here, we will use the PercentageFeatureSet argument to calculate these, and add to the Seurat object.
>
> Generate a help message for this command.
>
> ```
> ?PercentageFeatureSet
> ```
> {: .language-r}
>
> In the example at the bottom, they run the command like so to add a variable called `percent.mt` to the object, containing the mitochondrial rates.
>
> ```
> pbmc_small[["percent.mt"]] <- PercentageFeatureSet(object = pbmc_small, pattern = "^MT-")
> ```
> {: .output}
>
> The pattern argument here means that we sum up the percent of reads going to all genes starting with `MT-`, e.g. `MT-ND1` and `MT-CO1`.
>
> Let’s run this on our object, but replace `pbmc_small` with the name of the Seurat object you just made in the previous step.
>
> > ## Solution
> >
> > ```
> > seurat_object[["percent.mt"]] = PercentageFeatureSet(object = seurat_object, pattern="MT-")
> > ```
> > {: .language-r}
> {: .solution}
>
> ### QC metrics extraction
>
> Besides the mitochondrial rates we just calculated, the following metrics are also calculated automatically when we create the Seurat object:
>
> - [nCount_RNA](#nCount_RNA) - Number of total UMIs per cell
> - [nFeature_RNA](#nFeature_RNA)  -  Number of genes expressed per cell
>
> And then based on the code we already ran, we have:
>
> - [percent.mt] (#percent.mt) - Mitochondrial rate, aka % of reads in a cell that go to mitochondrial genes
>
> We can extract these metrics from the object by using the `$` operator.
>
> Let's extract each of these metrics into a series of new objects. The example below is for `nCount_RNA`. Replace the Seurat object name with the name of your object, and the name for `nFeature_RNA` and `percent.mt` as appropriate.
>
> ```
> nCount_RNA = your_seurat_object_name$nCount_RNA
> ```
> {: .language-r}
>
> > ## Solution
> > ```
> > nCount_RNA = seurat_object$nCount_RNA
> > nFeature_RNA = seurat_object$nFeature_RNA
> > percent.mt = seurat_object$percent.mt
> > ```
> > {: .language-r}
> {: .solution}
> {: .challenge}

## QC visualization

Plot `nCount_RNA` vs. `percent.mt`

```
plot(nCount_RNA, percent.mt, cex=0.1)
```
{: .language-r}

![percent_mt plot]({{ page.root }}/fig/sc-rna-1.png)

Plot `nCount_RNA` vs. `nFeature_RNA`

```
plot(nCount_RNA, nFeature_RNA, cex=0.1)
```
{: .language-r}


![nFeature_RNA plot]({{ page.root }}/fig/sc-rna-2.png)


We find that very few cells have mitochondrial rate more than 5%, and when they do they often have very low `nCount_RNA`.

We also find that `nCount_RNA` and `nFeature_RNA` are very correlated, even at very low and very high counts. 

Based on this, don’t think we need to apply any filters on `nCount_RNA` or `nFeature_RNA`. Let’s just remove cells with mitochondrial rate more than 5% (require `percent.mt < 5`).

How do we do this? Let’s look up the `subset` function, which is actually a method for a Seurat object (so look up `SeuratObject::subset` instead of `Seurat::subset`).

```
?SeuratObject::subset
```
{: .language-r}

Syntax from the help message:

```
subset(x = seurat_object_name, subset = insert_argument_here)
```
{: .output}

Example from the help message:

```
subset(pbmc_small, subset = MS4A1 > 4)
```
{: .output}

Replace `pbmc_small` here with the name of your Seurat object, and `MS4A1 > 4` with an expression to require `percent.mt` to be less than 5.

```
your_seurat_object_name = subset(x = your_seurat_object_name, subset = your_logical_expression)
```
{: .language-r}

## Data normalization, variable feature selection, scaling, and dimensional reduction (PCA)

Next, we need to run the following steps:

1. Normalize the data by total counts, so that cells with more coverage/RNA content can be made comparable to those with less. Also log-transform. This is done using the `NormalizeData` command in Seurat.
2. Choose a subset of genes with the most variability between cells, that we will use for downstream analysis. This is done using the `FindVariableFeatures` command in Seurat.
3. Z-score (scale) the expression within each gene, for all genes in the subset from step 2. This is done using the `ScaleData` command in Seurat.
4. Run dimensional reduction (PCA) on the matrix from step 3. This is done using the `RunPCA` command in Seurat.

> ## Step summary
> - `NormalizeData`
> - `FindVariableFeatures`
> - `ScaleData`
> - `RunPCA`
{: .testimonial}

And we will just use default arguments for all of these.

> ## Running commands on an object
> This is the syntax to run a set of commands like this on an object, and add the results to the existing object:
> 
> ```
> your_object_name = command1(object=your_object_name)
> your_object_name = command2(object=your_object_name)
> your_object_name = command3(object=your_object_name)
> your_object_name = command4(object=your_object_name)
> ```
> {: .language-r}
>
> In this lesson you will read the R package Seurat's' help menus and replace the example syntax above with the Seurat commands. Then replace the object name with the name of the object you have made.
{: .callout}



## Selecting number of principal components (PCs) to use downstream

Next, we need to do some exploration of the results of the PCA analysis (from `RunPCA`) to decide how many principal components to use for downstream analysis.

```
ElbowPlot(object = seurat_object)
DimHeatmap(object=seurat_object, dims=1:9, cells=500, balanced=TRUE)
DimHeatmap(object=seurat_object, dims=10:13, cells=500, balanced=TRUE, ncol=2)
```
{: .language-r}

Based on the elbow plot, it looks like it would make sense to take either the first 7 or the first 10 PCs.

Based on the heatmaps, if you know immune biology, you could also make an argument for including PCs 12 and 13 (these heatmaps contain markers of rare immune cell types).

We are going to proceed with using the first 10 PCs here, following what they did in the Seurat vignette.


## Run non-linear dimensional reduction (UMAP/tSNE)

For visualization purposes, we need to be able to have a two-dimensional representation of the data. Versus PCA, where here we have 10 dimensions that we are saying all capture an important percent of the variance.

This is where non-linear techniques like tSNE and UMAP come in. Here, we will use `UMAP`, with the principal component coordinates as input.

Let’s look at the help message for the command `RunUMAP`.


```
?RunUMAP
```
{: .language-r}

Scroll down to the example.

```
# Run UMAP map on first 5 PCs
pbmc_small <- RunUMAP(object = pbmc_small, dims = 1:5)

```
{: .output}

Here, replace `pbmc_small` with the name of your Seurat object name, and run on the first 10 PCs instead of the first 5.

```
# Write your command
```
{: .language-r}

There is also a command to plot the UMAP in the example.

From the example:

```
DimPlot(object = pbmc_small, reduction = 'umap')
```
{: .output}

> ## Write your own DimPlot command
> Replace the object name with our Seurat object name.
>
> > ## Solution
> > ```
> > DimPlot(object = seurat_object, reduction = 'umap')
> > ```
> > {: .language-r}
> {: .solution}
{: .challenge}

This gives the following plot:

![UMAP plot]({{ page.root }}/fig/sc-rna-3.png)

It looks like there are probably at least 3 major clusters in this data. But right now, we have not calculated those clusters yet. We will do so in the next step.

## Cluster the cells

Seurat applies a graph-based clustering approach, which means that clustering requires two steps. First, we run the `FindNeighbors` command, which means that we calculate the distance between each cell and all other cells, and obtain the "k" nearest neighbors for each cell from these distances. Then, we use this information to obtain the clusters using the `FindClusters` command, where the "resolution" parameter tunes the number of clusters to be reported (higher resolution = more clusters).

Just like the `UMAP` command, the `FindNeighbors` command also works on the principal components (PCA) results.

```
?FindNeighbors
```
{: .language-r}

Scroll down to the example, where we see the following:


```
pbmc_small <- FindNeighbors(pbmc_small, reduction = "pca", dims = 1:10)
```
{: .output}

We want to use the first 10 PCs as input like we did for RunUMAP, so this command is almost ready to use as-is.

Just replace `pbmc_small` with the name of your Seurat object.

Insert your command here


