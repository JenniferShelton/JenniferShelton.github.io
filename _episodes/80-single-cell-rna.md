---
title: "Single-cell RNA"
teaching: 3h
exercises: 10
questions:
- What is a command shell and why would I use one?
objectives:
- Explain how the shell relates to the keyboard, the screen, the operating system, and users’ programs.
- Explain when and why command-line interfaces should be used instead of graphical interfaces.
keypoints:
- Many bioinformatics tools can only process large data in the command line version not the GUI.
- The shell makes your work less boring (same set of tasks with a large number of files)"
- The shell makes your work less error-prone
- The shell makes your work more reproducible.
- Many bioinformatic tasks require large amounts of computing power
---

## About this tutorial

This Single-cell RNA (scRNA) workflow is based on material from the following link:

https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

10X Genomics (which makes a popular technology for single-cell sequencing), regularly publishes datasets to showcase the quality of their kits. This is commonly done on peripheral blood mononuclear cells (PBMCs), as these are a readily available source of human samples (blood) compared to other tissues. These samples also contain a set of easily distinguishable, well-annotated cell types (immune cells).

You can find more information on the dataset we will be using today here:

https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0

## Download and unpack the data

We will be using the filtered gene-cell matrix downloaded from the following link:

https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

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

## Start an R session

## Load libraries

Load all the libraries you will need for this tutorial using the `library` command. Today we will load `dplyr`, `Seurat`, `patchwork`. 


```R
library(dplyr)
library(Seurat)
library(patchwork)
```

## Read in counts and create a Seurat object.

Use the `Read10X` command to read in the downloaded counts, and `CreateSeuratObject` to create a Seurat object from this count matrix.

First, we need the path to the directory with the matrix.mtx, genes.tsv, and features.tsv files. What path is this?

```
filtered_gene_bc_matrices/hg19/barcodes.tsv
filtered_gene_bc_matrices/hg19/genes.tsv
filtered_gene_bc_matrices/hg19/matrix.mtx
```
{: .output}


Next, let’s look at the help message for Read10X.

```R
?Read10X
```

If we scroll down to the examples section, we get the following:

```
data_dir <- 'path/to/data/directory'
expression_matrix <- Read10X(data.dir = data_dir)

```
{: .output}

Let’s do something similar here, but replace 'path/to/data/directory' with the appropriate path.

```R
data_dir <- 'filtered_gene_bc_matrices/hg19/'
expression_matrix <- Read10X(data.dir = data_dir)
```
