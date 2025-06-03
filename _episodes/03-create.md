---
title: "Working With Files and Directories"
teaching: 31h
exercises: 20
questions:
- How can I create, copy, and delete files and directories?
- How can I edit files?
objectives:
- Delete, copy and move specified files and/or directories.
- Create files in that hierarchy using an editor or by copying and renaming existing files.
---

## Creating directories

We now know how to explore files and directories, but how do we create them in the first place?

In this episode we will learn about creating and moving files and directories, using our home directory.

First we see where we are and what we already have. We want to be in the home directory on the Desktop, which we can check using:

```bash
pwd
# if needed we can change to the home directory with
cd
```

Let’s create a new directory called workshop using the command `mkdir` (which has no output)

**mkdir** : make directory

```bash
mkdir workshop
ls -F ~/
mkdir workshop/output/
mkdir workshop/logs/
mkdir workshop/scripts/
mkdir workshop/mock_data/
mkdir workshop/mock_data/fasta/
```

**touch** : create new, empty files

```bash
touch genome.fa
```

**vi** : a plain text editor that is already installed on many linux computers. This is not the most 
user-friendly plain text editor. A more user-friendly option is [BBEdit](https://www.barebones.com/products/bbedit/). You may independtly install this on your laptop if you need an editor for future programming tasks. It is free and allows you edit
custom scripts that are stored on a remote cluster.

```bash
vi genome.fa
```

In `vi` type <kbd>i</kbd>. Then type out a few [FASTA](https://en.wikipedia.org/wiki/FASTA_format) records. Like the ones below. 
When you are finished editing type <kbd>esc</kbd>, then type `:wq` follwed by <kbd>return</kbd> to save and quit.

```
>chr1
AACCGCGCGTACGCGCGCGCGC
CTCTACNNNNNN
>chr2
TGTGTCTGAAANTGTCATGTCA
NNNTCNTGTGTGTGAAAAAAAA
NNTGCNA
>chr3
AACCGCGCGTATCTGAACGCGC
CTCTAAATCGAT
```

**cp** : copy file(s) to a new location

```bash
# makes a copy with the same filename
cp ~/genomes.fa ~/workshop/mock_data/fasta/
# makes a copy with a new filename (both .fa and .fasta and .fna are valid FASTA file suffixes)
cp ~/genomes.fa ~/genomes.fasta
# see new file
ls -thl ~/genomes.fa ~/workshop/mock_data/fasta/
```

**mv** : move files or directories (if they are moved in the same location then they are renamed)

```bash
mv ~/genomes.fasta ~/workshop/mock_data/fasta/
mv ~/genomes.fa ~/new_name_genomes.fa
```

**rm** : remove file (CAUTION: no undelete here !!)

```bash
rm ~/new_name_genomes.fa
```

**rmdir** : remove empty directory

```bash
mkdir ~/workshop/mock_data/wrong
rmdir ~/workshop/mock_data/wrong
```

## File permissions

**chmod** : command to set permissions of a file or directory

```bash
ls -thl ~/workshop/mock_data/fasta/genomes.fa
```

The command `ls -thl` lists the files in the current folder and displays them in the long listing format. While this may initially look complex, we can break this down in the following left to right order:

A set of ten permission flags
 - Link count (which is irrelevant to this course)
 - The owner of the file
 - The associated group
 - The size of the file
 - The data that the file was last modified
 - The name of the file
 - The permission flags are the important thing we want to look at here. We can further break these down into the following three basic

Permission Types:

 - Read - Which refers to a user’s capability to read the contents of the file.
 - Write - Which refer to a user’s capability to write or modify a file or directory.
 - Execute - Which affects a user’s capability to execute a file or view the contents of a directory.

Each of these permission types is listed in the `_rwxrwxrwx` section of the `ls` output. The first character marked by an underscore is the special permission flag that can vary. It shows things like whether the item is a directory.

The following set of three characters (`rwx)` is for the owner permissions.
The second set of three characters (`rwx`) is for the Group permissions.
The third set of three characters (`rwx`) is for the All Users permissions.

See all options at http://www.onlineconversion.com/html_chmod_calculator.html.
The value `7` is the file/directory can be read, written to and executed.
With a `5` the file/directory can be read and executed (but not written to!).
With a `0` the file/directory can't be read, written to or executed.

Your bioinformatics sequence and reference files should be read only.

```bash
chmod -R 550 ~/workshop/mock_data/fasta
echo "oops" > ~/workshop/mock_data/fasta/genomes.fa
```

Tokens and key files work something like passwords and should be read only by you (with no permissions for your group and other users).

```bash
touch ~/workshop/logs/gdc_token_example.txt
chmod  500 ~/workshop/logs/gdc_token_example.txt
```
