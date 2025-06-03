---
title: "Navigating and Editing Files"
teaching: 4h
exercises: 30m
questions:
- How can I move around on the computer/vm?
- How can I see what files and directories I have?
- How can I specify the location of a file or directory on the computer/vm?
- How can I specify the location of a file or directory in a bucket?
objectives:
- Translate an absolute path into a relative path and vice versa.
- Construct absolute and relative paths that identify specific files and directories.
- Use options and arguments to change the behaviour of a shell command.
- Demonstrate the use of tab completion and explain its advantages.
keypoints:
- The file system is responsible for managing information on the disk.
- Information is stored in files, which are stored in directories (folders).
- Directories can also store other directories, which then form a directory tree.
- The command `pwd` prints the user’s current working directory.
- The command `ls [path]` prints a listing of a specific file or directory; `ls` on its own lists the current working directory.
- The command `cd [path]` changes the current working directory.
- Most commands take options that begin with a single `-`.
- Directory names in a path are separated with `/` on Unix.
- Slash (`/`) on its own is the root directory of the whole file system.
- An absolute path specifies a location from the root of the file system.
- A relative path specifies a location starting from any location other than the root.
- A `~` indicates your home directory
- A `-` indicates the last directory that you were in
- Dot (`.`) on its own means ‘the current directory’; `..` means ‘the directory above the current one’.
---

## The Filesystem

The part of the operating system responsible for managing files and directories is called the file system. It organizes our data into files, which hold information, and directories (also called ‘folders’), which hold files or other directories.

Several commands are frequently used to create, inspect, rename, and delete files and directories. To start exploring them, we’ll go to our open shell window.

First, let’s find out where we are by running a command called pwd (which stands for ‘print working directory’). Directories are like places — at any time while we are using the shell, we are in exactly one place called our current working directory. Commands mostly read and write files in the current working directory, i.e. ‘here’, so knowing where you are before running a command is important. pwd shows you where you are:

```bash
$ pwd
```

Here, the response may be different on different computers. Often a session begins in the users home directory.

To understand what a file system is, let’s have a look at how the file system as a whole is organized. For the sake of this example, we’ll be illustrating a portion of the file system on our workshop VM. After this illustration, you’ll be learning commands to explore your own filesystem, which will be constructed in a similar way, but not be exactly identical.

On the workshop VM, part of the filesystem looks like this:

```bash
$ tree -L 2 /data/RNA/
/data/RNA/
├── bulk
│   ├── airway_raw_counts.csv.gz
│   └── airway_sample_metadata.csv
└── single_cell
    └── README
```
Another way to diagram the filesystem is like this. The root directory is always named `/`.

![Simple File system]({{ page.root }}/fig/‎file_system_diagram.png)

## Relative Paths

`/` : root directory

**Absolute Path** : a path that starts from the root of the file system. Any path that starts with `/` is an absolute path.



**Relative Path** : a path that starts from current location or any location other than the root. 

`.` : current working directory

```bash
cd /data/alignment/references/
ls ./GRCh38_1000genomes/
```

`..` : one level up directory, also known as the parent directory of the current directory

```bash
ls -F ../
ls -F ../combined
ls -a
```

`~` : user's home directory

```bash
ls -F ~/
```

`cd` : a command to change your current working directory

`-` : previous directory. The dash is interpreted as the last directory that the user was in.

```bash
cd -
```

**Note:** if no special characters are used UNIX assumes your path begins in the current working directory

> ## Create a Relative Path
> 
> Without changing directories create a **relative** path to list the contents of the `/data/alignment/combined`
> directory (showing a trailing slash to see which are directories). For the second part of the code challenge,
> what does the command `cd` without a directory name do?
>
> ```bash
> $ cd /data/alignment/references/GRCh38_1000genomes/
> ```
> > ## Solution
> >
> > ```bash
> > ls -F ../../combined/
> > ```
> > 
> > Second part: changes the current working directory to the home directory
> {: .solution}
{: .challenge}

## Software on the File System the PATH Variable

Commands like `ls` and (on our VM) `samtools` seem to exist as special words that the user can type to call a single version of a program. However, these programs are actual files on the file system that we can call because they are in one of the many locations that the shell knows to search when a command is executed.

How can we run samtools when we don’t see any program named 
samtools in our current working directory?

## Location of Samtools
```bash
# generate a samtools help menu
samtools
# show the absolute path to samtools
which samtools
```

> VM path for samtools
>```
>/home/student/miniconda3/envs/siw/bin/samtools
>```
>{: .output}

## Location of ls
```bash
# show the absolute path to ls
which ls
```

> VM path for ls
>```
>/usr/bin/ls
>```
>{: .output}

you or a systems administrator will probably install some bioinformatics programs that researchers use commonly
In this workshop those have been installed at `/home/student/miniconda3/envs/siw/bin` using a environment manager called 
conda. Ask your systems administrators to assist with software installation and/or tips 
for installing tools.

What if we want to know the [version](https://github.com/samtools/samtools/releases/) of samtools? 

```bash
samtools --version
```

>```
>samtools 1.20
>Using htslib 1.20
>Copyright (C) 2024 Genome Research Ltd.
>
>Samtools compilation details:
>    Features:       build=configure curses=yes 
>    CC:             /opt/conda/conda-bld/samtools_1720645213030/_build_env/bin/x86_64-conda-linu
>...
>```
>{: .output}

You may want to start with the most recent version of a tool or need to use a previous tool to match prior analysis runs. It can be useful to record the absolute path to bioinformatics tools 
in commands that you run for publication or intend to have to run again in a consistant fashion. It can also be useful to include the bioinformatics tool version in the path to the tool for clarity.

If you are just glancing at the alignment header to see what genome it was aligned to (e.g. GRCh38) then you don't need to be so explicit.

```bash
$ samtools view -H /data/alignment/combined/NA12878.dedup.bam
```

>```
>...
>@RG     ID:NA12878_TTGCCTAG-ACCACTTA_HCLHLDSXX_L001     PL:illumina     PM:Unknown      LB:NA12878      DS:GRCh38       SM:NA12878       CN:NYGenome     PU:HCLHLDSXX.1.TTGCCTAG
>@RG     ID:NA12878_TTGCCTAG-ACCACTTA_HCLHLDSXX_L002     PL:illumina     PM:Unknown      LB:NA12878      DS:GRCh38       SM:NA12878       CN:NYGenome     PU:HCLHLDSXX.2.TTGCCTAG
>@RG     ID:NA12878_TTGCCTAG-ACCACTTA_HCLHLDSXX_L003     PL:illumina     PM:Unknown      LB:NA12878      DS:GRCh38       SM:NA12878       CN:NYGenome     PU:HCLHLDSXX.3.TTGCCTAG
>...
>```
>{: .output}

## How the Shell Finds Programs
The PATH environment variables defines the shell's search path. 
In the shell a variable is defined without a starting dollar sign but when the value 
of the variable is retrived you add the `$` begining of the variable name. **Tips:** also wrap the variable name in curly braces `{}` so that the shell can clearly see the last character that belongs to the variable name. There cannot be a space on either side of the `=` sign.

```bash
# define a variable
$ project_name="LUAD"
# retrieve the value of the variable
echo ${LUAD}
# use export to define the variable for the shell session and for any programs called during the session
$ export project_name="LUAD"
```
When you run a command like `ls` or `samtools`, the shell:
Splits `$PATH` into components to get a list of directories
Unix uses `:` as a separator
The shell looks for the program in each directory in left-to-right.
The shell runs the first program with that name that it finds

`which` reported that samtools was in `/home/student/miniconda3/envs/siw/bin/`. This is the second directory listed in our `$PATH`.

```bash
$ echo $PATH
```

>```
>/home/student/bin:/home/student/miniconda3/envs/siw/bin:/home/student/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/snap/bin:/software/manta-1.6.0.centos6_x86_64/bin:/home/student/paragraph-v2.4a/bin:/home/student/gatk-4.6.0.0
>```
>{: .output}



> ## CLI typing hints
> - <kbd>Tab</kbd> : autocompletes paths (use this for speed and to avoid mistakes !!)
> - <kbd>↑</kbd>/<kbd>↓</kbd> arrow : moves through previous commands
> - <kbd>Ctrl</kbd><kbd>a</kbd> : goes to the beginning of a line
> - <kbd>Ctrl</kbd><kbd>e</kbd>: goes to the end of the line
> - short flags generally `-` followed by a single letter
> - long flags generally `--` followed by a word
> - flags are often called options in manuals (both terms are correct)
> - command/program will be used interchangeably (a whole line of code is also called a command)
> - To list your past commands: Type `history` in the command line
{: .testimonial}


## Buckets (not covered in workshop)

On your computer files are often stored "locally" on that computer in a directory. On the cloud permanent storage areas are called a “bucket.” The console that we are using is running on an ephemeral virtual machine (VM). We will copy files to our vm or read them from the bucket to use them. Any file we create or modify in our vm will be deleted when we turn off the vm. If your lab is working on the cloud then users will use a bucket to save files needed for analysis after the vm is stopped.

On google cloud the program `gcloud storage` allows you to run [`ls`](https://cloud.google.com/sdk/gcloud/reference/storage/ls) and [`cp`](https://cloud.google.com/sdk/gcloud/reference/storage/cp) commands to search and transfer files between VMs and your buckets.

> ## Example of a file in a bucket
>
> List a file in a bucket:
>
>```
> gcloud storage ls gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list
>```
>
> Copy a file from a bucket to your current working directory.
>```
> gcloud storage cp gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list .
>```
>{: .discussion}

