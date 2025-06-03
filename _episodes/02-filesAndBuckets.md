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
- A relative path specifies a location starting from the current location.
- A `~` indicates your home directory
- A `-` indicates the last directory that you were in
- Dot (`.`) on its own means ‘the current directory’; `..` means ‘the directory above the current one’.
---
## Open your terminal

To start we will open a terminal.
 
1. Go to the [link](https://pad.carpentries.org/Siw_atlanta) given to you at the workshop
2. Paste the notebook link next to your name into your browser
3. Select "Terminal" from the "JupyterLab" launcher (or blue button with a plus in the upper left corner)
4. After you have done this put up a green sticky not if you see a flashing box next to a `$`

What am I seeing: when the shell is first opened, you are presented with a **prompt**,
indicating that the shell is waiting for input.


## Background

Humans and computers commonly interact in many different ways, such as through a keyboard and mouse,
touch screen interfaces, or using speech recognition systems.
The most widely used way to interact with personal computers is called a
**graphical user interface** (GUI).
With a GUI, we give instructions by clicking a mouse and using menu-driven interactions.

While the visual aid of a GUI makes it intuitive to learn,
this way of delivering instructions to a computer scales very poorly.
Imagine the following task:
for a literature search, you have to copy the third line of one thousand text files in one thousand
different directories and paste it into a single file.
Using a GUI, you would not only be clicking at your desk for several hours,
but you could potentially also commit an error in the process of completing this repetitive task.
This is where we take advantage of the Unix shell.
The Unix shell is both a **command-line interface** (CLI) and a scripting language,
allowing such repetitive tasks to be done automatically and fast.
With the proper commands, the shell can repeat tasks with or without some modification
as many times as we want.
Using the shell, the task in the literature example can be accomplished in seconds.

## The Shell

The shell is a program where users can type commands.
With the shell, it's possible to invoke complicated programs like climate modeling software
or simple commands that create an empty directory with only one line of code.
The most popular Unix shell is Bash (the Bourne Again SHell ---
so-called because it's derived from a shell written by Stephen Bourne).
Bash is the default shell on most modern implementations of Unix and in most packages that provide
Unix-like tools for Windows.
Note that 'Git Bash' is a piece of software that enables Windows users to use a Bash like interface
when interacting with Git.

Using the shell will take some effort and some time to learn.
While a GUI presents you with choices to select, CLI choices are not automatically presented to you,
so you must learn a few commands like new vocabulary in a language you're studying.
However, unlike a spoken language, a small number of "words" (i.e. commands) gets you a long way,
and we'll cover those essential few today.

The grammar of a shell allows you to combine existing tools into powerful
pipelines and handle large volumes of data automatically. Sequences of
commands can be written into a *script*, improving the reproducibility of
workflows.

In addition, the command line is often the easiest way to interact with remote machines
and supercomputers.
Familiarity with the shell is near essential to run a variety of specialized tools and resources
including high-performance computing systems.
As clusters and cloud computing systems become more popular for scientific data crunching,
being able to interact with the shell is becoming a necessary skill.
We can build on the command-line skills covered here
to tackle a wide range of scientific questions and computational challenges.

## Read Evaluate Print Loop

When the shell is first opened, you are presented with a **prompt**,
indicating that the shell is waiting for input.

1) the shell presents a prompt (like `$`)
2) user types a command and presses the enter (or return) key
3) the computer reads it
4) the computer executes it and prints its output (if any)
Loop from step #4 back to step #1

## Reasons to learn about the shell

   - Many bioinformatics tools can only process large data in the command line version not the GUI.
   - The shell makes your work less boring (same set of tasks with a large number of files)
   - The shell makes your work less error-prone
   - The shell makes your work more reproducible.
   - Many bioinformatic tasks require large amounts of computing power

> ## Let's call some programs
>
> A command to find which user we are:
>
>```
>whoami
>```
>
>Or a command to find which shell we are using:
>
>```
>echo $SHELL
>```
>{: .discussion}

## The Filesystem




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

