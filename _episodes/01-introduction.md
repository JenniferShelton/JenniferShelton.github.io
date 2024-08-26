---
title: "Introduction"
teaching: 30
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

## The Shell
To start we will open a terminal.
 
1. Go to the [link](https://docs.google.com/spreadsheets/d/1t78ladvpt8t-uEaG1EeYRATiBIrelZMBJcKANbAGdHw/edit?usp=sharing) given to you at the workshop
2. Select "Terminal" from the "JupyterLab" launcher (or blue button with a plus in the upper left corner)
3. After you have done this put up a green sticky not if you see a flashing box next to a `$`

When the shell is first opened, you are presented with a **prompt**,
indicating that the shell is waiting for input.

```bash
$
```

The shell typically uses `$` as the prompt, but may use a different symbol.
In the examples for this lesson, we'll show the prompt as `$`.
Most importantly, *do not type the prompt* when typing commands.
Only type the command that follows the prompt.
This rule applies both in these lessons and in lessons from other sources.
Also note that after you type a command, you have to press the <kbd>Enter</kbd> key to execute it.

The prompt is followed by a **text cursor**, a character that indicates the position where your
typing will appear.
The cursor is usually a flashing or solid block, but it can also be an underscore or a pipe.
You may have seen it in a text editor program, for example.

Note that your prompt might look a little different. In particular, most popular shell
environments by default put your user name and the host name before the `$`. Such
a prompt might look like, e.g.:

```bash
student@workshop-1:~$ 
```

There are many ways for a user to interact with a computer. For example, we often use a Graphical User 
Interface (GUI). With a GUI we might roll a mouse to the logo of a folder and click or tap (on a touch 
screen) to show the content of that folder. In a Commandline Interface the user can do all of the same 
actions (e.g. show the content of a folder). On the Commandline the user passes commands to the 
computer as lines of text.

- What is a Read Evaluate Print Loop (REPL)?

1. the shell presents a prompt (like `$`)
2. user types a command and presses the <kbd>Enter</kbd> key
3. the computer reads it
4. the computer executes it and prints its output (if any)
5. loop from step #4 back to step #1

The most basic command is to call a program to perform its default action. For example, call the program `whoami` to 
return your username. 

```bash
$ whoami
```

You can also call a program and pass arguments to the program. For example, the `ls` command will list the contents of 
your current directory (directory is synonymous with folder). Any line that starts with `#` will not be executed. We can 
write comments to ourselves by starting the line with `#`. 

```bash
# call ls to list current directory
ls 

# pass one or more paths of files or directories as argument(s)
ls /bin
```

In a GUI you may customize your finder/file browser based on how you like to search. In general if you can do it on a GUI 
there is a way to use text to do it on the commandline. I like to see my most recently changed files first. I also like 
to see the date they were edited.

```bash
# call ls to list bin and show the most recently changed files first (with the `-t` option/flag)
ls -t /usr

# add the `-l` to show who owns the file, file size, and what date is was last edited
ls -t -l /usr
```

The basic syntax of a unix command is:

1. call the program
2. pass any flags/options
3. pass any "order dependent arguments"

## Getting help

`ls` has lots of other options. There are common ways to find out how to use a command and what options it 
accepts (depending on your environment). Today we will call the unix command and the use the flag `--help`.

```bash
$ ls --help
```

```
Usage: ls [OPTION]... [FILE]...
List information about the FILEs (the current directory by default).
Sort entries alphabetically if none of -cftuvSUX nor --sort is specified.

Mandatory arguments to long options are mandatory for short options too.
  -a, --all                  do not ignore entries starting with .
  -A, --almost-all           do not list implied . and ..
      --author               with -l, print the author of each file
  -b, --escape               print C-style escapes for nongraphic characters
```
{: .output}


Help menus show you the basic syntax of the command. Optional elements are shown in square brackets. Ellipses indicate that you can type more than one of the elements.

Help menus show both the long and short version of the flags. Use the short option when typing commands directly into the shell to minimize keystrokes and get your task done faster. Use the long option in scripts to provide clarity. It will be read many times and typed once.

> ## When `--help` does not work
> If `--help` does not show you a help menu there are other common ways to show help menus that you can try.
> 
> - call the program `man` and pas the name of the command that you are curious about as the argument (`man ls`). 
> Type `q` to exit out of this help screen.
> - some bioinformatics programs will show a help menu if you call the tool without any flags or arguments.
{: .callout}

> ## Command not found
> 
> If the shell can't find a program whose name is the command you typed, it
> will print an error message such as:
>
> ```bash
> $ ks
> ```
> > ## Solution
> >
> > ```
> > ks: command not found
> > ```
> > {: .output}
> > This might happen if the command was mis-typed or if the program corresponding to that command
> > is not installed. When your get an error message stay calm and give it a couple of read throughs. Error 
> > messages can seem akwardly worded at first but they can really help guide your debugging.
> {: .solution}
{: .challenge}


## The Cloud

There are a number of reasons why accessing a remote machine is invaluable to any scientists working 
with large datasets. In the early history of computing, working on a remote machine was standard 
practice - computers were bulky and expensive. Today we work on laptops or desktops that are more 
powerful than the sum of the world’s computing capacity 20 years ago, but many analyses (especially in 
genomics) are too large to run on these laptops/desktops. These analyses require larger machines, 
often several of them linked together, where remote access is the only practical solution. 

Schools and research organizations often link many computers into one High Performance Computing (HPC) 
cluster on or near the campus. Another model that is becoming common is to "rent" space on a 
cluster(s) owned by a large company (Amazon, Google, Microsoft, etc). In recent years, computational 
power has become a commodity and entire companies have been built around a business model that allows 
you to “rent” one or more linked computers for as long as you require, at lower cost than owning the 
cluster (depending on how often it is used vs idle, etc). This is the basic principle behind the 
cloud. You define your computational requirements and off you go.

The cloud is a part of our everyday life (e.g. using Amazon, Google, Netflix, or an ATM involves remote computing). The topic is fascinating, but this lesson says a few minutes or less so let’s get back to working on it for the workshop.

For this workshop starting a vm and setting up your working environment has been done for you. Going forward reach out to your organizations system administrators for your cluster for suggestions. To read more on your own here are lessons about working on the cloud and a local HPC. Additional lesson's that you can run from a remote computer:
- [Using a HPC Cluster](https://carpentries-incubator.github.io/hpc-intro/)
- [Cloud computing](https://datacarpentry.org/cloud-genomics)




> ## HUMAN GENOMIC DATA & SECURITY
> Note that if you are working with human genomics data there might be ethical and legal considerations that affect your choice of cloud resources to use. The terms of use, and/or the legislation under which you are handling the genomic data, might impose heightened information security measures for the computing environment in which you intend to process it. This is a too broad topic to discuss in detail here, but in general terms you should think through the technical and procedural measures needed to ensure that the confidentiality and integrity of the human data you work with is not breached. If there are laws that govern these issues in the jurisdiction in which you work, be sure that the cloud service provider you use can certify that they support the necessary measures. Also note that there might exist restrictions for use of cloud service providers that operate in other jurisdictions than your own, either by how the data was consented by the research subjects or by the jurisdiction under which you operate. Do consult the legal office of your institution for guidance when processing human genomic data.
{: .callout}

