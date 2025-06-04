---
title: "Repeating yourself"
teaching: 30
exercises: 20
questions:
- How can I combine existing commands to produce a desired output?
- How can I perform the same actions on many different files?
- How can I save and re-use commands?
objectives:
- Explain the advantage of linking commands with pipes and filters.
- Redirect a command’s output to a file.
- Write a loop that applies one or more commands separately to each file in a set of files.
- Trace the values taken on by a loop variable during execution of the loop.
- Explain the difference between a variable’s name and its value.
- Write a shell script that runs a command or series of commands for a fixed set of files.
- Run a shell script from the command line.
keypoints:
- wc counts lines, words, and characters in its inputs (see [SWC primer](https://swcarpentry.github.io/shell-novice/04-pipefilter.html)).
- cat displays the contents of its inputs.
- sort sorts its inputs (see [SWC primer](https://swcarpentry.github.io/shell-novice/04-pipefilter.html)).
- head displays the first 10 lines of its input by default without additional arguments.
- tail displays the last 10 lines of its input by default without additional arguments.
- command > [file] redirects a command’s output to a file (overwriting any existing content!!).
- command \>\> [file] appends a command’s output to a file.
- first <kbd>|</kbd> second is a pipeline. The output of the first command is used as the input to the second.
- The best way to use the shell is to use pipes to combine simple single-purpose programs (filters)
- A for loop repeats commands once for every thing in a list.
- Every for loop needs a variable to refer to the thing it is currently operating on.
- Use $name to expand a variable (i.e., get its value). ${name} can also be used.
- Do not use spaces, quotes, or wildcard characters such as * or ? in filenames, as it complicates variable expansion.
- Give files consistent names that are easy to match with wildcard patterns to make it easy to select them for looping.
- Use the <kbd>↑</kbd>/<kbd>↓</kbd> keys to scroll through previous commands to edit and repeat them.
- Use <kbd>Ctrl</kbd>+<kbd>R</kbd> to search through the previously entered commands (see [SWC primer](https://swcarpentry.github.io/shell-novice/06-script.html)).
- Use history to display recent commands, and !number to repeat a command by number (see [SWC primer](https://swcarpentry.github.io/shell-novice/06-script.html)).
- Save commands in files (usually called shell scripts) for re-use.
- bash filename runs the commands saved in a file.
- $@ refers to all of a shell script’s command-line arguments (see [SWC primer](https://swcarpentry.github.io/shell-novice/06-script.html)).
- $1, $2, etc., refer to the first command-line argument, the second command-line argument, etc.
- Place variables in quotes if the values might have spaces in them.
---

## Pipes and Filters

Now that we know a few basic commands, we can finally look at the shell’s most powerful feature: 
the ease with which it lets us combine existing programs in new ways. 
We’ll start with the directory `/data/alignment/references/GRCh38_1000genomes/` that contains `GRCh38` reference files. 
The `.fa` extension indicates that a files is in [FASTA format](https://zhanggroup.org/FASTA/), a simple text format that specifies the nucleic or amino acid sequences.

**head** : returns the first lines of a file (default number of lines is 10)

**tail** : returns the last lines of a file (default number of lines is 10)

```bash
$ cd /data/alignment/references/GRCh38_1000genomes/
$ head GRCh38_full_analysis_set_plus_decoy_hla.fa
$ tail GRCh38_full_analysis_set_plus_decoy_hla.fa
```

**grep** : program to find lines that match a pattern

**^** : regex (regular expression) which matches the first character

**|** : pipes output from the command on the left as input to the command on the right

The vertical bar, <kbd>|</kbd>, between the two commands is called a pipe. It tells the shell that we want to use the output of the command on the left as the input to the command on the right. Nothing prevents us from chaining pipes consecutively. We can for example send the output of `head` directly to `grep`, and then send the resulting output to `sort` (a command that sorts lines of text). This removes the need for any intermediate files.

This idea of linking programs together is why Unix has been so successful. Instead of creating enormous programs that try to do many different things, Unix programmers focus on creating lots of simple tools that each do one job well, and that work well with each other. This programming model is called ‘pipes and filters’. We’ve already seen pipes; a filter is a program like wc or sort that transforms a stream of input into a stream of output. Almost all of the standard Unix tools can work this way. Unless told to do otherwise, they read from standard input, do something with what they’ve read, and write to standard output.

The key is that any program that reads lines of text from standard input and writes lines of text to standard output can be combined with every other program that behaves this way as well. You can and should write your programs this way so that you and other people can put those programs into pipes to multiply their power.

```bash
$ # returns lines that contain GRCh38
$ cat GRCh38_full_analysis_set_plus_decoy_hla.fa | grep "GRCh38"
$ # returns lines that start with (^) the > symbol
$ cat GRCh38_full_analysis_set_plus_decoy_hla.fa | grep "^>"
$ # counts lines that start with (^) the > symbol
$ cat GRCh38_full_analysis_set_plus_decoy_hla.fa | grep -c "^>"
```

**\*** : glob (this translates to zero of more of any character)

```bash
$ ls /data/alignment/references/*/*fa
$ cat /data/alignment/references/*/*fa | grep -c "^>"
```

## Writing to files

**>** : redirects output to a file

The greater than symbol, <kbd>></kbd>, tells the shell to redirect the command’s output to a file instead of printing it to the screen. This command prints no screen output, because everything that wc would have printed has gone into the file lengths.txt instead. If the file doesn’t exist prior to issuing the command, the shell will create the file. **If the file exists already, it will be silently overwritten, which may lead to data loss.** Thus, redirect commands require caution.

```bash
$ cat /data/alignment/references/*/*fa | grep -c "^>" > ~/workshop/output/seq_counts.txt
$ cat ~/workshop/output/seq_counts.txt
```
**\>\>** : appends output to the end of a file

```bash
$ cat /data/alignment/references/*/*fa | grep -c "^>" >> ~/workshop/output/seq_counts.txt
$ cat ~/workshop/output/seq_counts.txt
```

## Loops

Loops are a programming construct which allow us to repeat a command or set of commands for each item in a list. As such they are key to productivity improvements through automation. Similar to wildcards and tab completion, using loops also reduces the amount of typing required (and hence reduces the number of typing mistakes).

Suppose we have dozens of genome sequence alignment files in [SAM, BAM or CRAM format](https://samtools.github.io/hts-specs/SAMv1.pdf). For this example, we’ll use the `/data/*/*cram` files which only have two example files, but the principles can be applied to many many more files at once.

The structure of these files is the same and that allows many tools (including QC tools to evaluate them). Let’s look at the files:

```bash
$ # the header describes the reference and read groups
$ samtools view -H /data/alignment/combined/NA12878.dedup.bam
$ # the first three alignments are printed to the screen to look at
$ samtools view /data/alignment/combined/NA12878.dedup.bam | head -3
```

We would like to get QC metrics for each alignment file. For each file, we would need to execute the command `samtools flagstat`. We’ll use a loop to solve this problem, but first let’s look at the general form of a loop, using the pseudo-code below:

```bash
# The word "for" indicates the start of a "For-loop" command
for thing in list_of_things 
#The word "do" indicates the start of job execution list
do 
    # Indentation within the loop is not required, but aids legibility
    operation_using/command $thing 
# The word "done" indicates the end of a loop
done
```

and we can apply this to our example like this:

```bash
for cram in /data/*/*cram; do
    echo ${cram}
    echo "working..."
    samtools flagstat ${cram}
    echo "done"
done
```

> ## Follow the Prompt
> The shell prompt changes from `$` to `>` and back again as we were typing in our loop. The second prompt, `>`, is
> different to remind us that we haven’t finished typing a complete command yet. A semicolon, `;`, can be used
> to separate two commands written on a single line.
{: .testimonial}

When the shell sees the keyword for, it knows to repeat a command (or group of commands) once for each item in a list. Each time the loop runs (called an iteration), an item in the list is assigned in sequence to the variable, and the commands inside the loop are executed, before moving on to the next item in the list. Inside the loop, we call for the variable’s value by putting `$` in front of it. The `$` tells the shell interpreter to treat the variable as a variable name and substitute its value in its place, rather than treat it as text or an external command.

In this example, the list is three filenames: `/data/cancer_genomics/COLO-829BL_1B.variantRegions.cram` and `/data/cancer_genomics/COLO-829_2B.variantRegions.cram`. Each time the loop iterates, we first use echo to print the value that the variable `${cram}` currently holds. This is not necessary for the result, but beneficial for us here to have an easier time to follow along. Next, we will run the flagstat command on the file currently referred to by `${cram}`. The first time through the loop, `${cram}` is `COLO-829BL_1B.variantRegions.cram`. The interpreter runs the command head on `COLO-829BL_1B.variantRegions.cram` and prints the QC metrics. For the second iteration, `${cram}` becomes `COLO-829_2B.variantRegions.cram`. This time, the shell runs head on `COLO-829_2B.variantRegions.cram` and prints the QC metrics. Since the list was only two items, the shell exits the for loop.

> ## Same Symbols, Different Meanings
> Here we see `>` being used as a shell prompt, whereas `>` is also used to redirect output. Similarly, `$` is
> used as a shell prompt, but, as we saw earlier, it is also used to ask the shell to get the value of a variable.
>
> If the shell prints `>` or `$` then it expects you to type something, and the symbol is a prompt.
>
> If you type `>` or `$` yourself, it is an instruction from you that the shell should redirect output or get
> the value of a variable.
{: .testimonial}

When using variables it is good practice to put the names into curly braces to clearly delimit the variable name: `$cram` is equivalent to `${cram}`, but is different from `${cr}am`. You may find this notation in other people’s programs.

We have called the variable in this loop cram in order to make its purpose clearer to human readers. The shell itself doesn’t care what the variable is called; if we wrote this loop as:

```bash
for x in /data/*/*cram; do
    echo ${x}
    echo "working..."
    samtools flagstat ${x}
    echo "done"
done
```

It would work exactly the same way. Don’t do this. Programs are only useful if people can understand them, so meaningless names (like x) or misleading names (like `temperature`) increase the odds that the program won’t do what its readers think it does.

In the above examples, the variables (thing, filename, x and temperature) could have been given any other name, as long as it is meaningful to both the person writing the code and the person reading it.

Note also that loops can be used for other things than filenames, like a list of numbers or a subset of data.

