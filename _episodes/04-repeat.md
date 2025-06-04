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

**>** : redirects output to a file

The greater than symbol, <kbd>></kbd>, tells the shell to redirect the command’s output to a file instead of printing it to the screen. This command prints no screen output, because everything that wc would have printed has gone into the file lengths.txt instead. If the file doesn’t exist prior to issuing the command, the shell will create the file. **If the file exists already, it will be silently overwritten, which may lead to data loss.** Thus, redirect commands require caution.

```bash
$ cat /data/alignment/references/*/*fa | grep -c "^>" > ~/workshop/output/seq_counts.txt
$ cat ~/workshop/output/seq_counts.txt
```
**>>** : appends output to the end of a file

```bash
$ cat /data/alignment/references/*/*fa | grep -c "^>" >> ~/workshop/output/seq_counts.txt
$ cat ~/workshop/output/seq_counts.txt
```



