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
- `wc` counts lines, words, and characters in its inputs (see [SWC primer](https://swcarpentry.github.io/shell-novice/04-pipefilter.html)).
- `cat` displays the contents of its inputs.
- `sort` sorts its inputs (see [SWC primer](https://swcarpentry.github.io/shell-novice/04-pipefilter.html)).
- `head` displays the first 10 lines of its input by default without additional arguments.
- `tail` displays the last 10 lines of its input by default without additional arguments.
- `command > [file]` redirects a command’s output to a file **(overwriting any existing content)**.
- `command >> [file]` appends a command’s output to a file.
- `[first] | [second]` is a pipeline: the output of the first command is used as the input to the second.
- The best way to use the shell is to use pipes to combine simple single-purpose programs (filters)
- A `for` loop repeats commands once for every thing in a list.
- Every for loop needs a variable to refer to the thing it is currently operating on.
- Use `$name` to expand a variable (i.e., get its value). `${name}` can also be used.
- Do not use spaces, quotes, or wildcard characters such as `*` or `?` in filenames, as it complicates variable expansion.
- Give files consistent names that are easy to match with wildcard patterns to make it easy to select them for looping.
- Use the <kbd>↑</kbd>/<kbd>↓</kbd> keys to scroll through previous commands to edit and repeat them.
- Use <kbd>Ctrl</kbd>+<kbd>R</kbd> to search through the previously entered commands (see [SWC primer](https://swcarpentry.github.io/shell-novice/06-script.html)).
- Use `history` to display recent commands, and `![number]` to repeat a command by number (see [SWC primer](https://swcarpentry.github.io/shell-novice/06-script.html)).
- Save commands in files (usually called shell scripts) for re-use.
- `bash [filename]` runs the commands saved in a file.
- `$@` refers to all of a shell script’s command-line arguments (see [SWC primer](https://swcarpentry.github.io/shell-novice/06-script.html)).
- `$1`, `$2`, etc., refer to the first command-line argument, the second command-line argument, etc.
- Place variables in quotes if the values might have spaces in them.

---

## Pipes and Filters

Now that we know a few basic commands, we can finally look at the shell’s most powerful feature: 
the ease with which it lets us combine existing programs in new ways. 
We’ll start with the directory shell-lesson-data/exercise-data/alkanes that contains six files 
describing some simple organic molecules. The .pdb extension indicates that these files are in 
Protein Data Bank format, a simple text format that specifies the type and position of each atom in the molecule.

```bash
$
```


