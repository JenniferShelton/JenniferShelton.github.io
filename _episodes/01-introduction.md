---
title: "Introduction"
teaching: 10
exercises: 10
questions:
- "What is a command shell and why would I use one?"
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



When the shell is first opened, you are presented with a **prompt**,
indicating that the shell is waiting for input.

```bash
$
```

~~~
echo "finished running"
~~~
{: .language-bash}

The shell typically uses `$ ` as the prompt, but may use a different symbol.
In the examples for this lesson, we'll show the prompt as `$ `.
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
jshelton@localhost $
```

There are many ways for a user to interact with a computer.
In a Commandline Interface the user passes commands to the computer as lines of text.

- What is a Read Evaluate Print Loop (REPL)?

1) the shell presents a prompt (like $)
2) user types a command and presses the enter (or return) key
3) the computer reads it
4) the computer executes it and prints its output (if any)
Loop from step #4 back to step #1

```bash
$ ls
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Command not found

If the shell can't find a program whose name is the command you typed, it
will print an error message such as:

```bash
$ ks
```

```output
ks: command not found
```

This might happen if the command was mis-typed or if the program corresponding to that command
is not installed.


::::::::::::::::::::::::::::::::::::::::::::::::::

