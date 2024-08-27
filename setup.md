---
title: Setup
---

## Open a new shell
FILL IN INSTRUCTIONS FOR CURRENT WORKSHOP

:::::::::::::::::::::::::::::::::::::::::  callout

## Where to type commands: How to open a new shell

The shell is a program that enables us to send commands to the computer and receive output.
It is also referred to as the terminal or command line.

Some computers include a default Unix Shell program.
The steps below describe some methods for identifying and opening
a Unix Shell program if you already have one installed.
There are also options for identifying and downloading a Unix Shell program,
a Linux/UNIX emulator, or a program to access a Unix Shell on a server.

If none of the options below address your circumstances,
try an online search for: Unix shell [your computer model] [your operating system].


::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::: solution

### Windows {#windows}

Computers with Windows operating systems do not automatically have a Unix Shell program
installed.
In this lesson, we encourage you to use an emulator included in [Git for Windows][install_shell],
which gives you access to both Bash shell commands and Git.

Once installed, you can open a terminal by running the program Git Bash from the Windows start
menu.

**For advanced users:**

As an alternative to Git for Windows you may wish to [Install the Windows Subsystem for Linux][wsl]
which gives access to a Bash shell command-line tool in Windows 10 and above.

Please note that commands in the Windows Subsystem for Linux (WSL) may differ slightly
from those shown in the lesson or presented in the workshop.

::::::::::::

:::::::::::: solution

### MacOS {#macos}

For a Mac computer running macOS Mojave or earlier releases, the default Unix Shell is Bash.
For a Mac computer running macOS Catalina or later releases, the default Unix Shell is Zsh.
Your default shell is available via the Terminal program within your Utilities folder.

To open Terminal, try one or both of the following:

- In Finder, select the Go menu, then select Utilities.
  Locate Terminal in the Utilities folder and open it.
- Use the Mac 'Spotlight' computer search function.
  Search for: `Terminal` and press <kbd>Return</kbd>.

To check if your machine is set up to use something other than Bash,
type `echo $SHELL` in your terminal window.

If your machine is set up to use something other than Bash,
you can run it by opening a terminal and typing `bash`.

[How to Use Terminal on a Mac][mac-terminal]

::::::::::::

:::::::::::: solution

### Linux {#linux}

The default Unix Shell for Linux operating systems is usually Bash.
On most versions of Linux, it is accessible by running the
[Gnome Terminal][gnome-terminal] or [KDE Konsole][kde-konsole] or [xterm],
which can be found via the applications menu or the search bar.
If your machine is set up to use something other than Bash,
you can run it by opening a terminal and typing `bash`.

::::::::::::

{% include links.md %}
