# Getting Started

General instructions on the use of these tutorials

## Summary

This document describes how to connect to the lynx server, how the directory structure should be set up, how commands can be copied and executed from tutorial instructions, how text files can be used and edited, and how files required for tutorials can be downloaded from the tutorial repository.

## Table of contents

* [Connecting to lynx](#connecting)
* [Directory structure](#directory_structure)
* [Copying and pasting commands](#copy_pasting)
* [Working with text files](#text_files)
* [Downloading tutorial files](#downloading)

<a name="connecting"></a>
## Connecting to lynx

To work on the lynx server, you will need a console program (see [Requirements](requirements.md)). The console program will be important to connect to lynx, to execute commands on it, and to transfer files between lynx and your own computer. The latter will be required occasionally, for example to look at results with a GUI program.

You will recognize console commands in the tutorials by monospace font, gray background, and an outline, like for example this command:

		pwd
		
To execute commands like the above, the "Enter" key always needs to be hit after writing the command. The command `pwd` stands for "_print working directory_", which is exactly what it does.

* To access lynx, use the following command:

		ssh YOURUSERNAME@10.153.166.2
		
	You will be asked for your password, which you will need to type to get access. User names and passwords will be distributed on the first course day.

* If you need to upload data to lynx, you can use the following command:

		scp -R YOURFILE YOURUSERNAME@10.153.166.2:PATH_TO_YOUR_FOLDER_ON_LYNX

* Downloading from lynx works very similar. However, it's best to first navigate to the directory on your own computer into which you would like to download the files. Then, use the following command for the download:

		scp -R YOURUSERNAME@10.153.166.2:PATH_TO_YOUR_FOLDER_ON_LYNX/YOURFILE .
		
	Note that instead of a file name in the above two commands (`YOURFILE`), you could also specify the name of a directory.

<a name="directory_structure"></a>
## Directory structure

After connecting to lynx with the `ssh` command described above, your starting point in the directory structure of the server is your home directory. Typing the `pwd` command should print the path to your home directory: `/home/YOURUSERNAME`.

* Make a new directory inside of your home directory for all of your analyses in this course:

		mkdir phylogenomics

* Then, navigate into that directory:

		cd  phylogenomics
		
* As each individual tutorial of this course will generate a number of result files, it may be best to use subdirectories for each tutorial. To make these subdirectories, use `mkdir` followed by the tutorial name:

		mkdir getting_started
		mkdir ml_phylogeny_inference
		mkdir bayesian_phylogeny_inference
		...


<a name="copy_pasting"></a>
## Copying and pasting commands

On lynx, it should be possible to execute all commands that are mentioned in the tutorial instructions exactly as they are specified. All one-line commands can be copied from the instructions and pasted into the console window, but care should be taken not to include whitespace symbols before the first letter or after the last letter of the command. Sometimes, it is not evident from the text selection that is used for copying, but whitespace nevertheless appears after the last letter copied to the command line (recognizable as a space between the last letter and the cursor); in that case, the whitespace on the command line should be deleted before executing the command.

For multi-line commands, however, it is safer to copy and paste each line individually, and hit the Enter key each time a line has been copied (again, whitespace before and after the first and last letters of the line should not be copied).
	
* Try to copy and paste this command into the console window, and then execute it:

		for n in {1..10}
		do
			echo ${n}
		done

	This should result in the numbers 1 to 10 being written line by line in the console window; if this is not the resulting output, something must have gone wrong with copy-pasting.

<a name="text_files"></a>
## Working with text files

In some cases, text files will need to be written or edited, and for this, a text editor of some sort will be required. For example, instead of copying many lines of commands into the console window one by one, these could also be copied into a text file, the file could be named, e.g., `script.sh`, and then all commands from that file could be executed jointly the command

		bash script.sh
		
To write and edit text files, either a command-line text editor or one with a graphical user interface could be used. There are many suitable options in the latter category. Programs like TextEdit (on Mac OS X) or Notepad (on Windows) could be used, but only if the default settings are changed so that plain text instead of "rich text" is written. Using Microsoft Word would guarantee trouble sooner rather than later. More convenient than any of these, however, are those that include syntax highlighting, such as [BBEdit](https://www.barebones.com/products/textwrangler/) for Mac OS X, [Notepad++](https://notepad-plus-plus.org) for Windows, or [Sublime Text](https://www.sublimetext.com), [Geany](https://www.geany.org), and [Atom](https://atom.io) that run on all three platforms.

If you want to avoid installing any of these programs, you could also use one of the text editors that are usually available within the console, like Vim or Emacs. No installation should be required for these tools, but their use may not be as intuitive as that of text editors with graphical user interfaces. Short tutorials are available online for [Vim](https://www.howtoforge.com/vim-basics) and [Emacs](https://www.digitalocean.com/community/tutorials/how-to-use-the-emacs-editor-in-linux).

* Use the text editor of your choice to write the following four lines:

		for n in {1..10}
		do
			echo ${n}
		done

* Save the text to a new file named `script.sh`, and then place this script in the directory for this tutorial on lynx (`/home/YOURUSERNAME/phylogenomics/getting_started`).

* Navigate to the same directory on lynx and then execute the commands in file `script.sh` with this command:

		bash script.sh
		
	This should again result in the numbers 1 to 10 being written on separate lines.

<a name="downloading"></a>
## Downloading tutorial files

Links to input files and scripts are included in all tutorials and to run a given tutorial, the linked files should be downloaded and placed in your working directory. Most often, these files are hosted in the same GitHub repository, and will not automatically download when the link is clicked (assuming that the tutorial instructions are opened in a web browser).

* To illustrate the point that linked files from the same GitHub repository are not directly downloaded when the link is clicked, click on this link to an alignment file from the tutorial on substitution model selection:
[`locus_0001.fas`](data/locus_0001.fas)

	You should see that clicking this link opens another webpage in which the alignment file is merely embedded. The reason for this is simply the way in which GitHub repositories are organized.
	
* To get to the alignment file itself rather than the webpage in which it is embedded, click the button labeled with "Raw" at the top right of the window, just above the embedded file.

	Depending on the file and the web browser used, this should either open the file in the browser or download it directly to your computer, or a dialog window should appear that allows you to select the download location. If the file opened in the browser, use the menu of the browser (e.g. the File menu and the option "Save As...") to download the file and save it in the tutorial directory. If it downloaded directly, it can most likely be found in the Download directory of your computer, and should then be moved into the tutorial directory.
	
	A simpler way to download files is to click on the download button to the right of the "Raw" button; this should open a dialog that allows one to select the download location.
	
	Instead of downloading the file through the browser (which would download it to your own computer), you could also just copy its URL, and then download it through the command line, for example with the `wget` command:
	
		wget https://raw.githubusercontent.com/mmatschiner/phylogenomics/refs/heads/main/getting_started/data/locus_0001.fas
	
* Make sure that the file has been placed in the expected directory, either with the file manager program or the command line. To do so on the command line, you could use this command:

		ls
		
	This command lists all files in the current directory, and if that list includes the name of the alignment file (`locus_0001.fas`), the download worked fine. Also make sure that the web browser or your file system has not automatically added another file extension such as `.txt` to the file name. Depending on the settings of your file manager program, these file extensions may not be shown but nevertheless be there, so unless you are sure that your file manager is not hiding any file extensions, it may be safer to check the file name with on the command line.
	
* As a final test that the file has been downloaded correctly, you could also have a look at the content of the file with this command:

		less locus_0001.fas
		
	This should present the file content in the console window. To return to the command line, type the letter `q`.
