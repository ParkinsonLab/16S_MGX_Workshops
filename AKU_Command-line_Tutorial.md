# Command Line Tutorial
## Introduction 
The tutorials will be conducted in a virtual machine (VM), a simulation of a computer environment that can be shared and distributed to others. This allows us to standardize the workshop, so that everyone is working within the same environment. As a result, the installation of all the packages you will need for the workshop has been done within these smaller computers. The only thing you will need is the actual software to run the VM. 
We are using Oracle VirtualBox, an open source software for running VMs. If you've never used Oracle VirtualBox, or used command-line, this tutorial will run you through the necessary tools and syntax we need to know to go forward. All of the essential commands for running the workshop are provided in the text, but this will help you understand what is being done. 

## Lab 1A - Starting the virtual machine
### Importing and Starting
First make sure you are running the latest version of Oracle VirtualBox (v7.1.2), which you can check at `Help > About VirtualBox`. You can import the 16S and MGX VMs by navigating to `File > Import Appliance`. Click the folder icon next to the "File" field and navigate to the virtual box. From there, you can just select `Next`, then `Finish`.  Importing might take a little while!

To open the VM, simply select `Start` at the toolbar on the top. Under settings, make sure the base memory is set to a number appropriate for your machine. You might also have to change the network setting to `Attached to: NAT`. Once the VM boots up, you will see a desktop very similar to your own! You are now within the VM environment, a computer within a computer. To close the VM, select `File > Close...` and then `Power off the machine`. Any files you have created will be saved. Keep your VM on, as we will now run a couple of additional steps to make the VM more usable. 

### Guest Additions
With your VM powered on, we can add Guest Additions to improve its useability. One improvement allows us to import files directly from your host computer to the VM, which you will need to import the workshop files for the MGX tutorials. 

First, navigate to `Devices > Insert Guest Additions CD image...`. You should see a disk pop up on the quickbar to the left - click on that disk and right click within the explorer window. Select `Open in Terminal`. This is our accesspoint to command-line, a command-line interface, which lets us talk directly to the computer. The rest of the tutorial will take place mainly in this interface. To tell your VM to perform functions, you can enter "commands". We'll get into the important commands you'll need to know soon, but for now just enter the following:

```bash
./autorun.sh
```

This tells the VM to run the `autorun.sh` script, which itself contains commands that will install Guest Additions to the VM (this will leave your computer unaffected!). It may ask you for a password - enter `mg2024`. Next, we can create a shared folder to move files into or out of the VM. In the same terminal window, type:

```bash
sudo adduser $USER vboxsf
```

When the password prompt shows up, again type in `mg2024` (you won't be able to see text appear but that is normal). 

Exit and re-start the VM to let these permissions take place. This command tells the VM to allow access to shared folders. Once the VM is back up, select `Devices > Shared Folders > Shared Folder Settings`. Click the Folder icon (with the green plus) to the right. Pick "Other" on the folder path and select a folder on your computer that you want to use to import files. You can leave the other two fields empty, and select `Auto-mount` and `Make permanent`. Now, try moving a file on your computer into the shared folder (this will appear with a hard drive icon in your file explorer sidebar) and opening it up in the VM. We're now all set to tackle command line!

## Lab 1B - Command Line
While you have already run two commands in the previous section (and some of you will likely be very familiar with command-line), we will go through a quick tutorial to get everyone up to speed on the operation of the tutorials. As mentioned before, command-line interfaces let us talk to the computer and tell it to do things. You'll notice that the UI is not very informative or intuitive at all. Using the terminal can be difficult for newcomers, but we'll try to start off slow. We'll go over the basic tools you'll need. 

The terminal is a text-based interface where you can enter commands, often line-by-line. Command-line is very exact - the command has to be in the correct syntax and you have to be in the right folder or it won't work. These are the two of the biggest sources of errors for those first learning command-line - either they are not in the correct folder or they have switched a parameter or missed a dash somewhere. A typical command might look like this:
```bash
# example command, do not run
wc -l file.txt
```

There are three components here. From left to right:
- `wc` is the command
- `-l` is a flag
- `file.txt` is the input

Sometimes you won't have a flag or input, but you will almost always have a command. The command tells the computer what action it is to perform. The flag modifies this behaviour in a specific way (each command has its own unique flags). The input is telling the computer what to act upon. In this example, `wc` is a command that counts words. The `-l` flag is telling it to count lines rather than an individual word count. Altogther, this command is telling the computer to return the number of lines in `file.txt`. 

Try navigating to the 16S workshop folder and opening a terminal in there. Now, type `ls`. You should see a list of files and folders - `ls` is likely the command you will use the most, as it tells you what the terminal sees right now. This ties into the previous point about being in the wrong folder. If you are running a command that fails, type `ls` to make sure you can see the file. If the file is in another folder, you can type `ls folder_a/` to look at all the files and folders within `folder_a/` You can also use `..` as a shortcut for "parent directory". Typing `ls ..` shows you the files and folders of the folder one level above yours. Some important flags for `ls` are `-lh`. This lists all files with a human-readable output of size, date, and other fields. 

Another command you will need to know is `cd`. This moves the terminal window to a different folder. You can use this to either enter folders you see currently (ex. `cd folder_a/`), or you can enter the "absolute path" of a folder that might take too long to navigate by relative location (ex. `cd /home/mg_user/Desktop/16S_Workshop/`). If you get confused about where you are, use `pwd` - this will give you your current folder name and all parent directories. Try navigating to the subdirectories you see with `ls` and navigating back. 

> Question 1 - What do you think `cd ..` does?
> 
> Answer - You will navigate to the parent directory

Similarly, `mv` moves files to a different location. This will remove the original copy. `mv file.txt folder_a/` moves `file.txt` into `folder_a/`. You can also use this command to rename files. If the second parameter ends in a file name, the original file will be renamed to that file name. `mv file.txt renamed.txt` will rename `file.txt` into `renamed.txt`. 

> Question 2- What does `mv file.txt folder_a/renamed.txt` do?
> 
> Answer - This both moves `file.txt` to `folder_a/`, and at the same time renames it to `renamed.txt`

If you want a file to move somewhere else without removing the original file, you can use `cp`. If you want to remove a file entirely, use `rm`. Be aware that there is no undo in command-line! Do not use `rm` unless you are:
1. Sure you will not use this file again
2. Sure the command is typed in correctly

Now that we know how to rename and delete files, we should address how to make or read them. To make an empty file - perhaps to later fill in with a text editor - you can use `touch`. Ex. `touch list.txt` will make an empty text file named "list". You can also make folders using `mkdir`.

> Question 3 - For practice, try making a file named "test.txt", making a new folder named "test_directory", and moving your text file into the folder. 
>
> Answer - You could do this with the following commands:
```bash
touch test.txt
mkdir test_directory
mv test.txt test_directory/
```

If you wanted to read a file, there are several ways to do so. I will list them and some use cases here:
- `cat` simply prints the entire file on the screen. This is simple and easy but will be impractical for larger files
- `less` opens an interface for scrolling through the file
- `head` displays the first few lines of the file. Useful for looking at the columns of a `.tsv` file for example.
- `tail` displays the last few lines of the file. Useful for reading what happened most recently in a log file. 
	- Both `head` and `tail` have the `-n` flag which lets you specify how many lines to be read. 

Change directory to the `chicken_samples/` folder and try all four methods on `metadata.tsv`. 

The last command we will go over is `grep`. `grep` looks through a file and returns each instance of a given pattern within that file. For example, `grep chicken list.txt` looks for every appearance of "chicken" in the file "list.txt". This can be useful especially when sorting through large tabulated data, to find relevant entries. Finally, this is not a command, but a useful piece of command-line grammar: the pipe `|` is used to direct output from one command to the input of another. `head -n 5 test.txt | wc` will find the word count of the first 5 lines of test.txt. The output of `head` is being used as the input for `wc`. 

> Question 4 - How would you use the pipe to find the number of lines that contain the word "bacteria" in a document called `taxonomy.txt`
> 
> Answer - `grep bacteria taxonomy.txt | wc -l`

That concludes the command-line tutorial! As mentioned before, this can be very overwhelming for newcomers. There will be commands not mentioned here that we use in the tutorials. Most of those are tool-specific, like `fastqc` being a command from the FASTQC toolset for sequence quality checking. Regardless, feel free to ask whenever something looks confusing!
