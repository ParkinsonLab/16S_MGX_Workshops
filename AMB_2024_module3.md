---
layout: workshop_main_2day
permalink: /AMB_2024_module3
title: AMB 2024
header1: Workshop Pages for Students
header2: Advanced Microbiome Analysis 2024: Module 3
image: /site_images/AMB_2024_v1.png
length: 2 days
---

# Module 3: Metagenomic functional annotation

This tutorial is part of the 2024 Canadian Bioinformatics Workshops Advanced Microbiome Analysis (St John's, NL, May 29-30).

**Author**: Robyn Wright

## Table of Contents

[Introduction](#Introduction)\
[3.1 Copy in the MMseqs data we will be using](#31-Copy-in-the-MMSeqs-data-we-will-be-using)\
[3.2 Run MMseqs](#32-Run-MMseqs)\
[3.3 Get MMseqs top hits](#33-Get-MMSeqs-top-hits)\
[3.4 Combine Kraken taxonomy and MMseqs function, normalise, and reformat](#34-Combine-Kraken-taxonomy-and-MMSeqs-function-normalise-and-reformat)\
[3.5 Visualise with JarrVis](#35-Visualise-with-JarrVis)\
[3.6 Run HUMAnN](#36-Run-HUMAnN)\
[3.7 Combine HUMAnN output](#37-Combine-HUMANn-output)\
[3.8 Visualise HUMAnN output](#38-Visualise-HUMANn-output)\
[3.9 AMR annotation of reads using CARD RGI](#39-AMR-annotation-of-reads-using-CARD-RGI)\
[3.10 Answers](#answers)\

## Introduction

The main goal of this tutorial is to introduce students to functional profiling of taxonomic reads using MMseqs and HUMAnN and to visualise the results of both the taxonomic and functional annotations. As in the taxonomic composition workshop (AMB module 1), we want to emphasize that there is not a one-size-fits-all pipeline for analyzing MGS data. Each of MMseqs and HUMAnN are typically used for general functional annotation, while CARD RGI is used specifically for AMR annotation of reads.

Each of the MMseqs, HUMAnN, and CARD RGI parts of this tutorial can stand alone, so feel free to start with whichever of these you feel will be most useful to you. Note that *you should run step 3.1 first regardless of which tool you want to start with.*

Throughout this module, there are some questions aimed to help your understanding of some of the key concepts. You'll find the answers at the bottom of this page, but no one will be marking them. 

### MMSeqs

[MMseqs2](https://github.com/soedinglab/MMseqs2) (**M**any-against-**M**any **seq**uence **s**earching) is a software suite to search and cluster huge protein and nucleotide sequence sets. We'll be using MMseqs to assign functions to our samples on a read-by-read basis by mapping them to the UniRef90 protein database, which allows us to link the function with the taxonomy that we've obtained from Kraken2 (although MMseqs can also be used for taxonomy assignment). MMseqs2 works by taking sequenced reads, translating them into protein and then mapping them against this protein database (in this case, [UniRef90](https://www.uniprot.org/help/uniref), a large protein database clustered at 90% identity).

### HUMAnN

[HUMAnN3](https://github.com/biobakery/humann) (**HMP** **U**nified **M**etabolic **An**alysis **N**etwork) is a tool for profiling the presence/absence and abundance of microbial pathways in a community from metagenomic (or metatranscriptomic) sequencing data. HUMAnN3 works by: (1) identifying the species in the samples using MetaPhlAn, (2) mapping these reads to pangenomes of the species using Bowtie2, and (3) aligning the reads that could not be mapped to the pangenomes to a protein database (usually UniRef50) with DIAMOND. 

## 3.1 Copy in the MMSeqs data we will be using

As we've done previously, we'll start by activating the conda environment and creating a directory to use:
```
conda activate functional
cd workspace
mkdir amb_module3
cd amb_module3
```

Then, we'll create symlinks to some of the files that we'll be using:
```
ln -s ~/CourseData/MIC_data/AMB_data/raw_data/ .
ln -s ~/CourseData/MIC_data/AMB_data/amb_module1/cat_reads/ .
ln -s ~/CourseData/MIC_data/AMB_data/scripts/ .
ln -s ~/CourseData/MIC_data/AMB_data/MMSeqs2_db/ .
ln -s ~/CourseData/MIC_data/AMB_data/mgs_metadata.txt .
```

You will notice that we'll be using the reads that we concatenated in the first module (but copied across from CourseData just incase anything went wrong!) as well as the metadata and the MMSeqs database information. 

## 3.2 Run MMseqs

Now, we'll start running MMseqs2. Note that these commands can actually all be combined for each sample, but so that we can see and understand what's going on, we're going to run each of them separately.

First, make a directory to store the output:
```
mkdir mmseqs_U90_out
```

Now, we'll use ```parallel``` to create databases for all of our sample files:
```
parallel -j 4 --progress 'mmseqs createdb {} mmseqs_U90_out/mmseqs-{/.}-queryDB' ::: cat_reads/*
```
This command creates an MMseqs database from the the input fastq file. The creation of this database is necessary for MMseqs as it vastly increases the speed at which translated DNA sequences can be mapped against a protein database.

Next, we'll actually run the searches with MMseqs:
```
parallel -j 1 --progress 'mmseqs search mmseqs_U90_out/mmseqs-{/.}-queryDB MMSeqs2_db/mmseqsUniref90DB mmseqs_U90_out/mmseqs-{/.}-resultDB tmp --db-load-mode 3 --threads 2 --max-seqs 25 -s 1 -a -e 1e-5' ::: cat_reads/*
```
This command is the real meat of the job file and runs the freshly created sample database against the provided UniRef90 protien database. There are a number of parameters in this command:
- ```--db-load-mode 3``` - This parameter tells MMseqs how to deal with loading the database into memory. For more information you can check out this page. However, setting this parameter to 3 helps when running MMseqs on a cluster environment.
- ```--threads``` - The number of processors we want MMseqs to use during the search
- ```--max-seqs 25``` - This indicates that we want MMseqs to output at maximum 25 hits for each sequence
- ```-s 1``` - This indicates the sensitivity that we want MMseqs to run at. Increasing this number will lower the speed at which MMseqs runs but will increase its sensitivity. For well-explored environments such as the human gut, a setting of 1 should suffice.
- ```-a``` - This indicates that we want our results to output backtraces for each sequence match. These are needed to convert the resulting MMseqs file into a usable file format.
- ```-e 1e-5``` - This indicates that we only want to keep matches that are below an E-value of 1e-5 (E-values are a measure of how well two sequences match one another, and the closer they are to zero, the better the match is).
- ```> /dev/null 2>&1``` - We could add this part to the end of the command if we wanted to run the command without having too much text printed to our screen.

**Got an error message??**

We actually unfortunately don't have enough memory on these servers to run this command. If you haven't yet got an error message, you can stop this command with ```ctrl```+```c```.

We'll just delete any files that we could have made if you ran that command, so that we don't confuse any further steps:
```
rm mmseqs_U90_out/*resultDB*
```
Note that you may get an error saying that there's no such file or directory. That's fine! You can't remove files that don't exist. 

Copy over the output that we *would* have got from this command if we could run it:
```
cp ~/CourseData/MIC_data/AMB_data/mmseqs_U90_out/*resultDB* mmseqs_U90_out/
```

And now run the final command that allows us to convert the resulting file from the MMseqs2 format into one that is more usable:
```
parallel -j 2 --progress 'mmseqs convertalis mmseqs_U90_out/mmseqs-{/.}-queryDB MMSeqs2_db/mmseqsUniref90DB mmseqs_U90_out/mmseqs-{/.}-resultDB mmseqs_U90_out/mmseqs-{/.}-s1.m8 --db-load-mode 2' ::: cat_reads/*
```
This command is similar and takes as input the query database we made from our first command, the UniRef90 database we searched against and the resulting file from our search command. It will output the files ```mmseqs_U90_out/mmseqs-*-s1.m8```.

Again, if we didn't want to print the output of this then we could add ```> /dev/null 2>&1``` to the end of the command.

**This command will take a few minutes to run, so it's a good time for a break if you'd like one!!**

Now, we'll move these ```*.m8``` files to a new folder:
```
mkdir mmseqs_m8_files
mv mmseqs_U90_out/*.m8 mmseqs_m8_files/
```

Let's take a quick look at one of the files we just moved into the directory mmseqs_m8_files using the less command:
```less mmseqs_m8_files/mmseqs-CSM7KOMH-s1.m8```

We you will see is a file in BLAST tabular format:

| Column Number       | Data Type     |
| :------------- | :----------: |
| 0 |  query sequence ID  |
| 1 | Subject (database) sequence ID |
| 2 | 	Percent Identity |
| 3 | Alignment Length |
| 4 | Number of gaps |
| 5 | 	Number of mismatches |
| 6 | Start on the query sequence |
| 7 | End on the query sequence |
| 8 | Start on the database sequence |
| 9 | 	End on the database sequence |
| 10 | 	E value - the expectation that this alignment is random given the length of the sequence and length of the database |
| 11 | bit score - the score of the alignment itself |

Use ```grep``` to look at the matches for sequence ```CB1APANXX170426:1:1101:10000:16732/1``` in ```mmseqs-CSM7KOMH-s1.m8```.\
**Question 1**: How many matches are there?\
**Question 2**: What is the highest percent identity match?\
**Question 3**: What is the lowest E-value match (smallest value = best match)?\

## 3.3 Get MMSeqs top hits

The next step we need to take is to get the name of the protein sequence that had the best alignment for each sequence read in our samples. We can achieve this by running the command:
```
mkdir mmseqs_U90_out_tophit
python3 scripts/pick_uniref_top_hit.py --unirefm8Dir mmseqs_m8_files --output_path mmseqs_U90_out_tophit
```

Now that we have the best protein sequence that matches best with each of the sequences in our samples we can begin creating our final data table containing the stratified abundance of each function in our samples. Because we (in the Langille lab) are continually updating some of these scripts, and the ones that we are using are ones that we have only recently developed, there are a few different steps that we'll take to get the files in the format that we want them to be in. Often as we develop things in bioinformatics, we'll try out a lot of different things, and as the protocols that we use things mature we can consolidate them into a single script. We have just about reached that point here, but we haven't consolidated things yet so we'll run these few steps.

### Create a file that maps from our Kraken output to the MMseqs output

While we don't have that many files and could just make a file by hand, it's good practice to make these files with scripts for when we're working on our own data and may have hundreds of samples!

We'll be using Python for this, and we can open up Python (version 3) by typing in ```python3``` and pressing enter. You should see some information about the Python version print out, and a new command prompt with ```>>>``` pop up. This shows us that we're now using Python rather than bash. 

Next, we'll import the package that we'll use:
```
import os
```
**Note**: You'll need to press enter after each line to make sure that it runs!

Now, we'll make a list of our samples:
```
samples = os.listdir('cat_reads/') #this command creates a list of the files in cat_reads/
samples = [s.split('.')[0] for s in samples] #and this command uses a for loop to get only the sample names - the part of the file names before the '.'
```

Print out the list to see what's in it:
```
print(samples)
```

**Question 4**: What's in the list?

Now we'll set up a few variables so that we don't need to keep typing them out:
```
kraken_path = '~/CourseData/MIC_data/AMB_data/kraken2_outraw_RefSeqCompleteV205/'
kraken_suffix = '_0.1_RefSeqCompleteV205.kraken.txt'
mmseqs_path = 'mmseqs_U90_out_tophit/mmseqs-'
mmseqs_suffix = '-s1.m8-parsed.txt'
mmseqs_m8_path = 'mmseqs_m8_files-'
mmseqs_m8_suffix = '-s1.m8'
```

And then we'll create a new file (```multi-sample-outfiles-w-m8.txt```) and loop through the samples adding them to the new file:
```
with open('multi-sample-outfiles-w-m8.txt', 'w') as f:
  for sample in samples:
    string = sample+'\t'
    string += kraken_path+sample+kraken_suffix+'\t'+'kraken2'+'\t'
    string += mmseqs_path+sample+mmseqs_suffix+'\t'+'uniref'+'\t'
    string += mmseqs_m8_path+sample+mmseqs_m8_suffix+'\n'
    s = f.write(string)
```
Make sure that you press enter to execute this command.
You'll see that for each sample, we're adding the sample name, the path to the kraken output for that file, the path to the mmseqs output for that file, and the path to the mmseqs m8 output for that file, each separated by a tab (```\t```).

Now we can quit Python again by entering ```quit()``` (and pressing enter).

If you take a look at this file that we just made, you should see all of your samples in there along with the outputs from Kraken and MMseqs.

## 3.4 Combine Kraken taxonomy and MMSeqs function, normalise, and reformat

Now that we have this master file we can pass this information into the helper script to add all of it together for our samples:
```
ln -s MMSeqs2_db/*.pbz2 .
python3 scripts/parse_TaxonomyFunction.py --multisample multi-sample-outfiles-w-m8.txt --outputf Workshop_strat_matrix_RPKM.txt --stratified Y --map2EC Y
```

The first link command will grab all the databases that contain information about the length of each gene in our UniRef protein database. This will be important to normalize the abundance of each functional classification.

The second command will generate a final stratified data matrix that shows the abundance of each EC number stratified by the different taxonomic classifications within the sample. This script also normalizes the abundances of each of these ECS into reads per kilobase per million (RPKM). This abundance metric is the number of reads that mapped to that EC number per million reads mapped within the sample, divided by the gene length of that EC. It's important that the abundances are normalized by gene length or there would be an unfair bias toward longer genes due to the higher chance of them being sequenced. We can also run the same command as above without the ```--stratified Y``` option to generate functional profiles that are not broken down by the contributing taxa.

This might take a while, or it might say "killed". Either way, we can copy the output if we need to, and take a look at it:
```
ln -s ~/CourseData/MIC_data/AMB_data/Workshop_strat_matrix_RPKM.txt .
less Workshop_strat_matrix_RPKM.txt
```

**Question 5**: Describe what's in the ```Workshop_strat_matrix_RPKM.txt``` file. 

### Reformat output

Now we're going to be running a few different scripts that will get the output in the format that we're going to need.

First, we'll activate the PICRUSt2 conda environment:
```
conda activate picrust2
```

Then we'll open up Python again and use this to split the first column that currently contains EC numbers and taxon information into two separate columns:
```
import pandas as pd #import the pandas package, and instead of needing to type pandas each time, we'll just call it pd
file = pd.read_csv('Workshop_strat_matrix_RPKM.txt', header=0, sep='\t') #open up the file as a dataframe
samples = [s for s in file.columns if s != 'function'] #get all of the sample names
file[['EC','taxon']] = file['function'].str.split('|',expand=True) #split the first column named 'function' based on the '|' symbol
file = file.drop('function', axis=1).loc[:, ['EC', 'taxon']+samples] #remove the original 'function' column and reorder the others to have 'EC' and 'taxon' before the sample names
file.to_csv('Workshop_strat_matrix_RPKM_split.txt', index=False, sep='\t') #save this as a new file
```
Make sure to quit Python again.

Now, we can use one of the PICRUSt2 scripts to add descriptions for the EC numbers to this file:
```
add_descriptions.py -i Workshop_strat_matrix_RPKM_split.txt -o Workshop_strat_matrix_RPKM_split_des.txt -m EC
```

**Question 6**: Look at this new file. Any ideas how you might check that these descriptions are correct?

Currently, this file just contains the taxon names, but we're going to need the full taxonomy information for all of these taxa. The database that we've used for classification of the reads with Kraken is from NCBI, so this is where we'll get the most up to date information on these taxa. We'll now be making a data dictionary - a mapping from the NCBI taxonomy ID's to the full taxonomy information. It can be useful to learn how to manipulate information like this from the NCBI database, but if you don't want to do this, you can just copy across this dictionary and skip ahead to the next part:
```
cp ~/CourseData/MIC_data/all_output/amb_module3/taxonomy.dict .
```

First, we'll make a new folder for the NCBI taxonomy information and change into it:
```
mkdir taxonomy
cd taxonomy
```

Now we'll download the taxonomy information and extract the information frm the tar file:
```
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
```
There are a couple of key files that we're going to be using from this: ```nodes.dmp``` and ```names.dmp``` - you can see information on everything that's in these files in ```readme.txt```, but what we're interested in is the taxid (taxonomy ID) and parent taxid from ```nodes.dmp```, and the taxid and name in the ```names.dmp``` file. What we'll be doing is going through each of these files and essentially formatting the data to be how we want to have it for each of these taxa, so that it'll be quick and easy for us to access and we won't need to open up these large files every time that we want to access it. 

Change out of the ```taxonomy``` directory: ```cd ..```, start up ```python``` again, and run this code:
```
import pickle

nodes_dict = {}
names_dict = {}
nodes_level_dict = {}
```
First we are importing the pickle package - this allows us to easily save python objects like dictionaries - and then we're setting up some empty dictionaries. In some programming languages, different types of brackets will be used to show different types of objects. For example, in python, regular brackets ```()``` will be used for giving arguments to functions, square brackets ```[]``` are used for lists or accessing information contained in lists, and curly brackets ```{}``` are used to create dictionaries.

Now, we'll get the information that we need from the ```nodes.dmp``` file:
```
for line in open('taxonomy/nodes.dmp', 'r'):  
  this_line = line.replace('\t', '').split('|') 
  node, parent, level = this_line[0], this_line[1], this_line[2] 
  nodes_dict[node] = parent
  nodes_level_dict[node] = level
```
Here, we were saying that for each line of the ```nodes.dmp``` file:
- remove the tabs and split the remaining string by the '|' symbol
- define the node, parent and level of this as the first, second and third field (note that python indexing starts at 0, not 1 like R and some other programming languages)
- add the parent taxid to the name of this taxid in the nodes_dict
- add which taxonomy rank/level this is at to the nodes_level_dict

Now we'll do similar for the ```names.dmp``` file:
```
for line in open('taxonomy/names.dmp', 'r'):
  this_line = line.replace('\t', '').split('|')
  taxid, name, name_type = this_line[0], this_line[1], this_line[3]
  if name_type == 'scientific name':
    names_dict[taxid] = name
```
So again we're saying that for each line of the ```names.dmp``` file:
- remove the tabs and split the remaining string by the '|' symbol
- define the taxid, name and type of name as the first, second and fourth fields
- if the name type is the scientific name (the two == signs means that we're asking it if it matches, rather than defining the name_type as 'scientific name'). Other options here would include synonym, common name, equivalent name, type material... but we only want the scientific name for each taxonomy ID
- if it is the scientific name, then we'll save the name to the taxid in the names_dict

Now we'll reformat these dictionaries to a new dictionary full_taxonomy to be usable for us and have everything in the same place:
```
full_taxonomy = {}
for taxid in nodes_dict:
  parent_id = nodes_dict[taxid]
  tax_list = [names_dict[taxid]]
  nodes_list = [taxid]
  level_list = [nodes_level_dict[taxid]]
  while parent_id in nodes_dict:
    if parent_id in ['1', '0']:
      break
    tax_list.append(names_dict[parent_id])
    nodes_list.append(parent_id)
    level_list.append(nodes_level_dict[parent_id])
    parent_id = nodes_dict[parent_id]
  full_taxonomy[taxid] = [nodes_list, tax_list, level_list]
```
So here we're looping through all of the taxid's in the nodes_dict and we're:
- getting the parent taxid (and calling it parent_id)
- starting a list called tax_list where the first thing in the list is the scientific name that we saved in names_dict - we'll be adding to this with the names of the parent taxid's
- starting a list called nodes_list where the first thing is this taxid - we'll be adding the parent taxid's
- starting a list called level_list where the first thing is the name is this level - we'll be adding the names of all levels to this. These would be e.g. genus, species, etc.
- starting a while loop - while loops will continue indefinitely until a certain condition is met. You always need to be a little careful with them, because it's easy to code something that would run forever, but here we're just saying that we want to continue running something while the parent taxid (parent_id) is still in the nodes_dict
- next, we're making sure the loop will stop by saying that if the parent taxid is either '1' or '0' (which means either 'root' or 'unclassified', respectively), then the loop has to break
- if the parent taxid is still in the nodes_dict, then we'll add the name of the parent to the names_dict
- we'll add the parent taxid to the nodes list
- we'll add the name of the level of the parent to the level_list
- and finally, we'll change the parent taxid to be the parent of the current parent (i.e., we'll move up the list)
- after we've broken out of the while loop, we'll save the new lists that we've made to the full_taxonomy dictionary

Now we're ready to save the new dictionary that we've made:
```
with open('taxonomy.dict', 'wb') as f:
  pickle.dump(full_taxonomy, f)
```

And we can quit python:
```
quit()
```

### Add taxonomy to our file and save it in full stratified format

Next, we're going to add the full taxonomy information that we just made (or copied across), and we'll resave the file in the full stratified format that JarrVis needs.

First, start up python again: ```python```

Now import the packages that we'll need and open up the file:
```
import pandas as pd
import pickle

file = pd.read_csv('Workshop_strat_matrix_RPKM_split_des.txt', header=0, sep='\t')
```

Open up the taxonomy dictionary that we made:
```
with open('taxonomy.dict', 'rb') as f:
  full_taxonomy = pickle.load(f)
```

Make a new dictionary that contains the prefixes for each taxonomy level:
```
tax_levels = {'phylum':'p__', 'class':'c__', 'order':'o__', 'family':'f__', 'genus':'g__', 'species':'s__'}
```

Now we'll go through each row of this file and replace the existing taxon field with the full taxonomy information:
```
for row in file.index.values:
  tax = file.loc[row, 'taxon']
  taxid = tax.split('taxid ')[1].replace(')', '')
  full_tax = full_taxonomy[taxid]
  this_tax = []
  for level in tax_levels:
    if level in full_tax[2]:
      index = full_tax[2].index(level)
      name = full_tax[1][index]
      this_tax.append(tax_levels[level]+name)
    else:
      this_tax.append(tax_levels[level]+'Unclassified')
  this_tax = ';'.join(this_tax)
  file.loc[row, 'taxon'] = this_tax
```
In this for loop, we are:
- getting the name of the taxon from the dataframe (when we're accessing information from a dataframe like this, we always give the row name first and then the column name)
- getting the taxid of this taxon by splitting the string based on the 'taxid ' (and then removing the additional bracket)
- getting all of the taxonomy information for this taxid from our dictionary
- creating an empty list called this_tax
- going through each of the levels in the tax_levels dictionary
- checking whether the level is in our full taxonomy, and if it is:
- getting the index (location) in the list that it is)
- getting the name at this level
- adding the prefix and the name to the this_tax list
- else (if it's not in the full taxonomy)
- just add that this taxon is 'Unclassified' at this level (along with the prefix)
- create a string from the this_tax list that is separated by semi-colons
- rename the dataframe cell with the full taxonomy information

Now we can save this file:
```
file.to_csv('Workshop_strat_matrix_RPKM_split_des_full_tax.txt', sep='\t', index=False)
```
If you like, to help understand what the previous lines of code did, take a look at the previous file (```Workshop_strat_matrix_RPKM_split_des.txt```) and this new one that we've just made (```Workshop_strat_matrix_RPKM_split_des_full_tax.txt```).

**Question 7**: Can you work out what was added?

Now we're going to change the format of this file so that rather than have all samples on the same line for each function/taxon, each function/taxon/sample will have it's own line:
```
new_file = []
samples = [s for s in file.columns if s not in ['EC', 'description', 'taxon']]
for row in file.index.values:
  for sample in samples:
    this_row = [sample, file.loc[row, 'taxon'], file.loc[row, 'EC']+': '+file.loc[row, 'description'], file.loc[row, sample]]
    new_file.append(this_row)
```
Here, we first made an empty list called new_file, and got the names of the samples.
Then, we looped through each row of the dataframe and:
- for each sample:
- we made a new row/list (this_row) that contained the sample name, the name of the taxon, the name of the EC number, the description of the EC number, and the abundance value of that function/taxon/sample
- added this new row/list (this_row) to the new_file list

You should be able to tell that here we're going to have many more rows than we did previously as we're going to be multiplying the number of rows that we did have by the number of samples. An easy way to check that everything is as we hope it will be is to look at both the previous dataframe and the length of the new list:
```
print(file)
print(len(new_file))
```

Now we can convert this new list (new_file) to a dataframe called jarrvis_file (with the column names that JarrVis is expecting) and save it:
```
jarrvis_file = pd.DataFrame(new_file, columns=['Sample', 'Genus', 'Gene', 'Contribution'])
jarrvis_file.to_csv('Workshop_stratified_for_jarrvis.txt', sep='\t', index=False)
```

And then finally, we can quit python:
```
quit()
```

## 3.5 Visualise with JarrVis

Now we're finally ready to run JarrVis to visualise our taxonomy-function information from MMseqs!

1. First off, download both ```Workshop_stratified_for_jarrvis.txt``` and ```mgs_metadata.txt``` (you should be able to right-click them and click 'Save as') to your laptop. 
2. Open up and modify the ```mgs_metadata.txt``` so that 'sample' says 'sample_id' (it's important that this matches exactly!)
3. If you already have RStudio (desktop) installed, go ahead and open it up. If not, download and install it from [here](https://posit.co/download/rstudio-desktop/), and then open it up. 
4. Open the console and type: ```install.packages('shiny')``` - it may ask you more questions, but you probably want to say yes/y to everything.
5. Type: ```library(shiny)```
6. Type: ```runGist("943ff5fdbd94815cc27f302d9f56ff0b")``` - you should see some things running and then a box will pop up!
7. In "Upload stratified output File (TSV)", click browse and choose the ```Workshop_stratified_for_jarrvis.txt``` file
8. In "Upload Sample Metadata File (TSV)", click browse and choose the ```mgs_metadata.txt``` file
9. Click "Select Metadata Categories" - you should see the ```Metadata Categories``` box populate itself with "disease_state" (if you had more metadata categories in your metadata file, these would come up here for you to choose from)
10. Choose the taxonomy level to collapse at - you can choose whatever you'd like, but I went with genera
11. Click "Select Taxa Categories" - you won't see anything change, but the next "Filter by Taxa" box now has all of the taxa. Click on "Select All" here.
12. Click "Select Function Categories" - this is the same as for taxa, so now click on "Filter by Function" and "Select All" again.
13. Now click on "Update the Gene Contribution Threshold from data". You can move up the bottom filter so that we don't plot all of the really low abundance stuff - I've found that to ~500 seems to be sensible for visualising everything.
14. Now scroll to the bottom of the window and click "Display plot".
15. Now scroll back up and you should see the functions in each group and be able to start exploring them. Try clicking on some of the links and seeing how they go back to the sample groupings.

**Note**: Sometimes weird things happen! If you get an error when you try to display the plot, try quitting RStudio (if it asks you if you want to save your workspace click never), re-opening it, and starting again from step 5 here. 

If you're struggling with getting this to run, you can download the html output from JarrVis [here](https://drive.usercontent.google.com/download?id=1yVOCmU6iE6MOoP3F6t6ZxgRmVr_j8UrW) (this is equivalent to saving rather than viewing the files in JarrVis).

**Question 8**: What can you see in the JarrVis output?

## 3.6 Run HUMAnN

Now we can run HUMAnN. Note that you should still run step 3.1 to set up the directories and copy across the data that we'll be using.

Hopefully you got to the part of Module 1 where MetaPhlAn was run.

If you didn't:
```
conda activate taxonomic
mkdir ~/workspace/amb_module1/metaphlan_out
mkdir ~/workspace/amb_module3/humann_out
cd ~/workspace/amb_module3
parallel -j 2 'metaphlan --input_type fastq {} --bowtie2db ~/CourseData/MIC_data/tools/CHOCOPhlAn_201901/ > ~/workspace/amb_module1/metaphlan_out/{/.}.mpa' ::: cat_reads/*.fastq
```

If you did, first change directory and make the HUMAnN output folder:
```
cd ~/workspace/amb_module3
mkdir humann_out
```

Now activate the environment and run HUMAnN:
```
conda activate biobakery3
parallel -j 1 --eta --progress 'humann -i {} -o humann_out/ --threads 4 --taxonomic-profile ~/workspace/amb_module1/metaphlan_out/{/.}.mpa --protein-database ~/CourseData/MIC_data/tools/humann_databases/uniref' ::: cat_reads/*.fastq
```
**Note**: This might take quite a while to run so you may want to use ```tmux``` for this. Check back in the previous module if you've forgotten how to open this. Remember to activate the ```biobakery3``` environment if you do start a new tmux session!

The options here are hopefully quite straightforward:
- ```-i``` - the input files
- ```-o``` - the output folder to use
- ```--threads``` - the number of threads to run with
- ```--taxonomic-profile``` - the location of the MetaPhlAn taxonomic profiles. Note that this is optional - if we didn't include this option, MetaPhlAn would be run again (but it would run again using MetaPhlAn4, which apparently works better but requires more memory to run than we have here)
- ```--protein-database``` - the location of the uniref database to use. This is also optional, but if we don't provide it and HUMAnN doesn't find it in the folder that it expects it to be, it will download this database again, which could take quite a long time to do. 

This command could take a long time to run!! If you just want to copy across the output so that you can continue with the other steps, you can:
```
cp -r ~/CourseData/MIC_data/all_output/amb_module3/humann_out/ .
```

## 3.7 Combine HUMANn output

As this finishes with each of the samples, you can take a look at the output files if you like. The key files will be ```*_genefamilies.tsv```, ```*_pathabundance.tsv``` and ```*_pathcoverage.tsv``` for each sample, and you can look in the ```*_tmp``` folders, too. 

After they've all finished, you can run some HUMAnN scripts to join the tables for all of the samples together into one. We run this separately for each file type:
```
humann_join_tables -i humann_out -o HMP_humann_genefamilies.tsv --file_name genefamilies
humann_join_tables -i humann_out -o HMP_humann_pathabundance.tsv --file_name pathabundance
humann_join_tables -i humann_out -o HMP_humann_pathcoverage.tsv --file_name pathcoverage
```

**Question 9**: Have a look at these files. Do you see any issues with them?

I've also run HUMAnN on the full samples (not sub-sampled), and you can download that like this:
```
wget -O HMP_humann_pathabundance_full.tsv https://drive.usercontent.google.com/download?id=1KVMz0JWC4Q3F70QFJ_ZoGMyJqlvwMKOP
```

As well as downloading a version that I modified to have the metadata as an additional line just under the header here:
```
wget -O HMP_humann_pathabundance_full_metadata.tsv https://drive.usercontent.google.com/download?id=1uxHb3-xLmyOD7Qgdtk-RMWmhW8OcffGG
```

**Question 10**: Does this look more reasonable to you?

## 3.8 Visualise HUMANn output

HUMAnN has some built-in functions for plotting this output, so we'll run some of the pathways like so:
```
mkdir humann_plots
parallel -j 1 'humann_barplot --input HMP_humann_pathabundance_full_metadata.tsv --focal-metadata disease_state --last-metadata disease_state --output humann_plots/{}.png --focal-feature {} --sort sum metadata --scaling logstack' ::: PWY-7221 PEPTIDOGLYCANSYN-PWY PWY-1269 PWY-6277 FUCCAT-PWY P164-PWY PWY-7013 HSERMETANA-PWY P163-PWY
```
**Note**: these pathways were chosen more-or-less at random to show a few different patterns. Feel free to look through the file and run this with some other pathways, *e.g.,*:
```
humann_barplot --input HMP_humann_pathabundance_full_metadata.tsv --focal-metadata disease_state --last-metadata disease_state --output humann_plots/CALVIN-PWY.png --focal-feature CALVIN-PWY --sort sum metadata --scaling logstack
humann_barplot --input HMP_humann_pathabundance_full_metadata.tsv --focal-metadata disease_state --last-metadata disease_state --output humann_plots/ANAEROFRUCAT-PWY.png --focal-feature ANAEROFRUCAT-PWY --sort sum metadata --scaling logstack
```

Now take a look through these plots so that you can see how the functions are distributed among the taxa, and how this varies for different functions.

**Question 11**: Is this the same as you thought it would be based on the MMseqs output?

## 3.9 AMR annotation of reads using CARD RGI

We can use the Comprehensive Antibiotic Resistance Database (CARD) Resistance Gene Identifier (RGI) on reads as well as the contigs from MAGs. If you made it to this part in Module 2, then you already have the data that we need. If not, go to Section 2.10 of Module 2 and follow the instructions for downloading the databases.

If you already have the data, just activate the conda environment again:
```
conda activate rgi
```

Make the output directory and then change to the ```card_data``` directory to run the RGI from:
```
mkdir card_out
cd ~/workspace/card_data/
```

Now run CARD RGI:
```
parallel -j 1 --link 'rgi bwt --read_one {1} --read_two {2} --aligner bowtie2 --output_file ~/workspace/amb_module3/card_out/{1/.}_CARD --threads 2 --local --clean' \
 ::: ~/workspace/amb_module3/raw_data/*_R1_subsampled.fastq.gz ::: ~/workspace/amb_module3/raw_data/*_R2_subsampled.fastq.gz
```
You'll notice that this tool wants the information about the forward and reverse reads separately, we've told it to use the Bowtie2 aligned, where to put the output files, to use the local database, and to remove temporary files. 

Now we can change to the directory with the CARD output files and have a look at some of them:
```
cd ~/workspace/amb_module3/card_out
less HSMA33J3_R1_subsampled.fastq_CARD.allele_mapping_data.txt
less HSM6XRQY_R1_subsampled.fastq_CARD.overall_mapping_stats.txt
less HSM7J4QT_R1_subsampled.fastq_CARD.reference_mapping_stats.txt
less HSM7J4QT_R1_subsampled.fastq_CARD.gene_mapping_data.txt
```
You can see some details about the output files [here](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst#rgi-bwt-tab-delimited-output-details).

For now, we'll just separate the ```*gene_mapping_data*``` files from the others as these are what we're most interested in:
```
mkdir gene_mapping
mv *gene_mapping_data* gene_mapping/
mkdir other_read_files
mv *subsampled* other_read_files
cd ..
```

If you open up one of these files in Excel (or similar), you can have a look at all of the information that is in here. You'll see that there is a lot of information about the quality of the matches to the CARD database, including the number of completely mapped reads in the sample, the mapped reads with the flanking sequence (sequence found surrounding the gene sequence in the database), the average percent of the gene that is covered (along with the length of the reference sequence - based on these and our 100bp reads, it's not very surprising that these percentages are relatively low!), etc. We're often most interested in the number of reads that are mapped, but we might also want to perform some filtering based on some of these other parameters, and we can find this information on a read-by-read basis in some of the other files. There's also information in this file about the AMR gene families that each gene belongs to, as well as the drug class that they give resistance to and the mechanism of resistance.

How exactly we want to filter these results would depend on the application, but you can see that in a [recent study that I contributed to](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01662-3), where we were looking for AMR genes in the biofilms growing on plastics that were in a river next to a wastewater treatment plant, we chose to: remove AMR genes with <10 reads per sample, normalise the AMR genes to the length of the reference AMR gene within the CARD database, normalise the number of reads within each sample to give reads per kilobase per million (RPKM), and then we removed the target point mutations for antibiotic resistance (*e.g.* gyrase and ribosomal mutations) due to the elevated background noise that can occur when including point mutation ARGs from metagenomic data.

For this tutorial, we won't worry too much about filtering the data because these are subsampled sequencing data files anyway, but we will combine the files together and look at the number of reads giving resistance to different drug classes.

Start up python again: ```python```

Import the pandas and os packages:
```
import pandas as pd
import os
import matplotlib.pyplot as plt
```
You might notice that we've used the ```os``` package a few times now - it contains a lot of basic functions for working with directories, listing the contents, making new ones, etc.

Set up an empty list for all of the results to go to, a dictionary for the drug class that each AMR gene gives resistance to, and get the samples within the ```gene_mapping``` folder:
```
all_df = []
gene_drug_class = {}
samples = os.listdir('card_out/gene_mapping')
```

Now, we'll loop through the samples, opening the files, getting only the gene abundance data, and adding the resulting dataframes to the all_df list:
```
for sample in samples:
  this_sample = pd.read_csv('card_out/gene_mapping/'+sample, index_col=0, header=0, sep='\t')
  sample_name = sample.split('_')[0]
  for row in this_sample.index.values:
    if row not in gene_drug_class:
      gene_drug_class[row] = this_sample.loc[row, 'Drug Class']
  this_sample = this_sample.loc[:, ['All Mapped Reads']]
  all_df.append(this_sample.rename(columns={'All Mapped Reads':sample_name}))
```

Now we'll combine all of the sample data frames (in the all_df list) into a single data frame:
```
amr_df = pd.concat(all_df, axis=1).fillna(value=0)
```

And then we'll rename these AMR genes based on the drug class that they confer resistance to:
```
amr_df = amr_df.rename(index=gene_drug_class)
```
Remember that you can always print things out by either just typing the name or them, or using the ```print()``` function: ```print(amr_df)```

And now, as we have a lot of duplicate rows, we'll group them all together:
```
amr_df = amr_df.groupby(by=amr_df.index, axis=0).sum()
```
Note that we're taking the sum of the reads here, but in many cases there might be reasons to use the mean/median/etc.

Usually we might want to normalise these counts by the number of reads in each sample. This isn't really so necessary when we started with the same number of reads in each sample (100,000 in each sample), but we can still do it anyway to convert to relative abundance (%):
```
amr_df = amr_df.divide(100000).multiply(100)
```

Now we'll just plot this as a heatmap:
```
fig = plt.figure(figsize=(20,20))
ax = plt.subplot(111)
ax.pcolor(amr_df, cmap='plasma', edgecolor='black')
```
I've chosen to color this using the plasma colormap (cmap), but you could choose any of the others [here](https://matplotlib.org/stable/users/explain/colors/colormaps.html), or there are ways to define your own.

Then we can add some other things to this, like x and y tick labels, and numbers into the cells of the heatmap, and then save the figure:
```
xt = plt.xticks([a+0.5 for a in range(len(amr_df.columns))], amr_df.columns.values, rotation=90)
yt = plt.yticks([a+0.5 for a in range(len(amr_df.index))], amr_df.index.values)
max_val = max(amr_df.max(axis=1).values)

for r in range(len(amr_df.index.values)):
  for c in range(len(amr_df.columns)):
    val = amr_df.iloc[r, c]
    if val < (max_val/2): fc = 'white'
    else: fc = 'black'
    tx = ax.text(c+0.5, r+0.5, str(round(val, 4)), ha='center', va='center', fontsize=10, color=fc)

plt.savefig('AMR_heatmap.png', dpi=600, bbox_inches='tight')
```

In the for loop here, we're getting the number that's in each cell of the data frame (looping through the rows and then the columns), then if the value is less than half of the maximum % relative abundance that we have, we're making the text colour (fc) white, and if it's not, then it'll be black. Then we're adding text to the right cell of the heatmap.

If you take a look at the heatmap, you'll see that just a few drug classes seem to be abundant, with lots of the others having very low percentage abundances. 

**Question 12**: Are there any patterns that you can notice about the samples that seem to have higher relative abundances of AMR genes than others?

## Answers

Use ```grep``` to look at the matches for sequence ```CB1APANXX170426:1:1101:10000:16732/1``` in ```mmseqs-CSM7KOMH-s1.m8```.\
**Question 1**: How many matches are there?\
Using ```grep -c "CB1APANXX170426:1:1101:10000:16732/1" mmseqs_m8_files/mmseqs-CSM7KOMH-s1.m8``` we can see that there are 25 matches for that sequence.

**Question 2**: What is the highest percent identity match?\
Using ```grep "CB1APANXX170426:1:1101:10000:16732/1" mmseqs_m8_files/mmseqs-CSM7KOMH-s1.m8``` to print the matches, we can look at the third column to see the % identity matches. The first two matches (UniRef90_A0A016DTE9 and UniRef90_J9H7P3) are both 1.00 (or 100%). 

**Question 3**: What is the lowest E-value match (smallest value = best match)?\
Using ```grep``` as above for question 2, and looking at the E-values in the second to last column, we can see that the lowest values are again for the top two sequences (5.517E-15). The closer to 0 the E-values are, the better. 

**Question 4**: What's in the list?\
Your output should look like this: ```['MSM9VZHR', 'HSMA33KE', 'HSM7J4QT', 'MSMB4LXW', 'MSM79HA3', 'CSM7KOMH', 'PSM7J18I', 'HSM6XRQY', 'HSMA33J3', 'CSM79HR8']```\
You can see that this is a list of sample names - you can tell that it's a list because it's enclosed in square brackets (```[]```). 

**Question 5**: Describe what's in the ```Workshop_strat_matrix_RPKM.txt``` file.\
We can see that this file contains a feature table containing samples as columns and stratified taxa/functions as rows. The cells each contain the abundance converted to RPKM (Reads Per Kilobase per Million). You can see that the functions (EC numbers) and taxa are separated by the "|" symbol. 

**Question 6**: Look at this new file. Any ideas how you might check that these descriptions are correct?\
It's always a good idea to check that scripts are doing what we expect of them! I'd usually check a few of the EC numbers/descriptions on the KEGG website, e.g. for [EC:6.1.1.7](https://www.genome.jp/dbget-bin/www_bget?ec:6.1.1.7) - I just found this page by googling the EC number. 

**Question 7**: Can you work out what was added?\
We added in the full ranked taxonomy here. As above, it's always a good idea to check some of these. We can use the [NCBI taxonomy website](https://www.ncbi.nlm.nih.gov/taxonomy) if we want to have a look. Try searching the taxon name that we had before and looking at the full taxonomy information.

**Question 8**: What can you see in the JarrVis output?\
In the left side of JarrVis, we can see our two sample groups, CD (Crohn's Disease) and nonIBD. Then we can see links to taxa - g__Bacteroides (the Bacteroides genus) has the largest number of connections, and it has connections to both sample groups. On the right, we can see the functions with their descriptions, and we can see how these link back to the taxa and samples that contribute them. Hovering over each of the grey lines shows us the RPKM are in this taxon/function grouping (note: it might take a second after hovering sometimes for the information to pop up!). 

**Question 9**: Have a look at these files. Do you see any issues with them?\
We don't have very much that's well classified. This is because we're using very sub-sampled samples (*vs* the millions of reads that is typical for metagenome samples), and for useful functional annotation we really need higher sequencing depth than for taxonomy.

**Question 10**: Does this look more reasonable to you?\
Yes! Now we have lots more functions present in our samples. 

**Question 11**: Is this the same as you thought it would be based on the MMseqs output?\
It's a bit difficult to tell now that we're looking at pathways and not EC numbers! If we wanted to find out, we could either look at the same functions for each, or we could look at the pathways that the EC numbers contribute to (or vice versa).

**Question 12**: Are there any patterns that you can notice about the samples that seem to have higher relative abundances of AMR genes than others?\
We can see that HSM7J4QT seems to be higher in aminoglycoside antibiotic resistance than other samples, HSM6XRQY and MSMB4LXW in tetracycline antibiotic resistance, and CSM79HR8 in macrolide/lincosamide (sometimes called MLS) antibiotic resistance. These are in the CD, nonIBD, nonIBD and CD groups, so based on this cursory look alone it doesn't look like there are massive differences in AMR gene abundance based on metadata grouping, but we haven't done any statistical testing or looked at the individual genes yet so can't really draw conclusions based on this!
