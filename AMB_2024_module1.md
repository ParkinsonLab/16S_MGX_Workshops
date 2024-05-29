---
layout: workshop_main_2day
permalink: /AMB_2024_module1
title: AMB 2024
header1: Workshop Pages for Students
header2: Advanced Microbiome Analysis 2024
image: /site_images/AMB_2024_v1.png
length: 2 days
---

# Module 1: Introduction to metagenomics and read‐based profiling

This tutorial is part of the 2024 Canadian Bioinformatics Workshops [Advanced Microbiome Analysis](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-workshop-2024-advanced) (St John's, NL, May 29-30).

Author: Ben Fisher


The goal of this tutorial is to familiarize students with processes by which we analyze metagenomic data. Shotgun Metagenomic Sequencing, sometimes called MGS or WGS, is capable of capturing any DNA extracted from a given sample (however, this does not necessarily mean it captures ALL of the DNA). With MGS reads, we must consider that there will be significant host contamination, so we must filter against host sequences in our pipeline. There are many tools that can work to produce similar results, and this is not a one-size-fits-all solution. This is instead a foray into popular tools and processes, and how to appropriately use them for our analyses.

Throughout this module, there are some questions aimed to help your understanding of some of the key concepts. You'll find the answers at the bottom of this page, but no one will be marking them.


## Part 1:
[Working on the Command line](#working-on-the-command-line)
* [A Crash Course in GNU Parallel](#a-crash-course-in-gnu-parallel)
* [Quality Control](#quality-control)
  * [Visualization with FastQC](#visualization-with-fastqc)
  * [Filtering with Kneaddata](#filtering-with-kneaddata)
* [Generating Taxonomic Profiles](#generating-taxonomic-profiles)
  * [Annotation with Kraken2/Bracken](#annotation-with-kraken2bracken)
  * [Annotation with MetaPhlAn](#annotation-with-metaphlan)
  * [Preparing our Files for Analysis](#preparing-our-files-for-analysis)
## Part 2:
[Working in RStudio](#working-in-rstudio)
- [Module 1: Introduction to metagenomics and read‐based profiling](#module-1-introduction-to-metagenomics-and-readbased-profiling)
  - [Part 1:](#part-1)
  - [Part 2:](#part-2)
- [Working on the Command Line](#working-on-the-command-line)
- [A Crash Course in GNU Parallel](#a-crash-course-in-gnu-parallel)
- [Quality Control](#quality-control)
  - [Visualization with FastQC](#visualization-with-fastqc)
  - [Filtering with KneadData](#filtering-with-kneaddata)
- [Generating Taxonomic Profiles](#generating-taxonomic-profiles)
  - [Annotation with Kraken2/Bracken](#annotation-with-kraken2bracken)
  - [Annotation with MetaPhlAn](#annotation-with-metaphlan)
  - [Preparing our files for Analysis](#preparing-our-files-for-analysis)
- [Working in RStudio](#working-in-rstudio)
- [Setup, Importing and Formatting Data](#setup-importing-and-formatting-data)
- [Remove Rare Taxa and Rarefy](#remove-rare-taxa-and-rarefy)
  - [Pruning](#pruning)
  - [Rarefy](#rarefy)
- [Alpha Diversity](#alpha-diversity)
- [Beta Diversity](#beta-diversity)
  - [Preparing the Data](#preparing-the-data)
  - [Plot the Data](#plot-the-data)
- [Visualization with Stacked Bar Charts](#visualization-with-stacked-bar-charts)
  - [Formatting the Data](#formatting-the-data)
  - [Selecting our Colours](#selecting-our-colours)
  - [Creating the Plot](#creating-the-plot)
- [Visualization with Heatmaps](#visualization-with-heatmaps)
  - [Formatting the data](#formatting-the-data-1)
  - [Annotate and Generate the Heatmap](#annotate-and-generate-the-heatmap)

# Working on the Command Line

Copy the files to your working directory.

```
cp -r ~/CourseData/MIC_data/AMB_data/raw_data/ .
```



# A Crash Course in GNU Parallel

First, activate our conda environment for this tutorial:

```
conda activate taxonomic
```

Sometimes in bioinformatics, the number of tasks you have to complete can get VERY large. Fortunately, there are several tools that can help us with this. One such tool is [GNU Parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html). This tool can simplify the way in which we approach large tasks, and as the name suggests, it can iterate though many tasks in _parallel_, i.e. concurrently. We can use a simple command to demonstrate this.

```
parallel 'echo {}' ::: a b c
```

With the command above, the program contained within the quotation marks `' '` is `echo`. This program is run 3 times, as there are 3 inputs listed after the `:::` characters. What happens if there are multiple lists of inputs? Try the following:

```
parallel 'echo {}' ::: a b c ::: 1 2 3
```

Here, we have demonstrated how `parallel` treats multiple inputs. It uses all combinations of one of each from `a b c` and `1 2 3`. But, what if we wanted to use 2 inputs that were sorted in a specific order? This is where the `--link` flag becomes particularly useful. Try the following:

```
parallel --link 'echo {}' ::: a b c ::: 1 2 3
```

In this case, the inputs are "linked", such that only one of each is used. If the lists are different lengths, `parallel` will go back to the beginning of the shortest list and continue to use it until the longest list is completed. You do not have to run the following command, as the output is provided to demonstrate this.

```
$ parallel --link 'echo {}' ::: light dark ::: red blue green
> light red
> dark blue
> light green
```
Notice how `light` appears a second time (on the third line of the output) to satisfy the length of the second list.

Another useful feature is specifying _which_ inputs we give parallel are to go _where_. This can be done intuitively by using multiple brackets `{ }` containing numbers corresponding to the list we are interested in. Again, you do not have to run the following command, as the output is provided to demonstrate this.

```
$ parallel --link 'echo {1} {3}; echo {2} {3}' ::: one red ::: two blue ::: fish
> one fish
> two fish
> red fish
> blue fish
```

Finally, a handy feature is that `parallel` accepts files as inputs. This is done slightly differently than before, as we need to use four colon characters `::::` instead of three. Parallel will then read each line of the file and treat its contents as a list. You can also mix this with the three-colon character lists `:::` you are already familiar with. Using the following code, create a test file and use `parallel` to run the `echo` program:

```
echo -e "A\nB\nC" > test.txt
parallel --link 'echo {2} {1}' :::: test.txt ::: 1 2 3
```

And with that, you're ready to use parallel for all of your bioinformatic needs! We will continue to use it throughout this tutorial and show some additional features along the way. There is also a cheat-sheet [here](https://www.gnu.org/software/parallel/parallel_cheat.pdf) for quick reference.

# Quality Control

Quality control is an important step in any pipeline. For this tutorial, we will use FastQC to inspect the quality of our samples.

## Visualization with FastQC

First, make your desired output directory (if it doesn't already exist). Then, run FastQC as follows:

```
fastqc -t 4 raw_data/*fastq.gz -o fastqc_out
```
Go to `http://##.uhn-hpc.ca/` (substituting ## for your student number) and navigate to your FastQC output directory. Click on the **.html** files to view the results for each sample. Let's look at the Per base sequence quality tab. This shows a boxplot representing the quality of the base call for each position of the read. In general, due to the inherent degradation of quality with increased sequence length, the quality scores will trend downward as the read gets longer. However, you may notice that for our samples this is not the case! This is because for the purpose of this tutorial, your raw data has already been trimmed.

More often, per base sequence quality will look like the following. The FastQC documentation provides examples of "[good](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html#M1)" and "[bad](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)" data. These examples are also shown below:

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/1ceec41e-f7ae-4763-889d-a24e95384d57" width="800">

**Question 1:** Which of the graphs does your data resemble more closely?\
**Question 2:** What can we do if data fails the Per Base Sequence Quality module?

Now, you may have also noticed that most of the samples fail the "Per Base Sequence Content" module of FastQC. Let's look at our visualization:

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/10887a71-fe6d-495a-9313-f4e11101d6ce" width="800">

This specific module plots out the proportion of each base position in a file, and raises a warning/error if there are large discrepancies in base proportions. In a given sequence, the lines for each base should run in parallel, indicating that the base calling represents proper nucleotide pairing. Additionally, the A and T lines may appear separate from the G and C lines, which is a consequence of the GC content of the sample. The relative amount of each base reflects the overall amount of the bases in the genome, and should not be imbalanced. When this is not the case, the sequence composition is **biased**. A common cause of this is the use of primers, which throws off our sequence content at the beginning of the read. Fortunately, although the module error stems from this bias, according to the [FastQC documentation for this module](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html#:~:text=Whilst%20this%20is%20a%20true%20technical%20bias%2C%20it%20isn%27t%20something%20which%20can%20be%20corrected%20by%20trimming%20and%20in%20most%20cases%20doesn%27t%20seem%20to%20adversely%20affect%20the%20downstream%20analysis.) it will not largely impact our downstream analysis.

The other modules are explained in depth in the [FastQC Module Help Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)


## Filtering with KneadData

KneadData is a tool which "wraps" several programs to create a pipeline that can be executed with one command. Remember that these reads have already been trimmed for you - this is in fact one of KneadData's functionalities. For this tutorial though, we will use KneadData to filter our reads for contaminant sequences against a human database.

Kneaddata outputs many files for each sample we provide it. These include:
* paired sequences which match our database;
* singletons which match our database;
* paired sequences that **do not** match our database;
* singletons that **do not** match our database;
* and some others.

**If you have not already**, we want to activate our `conda` environment where KneadData and our other tools for this tutorial are installed:

```
conda activate taxonomic
```

Then, run KneadData, using parallel, with the following command:

```
parallel -j 1 --eta --link 'kneaddata -i1 {1} -i2 {2} -o kneaddata_out -db ~/CourseData/MIC_data/tools/GRCh38_PhiX --bypass-trim --remove-intermediate-output' ::: raw_data/*R1_subsampled.fastq.gz ::: raw_data/*R2_subsampled.fastq.gz
```

While `kneaddata` is running, consider the following:
The GRCh38_PhiX database is made up of the human genome. 

**Question 3:** Of the four output files specified above, which should we choose for analyzing the microbiome?

You can check out all of the files that `kneaddata` has produced by listing the contents of the output directory (there is a lot!). Take note of how the files are differentiated from one another, and try to identify some of the files we are interested in. Once `kneaddata` is complete, we want to stitch our reads together into a single file. This is accomplished with a Perl script from our very own [Microbiome Helper](https://github.com/LangilleLab/microbiome_helper/blob/master/concat_paired_end.pl). For your convenience, it is already on your student instance.

with the Perl script, concatenate the reads into a single file:

```
perl ~/CourseData/MIC_data/AMB_data/scripts/concat_paired_end.pl -p 4 --no_R_match -o cat_reads kneaddata_out/*_paired_contam*.fastq
```
The script finds paired reads that match a given _regex_ and outputs the combined files.
* We first specify that our program is to be run with Perl, and then provide the path to the program.
* The `-p` flag specifies how many processes to run in parallel. The default is to do one process at a time, so using `-p 4` speeds things up.
* The `--no_R_match` option tells the script that our read pairs are differentiated by `*_1.fastq` instead of `*_R1.fastq`.
* The `-o` flag specifies the directory where we want the concatenated files to go.
* Our regex matches the paired reads that do not align to the human database from the KneadData output. This is because the reads that aren't "contaminants" actually align to the human genome, so what we are left with could contain microbial reads.
  * Consider that our files of interest are named something like `MSMB4LXW_R1_subsampled_kneaddata_GRCh38_PhiX_bowtie2_paired_contam_1.fastq`. If we want to match all of our paired contaminant files with a regex, we can specify the string unique to those filenames `_paired_contam`, and use wildcards `*` to fill the parts of the filename that will change between samples.

If the above does not work, you may need to install Perl:

```
conda install conda-forge::perl
```

If it still does not work or you already have Perl installed, you may get an error saying you require Parallel::ForkManager. Fix by executing the following inside your conda environment:

```
conda install bioconda::perl-parallel-forkmanager

```

# Generating Taxonomic Profiles

First, we should see how many reads are in our concatenated samples. Since `.fastq` files have 4 lines per read, we can divide the number of lines in the file by 4 to count the reads. Use the following command to check number of lines in the output files:
```
wc -l cat_reads/*
```

Woah! There's almost nothing left in most of these files! One even has zero reads! From this, we can infer that KneadData found the majority of our reads aligned to the human genome, leaving us with very few sequences to look for microbial reads. This is not entirely uncommon, however the fact that our input reads are actually subsets of _much_ larger samples exacerbated this effect.

For the purpose of this tutorial, we will instead continue with the **raw data** instead of our filtered data. **We are only doing this to demonstrate the tools for the purposes of this tutorial. This is NOT standard practice. In practice, you SHOULD NOT use unfiltered metagenomic data for taxonomic annotation.**

To continue, we will concatenate the raw data, then unzip it (the `;` lets you enter multiple command lines that will execute in series when you press enter).
```
perl ~/CourseData/MIC_data/AMB_data/scripts/concat_paired_end.pl -p 4 -o cat_reads_full raw_data/*.fastq.gz;
gunzip cat_reads_full/*.gz
```

Now let's check how many reads we are dealing with:
```
wc -l cat_reads_full/*
```
**Question 4:** How many reads are in each sample?


## Annotation with Kraken2/Bracken

Now that we have our reads of interest, we want to understand _what_ these reads are. To accomplish this, we use tools which **annotate** the reads based on different methods and databases. There are many tools which are capable of this, with varying degrees of speed and precision. For this tutorial, we will be using Kraken2 for fast exact k-mer matching against a database. 

Our lab has also investigated which parameters impact tool performance in [this Microbial Genomics paper](https://pubmed.ncbi.nlm.nih.gov/36867161/). One of the most important factors is the contents of the database, which should include as many taxa as possible to avoid the reads being assigned an incorrect taxonomic label. Generally, the bigger and more inclusive database, the better. However, due to the constraints of our cloud instance, we will be using a "Standard 8GB" index [provided by the Kraken2 developers](https://benlangmead.github.io/aws-indexes/k2). For your convenience, the database is already available on your instance.

First, you must create the appropriate output directories, or Kraken2 will not write any files. Using `parallel`, we will then run Kraken2 for our concatenated reads as shown below:

```
parallel -j 1 --eta 'kraken2 --db ~/CourseData/MIC_data/tools/k2_standard_08gb --output kraken2_outraw/{/.}.kraken --report kraken2_kreport/{/.}.kreport --confidence 0 {}' ::: cat_reads_full/*.fastq
```
This process can take some time. While this runs, let's learn about what our command is doing!
* We first specify our options for parallel, where:
  * the `-j 1` option specifies that we want to run two jobs concurrently;
  * the `--eta` option will count down the jobs are they are completed;
  * after the program contained in quotation marks, we specify our input files with `:::`, and use a regex to match all of the concatenated, unzipped `.fastq` files.
* We then describe how we want `kraken` to run:
  * by first specifying the location of the database with the `--db` option;
  * then specifying the output directory for the raw kraken annotated reads;
    * notice that we use a special form of the brackets here, `{/.}`, this is a special function of `parallel` that will remove both the file path and extension when substituting the input into our `kraken` command. This is useful when files are going into different directories, and when we want to change the extension.
  * similarly, we also specify the output of our "report" files with the `--report` option;
  * the `-confidence` option allows us to filter annotations below a certain threshold (a float between 0 and 1) into the unclassified node. We are using `0` because our samples are already subset, however this should generally be higher. See [our paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10132073/#:~:text=Comparison%20of%20key%20performance%20metrics%20at%20the%20read%20level) for more information.
  * and finally, we use the empty brackets `{}` for `parallel` to tell `kraken` what our desired input file is.


With Kraken2, we have annotated the reads in our sample with taxonomy information. If we want to use this to investigate diversity metrics, we need to find the abundances of taxa in our samples. This is done with Kraken2's companion tool, Bracken (Bayesian Reestimation of Abundance with KrakEN). 

Let's run Bracken on our Kraken2 outputs! First, make an output directory, then run the following:

```
parallel -j 2 --eta 'bracken -d ~/CourseData/MIC_data/tools/k2_standard_08gb -i {} -o bracken_out/{/.}.species.bracken -r 100 -l S -t 1' ::: kraken2_kreport/*.kreport
```
Some notes about this command:
* `-d` specifies the database we want to use. It should be the same database we used when we ran Kraken2;
* `-i` is our input file(s);
* `-o` is where and what we want the output files to be;
* `-r` is the read length to get all classifications for, the default is 100;
* `-l` is the taxonomic level at which we want to estimate abundances;
* `-t` is the number of reads required prior to abundance estimation to perform re-estimation

## Annotation with MetaPhlAn

Another tool that is commonly used for taxonomic annotation of metagenomic sequences is MetaPhlAn. This tool is different from Kraken2 in that it uses a database of marker genes, instead of a collection of genomes. We will use MetaPhlAn 3 for this tutorial, but a newer version (4) utilizing a larger database is available. First, let's make an output folder. Then, we can run the following:

```
parallel -j 2 --eta 'metaphlan --input_type fastq --bowtie2db ~/CourseData/MIC_data/tools/CHOCOPhlAn_201901 --no_map -o metaphlan_out/{/.}.mpa {}' ::: cat_reads/*.fastq
```


Similar to Kraken2, MetaPhlAn will output individual files for each sample. We can use a utility script from MetaPhlAn to merge our outputs into one table.

```
python ~/CourseData/MIC_data/AMB_data/scripts/merge_metaphlan_tables.py metaphlan_out/*.mpa > metaphlan_out/mpa_merged_table.txt
```

If we view this new combined table, we will see three key things: 
1. First, the output format is different to that of Kraken2, where the full taxonomic lineages are expanded line by line. 
2. Second, MetaPhlAn only outputs relative abundances for each taxonomic node, whereas Kraken2 (before re-analysis with Bracken) will output absolute numbers of reads assigned to each taxonomic node.
3. Third, the number of taxa that MetaPhlAn finds is _much_ smaller than Kraken2. This is partially due to us using a low confidence threshold with Kraken2, but this discrepancy between the two tools tends to hold true at higher confidence thresholds as well. [See our paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10132073/#:~:text=Classification%20of%20real%20metagenome%20datasets) for more info about how these tools perform compared to each other.


## Preparing our files for Analysis

After all of this, we are almost ready to create some profiles from our samples! The last step is to put everything into a format that R can easily handle. One standard format is the `.biom` format, which is a recognized standard for the [Earth Microbiome Project](http://www.earthmicrobiome.org/) and is supported by the [Genomics Standards Consortium](http://gensc.org/).

To transform our data, we will use **kraken-biom**, a tool which takes the `.kreport` files output by `bracken` and creates a new, combined file in the `.biom` format. This tool can also incorporate our sample metadata, and it is easier to merge this information now versus later.

First, we will have to copy the metadata file to our working directory:

```
cp ~/CourseData/MIC_data/AMB_data/amb_module1/mgs_metadata.tsv .
```

Let's have a look at this file, try reading it with `cat`

Now that we have our metadata, let's run `kraken-biom` and merge our samples to one organized output file:

```
kraken-biom kraken2_kreport/*bracken_species.kreport -m mgs_metadata.tsv -o mgs.biom --fmt json
```
* The inputs are a positional argument, and `kraken-biom` will accept lists. We can match all of our desired bracken-corrected `.kreport` files with the regex `kraken2_kreport/*bracken_species.kreport`
* the `-m` option is for specifying our metadata file. Note that `kraken-biom` is picky about the order of the metadata, so the entries in the list should be in the same order that your files are found. Typically this is in lexicographic order.
* the `-o` option specifies our output file.
* the `--fmt` option specifies that we want the output to be in `json` format, which is not the default behavior of the program. JSON is a text-based version of the BIOM format, as opposed to HDF5 which is a binary file format.

Now, with all of that, we should have our final `mgs.biom` file ready to go! You can check out what this file looks like (it's a single line so `head` will not work here). It can be cumbersome to look at, but there are patterns. Fortunately, R is much better at reading these files than we are!


# Working in RStudio

### Create the R Notebook
Using the menus, click `File > New File > R Notebook`, which will open an Untitled R markdown (Rmd) document. R notebooks are helpful in that you can run several lines, or **_chunks_** of code at the same time, and the results will appear within the document itself _(in the whitespace following the chunk)_.

The default R Notebook contains a header and some information you may find helpful. Try running the chunk containing `plot(cars)` to see what happens!

You do not need to preserve most of the information in the new, untitled document. Select all of its contents by click+dragging your cursor or entering the **`ctrl+a`** shortcut, and press **`backspace`** or **`delete`** to clear the document.

The **chunks** are distinguished by the grey shading. Everything between the first ` ```{r} ` and subsequent` ``` ` belongs to the chunk. Anything written in the white space surrounding the chunk is meant to be annotation. _Although you can run lines of code outside of the chunks, the chunks are useful for running multiple lines in series with one click._

#### Adding new chunks
To add a new chunk into your R notebook, either:
1. Navigate to `Code > Insert Chunk` from the toolbar, or,
2. Use the shortcut **`ctrl+alt+I`**


# Setup, Importing and Formatting Data

It would be best to create a new R markdown document for this section. You can then paste the following lines of code into the first chunk and clink "Run". If you are not using a markdown document, omit the first line.

```
knitr::opts_chunk$set(echo = TRUE)
library(phyloseq)
library(vegan)
library(ggplot2)
library(taxonomizr)

#setwd("~/workspace/")
data<-import_biom("~/CourseData/MIC_data/AMB_data/amb_module1/mgs.biom", header = TRUE)

```
We can see what the phyloseq object looks like by using:
```
View(data)
```
With this, we can see that it is a special object that should contain `otu_table`, `tax_table`, and `sam_data` objects within it. We can additionally inspect each of these elements with the `View()` function, specifying what we are looking for inside the phyloseq object. This command takes the form: `View(data@sam_data)`, where the string following the `@` sign is what we want to see. 
**Question 5:** What do each of these objects within `data` represent?

Next, we will do some housekeeping with the phyloseq object we have created with the `import-biom()` function. If we look at the `tax_table` with the following:
```
View(data@tax_table)
```
We get something like this:

| | Rank1 | Rank2 | Rank3 | Rank4 | Rank5 | Rank6 | Rank7 |
|---|---|---|---|---|---|---|---|
820 | k__Bacteria | p__Bacteroidota | c__Bacteroidia | o__Bacteroidales | f__Bacteroidaceae | g__Bacteroides | s__uniformis
2755405 | k__Bacteria | p__Bacteroidota | c__Bacteroidia | o__Bacteroidales | f__Bacteroidaceae | g__Bacteroides | s__sp. CACC 737
2528203 | k__Bacteria | p__Bacteroidota | c__Bacteroidia | o__Bacteroidales | f__Bacteroidaceae | g__Bacteroides | s__sp. A1C1
|...|...|...|...|...|...|...|...|

We can see that the table contains all of our information, but the labels need some work. So, we can modify the contents of the `tax_table` as follows. Let's create a new chunk in our document and start pasting this code to run:

```
#Rename the taxa names in the taxa table.
data@tax_table@.Data <- substring(data@tax_table@.Data, 4)
#Rename columns of taxa table.
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
bacteria<-subset_taxa(data, Kingdom == "Bacteria")
```
So, what's this doing to our data?
* The first command strips all the taxonomic names of their first 3 characters, which removes those pesky leading `k__` strings. We are directly modifying the vectors within the data frame with the `substr()` command. The parameters of this command tell R we only want the strings starting at the 4th character, which skips the letter and two underscores (i.e. `k__).
* The second command creates a character vector of the taxonomic levels to replace the "Ranks" in the table, and replaces those column names.
* The third command subsets all of the taxa in our data to only bacteria, by selecting only taxa with `"Bacteria"` as their `Kingdom` column entry in the `tax_table`. 

Now if you run `View(data@tax_table)`, the resulting output should look much cleaner!

And, a quick way to check if this worked is to execute the following command in your console. This looks for unique values in the "Kingdom" row of the taxonomy table, after coercing it to a data frame:
```
unique(as.data.frame(bacteria@tax_table@.Data)$Kingdom)
```
Which should only return 1 unique string:
```
[1] "Bacteria"
```


# Remove Rare Taxa and Rarefy

An important step in analyzing sequencing data is rarefaction. Rarefaction involves randomly subsetting samples so that all samples have **even sequencing depth**

<details>
  <summary><b>A quick note on rarefying and rarefaction: </b>(<i>click to expand</i>)</summary>
  <!-- have to be followed by an empty line! -->
<br>
<i>It is worth knowing that the practice of rarefaction has been <a href="https://pubmed.ncbi.nlm.nih.gov/24699258/">called into question</a> in the past. However, the current thinking is that "<a href="https://journals.asm.org/doi/10.1128/msphere.00355-23">rarefaction is the most robust approach to control for uneven sequencing effort when considered across a variety of alpha and beta diversity metrics</a>."</i><br>
<br>
</details>

Generally, it is a good idea to start by manually removing rare taxa. It is common to remove taxa that have less than 20 reads across all samples in our dataset. To do this, we will use the `prune_taxa()` command from phyloseq.

## Pruning

Create a new chunk. If we view the `otu_table` of `Bacteria`, we will see as we scroll through that there are many taxa that appear sparsely across the different samples.
```
View(bacteria@otu_table)
```
We are not particularly interested in these rare taxa, so a quick way to deal with them is to "prune" taxa from our samples that have less than "n" reads across all samples. We can first look at the taxa sums in the form of a histogram.

```
#View the histogram of taxa sums in our "Bacteria" dataset.
hist(taxa_sums(bacteria), breaks = 2000, xlim = c(0,1000), main = "Taxa Sums before Pruning")
```
With this, we see that most frequently, the sum of all reads for a given taxa is in the 0-20 bin. This means that there are lots of taxa with low abundances (less than 20 reads total) in our dataset. So, we can prune these rare taxa with a built-in phyloseq command called `prune_taxa`:

```
#Prune rare taxa from the dataset. This removes taxa that have less than 20 occurances across all samples.
bacteria<-prune_taxa(taxa_sums(bacteria)>=20,bacteria)
```
Some notes about this command:
* The first argument we give to `prune_taxa()` is the condition we want met for the taxa to 'pass' this filtering.
* We are looking for the sum of the taxa (found by `taxa_sums`) in our phyloseq object `bacteria` to be greater than 20.
* The second argument is the phyloseq object we are pruning, which will still be `bacteria`.
* The result of this command is overwrites our previous `bacteria` R object.

With that, we can re-visit our histogram and see what this pruning has done:

```
#View the histogram of taxa sums after we remove the rare taxa.
hist(taxa_sums(bacteria), breaks = 2000, xlim = c(0,1000), main = "Taxa Sums after Pruning")
```

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/76c9e1bd-2e97-47e6-8352-5abc2bd0a709" width="800">

From the scale of the y-axes, we can see that this has pruned many of the rare taxa. This can be verified by viewing the `otu_table` with `View(bacteria@otu_table)`.


## Rarefy
Create a new chunk. To rarefy our dataset, we must first visualize the _rarefaction curve_ of our samples using the **vegan** package. To do this, we need to create a dataframe that vegan can work with from our Phyloseq object `bacteria`.

```
#Visualize the rarefaction curve of our data.
rarecurve(as.data.frame(t(otu_table(bacteria))), step=50,cex=0.5,label=TRUE, ylim = c(1,150))
```
Some notes about this command:
* We apply a number of transformations to `bacteria` before `rarecurve` works on it:
  * the otu_table of `bacteria` is returned by the `otu_table(bacteria)` command;
  * the otu_table is then tranposed by `t()` such that the rows and columns are switched, because this is the format `rarecurve` expects;
  * this transposed otu_table is then coerced to a dataframe by `as.data.frame()` so that `rarecurve` can read it.
* The remainder of the `rarecurve` parameters control how the output is displayed.

Your rarefaction curve should look similar to the following:

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/f900f6f3-1d0a-44e2-8673-1920d728541c" width="800">

Looking at the rarefaction curve, we can see that the number of species for all of the samples eventually begins to plateau, which is a good sign! This tells us that we have reached a sequencing depth where more reads does not improve the number of taxa we find in our sample. However, with the labels, it can be difficult to see exactly what these sample sizes are, so the following code will print it out for us:

```
print(c("Minimum sample size:",min(sample_sums(bacteria)), "Maximum sample size:", max(sample_sums(bacteria))))
```

So, we know that at a minimum we must rarefy to the smaller number. However, considering some samples plateau much before this number, we should choose a smaller cutoff. There is not a strict method to choosing a cutoff, but we should remember we want to include abundant taxa and exclude rare taxa. With this in mind, it is acceptable to rarefy our samples to a size of `10000` reads. Use the following lines of code to achieve this:

```{r}
#Set seed for reproducibility. Rarefy to a sample size equal to or smaller than the minimum sample sum.
set.seed(711)
rarefied<-rarefy_even_depth(bacteria, rngseed = FALSE, sample.size = 10000, trimOTUs = TRUE)

#See that the samples are all the same size now.
rarecurve(as.data.frame(t(otu_table(rarefied))), step=50,cex=0.5,label=TRUE, ylim = c(1,150))
```
Some notes about these commands:
* We first set the seed to some number, in this case `711`. This is for **reproducibility**, since otherwise, `rarefy_even_depth()` would use a random seed by default. This means that if you ran the same code twice without setting the seed, you would get two different rarefied subsets.
* The `rarefy_even_depth()` command needs to know to **not** use a random seed (`rngseed = FALSE`), that our cutoff is 10000 (`sample.size = 10000`), and that we are rarefying the `bacteria_pruned` phyloseq object.
* The `trimOTUs = TRUE` parameter of `rarefy_even_depth()` means that if a taxa is subsampled to an abundance of 0 across all samples, that taxa is removed from the OTU table. Having taxa with 0 reads can mess things up later in the analysis. 

**Question 6:** Why do we prune rare taxa before rarefying?

Excellent! Now that our data is imported, formatted, and rarefied, we can finally look at some diversity metrics!

# Alpha Diversity

Alpha diversity is a metric that evaluates the different types of taxa in a given sample. Different alpha diversity methods typically use calculations based on different components of the sample, which are:

* _Richness:_ the number of taxa in a sample.
* _Evenness:_ the distribution of abundances of the taxa in a sample (i.e. similarities/differences in read quantity per taxa).

Some methods use only one of these components, some will use a combination of both, and others use neither. There are many indices to choose from, but some common ones are:
* Observed taxa: The number of different taxa (richness)
* Chao1 index: another measure of richness, more sensitive to rare taxa
* Shannon index: a combination of richness and evenness, with more weight to the richness component.
* Simpson index: a combination of richness and evenness, with more weight to the evenness component.
* Fisher's diversity index / Fisher's alpha was one of the first methods to calculate the relationship between the number of different taxa and the abundance of individual taxa

First, make a new chunk. Then, try running the following lines of code:

```{r}
#Plot alpha diversity with several metrics.
plot_richness(physeq = rarefied, x="diagnosis", color = "diagnosis", measures = c("Observed", "Shannon", "Fisher")) + geom_boxplot() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
```

Your plots will look something like the following:

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/2c9dd861-585f-4926-a80c-0e1bd0cac0cb" width="800">

We can see that the groups are not identical, and that the different indices yield different plots. As well, some indices are similar to each other (like observed taxa and Fisher's alpha). Since we see differences in both the Shannon and Simpson plots, we can say that there are differences in both richness and evenness between our CD and non-IDB sample groups.

Try adding or changing the measures to see how they compare to one another. Also, try changing value of "x" to different (categorical) metadata variables. 

**Question 7:** How can you use the `View()` command to see what metadata you can choose from?


# Beta Diversity

Beta diversity is another useful metric we use to compare microbiome compositions. Specifically, beta diversity is used to evaluate the differences between entire samples in a dataset. You can imagine that each sample is compared to one another, and the "distance" between these samples is calculated, such that we get a symmetric matrix of pairwise comparisons.

| | Sample A | Sample B | Sample C |
|---|:---:|:---:|:---:|
Sample A | 1 | A vs B | A vs C |
Sample B | B vs A | 1 | B vs C |
Sample C | C vs A | C vs B | 1 | 

To capture beta diversity, we use some metric to calculate these pairwise distances between samples. There are several to choose from, but the most common include:
* Bray-Curtis dissimilarity: takes into account presence/absence and abundance of taxa 
* Jaccard distance: takes into account only presence/absence of taxa
* Weighted UniFrac: accounts for phylogenetic distances between present and absent taxa, weighted by taxa abundance
* Unweighted UniFrac: accounts for phylogenetic distances between present and absent taxa

Once we decide on the metric(s) to use, we then have to consider how to ordinate our data. Ordination is the method by which we take our data points in multidimensional space and project them them to lower dimensional, i.e. 2D space. This allows us to more intuitively visualize our data and the differences between groups.

Common ordination methods include:
* Principal Coordinate Analysis (PCoA) uses eigenvalue decomposition of the distance matrix for normally distributed data
* Non-metric Multidimensional Scaling (NMDS) uses iterative rank-order calculations of the distance matrix for non-normally distributed data

Any combinations of the above metrics and methods can be used, depending on the type of data to be analyzed. For our dataset, we will be using the Bray-Curtis dissimilarity and PCoA method. Fortunately, all of these can easily be implemented in R.

## Preparing the Data

Create a new chunk. First, we will create a relative abundance table from our OTU table. This can be done by transforming our sample counts to a percentage (float) between 0 and 1 using the `transform_sample_counts()` function from phyloseq:

```{r}
percentages <- transform_sample_counts(rarefied, function(x) x*100 / sum(x))
```
Now, we have a new phyloseq object `percentages` in which the abundances are relative. With this object we can create an ordination by selecting our `method` and `distance` metric. Fortunately, the `ordinate` function from phyloseq can do this transformation in one step:

```
ordination<-ordinate(physeq = percentages, method = "PCoA", distance = "bray")
```

## Plot the Data

We are ready to plot our ordination! In this step, we have to specify what data we are using and what our ordination object is. Additionally, we can select our metadata group of interest with the `color` parameter.

```
plot_ordination(physeq = percentages, ordination = ordination, color="diagnosis") +
  geom_point(size=10, alpha=0.1) + geom_point(size=5) + stat_ellipse(type = "t", linetype = 2) + theme(text = element_text(size = 20)) +  ggtitle("Beta Diversity", subtitle = "Bray-Curtis dissimilarity")
```

Your plot should look something like the following:

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/200623fd-fc27-4999-a64d-c94ea351ce2c" width="800">

Now, try substituting the `color` parameter with some different categories from our metadata. Recall how to view the metadata. What does this do to the plot?

As well, you can change the parameters of the ordination like `distance` or `method`, then run the `plot_ordination()` function again. How does this change the plot?


# Visualization with Stacked Bar Charts

Another useful way to visualize the microbial composition of our samples is through the use of stacked bar charts. These charts break down the relative abundances of different taxa in our samples, informing us more about _who_ is making up the communities. Additionally, we can compare our samples side-by-side to see if there are differences in the taxa abundances between samples. We will start by using the `percentages` object that we created in the beta diversity step.

## Formatting the Data

Create a new chunk. We can visualize any taxonomic level with stacked bar charts by collapsing the taxa to their respective nodes at that level. This is done by the `tax_glom()` phyloseq function, which "agglomerates" all of the taxa to the level we specify. Let's start by looking at the different families in our samples.

```{r}
percentages_glom <- tax_glom(percentages, taxrank = 'Family')
```

We can verify that this worked by viewing the `tax_table`. 
**Question 8:** What happened to the data "below" the family level?

Next, we want to "melt" our phyloseq object `percentages_glom` into a new dataframe for constructing our plot with ggplot2. We can do this with the `psmelt()` function and create a new data frame like so:

```
percentages_df <- psmelt(percentages_glom)
```

If we were to keep the data as is, we would be including every Family in our visualization. However, many of the families in our dataframe appear at low abundances across our dataset. This will make our bars look very cluttered, which is contrary to our objective of creating intuitive visualizations. So, it is useful to further collapse all of the families that occur below a certain abundance threshold. Let's try this with an arbitrary value of 1 percent. We can do this by checking all of the `Abundance` values in `percentages_df`, and if the value is less than 1 (i.e. 1 percent) the `Family` value will be replaced with the character string `Family < 1% abundance`.

```
percentages_df$Family[percentages_df$Abundance < 1] <- "Family < 1% abundance"
```

Our `Family` data is categorical, so we have essentially just created a new "Family" category for simplicity. In fact, nearly all of the data in our `percentages_df` dataframe is categorical, which presents a problem. Ggplot2 will want to separate our data by factor levels, so we have to coerce our metadata column of interest to the factor class. We can do this easily with the `as.factor()` function:

```
percentages_df$Family <- as.factor(percentages_df$Family)
```

## Selecting our Colours

Additionally, we can choose the colours of our chart. By default, it will be grayscale, and that's no fun! The RColorBrewer package [provides many palettes](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html) to choose from for different types of data. We can use the "Spectral" palette, but have a problem: "Spectral" only has 11 colours, but we have more than 11 families! We can fix this by using the `colorRampPalette()` function from the "grDevices" package, which interpolates new colours with a given palette.

This command works by first using `brewer.pal(11,"Spectral")`, which returns a character vector of 11 colours in hexadecimal form. Then, we specify how many colours we want to "ramp" to by finding the number of unique factor levels, or Families, in our dataframe with `length(levels(percentages_df$Family))` which returns a numerical value. Finally, `colorRampPalette()(n)` creates a function that returns a new character vector of `n` hexadecimal colour codes, which we store in a new vector `colours`. Funny enough, we have exactly 11 factor levels within `Family`. However this method is useful when we want to look at more or less factors, like if we agglomerated the data to the Fhylum or Genus levels instead.

```
colors<- colorRampPalette(brewer.pal(11,"Spectral"))(length(levels(percentages_df$Family)))
```

Additionally, it will be helpful to see which samples come from each metadata group. Let's consider the `Diagnosis` column of our dataframe. We can differentiate between them by creating another character vector containing colours assigned to each group. To do this quickly, we can use an `ifelse()` function that checks whether the value of each cell in the `Diagnosis` column is equal to `CD`. If this is true, we set the colour to red. If it is false, we know that it must be equal to `nonIDB` instead, and set the colour to blue. R natively supports hundreds of named colours too, which you can view with `colours()`; feel free to pick your favourites! This works for any category that contains two different values.

```
c <- ifelse(percentages@sam_data$diagnosis == "CD", "red", "blue")
```
## Creating the Plot

Finally, we can plot our stacked bar chart! Using the `ggplot()` function, we will combine all of the data and objects we have prepared. We must specify our dataframe with the `data` parameter. As well, we must tell ggplot how we want the graph to look with `aes()`, or the aesthetic map. First, we have to specify what our axes are. Then, we must specify how the bars are to be "filled", or separated. And, we want to add a `theme()` which colours our X-axis labels by the `Diagnosis` group, using the character vector we made above.

```
relative_plot <- ggplot(data=percentages_df, aes(x=ID, y=Abundance, fill=Family))+
  ggtitle("Relative Abundance", subtitle = "Family Level")+
  xlab("Sample ID")+
  ylab("Abundance (%)")+
  geom_bar(aes(), stat="identity", position="stack")+
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, colour = c))

relative_plot
```

Your plot will look similar to the following:

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/0a633076-b8d9-4052-a319-d10977eec2e9" width="800">


Try agglomerating to different taxonomic levels. Doing so, you will have to change _a few_ lines of code and run them before building and viewing the `relative_plot`.

Try using the `high_dysbiosis` metadata group. **Question 9:** How do you have to change the `ifelse()` function to make this work?


# Visualization with Heatmaps



So far, we have looked at how the groups differ from each other at different levels by diversity between `Diagnosis` groups, by comparing the distances between samples, and looking at how they differ taxonomically. However, we have yet to comprehensively look at the taxa inside of our samples. Heatmaps are a great way to visualize this information, and you will see that it allows us to easily directly compare the taxa making up the samples.

## Formatting the data

With taxonomic annotation tools, it can sometimes be difficult to accurately resolve the taxonomy at species and strain levels. For our demonstration, we will again agglomerate the taxa in our dataset, but to the Genus level instead of Family. We will create a new object using use the same `tax_glom()` function as before, but this time start with our `rarefied` phyloseq object. This is because it has not been transformed to relative abundances.

Create a new chunk, then paste and run the following line.

```{r}
hm <- tax_glom(rarefied, taxrank = "Genus")
```

This still leaves us with lots of genera, which we can count with the following command: `length(get_taxa_unique(rarefied, taxonomic.rank = "Genus"))`. This returns 64! We should now consider how we want our heatmap to look. With that many taxa, it will become more difficult to interpret. So, let's represent only the top `n` taxa in our dataset. 

The following command will use the `taxa_sums()` function to first sum up how many reads we have for each genus. Then, the `sort()` function will sort the taxa - to do this we have to specify `TRUE` for the `decreasing` parameter to indicate we want to go from largest to smallest. This results in a named numeric vector, from which we want the names of these top 20 taxa. This is extracted with the `names()` function, which returns a character vector of taxa names. Finally, the `prune_taxa()` function takes this character vector and takes those taxa out of `hm`, and overwrites our original `hm` object.

```
hm<-prune_taxa(names(sort(taxa_sums(hm),decreasing=TRUE)[1:20]), hm)
```

Let's also change the column names of our matrix. Currently, they are the names of our bracken output files, but we can swap them with the values from `ID` stored in the phyloseq object.

```
#change column names to sample IDs
colnames(hm)<-rarefied@sam_data$ID
```

Now, we should format our reduced dataset into a matrix, so our `Heatmap()` function will be able to handle it. We only need the `otu_table` from `hm`, so we can extract it and then coerce it to the matrix class. Doing this converts all of the numbers in the dataframe to the numeric class.

```
hm<-as.matrix(otu_table(hm))
```

If we inspect our matrix with `View()`, we will see that the row names are not actually taxonomic names, but the taxID's. This is not particularly useful for us because taking the time to look up what each of these taxIDs are on the NCBI taxonomy browser is very inefficient. So, we should interchange these taxIDs with the corresponding genus name.

A quick note, the `tax_glom()` function works somewhat counterintuitively with the taxIDs. If several species get collapsed into one genera, the taxID of the first species in that genera gets assigned. So, the taxIDs we have do not correspond to a genus node in NCBI, but one of the original species. This is ok though, because we can still get the genus name from this taxID using the taxonomizr package from NCBI. 

The first time you use this package, you will build and store a SQL database that the functions will rely on. In our case, however, the database is already stored on your virtual instance, so this command will just tell R where to look for it.

```
prepareDatabase(sqlFile = 'accessionTaxa.sql', tmpDir = "~/CourseData/MIC_data/tools/", getAccessions = FALSE)
```

Now, using the `taxa_names()` function we can get the 20 top taxIDs from our `hm` object. This provides a character vector to the `getTaxonomy()` function, which will return the genus name for each taxID, and store them in a new matrix `genera` like so:

```
genera<-getTaxonomy(ids = taxa_names(hm), sqlFile = "accessionTaxa.sql", desiredTaxa = "genus")
```

With this, we can replace the taxIDs with their corresponding genera names. This works because we obtained the genera in the same order that the IDs are in the matrix.

```
rownames(hm)<-genera
```

Finally, we are going to transform our data by taking the base 2 logarithm of all of the abundance values in our matrix. This will allow us to discern more easily between the values that fall in the middle of the range of abundances.

```
hm<-log2(hm)
```

Now if we inspect `hm` with `View()`, there are some new values in our matrix. As well, you may notice that `-Inf` has appeared in several places. This is because any logarithm of 0 is undefined, so R substituted this value when doing the transformation. We can fix this by simply replacing `-Inf` with `0` in our matrix.

```
hm[hm=="-Inf"]<-0
```

## Annotate and Generate the Heatmap

We are nearly ready to build the heatmap, however we should prepare some additional information to include. If we want to show which samples belong to particular metadata groups, we need to make an annotation. The ComplexHeatmap package has [extensive documentation and examples](https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html) to help us add elements to our visualization, but fortunately, adding a simple annotation to our heatmap is a relatively straightforward process.

The annotation will be a bar that sits on top of the heatmap, and indicates which `diagnosis` group the sample belongs to. This should make it apparent if these samples are clustering together in our plot! To do this, we will use the `HeatmapAnnotation()` function. The first parameter we provide is what our categorical variable is and where to find it. You may notice that we are going back to the `sam_data` table of `rarefied`, because `rarefied` contains our sample metadata while the `hm` matrix does not. This works though, as we generated `hm` from `rarefied` above. The function will also assign random colours to indicate the groups, and will include a legend in our final heatmap because we specify the `show_legend` parameter.

```
ha <- HeatmapAnnotation(Diagnosis = rarefied@sam_data$diagnosis, show_legend = TRUE)
```

Another useful step in preparation of the heatmap is sorting the dendrogram that accompanies it. `Heatmap()` will do hierarchical clustering by default and sort the nodes of the resulting dendrogram. However, we can elect to use a different sorting algorithm, such as `dendsort()`. To do this, we have to generate our own hierarchical cluster of the data with `hclust()`, then apply `dendsort()` to the product and save it as an object.

```
col_dend = dendsort(hclust(dist(t(hm))))
```

One last preparation step! We should pick what colour gradient we want the heatmap to use. Again, we can use the RColorBrewer package, but this time using a palette that is good for continuous data, called a [sequential palette](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html#:~:text=types%20of%20palettes%20%3A-,Sequential%20palettes,-are%20suited%20to). We will then ramp the colours to a larger number so that we create a nice gradient. I like "RdPu" as a palette, but you could choose your own as well!

```
colours<- colorRampPalette(brewer.pal(9,"RdPu"))(256)
```

We're finally there! now we can put it all together with the `Heatmap()` function from ComplexHeatmap. In this command, we will use the parameters to specify our matrix containing the data, our custom dendrogram, and our annotation bar. You can also elect whether to show the row dendrogram by changing/removing the `show_row_dend` parameter. Let's plot it!

```
Heatmap(hm, cluster_columns = col_dend,
                show_row_dend = FALSE, col = colours,
                top_annotation = ha,
                row_dend_width = unit(30,"mm"),
                heatmap_legend_param = list(title = "Log(2)\nAbundance"))
```

And your resulting heatmap should look something like the following:

<img src="https://github.com/LangilleLab/microbiome_helper/assets/106988687/e0bbfd25-7cf3-4438-864f-1e611b830bc5" width="800">

Try changing your heatmap annotation to be another metadata variable. Does this cluster better or worse than `diagnosis`?
Try changing the command at near the beginning of this module to create a different heatmap. Maybe look at the least abundant taxa? Or the top 30?


**Question 1:** Which of the graphs does your data resemble more closely?
_The FastQC plots for our samples should look more like the "good" Illumina data. The quality scores do not decrease rapidly like in the "bad" example._

**Question 2:** What can we do if data fails the Per Base Sequence Quality module?
_See [this documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html#:~:text=Failure,is%20often%20employed.) from FastQC for an explanation. One possibility would be trimming the adapter sequence and the low quality part of the read._

**Question 3:** Of the four output files specified above, which should we choose for analyzing the microbiome?
_We should use the _`*paired_contam*.fastq` _output files. This is because they contain the reads that 1) do not align to the human genome; and 2) are not singletons. If the reads do not align to the human reference, there is a chance they are microbial in origin, and we can use them for the next step._

**Question 4:** How many reads are in each sample?
_Using the _`wc` _program, we can see that there are 20000 reads in each of the concatenated .fastq files. Since we know that .fastq files have one read for every four lines, we divide 20000 by 4 to get 5000 reads.

**Question 5:** What do each of these objects within `data` represent?
`otu_table()` _contains the counts for each taxa in our sample._
`tax_table()` _contains the full taxonomy, as given by Kraken2, for each taxa sample._
`otu_table()` _contains the metadata for our samples._

**Question 6:** Why do we prune rare taxa before rarefying?
_We prune the rare taxa first because rarefying randomly subsamples the data to our specified read depth. If these rare taxa are not removed beforehand, they may be included in the rarefied dataset. Why is this bad? Well, rarefying changes the absolute abundances of our samples. When we rarefy, we also round to an integer. So, if a rare taxa with 10 reads in the original sample is rarefied to 0.01 and gets rounded to 1, its abundance is inflated by 100x. We can prune to an arbitrary cutoff, or pick one that is 1/1000th of the mean read depth._

**Question 7:** How can you use the `View()` command to see what metadata you can choose from?
_Try _`View(rarefied@sam_data)` _to view the metadata table_

**Question 8:** What happened to the data "below" the family level?
_It gets collapsed, or "agglomerated", up to the family level, so if family A has genera B, C, and D, the sum of the B, C, and D abundances becomes the abundance of A. So it is not removed, instead it is combined. Interestingly, taxa that have _`NA` _at the family level (or above) are actually removed._

**Question 9:** How do you have to change the `ifelse()` function to make this work?
_Try_ `c <- ifelse(percentages@sam_data$high_dysbiosis=="FALSE", "red", "blue")`
