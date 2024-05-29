# Module 2: Metagenomic Assembly and Binning

This tutorial is part of the 2024 Canadian Bioinformatics Workshops [Advanced Microbiome Analysis](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-workshop-2024-advanced) (St John's, NL, May 29-30).

**Author**: Robyn Wright

## Table of Contents

[Introduction](#introduction)\
[2.1 Initial setup](#21-initial-setup)\
[2.2 Assembly of raw reads with MEGAHIT](#22-assembly-of-raw-reads-with-megahit)\
[2.3 Make an Anvi'o contigs databases](#23-make-an-anvio-contigs-databases)\
[2.4 Run HMMs to identify genes](#24-run-hmms-to-identify-single-copy-genes)\
[2.5 Map samples onto contigs with Bowtie2](#25-map-samples-onto-contigs-with-bowtie2)\
[2.6 Cluster contigs into bins](#26-cluster-contigs-into-bins)\
[2.7 Interactive refinement of bins](#27-interactive-refinement-of-bins)\
[2.8 Export the final bins into MAGs](#28-export-the-final-bins-into-mags)\
[2.9 Run CheckM](#29-run-checkm)\
[2.10 Run CARD RGI annotation](#210-run-card-rgi-annotation)\
[Other things that we do with MAGs](#other-things-that-we-frequently-do-with-mags)\
[Answers](#answers)

## Introduction

The main goal of this tutorial is to introduce students to the assembly of genomes from metagenomic reads (Metagenome Assembled Genomes/MAGs). There is not a one-size-fits-all pipeline for assembling MAGs. MAG assembly is incredibly computationally intensive, and the approach that we are taking here is therefore a very "bare bones" approach designed to demonstrate the main steps involved and give you some familiarity with the methods used. At the end of this tutorial we've provided a few other pipelines for MAG assembly that you may wish to look into if you are looking to assemble MAGs with your own metagenome data.

Throughout this module, there are some questions aimed to help your understanding of some of the key concepts. You'll find the answers at the bottom of this page, but no one will be marking them. 

### Anvi'o

[Anvi'o](https://anvio.org/) is an open-source, community-driven **an**alysis and **vi**sualization platform for microbial **'o**mics. It packages together many different tools used for genomics, metagenomics, metatranscriptomics, phylogenomics, etc. and has great interactive visualisations that can be used to help this. We are just touching the surface of what Anvi'o can do today, but the website has great tutorials and learning resources for all of its capabilities - I recommend browsing through to get some inspiration!

## 2.1 Initial setup

First of all, we need to create the directory that we'll be working from and change into it:
```
mkdir amb_module2
cd amb_module2
```

Then we'll create a symlink to the data that we're going to be using, and copy in the sample names and metadata:
```
ln -s ~/CourseData/MIC_data/AMB_data/mapped_matched_fastq/ .
ln -s ~/CourseData/MIC_data/AMB_data/raw_data/ .
cp ~/CourseData/MIC_data/AMB_data/sample_ids.txt .
cp ~/CourseData/MIC_data/AMB_data/mgs_metadata.txt .
```

Finally, we'll activate the conda environment that we're going to be using:
```
conda activate anvio-7
```

## 2.2 Assembly of raw reads with MEGAHIT

The first step in the assembly of MAGs is the assembly of the shorter reads from sequencing into contigs. A contig is a contiguous sequence assembled from a set of sequence fragments, and there are different bioinformatic methods for assembling our reads into contigs. For this part of the tutorial, we are continuing to use reads from the same samples that we used for taxonomic annotation, but we have sub-sampled these to contain reads from only a few species, to ensure that we have the read depth required for MAG assembly while keeping the files small enough to be able to run on these small AWS instances. These are all still short reads of approximately 100 base pairs each. 

We are using a tool called ```MEGAHIT``` for this assembly because it allows co-assembly of the reads from multiple samples. This means that rather than the assembly being performed separately for each sample, it is performed on all of the reads from all of the samples at the same time. In this tutorial, we will co-assemble all of the samples together, but in your own analyses you should think about what makes the most sense. Are you expecting microbial taxa to overlap between different samples? Would it make sense to find similar genomes in multiple samples? If you answered "yes" to those questions, then it might be worth thinking about co-assembly. If not, then it is probably worth assembling each sample separately. You may also want to consider assembling the samples from different treatments separately.

Prior to assembling your samples you would usually run quality checks and remove potentially contaminating sequences, but seeing as we already did that during the Taxonomic annotation tutorial, we don't need to repeat this step here.

First, we'll make a directory for the output to go into:
```
mkdir anvio
```

Now we'll run MEGAHIT on our samples. You may already know about a program called ```tmux``` (and if you were in the beginner workshop then you definitely should), but it's really useful for running programs that may take a while, or just keeping track of where you are. There are several tools that are pre-installed on most Linux systems that we can use to make sure that our program carries on running even if we get disconnected from the server. ```tmux``` is one of the most frequently used. To activate it, just type in ```tmux``` and press enter. It should take a second to start up, and then load up with a similar looking command prompt to previously, but with a coloured bar at the bottom of the screen.

To get out of this window again, press ```ctrl```+```b``` at the same time, and then d. You should see your original command prompt and something like
```
[detached (from session 1)]
```

We can actually use tmux to have multiple sessions, so to see a list of the active sessions, use:
```
tmux ls
```

We can rename the tmux session that we just created with this:
```
tmux rename-session -t 1 anvio
```
Note that we know it was session 1 because it said that we detached from session 1 when we exited it.

If we want to re-enter this window, we use:
```tmux attach-session -t anvio```

Now, we can run MEGAHIT inside this tmux session, but we will need to activate the conda environment again inside here:
```
conda activate anvio-7
```

Run MEGAHIT:
```
R1=$( ls mapped_matched_fastq/*_R1.fastq | tr '\n' ',' | sed 's/,$//' )
R2=$( ls mapped_matched_fastq/*_R2.fastq | tr '\n' ',' | sed 's/,$//' )
megahit -1 $R1 \
         -2 $R2 \
         --min-contig-len 1000 \
         --num-cpu-threads 4 \
         --presets meta-large \
         --memory 0.8 \
         -o anvio/megahit_out \
        --verbose
```
The arguments here are:
- ```-1``` - The forward reads
- ```-2``` - The reverse reads
- ```--min-contig-len``` - The minimum length in base pairs for contigs that we want to use - 1000 is about the minimum length that you will ever want to use, but sometimes people might increase this to e.g. 2500 or 5000
- ```--num-cpu-threads``` - The number of threads to use
- ```--presets``` - This is for a set of parameters needed within MEGAHIT, and meta-large is suggested for large and complex metagenomes
- ```--memory``` - The amount of available memory that we want to allow MEGAHIT to use - this means it will use up to 80% of the memory available. Keeping this below 100% just means that we would still be able to use the Amazon instance for other things, and could be important if you're sharing a server with other people that might also need to be carrying out some work!
- ```-o``` - The output folder name
- ```--verbose``` - MEGAHIT will print out what it is doing

MEGAHIT can take quite a long time to run, so it's up to you if you want to wait for it to run - if not, press ```ctrl```+```c``` to stop the run and you can copy across the output that we already ran:
```
mkdir anvio/megahit_out
cp ~/CourseData/MIC_data/AMB_data/final.contigs.fa anvio/megahit_out/
```

If you ran it yourself, a step that we'll often do is removing the intermediate contigs to save space:
```
rm -r anvio/megahit_out/*/intermediate_contigs
```

The main output at this point is a fasta file containing the contigs ```anvio/megahit_out/final.contigs.fa```. You can take a look at this with the ```less``` command if you like (remember to press ```q``` to exit this view), and we can also count the number of contigs that we have with ```grep -c ">" anvio/megahit_out/final.contigs.fa```.

**Question 1**: How many contigs are there in the ```anvio/megahit_out/final.contigs.fa``` file?

## 2.3 Make an Anvi'o contigs databases

First of all, we'll run a script to reformat our final contigs file from MEGAHIT. This ensures that they're in the right format for reading into Anvi'o in the next step:
```
anvi-script-reformat-fasta anvio/megahit_out/final.contigs.fa \
                               --simplify-names \
                               --min-len 2500 \
                               -o anvio/megahit_out/final.contigs.fixed.fa
```
Take a quick look at both ```anvio/megahit_out/final.contigs.fa``` and ```anvio/megahit_out/final.contigs.fixed.fa``` using the ```less``` or ```head``` commands to see what was changed! You can read more about the ```anvi-script-reformat-fasta``` command on [this page](https://anvio.org/help/main/programs/anvi-script-reformat-fasta/). You should be able to see that we're giving this command a few options (aside from obviously giving it the ```final.contigs.fa``` file from the previous step:
- ```--simplify-names``` - this is simply telling the program to simplify the names of the contigs in the file - you should have seen that in the original ```anvio/megahit_out/final.contigs.fa``` file, they had a description as well as a name, and this included information that Anvi'o doesn't need and may confuse it.
- ```--min-len``` - although we already used a minimum contig length of 1000 in the previous step, we're further filtering here to ensure that we include only the best contigs, as well as reducing the computational time for the next steps.
- ```-o``` - the name of the output file.

Now that we've prepared the file, we can read this into Anvi'o:
```
mkdir anvio/anvio_databases
anvi-gen-contigs-database -f anvio/megahit_out/final.contigs.fixed.fa \
                              -o anvio/anvio_databases/CONTIGS.db \
                              -n HMP2
```
You can read more about what ```anvi-gen-contigs-database``` is doing [here](https://anvio.org/help/7/programs/anvi-gen-contigs-database/), but the key parts of this are:
- ```f``` - the input file of fixed contigs
- ```-o``` - what the output should be saved as
- ```-n``` - the name we're giving this contigs database

Throughout the next steps, we're going to be adding information to this contigs database.

## 2.4 Run HMMs to identify single copy genes

The first thing that we're going to add to the contigs database is information on the genes within the contigs:
```
anvi-run-hmms -c anvio/anvio_databases/CONTIGS.db \
                              --num-threads 4
```
Because the only arguments we've given this are the contigs database (```-c```) and the number of threads to use (```--num-threads```), it will run the default Hidden Markov Models (HMMs) within Anvi'o. You can read more about HMMs in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2766791/) if you're interested, but in short, they are "statistical models that can be describe the evolution of observable events that depend on internal factors, which are not directly observable". HMMs are often used in computational biology for predicting the function of a gene or protein; we can build HMMs using multiple sequence alignments of genes or proteins known to carry out the same function, and then we can run these HMMs against our protein (amino acid) or gene (DNA) sequences to find other proteins/genes that are likely to carry out this function. You can read about the default HMMs in Anvi'o [here](https://anvio.org/help/7/artifacts/hmm-source/), but they include Bacteria_71 (71 single-copy core genes for the bacterial domain), Archaea_76 (76 single-copy core genes for the archaeal domain), and Protista_83 (83 single-copy core genes for the protists, within the eukaryotic domain).

If we had our own set of HMMs, or HMMs for specific genes of interest, we could also use these here.

This next step exports the sequences for all of the genes that we've identified using the HMMs to a file called ```anvio/anvio_databases/gene_calls.fa```:
```
anvi-get-sequences-for-gene-calls -c anvio/anvio_databases/CONTIGS.db \
                              -o anvio/anvio_databases/gene_calls.fa
```

**Question 2**: How many genes were identified?

## 2.5 Map samples onto contigs with Bowtie2

The next steps are going to generate abundance profiles of our contigs within the samples - these will help us to bin the contigs together into groups that have similar abundance profiles and are therefore likely to come from the same genome. 

The first step here is to build a Bowtie2 database/index of our fixed contigs:
```
bowtie2-build anvio/megahit_out/final.contigs.fixed.fa anvio/megahit_out/final.contigs.fixed
```
You can see that here we just have two positional arguments - the file name for the fasta file that we want to use, and the prefix to use for the Bowtie2 index If you look at the files in ```anvio/megahit_out/``` now, you should see six files with this prefix and the file extension ```.bt2```. These files are what Bowtie2 will use in the next steps for mapping the samples to the contigs index. 

For the next step, we've already made a file that contains the names of all of the samples that we'll be using (```sample_ids.txt```) - we can use this to loop through each of the samples using a for loop, performing the same actions on the files for each sample:
```
mkdir anvio/bam_files

for SAMPLE in `awk '{print $1}' sample_ids.txt`
do

    # 1. do the bowtie mapping to get the SAM file:
    bowtie2 --threads 4 \
            -x anvio/megahit_out/final.contigs.fixed \
            -1 "raw_data/"$SAMPLE"_R1_subsampled.fastq.gz" \
            -2 "raw_data/"$SAMPLE"_R2_subsampled.fastq.gz" \
            --no-unal \
            -S anvio/bam_files/$SAMPLE.sam

    # 2. convert the resulting SAM file to a BAM file:
    samtools view -F 4 -bS anvio/bam_files/$SAMPLE.sam > anvio/bam_files/$SAMPLE-RAW.bam
    
    # 3. sort and index the BAM file:
    anvi-init-bam anvio/bam_files/$SAMPLE-RAW.bam -o anvio/bam_files/$SAMPLE.bam
    
    # 4. remove the intermediate BAM file that is no longer needed:
    rm anvio/bam_files/$SAMPLE-RAW.bam

done
```
Hopefully you either remember about for loops from Beginner Module 4 or you already have some experience with them! Otherwise, I'll explain them briefly here.

For loops are an incredibly useful thing that we do in programming - to my knowledge, they exist in every programming language, and they allow us to repeat a section of code many times, with a different input each time. To explain what we are doing above, we can go through this step-by-step. The code is essentially saying that for each of the rows in the ```sample_ids.txt``` file, we want to: (1) run Bowtie2, (2) convert the output SAM file from Bowtie2 to a BAM file, (3) sort and index the BAM file, and (4) remove the intermediate BAM file. You can see that we give ```$SAMPLE``` in several places, as these parts will change for each of the sample names. 

In each of the steps within the for loop, we are giving several options:
1. We are running ```Bowtie2``` using the contigs index that we created ```-x```, forward and reverse reads for this sample, ```-1``` and ```-2```, and the ```--no-unal``` options means it won't output unaligned reads. We'll then get ```-S``` as the output SAM file. SAM (Sequence Alignment MAP) files are used to store information about sequences aligned to a reference and they have 11 mandatory fields corresponding to various metrics about how well the sequence aligns to the reference. You can read more about them [here](https://en.wikipedia.org/wiki/SAM_(file_format)). 
2. ```samtools``` is being used to convert the SAM file to a BAM file. You can see that this uses the ```view``` module of samtools, which would by default print the output to the Terminal, but as we've added ```> anvio/bam_files/$SAMPLE-RAW.bam```, this tells it to save the output to a file. A BAM (Binary Alignment MAP) file is a compressed, binary version of a SAM file (see more [here](https://en.wikipedia.org/wiki/Binary_Alignment_Map)). 
3. The running of ```anvi-init-bam``` is fairly straight-forward - we are just giving it the raw BAM file and it is outputting a sorted, reindexed version.
4. Then we just remove the intermediate BAM file. 

### Add coverage and detection statistics to the Anvi'o profile

Now that we've created information on the coverage and detection of our contigs within our samples (how abundant they are and whether they appear in every sample or not), we can make Anvi'o profiles with this information: 
```
mkdir anvio/anvio_databases/profiles

for SAMPLE in `awk '{print $1}' sample_ids.txt`
do

    anvi-profile -c anvio/anvio_databases/CONTIGS.db \
                 -i anvio/bam_files/$SAMPLE.bam \
                 --num-threads 4 \
                 -o anvio/anvio_databases/profiles/$SAMPLE
done
```
Note that we're again using a for loop to carry this out on each of our samples. The ```anvi-profile``` program quantifies coverage per nucleotide position within the contigs, and averages them per contig. It also calculates single-nucleotide, single-codon, and single-amino acid variants, as well as structural variants such as insertion and deletions and stores these data into appropriate tables.

### Merge sample profiles

The command above created contig profiles for each sample individually, so now we need to merge these profiles into one single profile that we can use for clustering the contigs into bins. We'll do this with the ```anvi-merge``` program:
```
anvi-merge -c anvio/anvio_databases/CONTIGS.db \
           -o anvio/anvio_databases/merged_profiles \
           anvio/anvio_databases/profiles/*/PROFILE.db
```

And then it's always a good idea to have a look at some basic stats about the contigs before we go further:
```
anvi-display-contigs-stats --report-as-text \
                           --output-file contigs_stats.txt \
                           anvio/anvio_databases/CONTIGS.db
```
If you take a look at the ```contigs_stats.txt``` file, you'll see (you can also see these explanations [here](https://anvio.org/help/main/programs/anvi-display-contigs-stats/):
- ```Total Length``` - the total number of nucleotides in your contigs
- ```Num Contigs``` - the number of contigs in your database
- ```Num Contigs > X kb``` - the number of contigs that are longer than X
- ```Longest Contig``` - the longest contig in your databases (in nucleotides)
- ```Shortest Contig``` - the shortest contig in your databases (in nucleotides), hopefully this is longer than the 2500 that we set earlier!
- ```Num Genes (prodigal)``` - the number of genes that Prodigal predicts are in your contigs
- ```L50```, ```L75```, ```L90``` - if you ordered the contigs in your database from longest to shortest, these stats describe the number of contigs you would need to go through before you had looked at a certain percent of a genome. For example, L50 describes the number of contigs you would have to go through before you reached 50 percent of the entire dataset
- ```N50```, ```N75```, ```N90``` - if you ordered the contigs in your database from longest to shortest, these stats describe the length of the contig you would be looking when you had looked at a certain percent of a genome. For example, N50 describes the length of contig you would be on when you reached 50 percent of the entire genome length
- The number of HMM hits in your database
- The number of genomes that Anvi’o predicts are in your samples, based on how many hits the single-copy core genes got in your database

As long as we're satisfied with all of this, then we can carry on to the clustering.

**Question 3**: How many contigs are there? Is this the same as what we started with? Why or why not?\
**Question 4**: What are the longest and shortest contigs? What do you think of this? 

## 2.6 Cluster contigs into bins

Next we will bin - or cluster - the contigs to create genome "bins", or MAGs. To do this, we typically use information on the coverage, or abundance, of contigs within samples (that we generated in step two) and the binning algorithms will usually identify patterns in nucleotide composition or k-mer frequencies in order to group together contigs that they think are likely to have originated from the same genome. We are using ```CONCOCT``` for this, but as with every step of this pipeline, there are multiple tools that are capable of carrying this out and may be better or worse in different circumstances.

Now we'll run the clustering:
```
anvi-cluster-contigs -c anvio/anvio_databases/CONTIGS.db \
                         -p anvio/anvio_databases/merged_profiles/PROFILE.db \
                         -C "merged_concoct_2500" \
                         --driver CONCOCT \
                         --length-threshold 2500 \
                         --num-threads 4 \
                         --just-do-it
```
Here you can see that we're giving a few options:
- ```-c``` - the contigs database
- ```-p``` - the profile database
- ```-C``` - the name for the resulting collection of bins to be given
- ```--driver``` - the clustering algorithm to be used. We could also have used metabat2, maxbin2, dastool, or binsanity. If you're running this pipeline on your own data then you will probably want to think about which of these algorithms you should use; dastool has the advantage that it can take the input of multiple binning algorithms to calculate an optimised set of bins. It's also worth noting that no binning algorithm is perfect, and you will likely need to check the bins and refine them manually after any of these are run
- ```--num-threads``` - the number of threads to use
- ```--just-do-it``` - ignore any warnings (like concoct being implemented experimentally into Anvi'o) and just run it anyway

Next, we'll summarise the bins that concoct has given us:
```
mkdir anvio/concoct_summary_2500
anvi-summarize -c anvio/anvio_databases/CONTIGS.db \
                   -p anvio/anvio_databases/merged_profiles/PROFILE.db \
                   -C "merged_concoct_2500" \
                   -o anvio/concoct_summary_2500/summary/
```

And we'll take a look at the bin summaries in the Terminal quickly just to check that everything looks OK:
```
less anvio/concoct_summary_2500/summary/bins_summary.txt
less anvio/concoct_summary_2500/summary/bins_across_samples/abundance.txt
less anvio/concoct_summary_2500/summary/bins_across_samples/mean_coverage.txt
```
Remember that you can exit the ```less``` view by pressing ```q```!

We can also look at this in a more interactive way in our workspace: http://##.uhn-hpc.ca/ (remember to replace the ## with your number!)
Go to ```amb_module2/anvio/concoct_summary_2500/```. Right click on ```summary``` > open in new tab. You can scroll through and explore what this says about our bins so far, and then go back to the other tab with ```amb_module2/anvio/concoct_summary_2500/```. Add ```/summary/bins_summary.txt``` to the URL bar and copy and paste the resulting page into a new Excel (or whatever you usually use for viewing spreadsheets) document - it will be useful to refer back to. If you're in Excel, you can easily get this into columns by going to the "Data" tab > click on "Text to columns" > check "Delimited" > Next > Check "Space" > Finish.

Looking through the bins (the rows), you should see that there are a number of columns giving some stats on each of the bins (clusters of contigs), including:
- ```total_length``` - the total number of nucleotides in this bin
- ```num_contigs``` - the total number of contigs in this bin
- ```N50``` - N50 is a metric widely used to assess the contiguity of an assembly, and it is defined as the sequence length of the shortest contig that covers at least 50% of the total assembly length
- ```GC_content``` - GC content (%) is the percentage of bases in the contigs that make up this bin that are either guanine (G) or cytosine (C)
- ```percent_completion``` and ```percent_redundancy``` - completeness and contamination (also called redundancy) are both important when we're assessing the quality of the MAGs that we've assembled. Both scores are calculated based on the presence of the ubiquitous, single-copy marker genes - ones that all bacteria are known to possess - and the completion is a prediction of how complete the genome is likely to be, so whether it possesses a copy of all of the marker genes that are expected. The contamination/redundancy is a measure of whether those marker genes that are expected to be present in only a single copy are duplicated. Typically speaking, a MAG that has >50% completeness and <10% contamination/redundancy is considered to be reasonable, although >90% completeness is desirable. There is a good explanation on completeness and redundancy [here](https://merenlab.org/2016/06/09/assessing-completion-and-contamination-of-MAGs/).

You'll see that we have quite a lot of bins that are 0% completion (also 0% reduncancy), and typically these are made up of a single contig. Now, we'll just be focussing on the bins that have completion >50% and we'll be trying to refine them so that their redundancy is <10%. 

**Question 5**: How many bins are there?\
**Question 6**: How many bins >50% completion are there?\
**Question 7**: What is the redundancy in these bins?\

## 2.7 Interactive refinement of bins

Now we're going to be doing the interactive refinement, and first of all we'll have a look at the bin that's 100% complete and 0% redundant - for me, this is called Bin_17, but I don't know if this naming is always consistent, so if that's a different bin for you then don't worry!

To do this, we'll run the Anvi'o refind command:
```
anvi-refine -c anvio/anvio_databases/CONTIGS.db \
            -p anvio/anvio_databases/merged_profiles/PROFILE.db \
            -C "merged_concoct_2500" \
            -b Bin_17 \
            --server-only \
            -P 8081
```
You'll see that we're telling Anvi'o the contigs database, profile and collection name, as well as the name of the bin we want to look at and: 
```
--server-only
-P 8081
```
Both of these parts are to do with Anvi'o being run on the Amazon instances rather than on our local computers. The ```--server-only``` part is telling it that we will want to create an SSH tunnel to the server, and then the ```-P``` port is telling it which port to use. This could be one of many ports, just like we are using port 8080 for accessing RStudio. 

Now the instructions are different depending on whether we're on Mac/Linux or Windows.

### Mac

Open a second Terminal window or tab and run:
```
ssh -L 8081:localhost:8081 -i BMB.pem ubuntu@##.uhn-hpc.ca
```
Where ## is your number. It should just look like you logged into the server in a new window.
Now go to http://localhost:8081/ in your browser. This should have an Anvi'o page loaded up.

### Windows

Open a new Putty window. 

Click on "Session" and under "Saved Sessions", click on the "Amazon node" and then click "Load".

Now go to "Connection" > "SSH" > "Tunnels". In "Source port" type in ```8081```. In "Destination" type in ```90.uhn-hpc.ca:8081```.

Click "Add" and then click "Open". You should see a new Putty window open. Now go to http://localhost:8081/ in your browser. This should have an Anvi'o page loaded up.

### Back to both

Now press the "Draw" button. You should see something that looks like a phylogenetic tree get drawn. Each branch of this is for one of the contigs that makes up the bin, and you can see information about their abundance in different samples in the rings. 

Now click on the "Bins" tab.

If we select the whole tree (click on it), we can see that this is 100% complete and 0% redundant. If we unselect (right click) some of the tree branches, this makes the completion and the length of the genome go down slightly, but it doesn't affect the redundancy.

Now that we've seen how this looks for a bin with 100% completion and 0% redundancy, we can leave this and try another.

Go back to the first Terminal window (where we ran ```anvi-refine```) and press ```ctrl```+```c``` to stop.

Now run this with another bin - choose the one that is 59% complete but has 54.9% redundancy. I'm using Bin_27 here, but as I said before, if this is a different bin for you then that's fine:
```
anvi-refine -c anvio/anvio_databases/CONTIGS.db \
            -p anvio/anvio_databases/merged_profiles/PROFILE.db \
            -C "merged_concoct_2500" \
            -b Bin_27 --server-only -P 8081
```
Refresh your Anvi'o browser window and you should see a different page now.

Again, click "Draw" and then go to the "Bins" tab. We ideally want to get the redundancy below 10% while keeping the completion above 50%. 

Start by selecting the whole tree, and then begin by deselecting (right clicking) branches that appear to have a different abundance profile than others. For example, there is one that seems to only be present in sample CSM7K0MH - deselecting this one only reduces the completion a little (to 56.3%), but reduces the redundancy to 15.5%.

You can deselect large areas to determine whether they'll have an impact on the redundancy, and then go through in more detail and deselect/reselect individual parts of this. Make sure you use the zoom +/- and zoom in and out to see the tree in more detail, as well as checking the abundance/GC content profiles of the parts that you've deselected.

Sometimes you might want to add another bin to keep some of the contigs that you don't think belong in the first one (this particularly might be the case if you had really high completion initially).

When you are happy with what you've selected, click on "Store refined bins in the database". You can now go back to the first console window and press ```ctrl```+```c```.

The other two bins that are >50% completion are already <10% redundancy, so it's not necessary to do anything to these, but you can have a look at them in the same way as we did with these two if you like. 

## 2.8 Export the final bins into MAGs

Now that we've refined our bins into MAGs, we can rename them:
```
anvi-rename-bins -c anvio/anvio_databases/CONTIGS.db \
                     -p anvio/anvio_databases/merged_profiles/PROFILE.db \
                     --collection-to-read merged_concoct_2500 \
                     --collection-to-write FINAL \
                     --call-MAGs \
                     --min-completion-for-MAG 50 \
                     --max-redundancy-for-MAG 10 \
                     --prefix HMP2 \
                     --report-file renaming_bins.txt
```
You can see that here we're defining a new collection called "FINAL", and we're saying that to rename these bins as MAGs, they should be >50% completion and <10% redundancy. You can look at the file that's created with the renamed bins if you like, to check that this did what you expected: ```renaming_bins.txt```.

Now we'll summarise this new collection that we've made:
```
anvi-summarize -c anvio/anvio_databases/CONTIGS.db \
                   -p anvio/anvio_databases/merged_profiles/PROFILE.db \
                   -C "FINAL" \
                   -o anvio/FINAL_summary/
```
And we can have a look at the bins that we've created. Take a look in the folder ```anvio/FINAL_summary/bin_by_bin``` using the ```ls``` command. You should see that some of these are called ```*MAG*``` while some are called ```*Bin*``` - this is because in the command above, we're only calling things a MAG if they're >50% completion and <10% redundancy.

Have a look at the summary file again. You can see that for most of them, the source is shown as ```concoct```, but the MAG that we refined is shown as ```anvi-refine```, and all of the MAGs are identified as bacteria.

Now have a look in one of those MAG folders (```anvio/FINAL_summary/bin_by_bin/```) - you'll see a lot of statistics, as well as a ```*-contigs.fa```. This is a fasta file of the contigs used for each MAG, and we can use it as their genome for further analyses. 

Now we'll create a copy of these in a new folder:
```
mkdir MAG_fasta
cp anvio/FINAL_summary/bin_by_bin/*MAG*/*contigs.fa MAG_fasta/
ls MAG_fasta/
```
Now that we have our refined MAGs, we can do anything that we like with them!

## 2.9 Run CheckM

Next, we're going to run CheckM. CheckM can give us information on the quality of genomes as well as assigning taxonomy to them.

First, we'll activate the conda environment:
```
conda activate checkm
```

And then we'll tell CheckM where to find the data that it needs to run:
```
export CHECKM_DATA_PATH=/home/ubuntu/CourseData/MIC_data/tools/CheckM_databases
```

Now we'll start the ```lineage_wf``` of CheckM. This will take a little while to run, so you could do it in a ```tmux``` session if you like:
```
checkm lineage_wf --reduced_tree -t 4 -x fa MAG_fasta MAGs-CHECKM-lineage -f MAGs-CHECKM.txt --tab_table
```
In this command, we're telling CheckM to run it's ```lineage_wf``` on all of the ```.fa``` files within the ```MAG_fasta``` folder. It will do this using 4 threads (```-t 4```) and will write the output to the ```MAGs-CHECKM-lineage``` folder. Internally, CheckM will: (1) insert each of the MAGs into a reference genome tree; (2) use the tool Prodigal to identify genes within each of the fasta bins, and identifies the lineage-specific marker genes for each of the genomes; (3) analyze whether the lineage-specific marker genes are found in each of the genomes; and (4) summarize the quality of each genome bin/MAG.

If you didn't start this in a ```tmux``` session, you may want to open a new tab and start with the next CARD RGI steps while this runs. 

Once it's finished running, you can take a look at the output from it, ```MAGs-CHECKM.txt``` (either with ```less``` or in the file viewer). 

You can see that CheckM tells us which set of markers it used for assessing genome quality, the number of genomes (reference genomes used to infer the lineage-specific marker set), markers (number of marker genes within the inferred lineage-specific marker set) and marker sets (number of co-located marker sets within the inferred lineage-specific marker set) that were used, 0-5+ shows the number of times each marker gene is identified, and then we have the completeness/contamination as determined from the presence/absence of marker genes and the expected collocalization of these genes, and the strain heterogeneity. 

The estimated strain heterogeneity is determined from the number of multi-copy marker pairs which exceed a specified amino acid identity threshold (default = 90%). High strain heterogeneity suggests the majority of reported contamination is from one or more closely related organisms (i.e. potentially the same species), while low strain heterogeneity suggests the majority of contamination is from more phylogenetically diverse sources 

Now you can run the next step:
```
checkm tree_qa MAGs-CHECKM-lineage -f MAGs-CHECKM-tax.txt --tab_table
```

And then take a look at ```MAGs-CHECKM-tax.txt``` (either with ```less``` or in the file viewer).

This output is a little more straight-forward to understand, and you can see that it gives the taxonomy down to the genus level for each of our MAGs. 

**Question 8**: Are these taxa what you would have expected based on the read-based taxonomy of the samples?

## 2.10 Run CARD RGI annotation

While general functional annotation can be done inside Anvi'o, one of the things that we frequently are interested in when we have MAGs is the presence of antibiotic/antimicrobial resistance (AMR) genes. We're often interested in finding out which MAGs these might be in because finding out whether these bacteria are also likely to be human pathgens is really relevant for how concerned we may be that they possess AMR genes. 

There are many bioinformatic tools that have been designed for this purpose, but one downfall of many tools is that they aren't well maintained after they're developed. However, the Comprehensive Antibiotic Resistance Database (CARD) Resistance Gene Identifier (RGI) has been shown to work well, and is well maintained with the database being updated regularly. 

Often when we run tools like this, we'll need to download and setup the database before we can use it, so that's what we'll do now.

First, activate the conda environment and make the directory for the database:
```
conda activate rgi
cd ..
mkdir card_data
cd card_data
```

Now we can download the data (following the instructions from [here](https://github.com/arpcard/rgi/blob/master/docs/rgi_load.rst):
```
wget https://card.mcmaster.ca/latest/data --no-check-certificate
```

Unzip it:
```
tar -xvf data ./card.json
```

Load the database and then find out what version the database is:
```
rgi load --card_json ./card.json --local
rgi database --version --local
```
We should see that this is version 3.2.9.
Note that the ```--version``` command is something that works with many tools. Try it with Python like this: ```python --version```.

Now we'll get another part of the database and unzip it:
```
wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
mkdir -p wildcard
tar -xjf wildcard_data.tar.bz2 -C wildcard
gunzip wildcard/*.gz
```

Now we'll get the annotations:
```
rgi card_annotation -i localDB/card.json > card_annotation.log 2>&1
rgi wildcard_annotation -i wildcard --card_json localDB/card.json -v 3.2.9 > wildcard_annotation.log 2>&1
```

And finally load the database:
```
rgi load \
  --card_json localDB/card.json \
  --debug --local \
  --card_annotation card_database_v3.2.9.fasta \
  --wildcard_annotation wildcard_database_v3.2.9.fasta \
  --wildcard_index wildcard/index-for-model-sequences.txt \
  --wildcard_version 3.2.9 \
  --amr_kmers wildcard/all_amr_61mers.txt \
  --kmer_database wildcard/61_kmer_db.json \
  --kmer_size 61
```

Now we're ready to run CARD RGI!

So we'll set up the directories:
```
cd ..
mkdir amb_module2/card_out/
cd card_data
```

And then we'll run CARD RGI using parallel:
```
parallel -j 1 'rgi main -i {} -o ~/workspace/amb_module2/card_out/{/.} -t contig -a DIAMOND -n 4 --include_loose --local --clean' ::: ~/workspace/amb_module2/MAG_fasta/*
```
Hopefully by now you can work out what most of these parts must be, and you can see a description of them for yourself by typing ```rgi main --help```, but there are some that might not be so obvious:
- ```-t``` - whether the data input is contigs or proteins
- ```-a``` - the alignment tool to use (DIAMOND or BLAST)
- ```--include-loose``` - that we want to include loose hits in addition to strict and perfect hits
- ```--local``` - that we want to use the local database (i.e., that we don't need to download the database)
- ```--clean``` - that we want to remove temporary files when we're done

Once this finishes running, we can change to the results folder and make a heatmap:
```
cd ~/workspace/amb_module2/
rgi heatmap --input card_out/ --output card_heatmap
```

If you take a look at the output file, ```card_heatmap-4.png```, then you can see that there are just a few AMR gene hits in each genome. In these heatmaps, yellow represents a perfect hit, teal represents a strict hit, and purple represents no hit. 

Try searching for some of these genes on the [CARD website](https://card.mcmaster.ca/browse) to find out some more information about them.

**Question 9**: How many different ARGs are there?

## Other things that we frequently do with MAGs

While I hope that you now have a reasonable idea of the steps that go into assembling MAGs, how to assess their quality, and some of the ways in which we might annotate their functions, this has really only scratched the surface of what they might be used for. 

If you still have time in the workshop and want to have a look, I've added some papers that I think are nice uses of MAGs as well as some other things that can be done with MAGs that you might want to explore.

**Papers**:
- [Bacterial ecology and evolution converge on seasonal and decadal scales](https://www.biorxiv.org/content/10.1101/2024.02.06.579087v1)
- [Nitrogen-fixing populations of Planctomycetes and Proteobacteria are abundant in surface ocean metagenomes](https://www.nature.com/articles/s41564-018-0176-9)
- [Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life](https://www.nature.com/articles/s41564-017-0012-7)
- [A genomic catalog of Earth’s microbiomes](https://doi.org/10.1038/s41587-020-0718-6)

**Other tools**:
- [GTDB-toolkit](https://github.com/Ecogenomics/GTDBTk) - this can be used for assigning taxonomy to MAGs and has a much larger reference database than CheckM, which is designed for generating quality metrics of MAGs
- [Other Anvi'o workflows and tutorials](https://anvio.org/learn/)
- [Programs available within Anvi'o](https://anvio.org/help/main/)
- [BlastKOALA and GhostKOALA](https://www.kegg.jp/ghostkoala/)
- [IslandViewer](https://www.pathogenomics.sfu.ca/islandviewer)
- [METABOLIC](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8)

## Answers

**Question 1**: How many contigs are there in the ```anvio/megahit_out/final.contigs.fa``` file?\
5,535 contigs (sequences in the file)
	
**Question 2**: How many genes were identified?\
20,213 identified (sequences in the file)
	
**Question 3**: How many contigs are there? Is this the same as what we started with? Why or why not?\
There are 2,185 contigs. This is less than we started with because we only kept contigs with >2500 bp. 
	
**Question 4**: What are the longest and shortest contigs? What do you think of this?\
171,697 and 2,501 bp, respectively. This seems pretty good!
	
**Question 5**: How many bins are there?\
34
	
**Question 6**: How many bins >50% completion are there?\
4
	
**Question 7**: What is the redundancy in these bins?\
0%, 4.23%, 54.93% and 7.04%
	
**Question 8**: Are these taxa what you would have expected based on the read-based taxonomy of the samples?\
Yes! These are similar to the abundant taxa in the read-based analyses. 
	
**Question 9**: How many different ARGs are there?\
There are 5 shown in the plot. 
	
**Question 10**: What are they?\
CfiA14, adeF, vanT gene in vanG cluster, qacJ, vanW gene in vanI cluster
