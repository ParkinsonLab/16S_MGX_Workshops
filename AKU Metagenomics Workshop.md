# Metagenomics Workshop
## Introduction
We'll now begin the metagenomics portion of the tutorials! The processing of these samples is much more lengthy and involved than with 16S data. While there are pipelines that will get you from start to end, workflows can vary greatly depending on your purposes and it is important to run through each step carefully to appreciate what is happening with these samples. 

These metagenomics sequences were publically available on the European Nucleotide Archive, and were taken from stool samples collected from athletes (you can find them [here](https://www.ebi.ac.uk/ena/browser/view/PRJNA472785) and [here](https://www.ebi.ac.uk/ena/browser/view/PRJNA305507)). These data were most recently used for [this study](https://doi.org/10.1186/s40168-023-01470-9) on how the gut microbiome can affect an athlete's physical performance. For this tutorial, 8 paired-end samples are provided.

Throughout this workshop, we will run into several instances where it is not feasible to complete a step within the time alloted. We will provide pre-computed files in these cases, and we will indicate where to find them. 

This workshop is adapted from the [2024 Canadian Bioinformatics Workshops Advanced Microbiome Analysis](https://bioinformaticsdotca.github.io/AMB_2024) and the [2023 Metatranscriptomics Practical Lab](https://github.com/bioinformaticsdotca/MIC_2023/blob/main/module6/metatranscriptomics_lab.md) distributed by Canadian Bioinformatics Workshops. 

## Lab 1 - Taxonomic Annotation
Our first tutorial will cover quality control and taxonomic annotation. We are dealing with untargeted, shotgun metagenomics data now, so we'll actually find some unwanted human DNA in these files. 

### Processing Raw Reads
For the metagenomics labs, we will be using a set of 8 gut microbiome samples collected from various athletes. We can examine the files by checking how many reads are in each sample, and looking at the metadata associated with each sample. 
```bash
ls athlete_samples
wc -l athlete_samples/*1.fastq

cat sample_metadata.tsv
```

Our first step will be to run fastqc, to examine the quality of these samples. We can open each fastqc report as an html file and assess our quality metrics. 
```bash
mkdir fastqc_out
fastqc -t 4 athlete_samples/*fastq -o fastqc_out
```

Open a couple of the resulting .html files for review. Our sequences are generally of good quality!  

>Question 1.1
>
>Look at the *Per base sequence content* of some of these reports. What do you notice?

As mentioned before, these samples were taken from the stool of athletes, and will contain some human data we don't want in our final analysis. We are also using parallel, a tool that helps us run a command with a series of different inputs. parallel is used in this case to feed kneaddata each pair of your fastq files without the need for you to type out this command for each file. Removing contaminant human sequences using kneaddata can be done as follows:
```bash
parallel -j 1 --eta --xapply 'kneaddata --input {1} --input {2} -o kneaddata_out -db databases/kneaddata --bypass-trim --remove-intermediate-output' ::: athlete_samples/*1.fastq ::: athlete_samples/*2.fastq
```

Now we combine relevant output reads into one fastq using a perl script from MicrobiomeHelper.
```bash
perl scripts/concat_paired_end.pl -p 4 --no_R_match -o cat_reads kneaddata_out/*_paired_contam*.fastq
```

Unfortunately, we can see that we have very few reads remaining. In this case we chose samples with a low number of reads for a small scale lab, and as such there might not have been as much bacterial biomass. Just for this tutorial, we will go ahead with the raw reads (not recommended for a normal experiment!).
```bash
# raw reads
perl scripts/concat_paired_end.pl -p 4 --no_R_match -o cat_reads_full athlete_samples/*.fastq
```

>Question 1.2
>
>What are the dangers of using these raw reads in a real experiment?

### Taxonomic Annotation
Using kraken2 to annotate these files:
```bash
mkdir kraken2_kreport
mkdir kraken2_outraw
parallel -j 1 --eta 'kraken2 --db databases/kraken2  --output kraken2_outraw/{/.}.kraken --report kraken2_kreport/{/.}.kreport --confidence 0 {}' ::: cat_reads_full/*.fastq
```

Using bracken to quantify our taxa from kraken2 outputs.
```bash
mkdir bracken_out
parallel -j 2 --eta 'bracken -d databases/kraken2 -i {} -o bracken_out/{/.}.species.bracken -r 100 -l S -t 1' ::: kraken2_kreport/*.kreport
```

Exports our kraken stuff to R.
```bash
# note we altered the sample names column to contain the suffix "_bracken_species"
kraken-biom kraken2_kreport/*bracken_species.kreport -m sample_metadata_2.tsv -o mgs.biom --fmt json
```

### Exploration of Data
Open up RStudio, and load in a few packages we'll need for this next section.
```r
library(phyloseq) # organizational tools for microbiome data
library(vegan) # computes our diversity analyses
library(ggplot2) # R graphing package
library(taxonomizr) 
library(tidyverse)
```

For those of you new to R, it is a programming language that is often used by computational biologists. A lot of bioinformatics tools, including those useful for microbiome data, can be accessed here. Unlike command line, you can write out multiple lines of commands and execute them whenever you wish. Feel free to ask about anything confusing, but some of the key ideas are here:
- `library()` loads packages that are not included in base R. We need to load these packages to access their commands.
- At times you will see `x <- y` or `x = y`. Both indicate that a variable is being assigned, the syntax is interchangable.
- You will also see the pipe (`%>%`) show up often. This is a useful bit of programming grammar that "pipelines" the output of one command into the next. 
	- For example, `colnames(df) %>% table()` tells R to take the output of `colnames(df)` (the column names of object `df`) as input for the command `table()` (which counts the unique entries of the input). This would be the same as typing `table(colnames(df))`, but appears much more organized, and makes it easier when linking more than two commands. 

We'll use phyloseq to organize our data. phyloseq comes with a function to import our newly created .biom file as a phyloseq object, with which we can do the rest of our analyses for this section. 
```r
data<-import_biom("mgs.biom", header = TRUE)
```

Check out our taxonomic annotations here.
```r 
View(data@tax_table@.Data) 
```

We can see that kraken2, along with many other tools, likes to prepend taxonomic names with an indicator of the level they are at (ex. "g__" for genus). Since phyloseq organizes these names into the correct hierarchy for us, we can remove these, and rename the columns of the tax.table.
```r
data@tax_table@.Data <- substring(data@tax_table@.Data, 4) # removing the prefix
colnames(data@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # renaming the columns
bacteria<-subset_taxa(data, Kingdom == "Bacteria") # taking a bacteria-only subset of our original data 
```

We can examine the difference between these two phyloseq objects by checking which kingdoms are present. The original data shows that, in addition to one Eukaryota taxon, we also had some Archaea and Viruses.
```r
# original data
data@tax_table@.Data[,1] %>% table() # table() is a useful function for tabulating the unique factors in a list

# new subset data
bacteria@tax_table@.Data[,1] %>% table()
```

> Question 1.3
> 
> How many different phyla of bacteria are in our data? 
> **HINT:** the function `unique()` can shorten a list to unique elements and `length()` will return the length of a list. 

### Diversity Analyses
One of the most common metrics with which to examine your microbiome data is through alpha and beta diversity measures. As a reminder, alpha diversity can be measured in several ways, but primarily concerns the within-sample diversity of your data. 
```r
plot_richness(physeq = bacteria, x = "Athlete_type", color = "Athlete_type", 
			  measures = c("Observed", "Shannon")) + 
  geom_boxplot() + # telling gpglot to make a boxplot
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + # rotating x-axis labels
  theme_classic() + # changing overall look of chart
  scale_colour_manual(values = c("#B5FF66", "#FC7B56", "#567AFC")) # changing colour of boxes
```

Here we have two alpha diversity indices: 
- *Observed*
	- Simply represents the number of observed species.
- *Shannon*
	- Takes into account both the number of observed species and how evenly distributed these species are.

We can also change the comparison we want to make. Try looking at `bacteria@sam_data` to see which other columns you could compare.

> Question 1.4
> 
> The ultramarathoner Shannon index is noticeably lower than its Observed index. What could this indicate for alpha diversity in our ultramarathoner samples?

Next we can look at beta diversity. Where alpha diversity examined the within-sample diversity, beta diversity computes a distance matrix between all samples to convey *between*-sample diversity. 

Just like alpha diversity, we can measure beta diversity in different ways. We will go with Bray-Curtis distances, which calculate your distances according to the abundance of your taxa. First we calculate these distances then we project them as a PCoA plot, which resembles a scatter plot.  Samples that are closer to each other on this chart are thus closer to each other in terms of the structure of their taxa. 
```r
ordination = ordinate(physeq = bacteria, method = "PCoA", distance = "bray") # calculating distance matrix

plot_ordination(physeq = bacteria, ordination = ordination, color = "Athlete_type") +
  geom_point(size=10, alpha=0.1) + geom_point(size=5) + # telling ggplot to plot a scatterplot
  theme(text = element_text(size = 20)) + # changing font size
  theme_bw() + # changing overall look of our plot
  scale_colour_manual(values = c("#B5FF66", "#FC7B56", "#567AFC")) # changing dot colours
```

Again you can try changing out "Athlete_type" for a different grouping. Notice how changing these points don't actually change how the dots are grouped. 

> Question 1.5
> 
> What does the clustering of these samples by athlete tell you about your samples?

### Visualizing Taxa
Finally we can visualize the relative abundances of these taxa in each sample. Phyloseq comes with a useful function `tax_glom()`, that lets you group up everything by a given taxonomic rank. We will visualize at the family level to reduce the number of features. 
```r
bacteria_family <- tax_glom(rarefied, taxrank = "Family")
bacteria_family = psmelt(bacteria_family) 
```

`pmelt()` changes the format of our data to one that is more friendly to ggplot2, a versatile graphing library available for R. Let's look at the format of this data now using `View(bacteria_family)`.
```r
bacteria_family$Family %>% unique %>% length
```

However we still have nearly 500 families in this data. Let's limit to those families that show up most often. We'll set an arbitrary threshold of 2500 counts.
```r
bacteria_family$Family[bacteria_family$Abundance < 2500] = "Other"
```

This leaves us with 25 families. Before we move onto plotting the taxa, we'll have to make a colour palette that works with that many taxa. One way to do this is with the use of the RColorBrewer package. Here we expand a pre-set palette from RColorBrewer into a larger palette with the stated number of colours. You can check using `display.brewer.all()` for other palettes
```r
family_palette = colorRampPalette(brewer.pal(12, "Set3"))(25)

# And finally, we can visualize using ggplot
ggplot(data=bacteria_family, aes(x = Sample, y = Abundance, fill = Family))+
  xlab("Sample ID") + # X-axis title
  theme_classic() + # overall look of our graph
  ylab("Abundance (%)")+ # Y-axis title
  geom_bar(stat="identity", position="fill") + # telling ggplot to make a bar chart
  scale_fill_manual(values = family_palette) + # setting the colours of the bars
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + # rotating x axis labels
  facet_wrap(~Athlete_type, scales = "free") # separating our bars by athlete type
```

Try to play around with different taxonomic levels and palettes! Remember to recreate the palette and the threshold of how many counts below which the taxa will be set as "Other". Keep in mind, a higher number of taxa makes it difficult to distinguish between colours in your final chart. 

### Answers
- **Question 1.1 - Look at the *Per base sequence content* of some of these reports. What do you notice?**
	- The start of the file seems to have imbalanced or biased. It's normal to have a bit of an imbalance at the beginning, especially for RNA-seq data (though this is DNA). If this was imbalanced throughout the whole sequence, it would be worrying. 
- **Question 1.2 - What are the dangers of using these raw reads in a real experiment?**
	- Reads that are actually human might be misidentified as bacterial reads, giving us false positive hits and complicating the assembly of our bacterial genomes later on. 
- **Question 1.3 - How many different phyla of bacteria are in our data?** 
	- 43 unique phyla. You can use this code to count this number of phyla - remember to change the index to the second column.
	- `bacteria@tax_table@.Data[,2] %>% unique() %>% length()`
	- OR
	- `length(get_taxa_unique(bacteria, taxonomic.rank = "Phylum"))`
- **Question 1.4 - The ultramarathoner Shannon index is noticeably lower than its Observed index. What could this indicate for alpha diversity in our ultramarathoner samples?**
	- Since Shannon index takes evenness into account and Observed does not, this indicates that the evenness of our ultramarathoner samples are low despite having a relatively high abundance of species. 
- **Question 1.5 - What does the clustering of these samples by athlete tell you about your samples?**
	- For the most part, the samples of each athlete type seem to cluster with each other, indicating that samples of each athlete type are closer to each other than other athletes. Note there is one cyclist sample that seemingly does not resemble other cyclists, and actually seems to be closer to the ultramarathoner samples.   

## Lab 2 - MAG Assembly and Annotation
One of the strengths of shotgun metagenomics analyses over marker gene analyses is the ability to examine the entire genomes of your microbes. To do this, we assemble reads into metagenome assembled genomes, or MAGs. This lab will take you through this process, and provide a couple of options for annotating these MAGs so that we can look at things at a functional level. 

### MAG Assembly
The first step is to arrange our reads into contigs, which are longer, contiguous segments of DNA. We'll use a tool called MEGAHIT to do this. In this tutorial, we are co-assembling our contigs, meaning the contigs will be formed from sequences from all samples rather than contigs being unique to an individual sample. While there are legitimate reasons to do this with real data, our athlete data is collected from public repositories and were taken from different people from different parts of the world. In reality, we would not want to co-assemble these reads, but for the sake of time and ease of analysis, we will do so here. 

>Question 2.1
>
>When would we want to co-assemble our contigs?

We'll use tmux to run megahit in the background, as it is a very time-consuming step. 
```bash
# Creating a new tmux environment
tmux new-session -A -s anvio
R1=$( ls athlete_samples/*1.fastq | tr '\n' ',' | sed 's/,$//' )
R2=$( ls athlete_samples/*2.fastq | tr '\n' ',' | sed 's/,$//' )
megahit -1 $R1 \
         -2 $R2 \
         --min-contig-len 1000 \ # minimum length of contigs
         --num-cpu-threads 4 \ # number of threads to use 
         --presets meta-large \ # an option for assembling large and complex MAGs
         --memory 0.8 \ # this means megahit will use 80% of the memory available
         -o anvio/megahit_out \ # output folder
        --verbose # verbose is a common option in a lot of functions to print out additional information about what the command is actually doing - useful for debugging!
        
# remember that to exit tmux, hit CTRL+B and then D
```

In the meantime, we have pre-computed the output files for you to move onto further steps. They are stored at `anvio/megahit_out_precomputed`. Use `grep -c ">" anvio/megahit_out_precomputed/final.contigs.fa | wc -l` to check how many contigs are in our final assembly. Don't forget the quotation marks around the `>` or you will overwrite the file!

You can check on your progress at any time using `tmux attach`. Because the command was sent with the `--verbose` option, you can read the output to see what MEGAHIT is doing. 

Next we'll run a quick formatting step on our final.contigs.fa before integrating it into a database. This database will serve as a jumping off point for the next few anvio tools, and stores information like the position of ORFs and any annotations. 
```bash
# Formatting 
anvi-script-reformat-fasta anvio/megahit_out_precomputed/final.contigs.fa \
	                       --simplify-names \ # reducing descriptive names
	                       --min-len 2000 \ # further filtering our contigs to this length
	                       -o anvio/final.contigs.fixed.fa # our output file

# Building our contig database 
mkdir anvio/anvio_databases
anvi-gen-contigs-database -f anvio/final.contigs.fixed.fa \ # input
                              -o anvio/anvio_databases/contigs.db \ # output
                              -n athlete_db # name of the database
```

You'll note that the output for this step mentions Prodigal. Prodigal is a tool for recognizing prokaryotic genes, and anvio is using it here to identify ORFs. 

>Question 2.2
>
>How many contigs are in our database? 

### Identifying Contigs
Now that we have our contigs assembled and organized into a database, we can start to identify genes are here. First we will use hidden Markov models (HMMs) that come with anvio. HMMs are probabilistic models that are used often in computational biology, with a broad variety of applications. For our purposes, HMMs can be used to determine the probability that a particular gene exists within our contigs. 
```bash
anvi-run-hmms -c anvio/anvio_databases/contigs.db \ # testing anvio HMMs on our contigs
                              --num-threads 4
```

We can then export the results of these HMMs into a separate fasta file.
```bash
anvi-get-sequences-for-gene-calls -c anvio/anvio_databases/contigs.db \
                              -o anvio/gene_calls.fa
```

>Question 2.3
>
>How many genes did anvio identify?
>**HINT:** we can use the same command we used earlier to find the number of contigs in our `final_contigs.fa` file.

### Quantifying Contigs
The next big step of a metagenomics workflow is binning. This is the act of sorting our contigs into genomes, or "bins" that represent a taxon. Before we get to the actual binning step, we have to count how abundant these contigs are among our samples. We will use Bowtie2 to do this, a tool that groups contigs with similar abundances into a bin. This first step creates a database for our formatted contigs. 
```bash
mkdir anvio/bowtie2_output
bowtie2-build anvio/final.contigs.fixed.fa anvio/bowtie2_output/final.contigs.fixed
# the two parameters here are the input and the prefix and filepath for outputs
```

Now that these bins are sorted, we can run this loop. It looks complicated but it basically performs 4 functions for each of our samples:
1. Performs bowtie mapping to get a SAM file
2. Turns that SAM file into a BAM file, which we want for downstream tools
3. Sort and index that BAM file to make it easier to search
4. Remove intermediate files
```bash
mkdir anvio/bam_files

for SAMPLE in `awk '{print $1}' sample_list.txt`
do

    # 1. do the bowtie mapping to get the SAM file:
    bowtie2 --threads 4 \
            -x anvio/bowtie2_output/final.contigs.fixed \
            -1 "athlete_samples/"$SAMPLE"_1.fastq" \
            -2 "athlete_samples/"$SAMPLE"_2.fastq" \
            --no-unal \
            -S anvio/bam_files/$SAMPLE.sam

    # 2. convert the resulting SAM file to a BAM file:
    samtools view -F 4 -bS anvio/bam_files/$SAMPLE.sam > anvio/bam_files/$SAMPLE-RAW.bam
    
    # 3. sort and index the BAM file:
    anvi-init-bam anvio/bam_files/$SAMPLE-RAW.bam -o anvio/bam_files/$SAMPLE.bam
    
    # 4. remove the intermediate BAM file that is no longer needed:
    rm -f anvio/bam_files/$SAMPLE-RAW.bam

done
```

To organize the information we've added to our contigs, we can make anvio profiles, then merge them into one profile for binning.
```bash
# Making profiles
mkdir anvio/profiles
for SAMPLE in `awk '{print $1}' sample_list.txt`
do

    anvi-profile -c anvio/anvio_databases/contigs.db \
                 -i anvio/bam_files/$SAMPLE.bam \
                 --num-threads 4 \
                 -o anvio/profiles/$SAMPLE
done

# Merging profiles
anvi-merge -c anvio/anvio_databases/contigs.db \
           -o anvio/profiles/merged_profiles \
           anvio/profiles/*/PROFILE.db
```

As a checkpoint, we can look at some key statistics of the contigs so far.
```bash
anvi-display-contigs-stats --report-as-text \
                           --output-file contigs_stats.txt \
                           anvio/anvio_databases/contigs.db
```

Remember we can read the resulting file with `cat contigs_stats.txt`. While a lot of these fields are self-explanatory, there are a few that are less so:
- L50, L75, L90: The number of contigs you would need to go through before reaching X% (ex. L50 is 50%) of your database
- N50, N75, N90: Like above, but if your contigs were sorted by length (longest to shortest) and you wanted to know the length of the contig at X% (ex. the contig at 75% would be 3039 nts long)
- Bacteria_71: HMMs for 71 single-copy core bacterial genes 
- Archaea_76: HMMs for 76 single-copy core archaeal genes 
- Protista_83: HMMs for 83 single-copy core protist genes
- Ribosomal HMMs
- Number of genomes detected, as predicted by anvio (we have 17 bacterial genomes)

>Question 2.4
>
>How many contigs would we need to look through before getting 90% of the way through the database?

### Binning
Finally, we can cluster our contigs into bins. Anvio uses CONCOCT for this, which takes into account both sequence composition and abundance across samples (which we get from bowtie!). This step also takes a while, so we will provide the pre-computed `PROFILE.db` at `anvio/precomputed_profiles/PROFILE.db`.
```bash
# do not run!
	anvi-cluster-contigs -c anvio/anvio_databases/contigs.db \ # contig database
                         -p anvio/profiles/merged_profiles/PROFILE.db \ # profile
                         -C "merged_concoct_2000" \ # name of the bins
                         --driver CONCOCT \ # clustering tool
                         --length-threshold 2000 \ # minimum contig length
                         --num-threads 4 \ # number of threads
                         --just-do-it # ignoring warnings
```

We want to examine the quality of these bins before moving ahead, so let's generate summaries of our CONCOCT output.
```bash
mkdir anvio/concoct_summary_2000
anvi-summarize -c anvio/anvio_databases/contigs.db \ # contig database
                   -p anvio/precomputed_profiles/PROFILE.db \ # profile database
                   -C "merged_concoct_2000" \ # bin names
                   -o anvio/concoct_summary_2000/summary/ # output folder
```

To open some of these summaries, recall we can use `less`
```bash
less anvio/concoct_summary_2000/summary/bins_summary.txt # overview of bin metrics
less anvio/concoct_summary_2000/summary/bins_across_samples/abundance.txt # bin abundances over samples
```

Bins are usually assessed by two key properties: Completeness and Redundancy. Completeness refers to how much coverage of the expected genome has been achieved in your bin. Redundancy refers to how many times we see duplicates of a gene that is expected to only occur once. You'll notice that many of these bins have 0% completeness and redundancy - this is probably a single contig that was unable to be sorted elsewhere. 

Normally we would move onto a bin refinement step - this allows us to try multiple binning methods to find the optimal makeup of our bins, ideally improving completeness and redundancy. For the sake of time we will skip this, and use our original CONCOCT bins. But this step might look like:
```bash
# do not run!
anvi-refine -c anvio/anvio_databases/contigs.db \
            -p anvio/profiles/merged_profiles/PROFILE.db \
            -C "merged_concoct_2000" \
            -b Bin_17 

```

With our bins formed, we want to export everything to MAGs. We'll filter our bins first, to have at least 50% completeness and at most 10% redundancy. 
```bash
anvi-rename-bins -c anvio/anvio_databases/contigs.db \ # contig database
                     -p anvio/profiles/merged_profiles/PROFILE.db \ # profile database
                     --collection-to-read merged_concoct_2000 \ # concoct bins
                     --collection-to-write athlete_mags \ # new name of bins
                     --call-MAGs \ # renaming bins to mags
                     --min-completion-for-MAG 50 \ # completion threshold
                     --max-redundancy-for-MAG 10 \ # redundancy threshold
                     --prefix athlete_db \ # name of our contigs
                     --report-file renaming_bins.txt # output file
```

You might want to implement stricter parameters for your own data if you expect lower quality data, but it is not recommended to allow bins with less than 50% completeness or over 10% redundancy. Bins that do not meet these requirements should either be refined with a different binning method, or removed from further analysis. 

Now we'll summarize these MAGs like we did before. 
```bash
anvi-summarize -c anvio/anvio_databases/contigs.db \
                   -p anvio/profiles/merged_profiles/PROFILE.db \
                   -C "athlete_mags" \
                   -o anvio/final_mags_summary/
```

And to examine the summaries:
```bash
ls anvio/final_mags_summary/bin_by_bin/ # bins and mags
cat anvio/final_mags_summary/bins_summary.txt # summary of both bins and mags
```

We can see that our original bins are included here, along with our newly renamed MAGs. We only want the MAGs from this point on, so let us move the desired fastas to a new folder.
```bash
mkdir MAG_fasta
cp anvio/final_mags_summary/bin_by_bin/*MAG*/*contigs.fa MAG_fasta/
ls MAG_fasta/
```

And now we have our MAGs! This was a long and involved process, but this is where your workflow really opens up; you can annotate your contigs with any number of databases depending on your particular aims. We can use more HMMs, sequence search tools, or any number of methods for associating our contigs with biological meaning. The first thing we will do is use CheckM to assess quality of genomes and assign taxonomy.

### Annotating Bins with Taxonomy
We want to point CheckM to its database, which we can do with this command:
```bash
export CHECKM_DATA_PATH=/scratch/j/jparkin/rchieu/pakistan_workshop/databases/checkm

# Then, we run the lineage workflow command
checkm lineage_wf --reduced_tree -t 4 -x fa MAG_fasta MAGs-CHECKM-lineage -f MAGs-CHECKM.txt --tab_table
```
 
 The lineage workflow command will go through a few steps with our MAGs, producing output that records which lineages each MAG has been assigned and how that lineage was determined. The columns include:
 - genomes: the number of genomes used to identify this lineage
 - markers: the number of marker genes within a lineage used 
 - marker sets: number of marker sets (groups of markers associated with a specific lineage) within the lineage's marker set used
 - The numerical columns are used to assess redundancy, and represent the number of repeat marker genes found. Ex. an 8 under column 3 means 8 marker genes were found 3 times.
 - Completeness and redundancy
 - Strain heterogeneity, which attempts to explain where your redundancy could come from. High strain heterogeneity might suggest those additional marker gene copies come from more divergent species. Low strain heterogeneity would suggest the opposite. 

>Question 2.5
>
>What does a strain heterogeneity of 50 indicate?

You'll also notice that many of these lineages are at quite a high level. We can run `tree_qa` to assign more specific taxonomy. 
```bash
checkm tree_qa MAGs-CHECKM-lineage -f MAGs-CHECKM-tax.txt --tab_table
```

We can see that a lot of these bins have acquired genus-level annotations! Not all of these MAGs are at a lower level of taxonomy. This is expected, and just indicates that there is not enough information to be more specific. This is a good stopping point for now - in our final lab, we will cover additional methods of functional annotation and explore ways of visualizing this data.

### Answers
- **Question 2.1 - When would we want to co-assemble our contigs?**
	- Co-assembly is valid when your samples come from very similar sources and went through the same DNA extraction/ sequencing steps. For example, gut samples from laboratory mice. 
- **Question 2.2 - How many contigs are in our database?** 
	- According to the contigs DB output, we have 12,802 contigs in total. 
- **Question 2.3 - How many genes did anvio identify?**
	- Using `grep '>' anvio/gene_calls.fa | wc -l`, we see that we've found 66779 genes in our data. 
- **Question 2.4 - How many contigs would we need to look through before getting 90% of the way through the database?**
	- According to `contigs_stats.txt`, our L90 is 9917.
- **Question 2.5 - What does a strain heterogeneity of 50 indicate?**
	- A middling strain heterogenity might suggest that the contaminant marker gene copies come from a group of moderately divergent species.


## Lab 3 - Functional Annotation and Visualization
Now that we have MAGs, we can move onto functional annotation! We'll be using tools to link our sequences with genes and functions, and then covering a few ways you might want to visualize this data. In your own workflow, you might want to run these tools on your contigs, but we will run them on our raw sequences as our contigs are not currently split up between samples. Before we move onto more specific databases of function, it is a good idea to start with gene annotations. We'll do this with a tool called MMSeqs2. 

### MMSeqs2
MMSeqs (Many-against-Many sequence searching) is a tool we will use to search through our reads and annotate them using the UniRef database. UniRef is a popular, well-maintained resource for protein sequences and functional information. We'll be using a smaller version of UniRef, UniClust30, which groups sequences together at 30% sequence identity. 

Our first step will be to create databases for our samples.
```bash
mkdir mmseqs_out
parallel -j 4 --progress 'mmseqs createdb {} mmseqs_out/mmseqs-{/.}-queryDB' ::: cat_reads_full/*

cd databases/uniref
mmseqs createdb demo_uniref.fasta demo_uniref
mmseqs createindex demo_uniref tmp
cd ../..
```

Then we will use MMSeqs2 to look though our database. You can attempt to run this - if it takes too long or you run into issues, skip past the next step and use the pre-computed data at 
`mmseqs_out/mmseqs_m8s/precomputed/`
```bash
parallel -j 1 --progress 'mmseqs search mmseqs_out/mmseqs-{/.}-queryDB databases/uniref/demo_uniref mmseqs_out/mmseqs-{/.}-resultDB databases/uniref/tmp --db-load-mode 3 --threads 2 --max-seqs 25 -s 1 -a -e 1e-5' ::: cat_reads_full/*
```

The parameters here are as follows:
- query sequences
- target sequences
- output
- index for target sequences
- `--db-load-mode 3`  technical option that tells mmseq how to load your databases
- `--threads 2` number of threads used
- `--max-seqs 25` maximum hits per sequence
- `-s 1` sensitivity - higher values will result in more hits but longer runtime
- `-a` a formatting option that makes reading the filein the next step easier
- ``-e 1e-5`` keep a threshold of below 1e-5 e-value 

Now we format the result files into one that is nicer to work with.
```bash
parallel -j 2 --progress 'mmseqs convertalis mmseqs_out/mmseqs-{/.}-queryDB databases/uniref/demo_uniref mmseqs_out/mmseqs-{/.}-resultDB mmseqs_out/mmseqs-{/.}-s1.m8 --db-load-mode 2' ::: cat_reads_full/*

# moving them to a new folder
mkdir mmseqs_m8s
mv mmseqs_out/*.m8 mmseqs_m8s
```

Looking at these files, we can see a table with the original sample ID and read number next to a UniRef ID. The rest of the columns correspond to:

Column Number | Data Type
-|-
0|query sequence ID
1|Subject (database) sequence ID
2|Percent Identity
3|Alignment Length
4|Number of gaps
5|Number of mismatches
6|Start on the query sequence
7|End on the query sequence
8|Start on the database sequence
9|End on the database sequence
10|E value - the expectation that this alignment is random given the length of the sequence and length of the database
11|bit score - the score of the alignment itself

Notice that some of the reads have multiple, even many UniRef hits. We just want the top hit, so we'll run a script to pick those out for us.
```bash
mkdir mmseqs_m8s/top
python scripts/pick_uniref_top_hit.py --unirefm8Dir mmseqs_m8s --output_path mmseqs_m8s/top
```

These new files no longer have multiple hits for each entry, and also removed the other columns. The UniRef IDs are useful in getting additional annotations as well. Let's look at one sample in particular and visualize which gene ontology (GO) terms we can find within.

### GO Term Map
We'll use the sample SRR2992927. We used the [UniProt ID Mapper](https://www.uniprot.org/id-mapping) to retrieve GO terms, which should act as more interesting units of function than genes. 
```R
library(readr)
library(stringr)
library(data.table)
library(tidyverse)
library(readxl)
library(rrvgo)
library(org.EcK12.eg.db)

# our mmseq top hit file
go_SRR2992927 <- read_delim("mmseqs_m8s/top/mmseqs-SRR2992927-s1.m8-parsed.txt", 
    delim = "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE)
# Renaming columns and making one that maps to idmapper
colnames(test_sample_SRR2992927) = c("Read", "UniRef")
test_sample_SRR2992927$UniRef_Code = test_sample_SRR2992927$UniRef %>% substr(10, 15)

# if you wanted to run the id mapper for yourself, you could export the IDs using this:
# test_sample_SRR2992927$UniRef_Code %>% as.list %>% fwrite("uniprot.txt")
# then upload the file to the link above

# output from uniprot id mapper
idmapping_SRR2992927 <- read_delim("precomputed_files/idmapping_SRR2992927.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
```

Now we're going to do some data manipulation to fit this data for our visualizing packages, `rrvgo`.
```R
# combining the m8 file and idmapper results
test_sample_SRR2992927$GO_Bio = idmapping_SRR2992927$`Gene Ontology (biological process)`[match(test_sample_SRR2992927$UniRef_Code, idmapping_SRR2992927$From)]
# we are just using the go terms that fall under biological process, as it makes the most sense for bacteria

# Now we'll make a new data frame of just Reads:Go_Term
ts_df = test_sample_SRR2992927[, c("Read", "GO_Bio")]
ts_df$GO_Bio = str_split(ts_df$GO_Bio, "; ")
ts_df = unnest(ts_df, cols = GO_Bio) %>%  na.omit()

# slicing off the term description, won't need it
temp = ts_df$GO_Bio %>% 
    str_match("\\[\\s*(.*?)\\s*\\]") 
ts_df$GO_Bio = temp[,2]
```

Our new data frame may be minimal, but we can take all the parameters we need from it. `rrvgo` requires the GO terms, a genome annotation with which it will look at relationships between GO terms, and the abundance of the GO terms in your data. For the second parameter, we will use the *E. coli* genome annotation, as there is not a microbiome wide file we can use. 
```R
# And finally put this all into a treemap as visualized by rrvgo
treemap_df = ts_df$GO_Bio %>% calculateSimMatrix(orgdb = "org.EcK12.eg.db", ont = "BP", method = "Resnik") %>%
  reduceSimMatrix(table(ts_df$GO_Bio), 
                  threshold=0.85, orgdb="org.EcK12.eg.db")
treemapPlot(treemap_df, size = "size")
```

This last function sorts our GO terms into "parent" terms based on how similar they are, then takes into consideration how many times that term pops up. The final output creates a map of these parent terms ranked by abundance! 

>Question 3.1
>
>Which collection of GO terms was most common in our data?

### Using the Resistance Gene Identifier (RGI)
The second tool we will use is the Resistance Gene Identifier, or RGI. RGI is used to find genes related to antibiotic resistance in our sequences - this can be helpful especially for clinical data. The RGI relies on CARD, the Comprehensive Antibiotic Resistance Database.
```bash
# First we load the database.
rgi load --card_json databases/card/card.json --local

# Then we grab the annotations from the database
rgi card_annotation -i localDB/card.json > card_annotation.log 2>&1
rgi wildcard_annotation -i databases/card/wildcard --card_json localDB/card.json -v 3.3.0 > wildcard_annotation.log 2>&1

# Moving the annotations to the localDB folder for organization
mv card_* localDB/
mv wildcard_* localDB/

# Now loading the whole database
rgi load \
  --card_json databases/card/card.json \
  --debug --local \
  --card_annotation localDB/card_database_v3.3.0.fasta \
  --wildcard_annotation localDB/wildcard_database_v3.3.0.fasta \
  --wildcard_index databases/card/wildcard/index-for-model-sequences.txt \
  --wildcard_version 3.2.9 \
  --amr_kmers databases/card/wildcard/all_amr_61mers.txt \
  --kmer_database databases/card/wildcard/61_kmer_db.json \
  --kmer_size 61
```

Now that the CARD database has been properly loaded into RGI, we can finally run the tool itself. There are two commands we're running here, `main` and `bwt`. `main` acts on a genomic level, while `bwt` analyzes metagenomics reads specifically. 
```bash
mkdir card_out_main
parallel -j 1 'rgi main -i {} -o card_out_main/{/.} -t contig -a DIAMOND -n 4 --include_loose --local --clean' ::: MAG_fasta/*

mkdir card_out_bwt
parallel -j 1 --xapply 'rgi bwt --read_one {1} --read_two {2} --aligner bowtie2 --output_file card_out_bwt/{1/.}_CARD --threads 2 --local --clean' \
 ::: athlete_samples/*_1.fastq ::: athlete_samples/*_2.fastq
```

To visualize the `main` results:
```bash
rgi heatmap --input card_out/ --output card_heatmap
```

The tool helpfully outputs instructions on reading this heatmap:
>*Output file card_heatmap-12: Yellow represents a perfect hit, teal represents a strict hit, purple represents no hit.*

>Question 3.2
>
>Which antibiotic resistance gene occurs most commonly among our MAGs?

The BWT output will take a bit more work to visualize. We'll be doing some data manipulation first with a python script, then in R.
```bash
mkdir bwt_maps
mv card_out_bwt/*gene_mapping* bwt_maps/
python scripts/bwt_formatter.py
```

This outputs a file, all_bwt.npy, that contains all the CARD hits by sample. We can then import this file into R so that we can use ggplot again.
```R
library(stringr)
library(tidyverse)
library(ggplot2)
library(readr)

all_bwt <- read_csv("all_bwt.csv")
colnames(all_bwt)[1] = "AMR_Gene" # python did not export the first column name

# reading in our metadata for later
sample_metadata <- read_delim("sample_metadata.tsv", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
```

This time, we will make a heatmap. Start by formatting our data frame into the correct form for ggplot.
```R
all_bwt = all_bwt %>% pivot_longer(cols = - AMR_Gene, names_to = "Sample", values_to = "Hits")

# some of these names are really long. let's cut those off
all_bwt$Short_Names = paste0(substr(all_bwt$AMR_Gene, start = 1, stop = 35), ifelse(sapply(all_bwt$AMR_Gene, nchar) > 35, "...", ""))

# now adding in our metadata (ex. athlete type)
all_bwt = cbind(all_bwt, sample_metadata[match(all_bwt$Sample, sample_metadata$Sample_name), c("Athlete_type", "Sex", "Diet")])
```

And finally visualizing! 
```R
ggplot(all_bwt, aes(Sample, Short_Names, fill= log2(Hits))) + 
  geom_tile() + # telling ggplot to make a heatmap
  scale_fill_viridis_c(option = "mako") + # colour scheme
  theme_classic() + # changing overall design of plot
  facet_grid(~Athlete_type, scales = "free_x", space = "free") # sorting heatmap by athlete type
```

>Question 3.3
>
>How do these results compare to the `rgi main` results?

### Cytoscape
Depending on your goals, you might want to visualize how your genes lie on a particular pathway. For example, you might be interested in how your microbes digest amino acids, so you could map the abundance of related genes onto a tryptophan metabolism pathway. This is where Cytoscape comes in. Cytoscape is a program that allows for visualization of networks, and we can use it to illustrate biological pathways using microbiome data. An important resource for viewing these pathways is the Kyoto Encyclopedia of Genes and Genomes (KEGG), which contains information on pathways and the enzymes contained within. 

For this portion of the workshop, we will focus on operating Cytoscape itself. We've provided metatranscriptomics data that has been prepared for this program. This data comes from stool samples taken from mice. For metagenomics data, you could use the [eggNOG mapper](http://eggnog-mapper.embl.de/) to find enzymes, and the [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html) to map these to pathways. We'll also be using the Cytoscape plugins `enhancedGraphics` and `KEGGScape` for this. 

To import an XML into KEGG
- Select `File` -> `Import` -> `Network` -> `File...`
- Select the XML file, `ec00010.xml` or `ec00500.xml` and click `Open`
- Check `Import pathway details from KEGG Database` box then select `OK`

Then, we can import some of the metadata associated with this data. 
- Navigate to the Node Table panel (bottom right) and click on the Import Table from file icon (arrow pointing into document)
- Select the `mouse1_cytoscape.txt` file and click `Open`
- Change the `Import Data as` from `shared name` to `KEGG_NODE_LABEL`
- Click OK

To change the appearance of your network:
- In the left `Control Panel` select the `Style` tab
- Check the `Lock node width and height` box
- Click the left-most box by the `Size` panel and change the default node size to 20.0
- Click the blank box immediately to the right of the box you clicked to change the default size, change the `Column` field to `RPKM` and the `Mapping Type` field to `Continuous Mapping`
- Click the left-most box by the `Image/Chart 1` panel, switch to the `Charts` tab, Click the doughnut ring icon, and press the `>>` "add all" button between the two column fields before clicking apply (make sure to remove overall RPKM from the fields that are added to the doughnut ring)
- If you do not see the `Image/Chart 1` panel, select `Properties` -> `Paint` -> `Custom Paint 1` -> `Image/Chart 1` from the to left corner of the control panel
- To improve the visualization you can modify colour properties under `Image/Chart 1` -> `Charts` -> `Options`, or modify other properties such as Label Font Size, Label Position, Fill Color, Node location, and edge properties**

And you're done! Try visualizing the other `.xml` file included in this workshop. Cytoscape networks are useful for showing complicated relationships in your data overlaid upon a pathway of interest. For our athlete data, it might be interesting to visualize differences in metabolic pathways between the different athletes. If you have any questions on how to extend the tools you've learned today or how they might change on real datasets, feel free to ask!

### Answers
- **Question 3.1 - Which collection of GO terms was most common in our data?**
	- 'de novo' pyrimidine nucleobase biosynthetic process. This is just the creation of pyrimidines. 
- **Question 3.2 - Which antibiotic resistance gene occurs most commonly among our MAGs?**
	- vanT gene in vanG cluster appears to be most common. vanT is a membrane-bound protein involved in vancomycin resistance.
- **Question 3.3 - How do these results compare to the `rgi main` results?**
	- They're quite different! `rgi main` and `rgi bwt` use two different methods to analyze two different things. Since `rgi bwt` was run on our raw reads and `rgi main` on our MAGs, there are likely a lot of reads parsed by `bwt` that were not seen by `main`. 

## License
Workshop materials are licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.


