# 16S Workshop
## Introduction
Welcome to the data analysis portion of the workshops! This tutorial will run through the steps of processing 16S rRNA microbiome data, from raw reads to useable files that you will be able to perform analyses on. 16S studies are a great entry point to the microbiome, as samples are relatively cheap to sequence and the data processing and analysis is often uncomplicated. 

We'll be using 10 samples sequenced from chicken cecal content. These samples were originally studied for differences in antibiotic usage, and were sequenced using the V4/V5 hypervariable region of the 16S rRNA gene. You can read more about this study [here](https://doi.org/10.1186/s42523-021-00157-6). 

This workshop was adapted from: [Processing Raw 16S rRNA Data (QIIME 2 Basics)](https://github.com/ParkinsonLab/16S_rRNA_analysis_tutorials/wiki/Processing-Raw-16S-rRNA-Data-(QIIME-2-Basics)) and the [2024 Canadian Bioinformatics Workshops Beginner Microbiome Analysis](https://bioinformaticsdotca.github.io/BMB_2024_module1). 

## Lab 2 - QIIME Part 1
The next two workshops will focus on processing 16S samples using QIIME 2. QIIME is a platform for microbiome data that contains all the tools for processing 16S reads, and can run on the command line. 

As mentioned in the introduction, we'll be working with 10 chicken gut microbiome samples with this workshop. Let's list them in their directory.
```bash
ls chicken_samples/
```

You can see here 20 fastq samples - R1 and R2 denoting the two pairs for one gut sample. There's also a metadata file, which you can view:
```bash
less chicken_samples/metadata.tsv

# Now let's move that out of the sample folder so that it doesn't get caught up in later steps
mv chicken_samples/metadata.tsv .
```

Our metadata file is formatted to work with the QIIME2 workflow so now additional commands necessary, and we can begin to process our sequences. The first step will be to perform a quality check to see whether there are any problems with our samples. fastqc is a tool designed to do this. 
```bash
# First we make a folder for the output
mkdir fastqc_out

# Then we run the command on all our files
fastqc -t 4 chicken_samples/*.fastq.gz -o fastqc_out/
```

Fastqc outputs `.html` files that can be opened and viewed on a web browser. Open one up and you will see a variety of metrics that can be used to assess the quality of your reads. For 16S data, we are primarily concerned with "Per base sequence quality".

>Question 2.1
>
>Some of the other sections of the fastqc report, like "Per base sequence content" shows up as a warning. Why are these sections not as important for our data?

### Importing into QIIME
QIIME stores your data in a `.qza` file, which keeps track of your files as they run through each step. To create this file, we can import our raw reads into QIIME.
```bash
mkdir qza_files
# This command imports our paired end fastqs qza file
qiime tools import \
   --type SampleData[PairedEndSequencesWithQuality] \ 
   --input-path chicken_samples/ \
   --output-path qza_files/raw_chicken.qza \
   --input-format CasavaOneEightSingleLanePerSampleDirFmt
# This last confusing looking parameter describes the format of our samples as paired end `fastq.gz` files with _R1_ and _R2_ as the forward and reverse reads

# Make directory for visualization files
mkdir qzv_files
```
 
 Throughout this tutorial we will also be producing visualizations of our data in the form of`.qzv` files. Both `.qza` and `qzv` files can be viewed with [QIIME2 View](https://view.qiime2.org/). We encourage you to view each output to see how your data is changing as it goes along?

Another important quality control step is to remove the adapters that were sequenced as part of the 16S sequencing protocol. This can be done by using cutadapt. 
```bash
# Trim reads
qiime cutadapt trim-paired \
   --i-demultiplexed-sequences qza_files/raw_chicken.qza \
   --p-front-f GTGCCAGCMGCCGCGGTAA \ # adapter 1
   --p-front-r CCGYCAATTYMTTTRAGTTT \ # adapter 2
   --p-discard-untrimmed \ # removing the untrimmed sequences from qza
   --p-no-indels \ # removing detected indels
   --o-trimmed-sequences qza_files/trimmed_chicken.qza

# We are telling cutadapt to remove primers used to sequence the 16S V4-V5 region of the bacterial genome. 

# Now let's summarize those reads
qiime demux summarize --i-data qza_files/trimmed_chicken.qza --o-visualization qzv_files/trimmed_chicken.qzv
```

### Denoising Reads
To obtain our ASVs, we can use deblur, which denoises our reads by removing low quality reads, improving ASV assignment. QIIME also included another tool, DADA2, which can support paired-end reads. Because deblur does not, we will join our reads together using VSEARCH.
```bash
qiime vsearch join-pairs --i-demultiplexed-seqs qza_files/trimmed_chicken.qza \
                         --o-joined-sequences qza_files/trimmed_chicken_joined.qza

# Summarize and visualize the joined qza file
qiime demux summarize --i-data qza_files/trimmed_chicken_joined.qza --o-visualization qzv_files/trimmed_chicken_joined.qzv
```

Then, we filter out the low quality reads.
```bash
qiime quality-filter q-score --i-demux qza_files/trimmed_chicken_joined.qza \
                                    --o-filter-stats qza_files/filter_stats.qza \
                                    --o-filtered-sequences qza_files/trimmed_chicken_jf.qza

# Summarize and visualize the joined and filtered qza file
qiime demux summarize --i-data qza_files/trimmed_chicken_jf.qza --o-visualization qzv_files/trimmed_chicken_jf.qzv
```

The visualization here is important for picking a trimming length for our last deblur step. Take a look at the visualization, on the interactive quality plot tab. 
It is ideal to keep the quality scores of these reads above 20-30. Accordingly, we can set our trim length at 160. You can highlight portions of the plot to zoom in, and see in this case that the lower quality portions of the reads tend to end here. 

A couple of guidelines to keep in mind for your own data:
- Even if our quality scores were perfect, we want a minimum of 10 for this value
- It is ideal to keep the quality scores of these reads above 20-30. 
- Keep in mind we do not want a trim length longer than our reads
```bash
# Run deblur
qiime deblur denoise-16S --i-demultiplexed-seqs qza_files/trimmed_chicken_jf.qza \
                         --p-trim-length 160 \ 
                         --p-sample-stats \ # print out some stats
                         --p-jobs-to-start 4 \ # num threads
                         --p-min-reads 1 \ # remove features without at least this number of reads among all samples
                         --output-dir deblur_output
```

And once again we will visualize these results!
```bash
# summarize output table.qza and create visualization artifact
qiime feature-table summarize --i-table deblur_output/table.qza --o-visualization qzv_files/deblur_table.qzv
qiime feature-table summarize \
   --i-table deblur_output/table.qza \
   --o-visualization qzv_files/deblur_table.qzv
```

If there does not appear to be enough reads left in this output, try going back and changing the value of `X`. If you run into this problem with your experimental data, it may be a sign of low-quality sequences that do not match to reference. 

You have ASVs now! Keep the output handy, as you will need to refer back to it to get parameters for a few of the upcoming commands. 

>Question 2.2
>
>Take a look at the `deblur_table.qzv` visualization. What is the average frequency of our ASVs?

### Identifying Taxa
Naturally with a new set of ASVs, you might want to find out what they are! To assign taxonomy to these files, we will search with a Naive Bayes classifier traned against the GreenGenes database.

```bash
# output folder
mkdir taxa

qiime feature-classifier classify-sklearn --i-reads deblur_output/representative_sequences.qza \
                                          --i-classifier databases/greengenes/gg-13-8-99-nb-classifier.qza\
                                          --p-n-jobs 4 \
                                          --o-classification taxa/taxonomy.qza

# As with all QZA files, you can export the output file to take a look at the classifications and confidence scores
# Note: it will output a TSV file not a QZV file into a newly made directory taxa. 
qiime tools export --input-path taxa/taxonomy.qza --output-path taxa
# Viewing the output
less taxa/taxonomy.tsv
```

>Question 2.3
>
>How many ASVs do we have?

That takes us to the end of this portion of the workshop. The next lab will explore ways to analyze and visualize these ASVs, but for now, look through your list of taxa. Are there any that seem unexpected? Which taxa seem to pop up more often?

### Answers
- **Question 2.1 - Some of the other sections of the fastqc report, like "Per base sequence content" shows up as a warning. Why are these sections not as important for our data?**
	- fastqc isn't just built for 16S data, it is meant to work with most sequencing data. These other sections aren't as important for 16S samples because 16S sequencing amplifies specific regions of the genome, and as a result can look "wrong" to fastqc.
- **Question 2.2 - Take a look at the `deblur_table.qzv` visualization. What is the average frequency of our ASVs?**
	- ~32.5. You can find this value under the Frequency per feature table.
- **Question 2.3 - How many ASVs do we have?**
	- Using `wc -l`, we see that there are 151 ASVs

## Lab 3 - QIIME Part 2 
With our set of ASVs, we can perform a wide variety of analyses to understand our data. QIIME2 also includes tools for this! We will be going over a few of these tools and their visualizations in this part of the workshop.

### Phylogenetic Tree
One such visualization that can depict the makeup of your community in an interesting way is a phylogenetic tree. QIIME includes a few tools for doing this. We'll begin by making a *de novo* multi-sequence alignment of our ASVs, using a tool called MAFFT. 
```bash
# output folder for tree files
mkdir tree_out

#run mafft
qiime alignment mafft --i-sequences deblur_output/representative_sequences.qza \ # input
                      --p-n-threads 4 \ # number of threads
                      --o-alignment tree_out/rep_seqs_aligned.qza # output
```

Then we filter out unconserved columns or columns with many gaps.
```bash
qiime alignment mask --i-alignment tree_out/rep_seqs_aligned.qza \
                     --o-masked-alignment tree_out/rep_seqs_aligned_masked.qza
```

And lastly we will use FastTree to construct the phylogenetic tree, 
```bash
qiime phylogeny fasttree --i-alignment tree_out/rep_seqs_aligned_masked.qza \
                         --p-n-threads 4 \
                         --o-tree tree_out/phylogenetic_tree.qza

# and we will run this command to add a root at the midway point between the two tips with the largest distance
qiime phylogeny midpoint-root --i-tree tree_out/phylogenetic_tree.qza \
                              --o-rooted-tree tree_out/rooted_phylogenetic_tree.qza
```

QIIME does not have any native tools for tree viewing, but we can view this online [here](https://itol.embl.de/upload.cgi)! Upload your `rooted_phylogenetic_tree.qza`. If that does not work, try the command below/
```bash
qiime tools export \
  --input-path tree_out/rooted_phylogenetic_tree.qza \
  --output-path tree_out
```

### Taxonomic Composition
A common quality control step for 16S data is rarefaction, where we subsample our data to find whether we have sufficiently sequenced our samples. Rarefaction can tell us whether our data is a good representation of the community it comes from. Be aware that there is [some debate](https://academic.oup.com/bioinformatics/article/38/9/2389/6536959) over the use of rarefaction as a normalization step in microbiome data, and as always, make sure you are applying the right steps for *your* data.

We can perform rarefaction with the following command. The value for `--p-max-depth` is the highest sampling depth across all your samples according to `qzv_files/deblur_table.qzv`
```bash
qiime diversity alpha-rarefaction --i-table deblur_output/table.qza \
                                  --p-max-depth 582 \
                                  --p-steps 20 \
                                  --i-phylogeny tree_out/rooted_phylogenetic_tree.qza \
                                  --o-visualization qzv_files/rarefaction_curves.qzv
```

Check out the file on the QIIME viewer. Do the samples look as expected? Try changing the metric to see whether that changes. 

>Question 3.1
>
>How do we know if a sample properly captures the original community?

### Diversity Metrics
Some of the core analyses for our data include measuring the alpha and beta diversity of our samples. Alpha diversity concerns within-sample diversity, while beta diversity measures between-sample diversity. You might consider alpha diversity as assessing the structure of your sample, asking how many taxa there are and how evenly they are distributed. Beta diversity is an overall measure of your dataset, and how they relate to each other. For this tutorial, we will focus on alpha diversity, but we will cover beta diversity in the metagenomics tutorial!

We can generate output for both diversity metrics using this command. `--p-sampling-depth` this time comes from the lowest sampling depth among our samples.
```bash
qiime diversity core-metrics-phylogenetic --i-table deblur_output/table.qza \
                                          --i-phylogeny tree_out/rooted_phylogenetic_tree.qza \
                                          --p-sampling-depth 297 \
                                          --m-metadata-file metadata.tsv \
                                          --p-n-jobs-or-threads 4 \
                                          --output-dir diversity
```

Alpha diversity is commonly visualized with boxplots:
```bash
qiime diversity alpha-group-significance --i-alpha-diversity diversity/shannon_vector.qza \
                                         --m-metadata-file metadata.tsv \
                                         --o-visualization diversity/shannon_compare_groups.qzv
```

>Question 3.2
>
>Are any of our conditions associated with higher alpha diversity?

### Differential Abundance
With our last analysis, we will perform a differential abundance analysis on our data to find which taxa are more or less abundant in one condition over another. These results can give us clues as to how microbes respond to a stimuli or which microbes might be associated with a condition. There are many tools for differential abundance available;  QIIME implements a tool called [ANCOM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/) that treats our data as compositions (proportions of a whole, rather than counts). As a result, we have to add pseudocounts to our reads, so that none of the ASVs have a value of 0. This function will add 1 to our counts.
```bash
qiime composition add-pseudocount --i-table deblur_output/table.qza \
                                  --o-composition-table deblur_output/table_pseudocount.qza
```

Now we run the command
```bash
qiime composition ancom --i-table deblur_output/table_pseudocount.qza \
                        --m-metadata-file metadata.tsv \
                        --m-metadata-column treatment \ 
                        --output-dir ancom_output
```

The ANCOM visualization outputs one figure and two tables:
- The figure is a volcano plot where high absolute value on the x-axis indicates degree of differential abundance, while the y-axis indicates `W`, the number of times ANCOM has determined that ASV to be significantly differentially abundant. +ve clr means higher in control samples.
- The first chart shows the ASVs that have sufficiently high `W` to be considered differentially abundant
- The second chart shows the abundances of those ASVs in each group over quartiles.

>Question 3.3
>
>Which taxon is significantly more abundant in control samples?
>**HINT:** Use `grep`!

And that will be it for the QIIME tutorials! Try looking at some of the differentially abundant taxa for both conditions. ANCOM is a fairly stringent tool, so there may be fewer hits, but we can be more confident that these are not false positives. 

### Answers
- **Question 3.1 - How do we know if a sample properly captures the original community?**
	- If the curve plateaus, it indicates that we have reached a "saturation" point in representing taxa in the original community. A curve that does not reach this point may only capture a portion of that community. 
- **Question 3.2 - Are any of our conditions associated with higher alpha diversity?**
	- While it is not statistically significant, the antibiotic samples surprisingly seem to have higher alpha diversity. 
- **Question 3.3 - Which taxon is significantly more abundant in control samples?**
	- k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae 
	- `grep 014653c29367d6fb822868254b756bd5 taxa/taxonomy.tsv`

## License
Workshop materials are licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

