# Module 4: Advanced microbiome statistics

This tutorial is part of the 2024 Canadian Bioinformatics Workshops [Advanced Microbiome Analysis](https://github.com/LangilleLab/microbiome_helper/wiki/CBW-workshop-2024-advanced) (St John's, NL, May 29-30).

Author: <a href="benfisher@dal.ca">Ben Fisher</a></p>



## Overview

The goal of this module is to help familiarize you with statistical methods, identifying important features, and machine learning for microbiome data in R. While we use specific tests for the following tutorial, they will not be applicable across all types of data. This tutorial aims to help you understand why these tests are chosen, as it is important to understand the meaning of the statistical tests you run, and demonstrate a rough framework with which you can conduct your own analyses. This will not be a definitive guide to statistics and machine learning in the microbial ecology field, however it should serve as a good starting point to inform your future endeavors.


## Contents

* [Statistics for Microbial Diversity](#statistics-for-microbial-diversity)
* [Alpha Diversity](#alpha-diversity)
  * [Alpha Diversity for Two Groups](#alpha-diversity-for-two-groups)
  * [Alpha Diversity for Three or More Groups](#alpha-diversity-for-three-or-more-groups)
* [Post-hoc Testing and Multiple Test Correction](#post-hoc-testing-and-multiple-test-correction)
* [Beta Diversity](#beta-diversity)
  * [Beta Diversity with Two Groups](#beta-diversity-with-two-groups)
  * [Beta Diversity with Three or More Samples](#beta-diversity-with-three-or-more-samples)
* [Covariates](#covariates)
  * [One-Way ANCOVA](#one-way-ancova)
  * [Two-Way ANCOVA](#two-way-ancova)
* [Visualizing Taxa Contributions](#visualizing-taxa-contributions)
* [Random Forest Classification](#random-forest-classification)
* [Random Forest Regression](#random-forest-regression)


# Statistics for Microbial Diversity

We have seen previously how to use different methods of visualizing diversity in microbiome analyses. Visualizations are an important tool in presenting our findings, however to quantitatively evaluate the differences we observe, we must employ statistics. The following will be an introduction to statistical tests that may be used following measures of alpha and beta diversity in microbiome samples.

We should begin by logging into the AWS instance and creating a new folder called `amb_module4`. This is where we will be working for the remainder of the tutorial.

To begin, import the following packages, ideally by creating an R Markdown Document and pasting the following in the setup chunk. If you are not using a markdown document, you can omit the `knittr::` line.

```
knitr::opts_chunk$set(echo = TRUE)

library(phyloseq)
library(stats)
library(FSA)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(rstatix)
library(randomForest)
library(rfUtilities)
library(dplyr)
library(car)
library(caret)
library(Maaslin2, lib.loc="/home/ubuntu/CourseData/MIC_data/.conda/envs/r-env/lib/R/library")
library(ROCR)
library(pdp)
library(RColorBrewer)
```

Next, we should set our working directory to the `amb_module4` folder using `setwd()` or the dialog box under the "session" tab.


# Alpha Diversity

As you may recall, alpha diversity evaluates the different types of taxa in a given sample and may consider parameters such as but not limited to  richness, evenness, or a combination of both. We will begin by importing our data. Instead of beginning from the `.biom` file as in Module 1, we have provided an RDS (R Data Serialization) file. Once this is imported, you will have a phyloseq object that has already been formatted and rarefied, so you can go straight to diversity.

```
data<-readRDS("~/CourseData/MIC_data/AMB_data/amb_module4/bacteria.Rds")
```

Now that the data is imported, you should see that it is a "Formal class phyloseq" object in the environment. Similar to before, we will visualize the alpha diversity plots to see whether it _looks_ like there are differences between our groups.

## Alpha Diversity for Two Groups

We can do this for any of our metadata categories that have **two** factor levels, such as `high_dysbiosis` or `diagnosis`. For this example, the code will calculate diversity between `high_dysbiosis` factors.

```
plot_richness(physeq = data, x="high_dysbiosis", color = "high_dysbiosis", measures = c("Observed", "Chao1", "Shannon")) + geom_boxplot() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
```

If you would like to substitute some different diversity metrics in this command, go right ahead! Some of the supported measures are: `"Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"`. Just know that if you want to do statistics for these measures, you will need to change some of the code in the following section.

Just by visual inspection, we can see that there is some separation of the means in our boxplots between our `high_dysbiosis` groups. However, to be sure, we should conduct a statistical test.

The Mann-Whitney U Test, also called Wilcoxon Rank Sum Test, is a non-parametric statistical test that is used to compare two groups, as we have here in our dataset. Being a non-parametric test, it does not assume normalcy in our data. The test assesses if the two groups are from the same population, i.e. whether they are equal. You can think of it as an extension of the unpaired Student's t-test, which does assume that the data follows a normal distribution. 

First, calculate the alpha diversity of our groups with `estimate_richness()` using the metrics we visualized above:

```
alpha_div<-estimate_richness(data, measures = c("Observed","Chao1","Shannon"))
```

If you use `View()` to inspect the `alpha_div` object, you will see a data frame with 10 rows and 5 columns that looks similar to the following:


	
| | Observed | Chao1 | se.chao1 | Shannon |
| --- | --- | --- | --- | --- |
| CSM79HR8_bracken_species | 88 | 91.27273 | 3.0864882 | 2.154705 |
| CSM7KOMH_bracken_species | 46 | 50.66667 | 4.4888950 | 1.396005 |
| HSM6XRQY_bracken_species | 62 | 65.11111 | 3.1004887 | 1.970721 |
| ... | ... | ... | ... | ... |

Each of the columns corresponds to the alpha diversity metric for each of the samples in our dataset. Chao1 also produces an `se.chao1` output column, but this is not necessary for the next steps. Remember, these are the numbers which are used "behind the scenes" to create the boxplot that was previously generated

We can conduct a Mann-Whitney U test on each of our alpha diversity metrics separately. We want to do this by comparing the `diagnosis` groups which our samples belong to, so we must append those groups to our `alpha_diversity` table using our sample metadata (contained in the phyloseq object):

```
alpha_div$diagnosis<-as.factor(data@sam_data$diagnosis)
```

Then, we can do our Mann-Whitney U test for each of the output columns in `alpha_diversity` that we are interested in:

```
mwu_observed<-wilcox.test(Observed ~ diagnosis, data = alpha_div)
mwu_chao1<-wilcox.test(Chao1 ~ diagnosis, data = alpha_div)
mwu_shannon<-wilcox.test(Shannon ~ diagnosis, data = alpha_div)
```

Now, we can view each of the test results separately by just executing the list names in our console (i.e. by typing `mwu_observed` and pressing **`enter`**. Or, we can summarize the results for viewing simplicity. There are several ways to do this, but one quick method is by doing the following:

```
mwu_summary<-cbind(mwu_observed,mwu_chao1,mwu_shannon)
```

**Question 1:** What does the `cbind()` function do? Why would a similar function, `rbind()`, not work here?

If you view the `mwu_summary` list, you will see that although it is not a perfect solution, we can see the p-values for each of the Mann-Whitney U tests that we ran for our chosen alpha diversity metrics. 

**Question 2** Think back to the boxplots. Are the results of our statistical tests expected?

## Alpha Diversity for Three or More Groups

Often, we may want to compare more than two groups. Since the Mann-Whiteny U test is similar to the unpaired t-test, it is only designed for making comparisons between two populations. If we want to go beyond this, we have to use a different test. Instead, we are going to employ the Kruskal-Wallis test, which is a further extension of the Mann-Whiteny U test for three or more groups. Again, we want to find out if the sample groups are representative of different populations, and therefore different from each other. Fortunately, we can apply a similar methodology to before with some additions.

First, to generate our three groups from the existing data, use the following command to create a vector which we coerce to the factor class, then append as a new column to the sample metadata. This assigns the groups in a specific order for the purposes of this tutorial:

```
group<-as.factor(c(3, 2, 2, 1, 3, 1, 3, 2, 1, 3))
data@sam_data$group <- group
```

And again, we should visualize our data. We will use the same alpha diversity boxplots as before, instead with the `group` variable as our factor. **If you run this command without changing something, it will not work. Use the error message to find the problem and fix it!**


```
plot_richness(physeq = data, x="group", color = "gruop", measures = c("Observed", "Chao1", "Shannon")) + geom_boxplot() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
```

Once you have your plot, take a look at our data. Does it look like there are differences between the groups? To find out, we will again use the `estimate_richness()` function to create a dataframe for our statistical tests, then append the `group` factors to the result.


```
alpha_div<-estimate_richness(data, measures = c("Observed","Chao1","Shannon"))
alpha_div$group<-group
```

**Question 3:** Ignoring the fact that we are appending `diagnosis` or `group` to our `alpha_div` dataframes, would this step produce the same dataframe that we generated for our Mann-Whitney U test? Look at the code to find out.

Now, we can conduct our Kruskal-Wallis tests. We will do this separately, and you can see the individual outputs as they are stored in our lists. We can also summarize the results quickly using the `cbind()` function, since our three output lists have the same columns.

```
kw_observed<-kruskal.test(Observed ~ group, data = alpha_div)
kw_chao1<-kruskal.test(Chao1 ~ group, data = alpha_div)
kw_shannon<-kruskal.test(Shannon ~ group, data = alpha_div)
kw_summary<-cbind(kw_observed,kw_chao1,kw_shannon)
```

Let's see what the summary says. Look at the p-values once you inspect the `kw_summary` dataframe with the `View()` command. You should see that the p-value is below 0.05 for all three diversity metrics in our summary table.

### Post-hoc Testing and Multiple Test Correction

What conclusions can we make about the three groups if the p-value is less than 0.05 for the Kruskal-Wallis test? We can only go so far as to say that there are differences between the groups. However, we don't know exactly _which_ groups are different. And, because we are comparing more than 2 groups, we must incorporate [multiple test correction](https://geneviatechnologies.com/blog/what-is-multiple-testing-correction/) as well. Because the Kruskal-Wallis test uses a rank-sum approach, it is appropriate to use a post-hoc test which preserves these ranks, such as Dunn's test (but an alternative could be the Conover-Iman test).

Instead of evaluating whether the groups are different, Dunn's test uses pairwise comparisons to evaluate which groups might be significantly different from each other. And, we will use the Benjamini–Hochberg multiple test correction method to mitigate the false-discovery rate that is inherent to a series of pairwise tests.

We will begin by conducting Dunn's test for each of the diversity metrics we are looking at:

```
d_observed<-dunn_test(data = alpha_div, Observed ~ group, p.adjust.method = "BH")
d_chao1<-dunn_test(data = alpha_div, Chao1 ~ group, p.adjust.method = "BH")
d_shannon<-dunn_test(data = alpha_div, Shannon ~ group, p.adjust.method = "BH")
```

And again, we can conveniently create a summary dataframe. Why might we be using `rbind()` instead of `cbind()` to do this?

```
dunn_summary<-rbind(d_observed, d_chao1, d_shannon)
```

Finally, we can inspect the `dunn_summary` table. Take note that there are p-values and adjusted p-values; we are interested in the adjusted p-values, because these have been corrected with the Benjamini–Hochberg method.

**Question 4:** Which groups, if any, are significantly different from one another? Is this the same for more than one diversity metric?

**Question 5:** Why do we not have to do a post-hoc test for the Mann-Whitney U test, but follow up our Kruskal-Wallis test with the Dunn test?

# Beta Diversity

Recall that beta diversity measures differences in composition between samples. 

**IF YOU ARE STARTING FROM THIS POINT AND HAVE NOT DONE THE ALPHA DIVERSITY STEP,** you need import and append the data first. Otherwise, skip this step.

```
data<-readRDS("~/CourseData/MIC_data/AMB_data/amb_module4/bacteria.Rds")
group<-as.factor(c(3, 2, 2, 1, 3, 1, 3, 2, 1, 3))
data@sam_data$group <- group
```

We will transform our absolute abundance data into relative abundance with the following function, which will normalize the entire sample to fraction of 1.

```
percentages <- transform_sample_counts(data, function(x) x*100 / sum(x))
```

## Beta Diversity with Two Groups

Next, we will create an ordination and plot it to first visualize the differences between our groups. This will look the same as the plot you generated in module 1.

```
ordination<-ordinate(physeq = percentages, method = "PCoA", distance = "bray")
plot_ordination(physeq = percentages, ordination = ordination, color="diagnosis")
```

Now, we can visually inspect our plot. From this, it is difficult to tell whether there is much separation between groups. Truthfully, visual inspection of beta diversity can be misleading, and groups that overlap a lot can be significantly different from each other. So, to do a statistical test, we must first create the matrix of bray-curtis dissimilarity between samples.

```
bray_dist<-distance(data, method = "bray")
```

Now that we have the raw data, we can do a permutational ANOVA (PERMANOVA) to see if there are significant differences. The PERMANOVA is used to analyze multivariate data that uses a dissimilarity metric. It is a consensus-based approach that uses the whole dataset and random subsets of the data to calculate a test statistic. In R, this can be done with the `adonis2()` function from the vegan package. PERMANOVA and Adonis tend to be used interchangeably to talk about the same test. To run our Adonis, use the following code:

```
adonis<-adonis2(bray_dist ~ sample_data(data)$diagnosis, permutations = 101)
```

Now, you can view the results of our test by inspecting the `adonis` object with the `View()` function or the environment browser. Are there significant differences between the groups?

**Question 6:** If we were to find a significant difference, would we need to use a post-hoc test?


## Beta Diversity with Three or More Samples

Again, let's plot our visualization. We can use the same `ordination` that we created above, and just substitute `"group"` for the `colour` parameter in the following command. You will notice that the ordination of the points in space is the same, but the colors indicating our groups are different.

```
plot_ordination(physeq = percentages, ordination = ordination, color="group")
```

Like above, we will construct our matrix of bray-curtis dissimilarity, and then conduct our PERMANOVA to compare the three groups.

```
bray_dist<-distance(data, method = "bray")
adonis<-adonis2(bray_dist ~ sample_data(data)$group, permutations = 100)
```

Now, unlike when we calculated alpha diversity, we do not see significant differences between the groups when measured by bray-curtis dissimilarity. So, while the taxa are different in richness and evenness, the composition of the samples is comparable.

Normally, we would not need to continue to a post-hoc test. However, for the sake of demonstration, we will do pairwise comparisons between the groups. This can be accomplished with the pairwiseAdonis package, as shown below:

```
pairwise<-pairwise.adonis(bray_dist, sample_data(data)$group)
```

Now, inspect the `pairwise` object. Are there any significant differences? Is this expected?


# Covariates

In biological studies, covariates are continuous numerical factors that are not implicitly part of the study design, but may contribute to differences we observe between groups. In microbiome studies, common covariates would be Age, BMI, and other numerical factors that have previously been shown to contribute to differences we observe in our measured variable of interest.

First, we will download our data containing covariates. The HMP metadata does not contain BMI from all samples, so we extrapolated the values between diagnosis groups for the sake of this tutorial. **This is not actual BMI data provided by the HMP, it has been added and manipulated for demonstration purposes.** We can use the base r `download_file()` function to get the .csv file from pastebin:

```{r}
download.file("https://pastebin.com/raw/v6z0WMCK",destfile = "~/workspace/amb_module4/covariates.csv")
```

Now, the file is stored on our instance, but we still need to read it into R. We will also ensure that our factors are correctly read as so. Let's do that with the following:

```
covariates<-read.csv("covariates.csv", header = TRUE, row.names = 1)
covariates$sex <- as.factor(covariates$sex)
```

To take a first look at our covariates, we can see how dysbiosis score might differ between `sex` factor levels. We will do this with boxplots. Additionally, we will use the `facet_grid()` "geom" from ggplot2 to separate our plots into different facets. This way, we can see if there are any differences between sexes for each diagnosis group.

```
box<-ggplot(covariates, aes(x=sex, y=dysbiosis_score, color=diagnosis)) + geom_boxplot() + facet_grid(diagnosis ~ .)
```

And view the resulting plot with:

```
box
```

Another intuitive way to look at our data is with a scatter plot. We will add regression lines to our plots to see if they are collinear or not. If we do not see multicollinearity, it may indicate that there is an effect of our covariates on our `dysbiosis_score` measurement.

```
scatter<-ggplot(covariates, aes(x=BMI, y=dysbiosis_score, color=diagnosis)) + geom_point() + geom_smooth(method = "lm", se = FALSE, fullrange = TRUE)  + facet_grid(sex ~ .)
```
And view the resulting plot with:

```
scatter
```

The following is not necessary, but is a demonstration of the `facet_rid()` geom. We can additionally separate our boxplots by an additional factor, `high_dysbiosis_score`. This would not be a practical covariate as it is a categorical interpretation of our measured variable, `dysbiosis_score`, but it is interesting to see as a summary.


```
scatter2<-ggplot(covariates, aes(x=BMI, y=dysbiosis_score, color=diagnosis)) + geom_point() + geom_smooth(method = "lm", se = FALSE, fullrange = TRUE)  + facet_grid(sex ~ high_dysbiosis_score)
scatter2
```

To assess the effect of our treatment groups and covariates on our response variable, we will do an analysis of variance, or ANOVA. First, for demonstration, we will do a two-way ANOVA for our data with the following. Our groups will be `sex` and `diagnosis`, and our response variable is `dysbiosis_score`.

```
aov <- aov(dysbiosis_score ~ sex + diagnosis, data = covariates)
summary(aov)
```

We can see that there are in fact significant differences. However, it is not clear what contribution our covariates are making to these differences. Additionally, we cannot use continuous variables as treatment levels with the ANOVA. So, we must run a variation of the ANOVA called the ANCOVA, or analysis of co-variance, with which our numerical data can be used. The setup is similar to an anova, but we will construct a linear regression and give it to the `Anova()` function.

## One-way ANCOVA

We will begin by conducting a one-way ANCOVA, to see how `BMI` contributes to the variance. It is a one-way ANCOVA because we are only using one "treatment" category, `sex`.

```
ancova <- lm(dysbiosis_score ~ BMI + sex, data = covariates)
Anova(ancova, type="II")
```

Here, we've output the results of our ANCOVA. We should see a significant contribution of `BMI` to our observed variance in `dysbiosis_score`. 

The ANCOVA actually conduct multiple tests to generate its test statistic. So, we must incorporate multiple test correction. We can do this by using pipes from the `dplyr` package to supply our data to `emmeans_test()` for Benjamini–Hochberg multiple test correction.

```
pwc <- covariates %>% 
  emmeans_test(dysbiosis_score ~ sex, covariate = BMI,
    p.adjust.method = "BH"
    )
```

Now, let's view the `pwc` object that stores our results. You should see that the adjusted p value is now above our significance threshold of 0.05. This can be disappointing, but highlights the importance of these corrections in statistical analyses so we are not mislead by our tests.

## Two-way ANCOVA

We can also conduct a two-way ANCOVA with our `BMI` covariate. This follows the same procedure as before, but we must specify our second categorical factor, `diagnosis`, in this model. Go ahead and run the chunk below to generate the results.

```
#Compute 2-way ANCOVA
ancova <- lm(dysbiosis_score ~ BMI + sex * diagnosis, data = covariates)
Anova(ancova, type="II")

#Compute pairwise comparisons between treatment groups, using multiple test correction
pwc <- covariates %>% 
  group_by(diagnosis) %>%
  emmeans_test(dysbiosis_score ~ sex, covariate = BMI, p.adjust.method = "BH")
```

Again, we should see that the adjusted p value shows no significance between categories with our covariate. In practice, you would not want to add additional factor levels if your design with less factors and/or covariates does not produce significance. This is just to demonstrate use cases of the ANCOVA. If you DID see significant differences following the multiple test correction, you would want to investigate the simple main effect.

# Visualizing Taxa Contributions

You may remember using Maaslin2 for differential abundance in our previous modules. Another useful feature of Maaslin2 is its ability to calculate and visualize the contributions of individual taxa to our response variable. To being, import our `rf_data.tsv` file.

```
data<-read.delim("~/CourseData/MIC_data/AMB_data/amb_module4/rf_data.tsv", header = TRUE, sep = "\t", row.names = 1)
```

This file contains metadata that we want to extract. So, we will create a metadata object that is a subset of the dataframe, then remove the metadata from our original data. We can accomplish this by specifying which columns to extract and subsequently remove.

```
input_metadata <- data[,1:3]
data<-data[,-1:-3]
```

And with that, we're ready to go! Use the following to run Maaslin2 and generate our visualizations.

```
fit_data = Maaslin2(input_data     = data, 
                    input_metadata = input_metadata, 
                    min_prevalence = 0,
                    normalization  = "NONE",
                    output         = "~/workspace/amb_module4/Maaslin_output", 
                    fixed_effects  = c("diagnosis", "dysbiosis_score"),
                    reference      = c("diagnosis,nonIBD"))
```

Maaslin2 will create several outputs. To view these, connect to your directory in the browser with `http://##.uhn-upc.ca`, and navigate to the output folder. You will be able to view the PDF files in your browser, and the TSV files will download to your machine when you click them. One of the output files is called `significant_results.tsv`. This contains the features that Maaslin identifies as being significant to accurately predicting our `dysbiosis_score`. Look for these features in the boxplots and scatterplots to see the magnitude of their contribution, and view the heatmap to see whether they contribute to high or low scores, as well as which diagnosis they may be positively or negatively assosciated with.


# Random Forest Classification

We can construct a random forest model to attempt to classify our data into categories based on features of our dataset. In this case, we want to see if by using bacteria abundances, we can predict which `diagnosis` our data would fall into.

To begin, import our `rf_data.tsv` file:

```{r}
data<-read.delim("~/CourseData/MIC_data/AMB_data/amb_module4/rf_data.tsv", header = TRUE, row.names = 1)
```

If you view the `data` dataframe, you will see that we have three metadata columns. Since we are actually only interested in how the features, and not the metadata, may inform our `diagnosis` factor, we need to remove the other columns. One way to do this is by sub-setting the data frame and keeping only our columns of interest. The first value in the character vector in our `data[]` function corresponds to the index of the column we will keep.

```
#subset the data to the metadata variable we are interested in
#1 = dysbiosis_score
#2 = high_dysbiosis
#3 = diagnosis
data<-data[,c(2,4:ncol(data))]
```

We want to ensure all of our data is numeric, and a quick way to do this is by coercing it to a matrix and changing it back. The two commands below will accomplish this, however you should try to see if you can do it in one line by "nesting" the functions.

```
mx<-data.matrix(data)
df<-as.data.frame(mx)
```

We also want to make sure that our metadata column is in the `factor` form, as the matrix transformation turned them into binary values as the integers 0 and 1. Change the class like so:

```
df[,1]<-as.factor(df[,1])
```

Now, if we want to validate our random forest classifier, we have to create a "training" and "test" dataset. To do this, we will create an index of random integers where 70% will = 1 and 30% will = 2. Then, we will append the `ind` object to our data and use it to create our two subsets.

```
set.seed(222)
ind <- sample(2, nrow(df), replace = TRUE, prob = c(0.7, 0.3))
train <- df[ind==1,]
test <- df[ind==2,]
```

Using our training data, we can "train" the classifier to find features that correspond with a high dysbiosis score (as our categorical variable).

```
rf <- randomForest(train[,2:ncol(train)], train[,"high_dysbiosis"], ntree = 501, proximity=TRUE, importance = TRUE)
```

Now, we can use our model to predict whether our data will fall into the `high_dysbiosis` category or not. We will predict using both the training and test data.

```
p1 <- predict(rf, train)
p2 <- predict(rf, test)
```

If you were to view the objects resulting from the above commands, you are given a list (in factor form) of each of the samples and their boolean value (0 or 1) corresponding to FALSE or TRUE for being classified as `high_dysbiosis`. We can use the following commands to generate a confusion matrix for our training and testing datasets, which indicate how well the model performed on our data.

```
#what is the model accuracy?
confusionMatrix(p1, train$high_dysbiosis)
#what is the accuracy with the test data?
confusionMatrix(p2, test$high_dysbiosis)
```

**Question 7:** Does the model accuracy differ between the `train` and `test` data? What does the accuracy value for each mean?


A useful step in evaluating our model is seeing which features contribute the most to the classification. The object resulting from `randomForest` will contain `importance` data. If we view this, we can sort by the parameters to see the most and least important factors.

```
View(rf$importance)
```

To more easily compare the features, let's create a visualization:

```
varImpPlot(rf, sort = TRUE, n.var = 15)
```

Using these features, we can gain a sense of how the model works. Check out the partial plots for some of these features to see how the abundance (x-axis) informs the `high_dysbiosis` category (y-axis). If you do this for several features, you will begin to see how the random forest classifier considers each of the feature values for a given sample. The combination of all of these partial plots is a representation of the entire model.

```
pd<-partial(rf, pred.var = "Negativicutes")
plotPartial(pd)
```

We can permutationally evaluate the significance of our model using functions from the rfUtilities package. This command will take some time to run as it is randomly subsetting the data to calculate the test statistic.

```
sig<-rf.significance(rf, xdata = train[,2:ncol(train)], ntree = 501, nperm = 11)
sig
```

We can additionally use cross-validation on our model, which is another permutational approach.

```
cv<-rf.crossValidation(rf, xdata = train[,2:ncol(train)], ntree = 501, n = 11, p=0.1)
cv
```

And, to generate a summary of our model, we can use the `rpart.plot()` function from rpart.plot. This may not be installed on your system, so if you get an error with the `library()` function importing this package, remove the hashtag (#) from the `install.packages()` line and run it again.

```
#install.packages("rpart.plot")
library(rpart.plot)
tree <- rpart(high_dysbiosis ~., data = data, method = "class")
rpart.plot(tree)
```

Another common plot you will likely come across is the Reciever Operator Curve, or ROC plot. This plot shows the performance of the model at all classification thresholds by plotting true positives vs false positives. This is therefore an extension of the confusion matrices we generated before. For a model to operate well, you want the area under the curve (AUC) to be higher than 0.5, which represents out-of-bag error, or random chance.

To do this, we will use functions from the caret and ROCR packages.

```
pred1<-predict(rf, type = "prob")
pred1[,2]
perf<-prediction(pred1[,2], test$high_dysbiosis)
auc = performance(perf, "auc")
pred2 = performance(perf, "tpr","fpr")
plot(pred2,main="ROC Curve for Random Forest",col=2,lwd=2)
```

And with that, you can see that our model does in fact seem to perform better than random!


# Random Forest Regression

Another approach to random forest is using it to predict a continuous variable. This is achieved by regression instead of classification. We just completed a model for finding whether our samples would fit into one of two categories. Now, we want to actually predict the value of our response variable. This is a powerful application of random forests which will allow us to attempt to predict the `dysbiosis_score` from the feature abundances in our samples. Since much of this section is a repeat of the previous, less commentary will follow.

We want to use the following chunk to import and format our dataframe for our model. This time, we are electing to keep the `dysbiosis_score` column as our response variable for regression.

```{r}
#import the data for random forest 
data <- read.delim("~/CourseData/MIC_data/AMB_data/amb_module4/rf_data.tsv", header=TRUE)

#subset the data to the metadata variable we are interested in
data<-data[,c(1,4:ncol(data))]
mx<-data.matrix(data)
df<-as.data.frame(mx)


#change the class of our metadata from character to factor.
#df[,1]<-as.factor(df[,1])

#Separate our dataset into 70% training data and 30% validation data.
set.seed(222)
ind <- sample(2, nrow(df), replace = TRUE, prob = c(0.7, 0.3))
train <- df[ind==1,]
test <- df[ind==2,]
```

Now, we will use our training data to construct the model. Similar to before, we are specifying our "x" (feature abundances) and "y" (dysbiosis score) parameters.

```
rf <- randomForest(train[,2:ncol(train)], train$dysbiosis_score, ntree = 501, proximity=TRUE, importance = TRUE)
```


Since we are not using a categorical response variable, we cannot use a confusion matrix. Instead, we could look at the model's predictive abilities:
```
p1 <- predict(rf, train)
p2 <- predict(rf, test)
```

Which will return the sample name and dysbiosis scores. Try to find a couple of samples in `data` that are also in `p2` (which is our predictions for the `test` subset). Are they similar?


Again, we can see which features are important to our model's performance. If the output is small, you may need to enlarge your viewing window.

```
varImpPlot(rf, sort = TRUE, n.var = 15)
```


Let's create some visualizations for our model. These will be different than before because we don't have categories, but a continuous range of response variables. This means that we cannot use ROC plots. Instead, let's see if we can view some partial importance plots with the pdp package. This may not be installed on your system, so if you get an error with the `library()` function importing this package, remove the hashtag (#) from the `install.packages()` line and run it again.

```
#install.packages("pdp")
library(pdp)
```

Using the importance plot you previously generated, pick a couple of the features to substitute into the following command. The first `partial()` function may take a moment to run:

```
pd <- partial(rf, pred.var = c("Alistipes", "Rikenellaceae"))
pdp1 <- plotPartial(pd, pred.var = c("Alistipes"))

rwb <- colorRampPalette(c("red", "white", "blue"))
pdp2 <- plotPartial(pd, contour = TRUE, col.regions = rwb)
```

Now, we have created two partial plots for viewing. Using the gridExtra package, we can view them simultaneously:

```
library(gridExtra)
grid.arrange(pdp1, pdp2, ncol = 2)
```

The resulting partial dependance plots show how our model will classify a sample based on the abundances of these two taxa. The two plots are just using separate colour schemes to demonstrate how this may change your interpretability of the visualization. You can imagine then how a sample could be represented as a point in this two-dimensional plot. The dysbiosis score, indicated by the colour in the legend, would change depending on its coordinates. Remember that these coordinates correspond to abundances of your features.

So, you can begin to conceptualize how the model is working as a whole. It places a point in n-dimensional space, where n is the number of features, and estimates the `dysbiosis_score` based on its position. To simplify, you can imagine constructing these partial dependance plots for all features in the sample, and by finding its position in each you could estimate the score.

Try constructing another plot using different features. How does the result change?




Again, we can test the significance of our model and cross validate with the following:

```
#test model significance
sig<-rf.significance(rf, xdata = train[,2:ncol(train)], ntree = 501, nperm = 11)
sig

#cross-validate the model
cv<-rf.crossValidation(rf, xdata = train[,2:ncol(train)], ntree = 501, n = 11, p=0.1)
cv
```

**Question 1:** What does the `cbind()` function do? Why would a similar function, `rbind()`, not work here?
_Using `help("cbind)`, we can view the manual for this function. It combines given sequence of vector, matrix or data-frame arguments into one output. It matches by columns, and our MWU test outputs all have the same columns, so this creates a merged table. We could not use `rbind()` because the MWU outputs do not share rows._

**Question 2** Think back to the boxplots. Are the results of our statistical tests expected?
_If you used the `high_dysbiosis` factor, we might expect to see significant differences when considering the boxplot. The means had some separation._

**Question 3:** Ignoring the fact that we are appending `diagnosis` or `group` to our `alpha_div` dataframes, would this step produce the same dataframe that we generated for our Mann-Whitney U test? Look at the code to find out.
_Yes, the code we use to generate `alpha_div` is the same, and we use the same data to begin, so the same dataframe would be produced **before** we append the factors._

**Question 4:** Which groups, if any, are significantly different from one another? Is this the same for more than one diversity metric?
_Group 1 should be significantly different from Group 2, but not 3. This is consistent across the three diversity metrics we used. Groups 2 and 3 are not significantly different from each other. Sometimes, different methods will produce different results for the same group comparisons._

**Question 5:** Why do we not have to do a post-hoc test for the Mann-Whitney U test, but follow up our Kruskal-Wallis test with the Dunn test?
_The MWU test does not do multiple comparisons. We need to use a post hoc test with multiple test correction if we are comparing more than two groups. For the KW test, we do 3 pairwise comparisons to calculate the test statistic, and therefore need a post-hoc test to see which groups are significantly different from one another._

**Question 6:** If we were to find a significant difference, would we need to use a post-hoc test?
_No, because we are only doing one comparison between our groups._

**Question 7:** Does the model accuracy differ between the `train` and `test` data? What does the accuracy value for each mean?
_Yes. The model accuracy should be close to 1 (100%) for the training data because it is the same data that produced the classifer. For the test dataset, we expect it will be lower. The accuracy value is how well our model is able to accurately predict our response variable._



