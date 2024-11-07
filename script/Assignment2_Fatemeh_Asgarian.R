## Introduction ------
# Pinaceae (Pine family) is the largest family of conifers globally, making up 11 genera (Ran et al., 2018). There is a lot of ecological and economic importance in these trees and so identifying the species and conducting phylogenetic analysis is vital for their conservation. Morphological studies are often not sufficient for correctly identifying members of the Pinaceae family (Ran et al., 2018). This increases the importance of using genetic markers and building classifiers, which is the objective of this assignment.
# The genes maturase K (matK) and ribulose-1,5-biphosphate carboxylase (rbcL) are both chloroplast genes often used as biomarkers for species identification of plants (Yong et al., 2024) I consequently decided to build a classifier to classify matK and rbcL for the Pinaceae family. For this assignment, I plan to use Random Forest for machine learning and classification using nucleotide proportions and repeat it using dinucleotide frequencies to compare the impact. I will then compare Random Forest to Support Vector Machine (SVM), another classification method, using the data with nucleotide proportions. I hypothesise that both methods will be equally accurate with classification because they are both popular for handling large datasets with many features and I further hypothesise that neither will be 100% successful due to both genes being chloroplast DNA.



## Packages used ------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(Biostrings)
# install.packages("rentrez")
library(rentrez)
# install.packages("seqinr")
library(seqinr)
library(ggplot2)
# install.packages("e1071")
library(e1071)
library(randomForest)
library(styler)

## CODE Section 1 - DATA PREPARATION ----

# I am retrieving my matK sequences from NCBI and I wish to know how long the sequence lengths on my dataset generally are so I can decide on what sequence length range is suitable for me to filter my data - I begin by doing a search and then fetching my data as a FASTA file

# PinaceaematK <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND matK[GENE]", use_history = T)
# PinaceaematK #We get 2714 hits on October 23
# matKmaxHits <- PinaceaematK$count
# Pinaceae_matK <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND matK[GENE]", retmax = matKmaxHits, use_history = T)
# Pinaceae_matK
# length(Pinaceae_matK$ids)
# Pinaceae_matK$web_history

# for (seq_start in seq(1, matKmaxHits, 500)) {
#   matK_fetch <- entrez_fetch(
#     db = "nuccore", web_history = Pinaceae_matK$web_history,
#     rettype = "fasta", retmax = 500, retstart = seq_start - 1
#   )
#   cat(matK_fetch, file = "../data/matK_fetch.fasta", append = TRUE)
#   cat(seq_start + 499, "sequences downloaded\r")
# }

# The next step I take is create a boxplot of the sequence length of all the sequences below 5000 nucleotides from my retrieved data because I am looking for matK genes, not whole chloroplast genomes

BoxplotFastaM <- "../data/matK_fetch.fasta"
SequencesBoxplotM <- readDNAStringSet(BoxplotFastaM)
SequenceLengthM <- width(SequencesBoxplotM)
dfSequenceLengthM <- data.frame(length = SequenceLengthM)
Filtered_dfSequenceLengthM <- subset(dfSequenceLengthM, length < 5000)

ggplot(Filtered_dfSequenceLengthM, aes(y = length)) +
  geom_boxplot(fill = "yellow", color = "black") +
  labs(
    title = "Boxplot of Sequence Lengths for matK",
    y = "Length (bp)",
    x = ""
  ) +
  theme_minimal()

# Judging by our boxplot, a 700-1700bp range seems to be a safe bet to obtain reasonable sequences of appropriate length - we redo our search and then fetch the sequence as a FASTA file

# TrimmedPinaceaematK <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND matK[GENE] AND 700:1700[SLEN]", use_history = T)
# TrimmedPinaceaematK #We get 1704 hits on October 23
# TrimmedmatKmaxHits <- TrimmedPinaceaematK$count
# Trimmed_Pinaceae_matK <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND matK[GENE] AND 700:1700[SLEN]", retmax = TrimmedmatKmaxHits, use_history = T)
# length(Trimmed_Pinaceae_matK$ids)
# Trimmed_Pinaceae_matK$web_history

# for (seq_start in seq(1, TrimmedmatKmaxHits, 500)) {
#   Trimmed_matK_fetch <- entrez_fetch(
#     db = "nuccore", web_history = Trimmed_Pinaceae_matK$web_history,
#     rettype = "fasta", retmax = 500, retstart = seq_start - 1
#   )
#   cat(Trimmed_matK_fetch, file = "../data/Trimmed_matK_fetch.fasta", append = TRUE)
#   cat(seq_start + 499, "sequences downloaded\r")
# }

# We can observe our retrieved data to see our data type and get a preview of the sequences
# class(Trimmed_matK_fetch)
# head(Trimmed_matK_fetch)

# We can write this into our working directory if we wish to
# write(Trimmed_matK_fetch, "Trimmed_matK_fetch.fasta", sep = "\n")

# The data must be written again in DNAStringSet format using readDNAStringSet so we can use the data to create a data frame with our titles and sequences - we can check to see if the data is in the right format with the right organisms and number of sequences

Final_matK <- "../data/Trimmed_matK_fetch.fasta"
FinalmatKSet <- readDNAStringSet(Final_matK)
class(FinalmatKSet)
names(FinalmatKSet)
length(FinalmatKSet)

# We then convert our dataset into a data frame and title our columns for better clarity of what it is that we are looking at
dfmatK <- data.frame(Title = names(FinalmatKSet), Sequence = paste(FinalmatKSet))
# View(dfmatK)

# We then create new columns, so you can see the species names and unique IDs derived from the Title and the gene name for the dataset
dfmatK$Species_Name <- word(dfmatK$Title, 2L, 3L)
dfmatK$ID <- word(dfmatK$Title, 1L)
dfmatK$Gene <- "matK"
# The columns are rearranged increase the ease of viewing the data
dfmatK <- dfmatK[, c("ID", "Title", "Gene", "Species_Name", "Sequence")]
# View(dfmatK)

# Next, I wish to "clean" my obtained sequences, meaning I wish to remove any potential sequences starting or ending with Ns and those with trailing gaps because we are not looking for aligned sequences but rather available nucleotides and k-mers in our code. Those sequences with internal Ns that exceed 2% of the code are removed as well - 2% is used as the cutoff point here to select for high quality sequences and exclude ambiguous sequences

dfmatKCleaned <- dfmatK %>%
  filter(!is.na(Sequence)) %>%
  mutate(Sequence2 = Sequence) %>%
  mutate(Sequence2 = str_remove(Sequence2, "^[-N]+")) %>%
  mutate(Sequence2 = str_remove(Sequence2, "[-N]+$")) %>%
  mutate(Sequence2 = str_remove_all(Sequence2, "-+")) %>%
  filter(str_count(Sequence2, "N") <= (0.02 * str_count(Sequence2)))
# view(dfmatKCleaned)
summary(nchar(dfmatKCleaned$Sequence2))

# We can compare our original sequence data with our new sequence data after cleaning it
dfmatKCleanedComparison <- cbind(dfmatKCleaned$Sequence, dfmatKCleaned$Sequence2)
view(dfmatKCleanedComparison)

# I wish to confirm that the dataset I have prepared is correct so I count the unique number of species in my data
length(unique(dfmatKCleaned$Species_Name))

# Creating a subset of the data by grouping it by species name and sampling only one of each to count rows
dfmatK_Subset <- dfmatK %>%
  group_by(Species_Name) %>%
  sample_n(1)
# The code above should give us 234 species if our dataset is prepared correctly, we confirm this by ensuring the code below gives us "TRUE"
all.equal(length(unique(dfmatKCleaned$Species_Name)), nrow(dfmatK_Subset))


# This procedure must be repeated for the second gene being observed and predicted in the code, the rbcL gene - we first identify the number of hits through an NCBI search and use it to fetch our data in FASTA format

# PinaceaerbcL <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND rbcL[GENE]", use_history = T)
# PinaceaerbcL # We get 2910 hits on October 23
# rbcLmaxHits <- PinaceaerbcL$count
# Pinaceae_rbcL <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND rbcL[GENE]", retmax = rbcLmaxHits, use_history = T)
# Pinaceae_rbcL
# length(Pinaceae_rbcL$ids)
# Pinaceae_rbcL$web_history


# for (seq_start in seq(1, rbcLmaxHits, 500)) {
#   rbcL_fetch <- entrez_fetch(
#     db = "nuccore", web_history = Pinaceae_rbcL$web_history,
#     rettype = "fasta", retmax = 500, retstart = seq_start - 1
#   )
#   cat(rbcL_fetch, file = "../data/rbcL_fetch.fasta", append = TRUE)
#   cat(seq_start + 499, "sequences downloaded\r")
# }

# The next step I take is create a boxplot of the sequence length of all the sequences below 5000 nucleotides from my retrieved data because I am looking for rbcL genes, not whole chloroplast genomes

BoxplotFastaR <- "../data/rbcL_fetch.fasta"
SequencesBoxplotR <- readDNAStringSet(BoxplotFastaR)
SequenceLengthR <- width(SequencesBoxplotR)
dfSequenceLengthR <- data.frame(length = SequenceLengthR)
Filtered_dfSequenceLengthR <- subset(dfSequenceLengthR, length < 5000)

ggplot(Filtered_dfSequenceLengthR, aes(y = length)) +
  geom_boxplot(fill = "red", color = "black") +
  labs(
    title = "Boxplot of Sequence Lengths for rbcL",
    y = "Length (bp)",
    x = ""
  ) +
  theme_minimal()

# We wish to search for sequences of rbcL within the 500-1500 range after looking at our boxplot as a safe bet of obtaining sequences of reasonable length and then we fetch the data from NCBI in FASTA format

# TrimmedPinaceaerbcL <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND rbcL[GENE] AND 500:1500[SLEN]", use_history = T)
# TrimmedPinaceaerbcL # We get 2001 hits on October 23
# TrimmedrbcLmaxHits <- TrimmedPinaceaerbcL$count
# Trimmed_Pinaceae_rbcL <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND rbcL[GENE] AND 500:1500[SLEN]", retmax = TrimmedrbcLmaxHits, use_history = T)
# Trimmed_Pinaceae_rbcL
# length(Trimmed_Pinaceae_rbcL$ids)
# Trimmed_Pinaceae_rbcL$web_history


# for (seq_start in seq(1, TrimmedrbcLmaxHits, 500)) {
#   Trimmed_rbcL_fetch <- entrez_fetch(
#     db = "nuccore", web_history = Trimmed_Pinaceae_rbcL$web_history,
#     rettype = "fasta", retmax = 500, retstart = seq_start - 1
#   )
#   cat(Trimmed_rbcL_fetch, file = "../data/Trimmed_rbcL_fetch.fasta", append = TRUE)
#   cat(seq_start + 499, "sequences downloaded\r")
# }

# We can observe our retrieved data to see our data type and get a preview of the sequences
# class(Trimmed_rbcL_fetch)
# head(Trimmed_rbcL_fetch)

# We can write this into our working directory if we wish to
# write(Trimmed_rbcL_fetch, "Trimmed_rbcL_fetch.fasta", sep = "\n")

# The data must be written again in DNAStringSet format using readDNAStringSet so we can use the data to create a data frame with our titles and sequences - we can check to see if the data is in the right format with the right organisms and number of sequences
Final_rbcL <- "../data/Trimmed_rbcL_fetch.fasta"
FinalrbcLSet <- readDNAStringSet(Final_rbcL)
class(FinalrbcLSet)
names(FinalrbcLSet)
length(FinalrbcLSet)

# We then convert our dataset into a data frame and title our columns for better clarity of what it is that we are looking at
dfrbcL <- data.frame(Title = names(FinalrbcLSet), Sequence = paste(FinalrbcLSet))
# View(dfrbcL)

# We then create new columns, so you can see the species names and unique IDs derived from the Title and the gene name for the dataset
dfrbcL$Species_Name <- word(dfrbcL$Title, 2L, 3L)
dfrbcL$ID <- word(dfrbcL$Title, 1L)
dfrbcL$Gene <- "rbcL"
# We rearrange the columns to observe our data with greater ease
dfrbcL <- dfrbcL[, c("ID", "Title", "Gene", "Species_Name", "Sequence")]
# View(dfrbcL)

# I "clean" the rbcL data with the same criteria I applied to the matK data
dfrbcLCleaned <- dfrbcL %>%
  filter(!is.na(Sequence)) %>%
  mutate(Sequence2 = Sequence) %>%
  mutate(Sequence2 = str_remove(Sequence2, "^[-N]+")) %>%
  mutate(Sequence2 = str_remove(Sequence2, "[-N]+$")) %>%
  mutate(Sequence2 = str_remove_all(Sequence2, "-+")) %>%
  filter(str_count(Sequence2, "N") <= (0.02 * str_count(Sequence2)))
# view(dfrbcLCleaned)
summary(nchar(dfrbcLCleaned$Sequence2))

# We can compare our original sequence data with our new sequence data after cleaning it to look for differences
dfrbcLCleanedComparison <- cbind(dfrbcLCleaned$Sequence, dfrbcLCleaned$Sequence2)
view(dfrbcLCleanedComparison)

# I wish to confirm that the dataset I have prepared is correct by looking at the unique species and creating a subset of my data with contains one of each species and checking to see if the numbers of species match
length(unique(dfrbcLCleaned$Species_Name))
dfrbcL_Subset <- dfrbcL %>%
  group_by(Species_Name) %>%
  sample_n(1)
# The code for the subset above should give us 246 species, we confirm this by checking if the following code gives us "TRUE"
all.equal(length(unique(dfrbcLCleaned$Species_Name)), nrow(dfrbcL_Subset))

# Combining our two datasets is important for the preparation of the classifier
dfCombined <- rbind(dfrbcLCleaned, dfmatKCleaned)
view(dfCombined)

# We can then run a few checks to see if the data was combined correctly
dim(dfCombined)
unique(dfCombined$Gene)
sum(is.na(dfCombined$Sequence2))
summary(str_count(dfCombined$Sequence2[dfCombined$Gene == "rbcL"]))
summary(str_count(dfCombined$Sequence2[dfCombined$Gene == "matK"]))

# Mean for rbcL is 845.7 vs 1172 for matK

# I will then remove the files that I do not need downstream
rm(dfrbcLCleanedComparison, dfmatKCleanedComparison)

## CODE Section 2 - DATA ANALYSIS ----

# We first need our sequence data in DNAStringSet format for our analysis
dfCombined <- as.data.frame(dfCombined)
dfCombined$Sequence2 <- DNAStringSet(dfCombined$Sequence2)
class(dfCombined$Sequence2)

# We then calculate the nucleotide frequencies and insert them onto our dataframe. letterFrequency is used because we retained 2% of the Ns in the data - nucleotide proportions are straightforward and can be an insightful first look into our classifiers
dfCombined <- cbind(dfCombined, as.data.frame(letterFrequency(dfCombined$Sequence2, letters = c("A", "C", "G", "T"))))

# We consider sequence length differences (remember the mean sequence length for rbcL is 845.7 vs 1172 for matK) so we calculate proportions of each nucleotide instead to have a fair comparison
dfCombined$Aprop <- (dfCombined$A) / (dfCombined$A + dfCombined$T + dfCombined$C + dfCombined$G)

dfCombined$Tprop <- (dfCombined$T) / (dfCombined$A + dfCombined$T + dfCombined$C + dfCombined$G)

dfCombined$Gprop <- (dfCombined$G) / (dfCombined$A + dfCombined$T + dfCombined$C + dfCombined$G)

dfCombined$Cprop <- (dfCombined$C) / (dfCombined$A + dfCombined$T + dfCombined$C + dfCombined$G)


# Here, we also add dinucleotide frequency into our dataset for future use later in the analysis - I wish to later compare the impact of using 2-mers instead and see if they are more informative for my classifier
dfCombined <- cbind(dfCombined, as.data.frame(dinucleotideFrequency(dfCombined$Sequence2, as.prob = TRUE)))

# We convert the DNAStringSet format back to character data so we can utilise tidyverse as columns in a tibble must be vectors to proceed
dfCombined$Sequence2 <- as.character(dfCombined$Sequence2)

# Seeing counts by gene name
table(dfCombined$Gene) # matK is 1692 and rbcL is 1990

# As our maximum sample size for these two genes is 1692 for the matK gene in dfCombined, we select for 425 individuals (about 25% of the total for the matK gene) to sample from each gene to serve as the validation data set

set.seed(303)
dfValidation <- dfCombined %>%
  group_by(Gene) %>%
  sample_n(425)

# We check to see if the correct number of sample size is available for both genes, expecting 425
table(dfValidation$Gene)

# Next we create a training dataset with different NCBI IDs from the ones picked for our validation set to ensure there's no overlap. Our maximum sample size will be 1267 genes (about 75% of the total for the matK gene)

set.seed(418)
dfTraining <- dfCombined %>%
  filter(!ID %in% dfValidation$ID) %>%
  group_by(Gene) %>%
  sample_n(1267)

# We check to see if the correct number of sample size is available for both genes, expecting 1267
table(dfTraining$Gene)

# We then build a classifier to separate the two genes using the columns assigned for our nucleotide proportions as the predictors and the genes as our response variable
Gene_Classifier <- randomForest::randomForest(x = dfTraining[, 11:14], y = as.factor(dfTraining$Gene), ntree = 500, importance = TRUE)
Gene_Classifier

# Here, 1265/1267 of the matK are correctly classified and 1261/1267 of rbcL are correctly classified. This is not 100% successful, but the out-of-the-bag error rate is very low at 0.32% meaning the classification predictions are correct 99.68% of the time.

# Additionally, we can visualise a variable importance plot for each predictor in our classifier (nucleotide proportions)
varImpPlot(Gene_Classifier)

# The results indicate that guanine (G) proportion has the highest MeanDecreaseAccuracy and MeanDecreaseGini - this indicates that proportion of G has the highest significance in our results

# We can further observe relative importance of each feature by
Gene_Classifier$importance

# We can also observe the number of times each row in the dataset was left out-of-the-bag and predicted
Gene_Classifier$oob.times

# We can look at out-of-the-bag error rate based on the trees preceding that iteration
Gene_Classifier$err.rate

# We can also view our confusion matrix separately if we wish to
Gene_Classifier$confusion

# We can look at the vote proportions for each row and see which rows there were disagreements on
Gene_Classifier$votes
# For example, row 55 here has 0.915254237 for matK 0.084745763 and rbcL.

# We must then test if our classifier works on sequences that we had excluded from our training data frame and included in our validation data frame (425 sequences for each gene)
predictValidation <- predict(Gene_Classifier, dfValidation[, c(3, 11:14)])
predictValidation
class(predictValidation)
length(predictValidation)

# We form our own confusion matrix for the validation data's observations and predicted data to see how well it worked
table(Observed = dfValidation$Gene, Predicted = predictValidation)

# We observe that 423/425 sequences were classified correctly for both matK and rbcL and whilst this is not a 100% success rate, it shows that our classifier functions with an extremely high success rate (99.5%)

# We can repeat this process using 2-mers (dinucleotide frequency) to assess its impact on our classifier using the following code (mostly adapted from our previous analysis with the nucleotide proportions - reminder that we already added dinucleotide frequencies to our data for analysis):

Gene_Classifier_2mers <- randomForest::randomForest(x = dfTraining[, 15:30], y = as.factor(dfTraining$Gene), ntree = 500, importance = TRUE)
Gene_Classifier_2mers

# This time, we observe a 100% success rate - 1267/1267 of our sequences are correctly identified with an out-of-the-bag error rate of 0%, we can test it with our validation data to see if this is replicated

predictValidation_2mers <- predict(Gene_Classifier_2mers, dfValidation[, c(3, 15:30)])
predictValidation_2mers
table(Observed = dfValidation$Gene, Predicted = predictValidation_2mers)

# It is true! 425/425 of our sequences are classified correctly for both matK and rbcL

# Moving on, to compare our randomForest analysis with SVM, we need to to run our data on SVM so we convert our gene column into a factor from a character
dfTraining$Gene <- as.factor(dfTraining$Gene)
class(dfTraining$Gene)

# We then must use our nucleotide proportions as predictors once again
SVM_Combined <- svm(
  Gene ~ Aprop + Tprop + Gprop + Cprop,
  data = dfTraining,
  kernel = "radial",
  cost = 1,
  gamma = 0.1
)
summary(SVM_Combined)

# We use the SVM dataset on our validation dataset after checking to see if all 850 sequences were used
SVM_PredictValidation <- predict(SVM_Combined, dfValidation)
SVM_PredictValidation
class(SVM_PredictValidation)
length(SVM_PredictValidation)
table(Observed = dfValidation$Gene, Predicted = SVM_PredictValidation)

# For SVM validation data, we observe that our matK sequences were all correctly classified but 20/425 of the rbcL sequences were misclassified as matK genes - I believed this to be something worth exploring personally. I decided to observe the G proportions (most significant in our data) of the misclassified rbcL sequences relative to other sequences that were correctly classified.
# To single out the 20 misclassified rbcL genes and look further into them, we create a separate column
dfValidation$Predicted_Gene <- SVM_PredictValidation

# We then filter the data to identify the 20 rows where the gene is rbcL but was predicted as matK
Misclassified_rbcL <- dfValidation %>%
  filter(Gene == "rbcL" & Predicted_Gene == "matK")
# view(Misclassified_rbcL)

# We must then add a column to indicate whether the rbcL sequence was misclassified or not using the ifelse function
dfValidation$Misclassified <- ifelse(dfValidation$Gene == "rbcL" & dfValidation$Predicted_Gene == "matK", "Misclassified", "Classified")

# Finally, we can create the boxplot for G proportions with the misclassified sequences labelled as points with a different colour

ggplot(dfValidation, aes(x = Gene, y = Gprop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Misclassified), width = 0.4, size = 2) +
  scale_color_manual(values = c("Classified" = "green", "Misclassified" = "red")) +
  labs(
    title = "Boxplot of Gprop for matK and rbcL Sequences",
    x = "Gene",
    y = "Gprop"
  ) +
  theme_minimal()

# We can clearly see that our misclassified sequences share similar G proportions as other matK sequences

# Returning to one the main objectives of the assignment, we wish to compare our Random Forest accuracy and SVM accuracy using the predicted data for validation and counting the correctly classified sequences
RandomForestAccuracy <- sum(predictValidation == dfValidation$Gene) / nrow(dfValidation)
SVMAccuracy <- sum(SVM_PredictValidation == dfValidation$Gene) / nrow(dfValidation)
print(paste("Random Forest Accuracy:", RandomForestAccuracy))
print(paste("SVM Accuracy:", SVMAccuracy))

# We can observe that using our validation dataset, Random Forest seems to have higher accuracy (0.995) than SVM (0.976)

## Discussion and Conclusion ----
# 

## References ----
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#using-ncbis-web-history-features
# Meyer, D. (2024). svm: Support Vector Machines (Version 1.7-16). In RDocumentation. https://www.rdocumentation.org/packages/e1071/versions/1.7-16/topics/svm
# Ran, J.-H., Shen, T.-T., Wu, H., Gong, X., & Wang, X.-Q. (2018). Phylogeny and evolutionary history of Pinaceae updated by transcriptomic analysis. Molecular Phylogenetics and Evolution, 129, 106–116. https://doi.org/10.1016/j.ympev.2018.08.011
# Rustam, Z., Sudarsono, E., & Sarwinda, D. (2019). Random-Forest (RF) and Support Vector Machine (SVM) Implementation for Analysis of Gene Expression Data in Chronic Kidney Disease (CKD). 9th Annual Basic Science International Conference 2019 (BASIC 2019), 546(5), 52066-. https://doi.org/10.1088/1757-899X/546/5/052066
# Yong, W. T. L., Mustafa, A. A., Derise, M. R., & Rodrigues, K. F. (2024). DNA barcoding using chloroplast matK and rbcL regions for the identification of bamboo species in Sabah. Advances in Bamboo Science (Online), 7, 100073-. https://doi.org/10.1016/j.bamboo.2024.100073

## Acknowledgements ----
# I would like to thank Moiz Ali Syed and Thomas Tekle from the class for proposing certain sources to look at on how I could obtain my data from NCBI as I was initially struggling with it and I ended up finding useful guidance as cited in my references (rentrez tutorial by David Winter). 
# I would also like to thank our TA, Brittany MacIntyre, for her guidance during class hours as well as office hours. She provided me with some essential pointers for my formatting.