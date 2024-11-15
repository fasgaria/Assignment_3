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
# install.packages("progress")
library(progress)
# install.packages("caret")
library(caret)
# install.packages("stringr")
library(stringr)



##1. CODE Section 1 - DATA PREPARATION ----

# Define the download_sequences function
download_sequences <- function(db, organism, gene, output_file, batch_size = 500) {
  library(rentrez)
  library(progress)  # Load progress bar library
  
  # Step 1: Perform the search and set up web history
  search_term <- paste0(organism, "[ORGN] AND ", gene, "[GENE]")
  search_result <- entrez_search(db = db, term = search_term, use_history = TRUE)
  total_hits <- search_result$count
  message("Total sequences found: ", total_hits)
  
  # Set up the progress bar
  pb <- progress_bar$new(
    format = "  Downloading [:bar] :percent in :elapsed",
    total = ceiling(total_hits / batch_size),
    clear = FALSE,
    width = 60
  )
  
  # Step 2: Loop through results in batches and write to the output file
  for (seq_start in seq(1, total_hits, by = batch_size)) {
    fetched_sequences <- entrez_fetch(
      db = db,
      web_history = search_result$web_history,
      rettype = "fasta",
      retmax = batch_size,
      retstart = seq_start - 1
    )
    
    # Write each batch to the file
    cat(fetched_sequences, file = output_file, append = TRUE)
    
    # Update the progress bar
    pb$tick()
  }
  
  message("\nDownload complete! Sequences saved to ", output_file)
}

# Running the download_sequences function for the matK gene.
download_sequences(
  db = "nuccore",
  organism = "Pinaceae",
  gene = "matK",
  output_file = "../data/new_matK_fetch.fasta"
)

# Running the download_sequences function for the rbcL gene.
download_sequences(
  db = "nuccore",
  organism = "Pinaceae",
  gene = "rbcL",
  output_file = "../data/new_rbcL_fetch.fasta"
)




# The next step I take is create a boxplot of the sequence length of all the sequences below 5000 nucleotides from my retrieved data because I am looking for matK genes, not whole chloroplast genomes

BoxplotFastaM <- "../data/new_matK_fetch.fasta"
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

TrimmedPinaceaematK <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND matK[GENE] AND 700:1700[SLEN]", use_history = T)
TrimmedPinaceaematK #We get 1704 hits on October 23
TrimmedmatKmaxHits <- TrimmedPinaceaematK$count
Trimmed_Pinaceae_matK <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND matK[GENE] AND 700:1700[SLEN]", retmax = TrimmedmatKmaxHits, use_history = T)
length(Trimmed_Pinaceae_matK$ids)
Trimmed_Pinaceae_matK$web_history

for (seq_start in seq(1, TrimmedmatKmaxHits, 500)) {
  Trimmed_matK_fetch <- entrez_fetch(
    db = "nuccore", web_history = Trimmed_Pinaceae_matK$web_history,
    rettype = "fasta", retmax = 500, retstart = seq_start - 1
  )
  cat(Trimmed_matK_fetch, file = "../data/Trimmed_matK_fetch.fasta", append = TRUE)
  cat(seq_start + 499, "sequences downloaded\r")
}

# We can observe our retrieved data to see our data type and get a preview of the sequences
class(Trimmed_matK_fetch)
head(Trimmed_matK_fetch)

# We can write this into our working directory if we wish to
write(Trimmed_matK_fetch, "Trimmed_matK_fetch.fasta", sep = "\n")

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
# Cleaning will happen downstream.




# The next step I take is create a boxplot of the sequence length of all the sequences below 5000 nucleotides from my retrieved data because I am looking for rbcL genes, not whole chloroplast genomes

BoxplotFastaR <- "../data/new_rbcL_fetch.fasta"
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

TrimmedPinaceaerbcL <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND rbcL[GENE] AND 500:1500[SLEN]", use_history = T)
TrimmedPinaceaerbcL # We get 2001 hits on October 23
TrimmedrbcLmaxHits <- TrimmedPinaceaerbcL$count
Trimmed_Pinaceae_rbcL <- entrez_search(db = "nuccore", term = "Pinaceae[ORGN] AND rbcL[GENE] AND 500:1500[SLEN]", retmax = TrimmedrbcLmaxHits, use_history = T)
Trimmed_Pinaceae_rbcL
length(Trimmed_Pinaceae_rbcL$ids)
Trimmed_Pinaceae_rbcL$web_history


for (seq_start in seq(1, TrimmedrbcLmaxHits, 500)) {
  Trimmed_rbcL_fetch <- entrez_fetch(
    db = "nuccore", web_history = Trimmed_Pinaceae_rbcL$web_history,
    rettype = "fasta", retmax = 500, retstart = seq_start - 1
  )
  cat(Trimmed_rbcL_fetch, file = "../data/Trimmed_rbcL_fetch.fasta", append = TRUE)
  cat(seq_start + 499, "sequences downloaded\r")
}

# We can observe our retrieved data to see our data type and get a preview of the sequences
class(Trimmed_rbcL_fetch)
head(Trimmed_rbcL_fetch)

# We can write this into our working directory if we wish to
write(Trimmed_rbcL_fetch, "Trimmed_rbcL_fetch.fasta", sep = "\n")

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



# Define the clean_sequences function.
clean_sequences <- function(data_frame, sequence_column) {
  library(stringr)
  library(dplyr)
  
  # Summary before cleaning
  cat("********** BEFORE CLEANING **********\n")
  cat("Number of sequences:", nrow(data_frame), "\n")
  total_Ns_before <- sum(str_count(data_frame[[sequence_column]], "N"))
  cat("Total number of Ns:", total_Ns_before, "\n")
  
  # Clean the sequences
  cleaned_data <- data_frame %>%
    filter(!is.na(!!sym(sequence_column))) %>%  # Remove NA sequences
    mutate(Sequence2 = !!sym(sequence_column)) %>%
    mutate(Sequence2 = str_remove(Sequence2, "^N+")) %>%  # Remove leading Ns
    mutate(Sequence2 = str_remove(Sequence2, "N+$")) %>%  # Remove trailing Ns
    filter(str_count(Sequence2, "N") <= (0.02 * str_count(Sequence2)))  # Keep sequences with <= 2% Ns
  
  # Summary after cleaning
  cat("\n********** AFTER CLEANING **********\n")
  cat("Number of sequences after cleaning:", nrow(cleaned_data), "\n")
  total_Ns_after <- sum(str_count(cleaned_data$Sequence2, "N"))
  cat("Total number of Ns after cleaning:", total_Ns_after, "\n")
  sequences_removed <- nrow(data_frame) - nrow(cleaned_data)
  cat("Number of sequences removed during cleaning:", sequences_removed, "\n")
  cat("Summary of sequence lengths after cleaning:\n")
  summary_after <- summary(nchar(cleaned_data$Sequence2))
  print(summary_after)
  
  return(cleaned_data)
}


# Cleaning the matK sequences using the new function.
new_dfmatKCleaned <- clean_sequences(dfmatK, "Sequence")

# Cleaning the rbcL sequences using the new function.
new_dfrbcLCleaned <- clean_sequences(dfrbcL, "Sequence")


# Combining the cleaned datasets
dfCombined <- rbind(new_dfrbcLCleaned, new_dfmatKCleaned)

# We can then run a few checks to see if the data was combined correctly. As seen previously in the summaries, the mean for rbcL is 845.7 vs 1172 for matK.
dim(dfCombined)
unique(dfCombined$Gene)
sum(is.na(dfCombined$Sequence2))

# I will then remove the files that I do not need downstream
rm(new_dfrbcLCleaned, new_dfmatKCleaned)







##2. CODE Section 2 - DATA ANALYSIS ----

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









## Random Forest Classifier ------

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











## Support Vector Machine Classifier ------

# Define the training data for SVM (1-mer features).
train_data_1mer <- dfTraining[, c("Gene", "Aprop", "Tprop", "Gprop", "Cprop")]
train_data_1mer$Gene <- as.factor(train_data_1mer$Gene)

# Define the validation data for SVM (1-mer features).
validation_data_1mer <- dfValidation[, c("Gene", "Aprop", "Tprop", "Gprop", "Cprop")]
validation_data_1mer$Gene <- as.factor(validation_data_1mer$Gene)

# Train an SVM model using caret (train()) with a radial kernel for 1-mer features.
set.seed(123)  # Setting a seed for reproducibility.
svm_1mer_caret <- train(
  Gene ~ Aprop + Tprop + Gprop + Cprop,
  data = train_data_1mer,
  method = "svmRadial",  # Radial kernel for SVM.
  trControl = trainControl(method = "cv", number = 10),  # 10-fold cross-validation.
  preProcess = c("center", "scale")  # Standardize the data.
)

# Use the trained model to predict the validation dataset (1-mer).
svm_1mer_predictions_caret <- predict(svm_1mer_caret, validation_data_1mer)
svm_1mer_accuracy_caret <- sum(svm_1mer_predictions_caret == validation_data_1mer$Gene) / nrow(validation_data_1mer)
print(paste("SVM Accuracy with 1-mer (train()):", svm_1mer_accuracy_caret))

# Define the training data for SVM (2-mer features).
train_data_2mer <- dfTraining[, c("Gene", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")]
train_data_2mer$Gene <- as.factor(train_data_2mer$Gene)

# Define the validation data for SVM (2-mer features).
validation_data_2mer <- dfValidation[, c("Gene", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")]
validation_data_2mer$Gene <- as.factor(validation_data_2mer$Gene)

# Train an SVM model using caret with a radial kernel for 2-mer features.
set.seed(123)  # Setting a seed for reproducibility.
svm_2mer_caret <- train(
  Gene ~ AA + AC + AG + AT + CA + CC + CG + CT + GA + GC + GG + GT + TA + TC + TG + TT,
  data = train_data_2mer,
  method = "svmRadial",  # Radial kernel for SVM.
  trControl = trainControl(method = "cv", number = 10),  # 10-fold cross-validation.
  preProcess = c("center", "scale")  # Standardize the data.
)

# Use the trained model to predict the validation dataset (2-mer).
svm_2mer_predictions_caret <- predict(svm_2mer_caret, validation_data_2mer)
svm_2mer_accuracy_caret <- sum(svm_2mer_predictions_caret == validation_data_2mer$Gene) / nrow(validation_data_2mer)
print(paste("SVM Accuracy with 2-mer (train()):", svm_2mer_accuracy_caret))




# Ensure predictions are added to dfValidation.
dfValidation$Predicted_Gene_1mer <- svm_1mer_predictions_caret

# Create the Misclassified column.
dfValidation$Misclassified <- ifelse(
  dfValidation$Gene == "rbcL" & dfValidation$Predicted_Gene_1mer == "matK",
  "Misclassified",
  "Classified"
)

# Create a boxplot for G proportions with the misclassified sequences labeled for 1-mer predictions.
ggplot(dfValidation, aes(x = Gene, y = Gprop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Misclassified), width = 0.4, size = 2) +
  scale_color_manual(values = c("Classified" = "green", "Misclassified" = "red")) +
  labs(
    title = "Boxplot of Gprop for matK and rbcL Sequences Based on SVM Predictions (1-mer Predictions)",
    x = "Gene",
    y = "Gprop"
  ) +
  theme_minimal()




## Model Comparison ------

# Returning to one of the main objectives of the assignment, we wish to compare our Random Forest accuracy and SVM accuracy using the predicted data for validation and counting the correctly classified sequences.

# Calculate Random Forest accuracy (1-mer).
RandomForestAccuracy <- sum(predictValidation == dfValidation$Gene) / nrow(dfValidation)

# Calculate Random Forest accuracy (2-mer).
RandomForestAccuracy_2mers <- sum(predictValidation_2mers == dfValidation$Gene) / nrow(dfValidation)

# Calculate SVM accuracy (1-mer).
SVMAccuracy_1mer <- sum(svm_1mer_predictions_caret == dfValidation$Gene) / nrow(dfValidation)

# Calculate SVM accuracy (2-mer).
SVMAccuracy_2mer <- sum(svm_2mer_predictions_caret == validation_data_2mer$Gene) / nrow(validation_data_2mer)

# Print the accuracy results for comparison.
print(paste("Random Forest Accuracy (1-mer):", RandomForestAccuracy))
print(paste("Random Forest Accuracy (2-mer):", RandomForestAccuracy_2mers))
print(paste("SVM Accuracy (1-mer, caret):", SVMAccuracy_1mer))
print(paste("SVM Accuracy (2-mer, caret):", SVMAccuracy_2mer))


# Observations (using the validation dataset):
# - Random Forest with 1-mer = an accuracy of 0.995 (99.5%).
# - Random Forest with 2-mer = an accuracy of 1 (100%).
# - SVM with 1-mer (caret) = an accuracy of 0.993 (99.29%).
# - SVM with 2-mer (caret) = an accuracy of 0.999 (99.88%).

#Random Forest with 2-mer achieved perfect classification accuracy, suggesting that including dinucleotide frequencies significantly improves performance. SVM (caret) also performed very well, with 2-mer features at 99.88%. 
#The results indicate that Random Forest may slightly outperform SVM for this dataset when using 2-mer features.

