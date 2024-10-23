# Project: Mycophen Root Year 1: ITS QC to create clean and rarefied data table
# Author: Peter Kennedy
# Date: 08/16/24

# load libraries
library(here)
library(dplyr)
library(vegan)
library(phyloseq)
library(tidyverse)
library(stringr) ## for manipulating character strings 

# set working directory
setwd("/Users/kennedyp/Library/CloudStorage/GoogleDrive-kennedyp@umn.edu/My Drive/Research_Active/MycoPhen/Year1_Root_ITS")

# Load in ASV table (csv, txt, tab delimited file, etc.)
data <- read.delim(here("/Users/kennedyp/Library/CloudStorage/GoogleDrive-kennedyp@umn.edu/My Drive/Research_Active/MycoPhen/Year1_Root_ITS/MycoPhen_Root_ITS_combined_sequences_taxa_bootstrap.txt"), sep="\t", header=TRUE)

####Step 1 - tidy up the raw data table

# Add a header to sequence column (currently column 1)
colnames(data)[1] <- "Sequence"

# Re-organize the columns to move "Sequence" column to last 
data2 <- data %>%
  select(-1, everything(), 1)

# To add an ASV_ID column, need to know number of rows 
dim(data2)

# Add a column called "asv_id" and place it as the first column in the data frame
data3 <- data2 %>%
  mutate(asv_id = paste0("ASV", sprintf("%04d", 1:6464))) %>%  # the sprintf helps keeping the zeroes before the numbers
  select(asv_id, everything())

# looking at the column names to make adjustments
colnames(data3)

# Edit column names for bootstrap values
colnames(data3)[81] <- "Kingdom_bootstrap"
colnames(data3)[82] <- "Phylum_bootstrap"
colnames(data3)[83] <- "Class_bootstrap"
colnames(data3)[84] <- "Order_bootstrap"
colnames(data3)[85] <- "Family_bootstrap"
colnames(data3)[86] <- "Genus_bootstrap"
colnames(data3)[87] <- "Species_bootstrap"

# Looking at the column names to check adjustments
colnames(data3)

## Will make additional adjustments to sample names but later when dealing with just sample columns

# Remove the "k__" from the Kingdom column in the taxonomic data:
data3$Kingdom <- gsub(".*k__","", data3$Kingdom) # removing characters before k__
# remove the "p__" from the Phylum column in the taxonomic data:
data3$Phylum <- gsub(".*p__","", data3$Phylum) # removing characters before p__
# remove the "c__" from the Class column in the taxonomic data:
data3$Class <- gsub(".*c__","", data3$Class) # removing characters before c__
# remove the "o__" from the Order column in the taxonomic data:
data3$Order <- gsub(".*o__","", data3$Order) # removing characters before o__
# remove the "f__" from the Family column in the taxonomic data:
data3$Family <- gsub(".*f__","", data3$Family) # removing characters before f__
# remove the "g__" from the Genus column in the taxonomic data:
data3$Genus <- gsub(".*g__","", data3$Genus) # removing characters before g__
# remove the "s__" from the species column in the taxonomic data:
data3$Species <- gsub(".*s__","", data3$Species) # removing characters before s__

#####Step 2 - Separate Sequence and Taxonomic Datasets in ASV table

# Before separating seq and tax data, sort to only keep Fungi

# List which kingdoms are present
table(data3$Kingdom)

# Now subset out kingdom to only keep Fungi
data4.fungi <-data3 %>%
  dplyr::filter(Kingdom=="Fungi")

# confirm Fungi are the only kingdom left
table(data4.fungi$Kingdom)

# ID which columns to keep that have sequence data
colnames(data4.fungi)

# now subset out columns with the sequence data only, but make sure to include asv_id
data4.fungi.seq <- data4.fungi[,c(1:73)] %>%
  tibble::column_to_rownames("asv_id")

# confirm the number of ASVs and samples:
dim(data4.fungi.seq)

# ID which columns to keep that have taxonomy data
colnames(data4.fungi)

# now subset out columns with the taxonomic data only, but make sure to include asv_id 
data4.fungi.tax <- data4.fungi[,c(1,74:88)] %>%
  tibble::column_to_rownames("asv_id")

#check dimensions
dim(data4.fungi.tax) #should match sequence matrix for rows

####Step 3 - Account for index bleed

## In this table, we don't have any positive or negative controls, so we have to guess with regards to possible index bleed effects. I am going to use the old advice of Lindahl et al. (2013) which suggested zeroing cells with 10 or less reads. It is arbitrary, but better than nothing.

# Converting values in any cell with a value of <=10 to zero 
data4.fungi.seq.thres <- data4.fungi.seq %>%
  mutate(across(everything(), ~ ifelse(. <= 10, 0, .)))

# Check if all values in each cell have counts above 10 (the outcome should say FALSE for all samples)
counts_above_10 <- data4.fungi.seq.thres %>%
  summarize(across(everything(), ~ sum(. > 10, na.rm = TRUE) == n()))

# Print result
print(counts_above_10)

# Eliminating any ASVs that now have 0 reads
data4.fungi.seq.thres <- data4.fungi.seq.thres %>%
  filter(rowSums(across(everything())) != 0) 

#Check new dimensions
dim(data4.fungi.seq.thres)
##Went down from 4605 to 2653 ASVs.

# Renaming trimmed and thresholded dataframe
data5.fungi.seq <- data4.fungi.seq.thres

#check total number of fungal reads
sum(data5.fungi.seq) 
##5,244,824 total reads

#####Step 4 - Assess sample sequence count totals

## Some samples can have very few total sequences and as such are not likely good to include in the final analyses. However, determining where to cut-off samples in terms of read totals is notably subjective.  

# Create a list of column sums
column_sums <- data5.fungi.seq %>%
  summarize(across(everything(), sum))
            
# Visualize the distribution of counts across all samples
column_sums <-as.numeric(unlist(column_sums))
hist(column_sums)

# Re-create a list of column sums
column_sums <- data5.fungi.seq %>%
  summarize(across(everything(), sum))

# Pivot the column sums object to a long list
column_sums_long <- column_sums %>%
  pivot_longer(cols = everything(), names_to = "Column_Name", values_to = "Sums")

# Print the ranked list of 50 lowest column sums
column_sums_ranked <- column_sums_long %>%
  mutate(rank = rank(Sums)) %>%
  arrange(Sums) %>%
  head(50)

print(column_sums_ranked)

# These samples have the lowest read counts
#1 X6CA6.S633.L001    57
#2 X4PI3.S654.L001   448 
#3 X8PI2.S672.L001   741 
#4 X8PI1.S671.L001  1013 
#5 X8CA1.S665.L001  1420 
#6 X6PI2.S631.L001  1977   
#7 X8PI6.S676.L001  2456   
#8 X8PI3.S673.L001  3069   
#9 X6PI3.S630.L001 10629 
#10 X7CA1.S626.L001 10839  

#Could drop first 3 and rarefy to 1013, or could drop first 8 and rarefy to 10629
#Looks like PI is the plot with lots of the low read totals. To keep as many PI plots as possible, I recommend rarefying to 1013.

# Set the threshold for column sums
threshold <- 1013

# Drop columns with sums less than the threshold 
data6.fungi.seq <- data5.fungi.seq %>%
  select(where(~ sum(.) >= threshold))

dim(data6.fungi.seq)
#That leaves 68 samples

####Step 4 - tidying up the taxonomy dataframe to match the final sequencing dataframe

# Convert taxonomy df row names to a column
df_tax <- rownames_to_column(data4.fungi.tax, var = "RowName")

# Convert sequence df row names to a column
df_seq <- rownames_to_column(data6.fungi.seq, var = "RowName")

# Perform the join based on df_seq
df_joined <- df_seq %>%
  inner_join(df_tax, by = "RowName")

# list column names
colnames(df_joined)

# change first column name to "asv_id"
colnames(df_joined)[1] <- "asv_id"

data5.fungi.tax <- df_joined[,c(1,71:85)] %>%
  tibble::column_to_rownames("asv_id")

data6.fungi.seq <- df_joined[,c(1:69)] %>%
  tibble::column_to_rownames("asv_id")

#####Step 5 - add FungalTraits data to the taxonomic dataset (based off code from Talia)

dim(data5.fungi.tax)
colnames(data5.fungi.tax)

## Load in FungalTraits database:
traits <- read.csv("/Users/kennedyp/Library/CloudStorage/GoogleDrive-kennedyp@umn.edu/My Drive/Kennedy Lab - Shared Folder/Sequencing - HTS/Adding FungalTraits Data/FungalTraitsDatabase.csv")

colnames(traits)

##rename GENUS to Genus
traits <-traits %>%
  dplyr::rename("Genus"= "GENUS")
names(traits)

###
names(data5.fungi.tax)
table(data5.fungi.tax$Genus)

##Note There are a few cases where GENUS is NA:

traitsOfInterest <- select(traits, Genus, primary_lifestyle, Secondary_lifestyle,
                           Endophytic_interaction_capability_template,
                           Plant_pathogenic_capacity_template, Decay_substrate_template,
                           Decay_type_template, Aquatic_habitat_template,
                           Animal_biotrophic_capacity_template, Specific_hosts,
                           Growth_form_template, Ectomycorrhiza_exploration_type_template,
                           Ectomycorrhiza_lineage_template) 
traitsOfInterest

class(data5.fungi.tax$Genus)
data5.fungi.tax$Genus<-as.factor(data5.fungi.tax$Genus)
class(data5.fungi.tax$Genus)
traitsOfInterest$Genus<-as.factor(traitsOfInterest$Genus)

# joining tax with traits database
data6.fungi.tax <- left_join(data5.fungi.tax, traitsOfInterest, by="Genus") 

## unfortunately this loses the ASV_IDs as the row names, so we need to add those back.

# Save the row names from data5.fungi.tax
row_names <- rownames(data5.fungi.tax)

# Apply the row names to data5.fungi.tax.trait
rownames(data6.fungi.tax) <- row_names

# viewing the join
View(data6.fungi.tax)

dim(data6.fungi.tax)
dim(data5.fungi.tax)

table(data6.fungi.tax$primary_lifestyle)
table(data6.fungi.tax$Kingdom) 

####Step 6 - Cleaning up sample names to make it easy to match with the mapping file 

## Will need to first modify the column names to retain only the sample information in the sequence file so that it can be matched with the names in the mapping file. The function below keeps everything in front of the second period in the column name. 

# Function to process each column name
process_column_name <- function(col_name) {
  split_string <- strsplit(col_name, "\\.")[[1]]
  output_string <- paste(split_string[1:2], collapse = ".")
  return(output_string)}

# Apply the function to all column names
new_colnames <- sapply(colnames(data6.fungi.seq), process_column_name)

# Assign the new column names to the dataframe
colnames(data6.fungi.seq) <- new_colnames

# Print the modified dataframe
print(data6.fungi.seq)

# Getting rid of the X and the remaining info after the period
colnames(data6.fungi.seq) <- gsub("^X|\\..*$", "", colnames(data6.fungi.seq))

####Step 7 - Create a rarefied version of the sequence dataset

# Transposing matrix to match format for vegan
data6.fungi.seq.t <- t(data6.fungi.seq)

# Generate rarefaction curves for all samples to see ASV accumulation
rarecurve(data6.fungi.seq.t, step = 20, sample = raremax, label = FALSE, col = "blue", cex = 0.6, xlab = "Sequence Reads", ylab = "Number of ASVs", main = "ASV Accumulation Curve")

# Determine the lowest number of reads in a sample
raremax <- min(rowSums(data6.fungi.seq.t))
raremax

# Define the number of counts to rarefy to (e.g., 1000)
rarefaction_depth <- 1013

# Perform rarefaction
data6.fungi.seq.rar <- rrarefy(data6.fungi.seq.t, sample = rarefaction_depth)

# Save the rarefied data frame as a new object 
data6.fungi.seq.rar <- as.data.frame(data6.fungi.seq.rar)

# Transposing matrix back to original format
data6.fungi.seq.rar.t <- t(data6.fungi.seq.rar)

# Save the transposed rarefied data frame as a new object 
data6.fungi.seq.rar.t <- as.data.frame(data6.fungi.seq.rar.t)

# Convert taxonomy df row names to a column
df_tax <- rownames_to_column(data6.fungi.tax, var = "RowName")

# Convert sequence df row names to a column
df_seq <- rownames_to_column(data6.fungi.seq.rar.t, var = "RowName")

# Perform the join based on df_seq
df_joined <- df_seq %>%
  inner_join(df_tax, by = "RowName")

# list column names
colnames(df_joined)

# change first column name to "asv_id"
colnames(df_joined)[1] <- "asv_id"

# Save the new data frame to a CSV file
write.csv(df_joined, "ITS_data6.fungi.seq.tax.rar.csv", row.names = FALSE)

##############

