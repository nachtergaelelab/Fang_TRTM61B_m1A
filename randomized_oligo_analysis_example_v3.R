### PURPOSE OF THIS SCRIPT
## Count occurrences of different nucleotides at randomized positions
## Last modified date: October 2, 2025
## Author: Richard Li

library(dplyr)

##############################################
# DEAL WITH INSERTIONS WHICH MESS WITH INDICES
##############################################

# Finds all the starting indices of insertions in a target sequence
find_insertions <- function(targ_seq){
  
  is_lowercase <- unlist(strsplit(targ_seq,'')) %in% letters
  insertion_indices <- list()
  insertion_lengths <- list()
  
  # Loop that identifies location of first insertion + length, then removes it
  # from is_lowercase, looping until no remaining insertions.
  while(TRUE %in% is_lowercase){
    
    # What is the start index of the first remaining insertion?
    current_start_index <- which(is_lowercase == TRUE)[1]
    
    # Set the end index of this insertion to the start index, and increment
    # index until we find end of insertion (no longer lowercase)
    current_end_index <- current_start_index
    while(is_lowercase[current_end_index + 1] == TRUE){
      current_end_index <- current_end_index + 1
      # If we reach the end of the whole string break and do not continue
      if(current_end_index == length(is_lowercase)) break
    }
    
    current_length <- current_end_index - current_start_index + 1
    
    # Append current index + length to vectors
    insertion_indices <- append(insertion_indices, current_start_index)
    insertion_lengths <- append(insertion_lengths, current_length)
    
    # Remove the current insertion from is_lowercase and start loop again
    is_lowercase <- is_lowercase[-c(current_start_index:current_end_index)]
  }
  
  # In this application, we only care about the starting indices for the
  # insertions
  return(unlist(insertion_indices))
}

# Find if target sequence insertions fall within focal indices
flag_insertions <- function(df, FOCAL_INDEX_START, FOCAL_INDEX_END){
  
  df$insertion_locations <- lapply(df$TargetSequence, find_insertions)
  # Is there an insertion in the focal region of the sequence
  df$insertion_focal <- ifelse(df$insertion_locations %in% FOCAL_INDEX_START:FOCAL_INDEX_END, TRUE, FALSE)
  
  return(df)
}


###########
# MAIN CODE
###########
#The input files for this is a sheet that has counted the occurrences of each unique sequence in the sequencing data 
#"TargetSequence" column contains actual sequences
#"Reads" column contains counts per sequence
data <- readxl::read_xlsx("path_to_file.xlsx", sheet = 1)

# Start and end indices indicate the indices with which we want to subset our
# target sequence
FOCAL_INDEX_START <- 114
FOCAL_INDEX_END <- 120

# Indicates index of the modification site, WITHIN THE ALREADY TRIMMED SUBSET
MOD_SITE <- 4
UPSTREAM_1 <- 1
UPSTREAM_2 <- 3
DOWNSTREAM_1 <- 5
DOWNSTREAM_2 <- 7

# Find if any insertions fall within the region of interest
data <- flag_insertions(data, FOCAL_INDEX_START, FOCAL_INDEX_END)

# Let's only consider sequences without insertions in the focal region
data2 <- data %>% filter(insertion_focal == FALSE)

# Remove all insertions so that all sequences are aligned in index
data2$TargetSequenceTrimmed <- substr(gsub("[a-z]", "", data2$TargetSequence),
                                      FOCAL_INDEX_START, FOCAL_INDEX_END)

# Remove all sequences that don't match the desired fixed sequence upstream or downstream - comment out the one you don't want
data2 <- data2 %>% filter(substr(data2$TargetSequenceTrimmed, start = UPSTREAM_1, stop = UPSTREAM_2) == "TCA")
#data2 <- data2 %>% filter(substr(data2$TargetSequenceTrimmed, start = DOWNSTREAM_1, stop = DOWNSTREAM_2) == "GAT")

# Add new column representing whether the modification site has been modified
# or not
unmodified <- data2 %>% filter(substr(data2$TargetSequenceTrimmed, start = MOD_SITE,
                                   stop = MOD_SITE) == "A")
modified <- data2 %>% filter(substr(data2$TargetSequenceTrimmed, start = MOD_SITE,
                                 stop = MOD_SITE) != "A")
unmodified$modified <- FALSE
modified$modified <- TRUE
data2 <- rbind(unmodified, modified)

# New version of data2 where all deletions (in focal region) are filtered out
data3 <- data2
data3$has_deletion <- grepl("-", data3$TargetSequenceTrimmed)
data3 <- data3 %>% filter(has_deletion == FALSE)

##########################################
# FIND NUCLEOTIDE RAW FREQUENCIES BY INDEX
# From ChatGPT, vetted by Richard
##########################################
#install.packages("tidyr")
# Split strings into individual characters and create a long format
library("tidyr")
freqs_by_site <- tibble(string = data3$TargetSequenceTrimmed, modified = data3$modified) %>%
  # Split each string into characters (produces a list)
  mutate(id = row_number(), 
         nucleotide = strsplit(string, "")) %>%
  # Unnest the list into rows
  unnest_longer(nucleotide) %>%
  # Label with index within string
  group_by(id) %>% mutate(position = row_number()) %>% ungroup() %>%
  
  # Count frequency of each character at each position
  count(position, nucleotide, modified, name = "frequency") %>%
  
  # Divide each frequency by group sum (position x modified status total
  # frequency)
  group_by(position, modified) %>%
  mutate(percent = frequency / sum(frequency)) %>%
  ungroup() 

#Write out tables of interest
#install.packages("writexl")
library(writexl)
write_xlsx(freqs_by_site, path = "example_output_data.xlsx")

##########################################################
# FIND FREQUENCIES OF DIFFERENT RANDOMIZATION PERMUTATIONS
##########################################################
# Compares the frequencies of randomized nucleotides surrounding modification
# site, with modification vs. without modification
# -- df: data frame including column TargetSequenceTrimmed
# -- MOD_SITE: index of modification site within TargetSequenceTrimmed
# -- locations: indices on which to summarize permutation frequencies, within
#               TargetSequenceTrimmed, as a numeric vector

randomization_frequencies <- function(df, MOD_SITE, locations){
  
  # Throw an error if MOD_SITE is within locations
  if(MOD_SITE %in% locations) stop("Error: Index of modification site is
                                   within set of indices to permute on")
  
  # Add a column representing the indices assessed (permutation group)
  df$permutation_group <- paste(locations, collapse = "|")
  
  # Function that takes in a target sequence string and masks nucleotides
  # not at locations of interest
  # Returns '-' when masked, X for MOD_SITE, and nucleotide value at locations
  # of interest
  mask_targ_seq <- function(sequence, MOD_SITE, locations){
    # Split string into component characters
    sequence_split <- unlist(strsplit(sequence, split = ""))
    # Replace all indices within locations with "_", can't use dash since
    # used for deletions
    sequence_split[-locations] <- "_"
    # Replace modification site with "X"
    sequence_split[MOD_SITE] <- "X"
    # Rejoin vector of characters into a single string
    sequence_masked <- paste(sequence_split, collapse = "")
    return(sequence_masked)
  }
  
  # Use masking function to create a column that indicates the unique
  # permutation on which to summarize frequencies
  df$permutation <- lapply(df$TargetSequenceTrimmed, mask_targ_seq,
                           MOD_SITE, locations)
  
  # Count the number of reads that fall into each unique permutation at
  # specified locations. Grouping by permutation_group is not actually
  # meaningful, we just want to include it in output
  frequency_table <- df %>% group_by(permutation_group, permutation, modified) %>%
    summarize(reads = sum(Reads)) %>% ungroup() %>%
    
    # Get number of reads as a percentage of total reads across all permutations, within modification status
    group_by(modified) %>% mutate(percent_across_permutations = reads / sum(reads)) %>% ungroup() %>%
    
    # Get number of reads as a percentage of total reads within permutation, across modification status
    group_by(permutation) %>% mutate(percent_for_sequence = reads / sum(reads)) %>% ungroup()
  
  return(frequency_table)
}

permutation_set_1 <- randomization_frequencies(df = data3, MOD_SITE = MOD_SITE, locations = c(7))
permutation_set_2 <- randomization_frequencies(df = data3, MOD_SITE = MOD_SITE, locations = c(5,7))
permutation_set_3 <- randomization_frequencies(df = data3, MOD_SITE = MOD_SITE, locations = c(5,6,7))

# Can store all permutation analyses in a single table, permutations_group indicates which analysis
permutations_table <- rbind(permutation_set_1, permutation_set_2, permutation_set_3)




