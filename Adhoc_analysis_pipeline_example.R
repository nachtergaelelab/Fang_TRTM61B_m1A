### PURPOSE OF THIS SCRIPT
## Ad-hoc cutoff based hit identification pipeline
## Last modified date: October 29, 2025
## Author: Isaac Vock

# Load dependencies ------------------------------------------------------------

##### Packages #####

library(data.table)
library(dplyr)
library(arrow)
library(tidyr)
library(broom)
library(stringr)
library(dtplyr)
#install.packages('R.utils')
library(R.utils)


# Path to mutational data table: 
# Columns in example table: sample = replicate name,rname = transcript,gloc = exact SNT position,n = number of A mutations,trials = total reads at position,DART = sample used for DART or not,Enzyme = NE, LE, or HE treatment
mutdata_path <- "path_to_mutationcounts_file.csv.gz"

##### Functions #####

### Workflow

# For uncertainty quantification
beta_se <- function(n, t){
  
  alpha <- n + 1
  beta <- (t - n) + 1
  
  num <- alpha * beta
  den <- ((alpha + beta)^2) * (alpha + beta + 1)
  
  return(
    sqrt(
      num / den
    )
  )
  
}

#' @param lazy_cU Table of mutation rates. Called "lazy" because by default
#' this script lazily loads it using a fancy package called dtplyr
#' @param high_mutrate_cutoff There needs to be at least `n_hits` replicates with
#' a mutation rate greater than `high_mutrate_cutoff` to get called a hit
#' @param low_mutrate_cutoff All "NE" treatment replicates need to have a mutation
#' rate below this to get called a hit
#' @param treatment Do you want to analyze the "HE" or "LE" treatment data?
#' @param min_trials All "NE" replicates and at least `n_hits` treated replicates
#' need to have this many reads at a site for it to be a hit
#' @param min_high_n At least `n_hits` treated replicates need to have this many
#' mutations at the site for it to be a hit
#' @param DART If `TRUE`, only DART samples are analyzed; if `FALSE`, only non-DART
#' samples are considered; set this to anything else (I arbitrarily chose the number
#' `2` in this script) to use all samples
#' @param n_hits Number of replicates for which above filters need to be passed
#' for a site to be called a hit 
id_adhoc_hits <- function(lazy_cU,
                          high_mutrate_cutoff = 0.01,
                          low_mutrate_cutoff = 0.01,
                          treatment = c("HE", "LE"),
                          min_trials = 1,
                          min_high_n = 1,
                          DART = TRUE,
                          n_hits = 1){
  
  treatment <- match.arg(treatment)
  
  if(DART == TRUE){
    lazy_cU <- lazy_cU %>%
      dplyr::filter(DART)
  }else if(DART == FALSE){
    lazy_cU <- lazy_cU %>%
      dplyr::filter(!DART)
  }
  
  adhoc_fit <- lazy_cU %>%
    dplyr::group_by(
      rname, gloc
    ) %>%
    dplyr::filter(
      (sum( (n / trials > high_mutrate_cutoff) & Enzyme == treatment & trials > min_trials & n > min_high_n & Enzyme != "NE") >= n_hits) &
        all(Enzyme != "NE" | (n / trials < low_mutrate_cutoff & Enzyme == "NE" & trials > min_trials))
    ) %>%
    dplyr::group_by(
      rname, gloc, Enzyme
    ) %>%
    # For case where there are multiple samples, just get maximum mutation rate in each
    dplyr::filter(
      n / trials == max(n / trials)
    ) %>%
    dplyr::filter(
      row_number() == 1
    ) %>%
    dplyr::group_by(
      rname, gloc
    ) %>%
    dplyr::summarise(
      treatment_mutrate = n[Enzyme == treatment] / trials[Enzyme == treatment],
      NE_mutrate = n[Enzyme == "NE"] / trials[Enzyme == "NE"],
      treatment_trials = trials[Enzyme == treatment],
      NE_trials = trials[Enzyme == "NE"],
      treatment_mutcount = n[Enzyme == treatment],
      NE_mutcount = n[Enzyme == "NE"]
    ) %>%
    as_tibble()
  
  adhoc_fit <- adhoc_fit %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      treatment_se = beta_se(
        treatment_mutcount, treatment_trials
      ),
      NE_se = beta_se(
        NE_mutcount, NE_trials
      ),
      total_se = sqrt(treatment_se^2 + NE_se^2),
      difference = treatment_mutrate - NE_mutrate,
      stat = difference / total_se,
      # Calculate log p-value to avoid numerical underflow
      log_pval = log(2) + pnorm(-abs(stat), lower.tail = TRUE, log.p = TRUE),
      # Cap log_pval at 0 (p-value at 1)
      log_pval = pmin(log_pval, 0),
      # Convert to -log10 space for better resolution
      neglog10_pval = -log_pval / log(10),
      # For standard p-value column (with floor to avoid exact 0)
      pval = pmax(exp(log_pval), .Machine$double.xmin)
    ) %>%
    dplyr::arrange(
      log_pval  # Sort by log p-value for accuracy
    ) %>%
    dplyr::mutate(
      # Perform BH correction in log-space to preserve resolution
      # BH correction: padj = pval * n / rank
      n_tests = n(),
      rank = row_number(),
      log_padj = log_pval + log(n_tests) - log(rank),
      # Cap at 0 (adjusted p-value at 1)
      log_padj = pmin(log_padj, 0),
      # Enforce monotonicity (adjusted p-values should increase with rank)
      log_padj = cummax(log_padj),
      # Convert to -log10 space
      neglog10_padj = -log_padj / log(10),
      # Standard adjusted p-value column
      padj = pmax(exp(log_padj), .Machine$double.xmin)
    ) %>%
    dplyr::select(-n_tests, -rank) %>%
    dplyr::select(rname, gloc, difference, pval, padj, everything()) %>%
    dplyr::arrange(
      pval
    )
  
  return(adhoc_fit)
}


# Identify hits ----------------------------------------------------------------

##### Load data #####

cU <- fread(mutdata_path)

lazy_cU <- lazy_dt(as_tibble(cU))


### Unit test

test_df <- tibble(
  sample = c(),
  rname = "chr1",
  gloc = 1,
  n = c(200, 200, 1,
        1, 1, 2),
  trials = 1000,
  DART = c(TRUE, TRUE, FALSE,
           TRUE, TRUE, FALSE),
  Enzyme = rep(c("HE", "NE"), each = 3)
)

id_adhoc_hits(
  test_df,
  treatment = "HE",
  n_hits = 2,
  DART = 2
)


### Low enzyme treatment, 2 or more replicates across all samples

adhoc_fit_LE <- id_adhoc_hits(
  lazy_cU,
  treatment = "LE",
  n_hits = 2,
  DART = 2 # silly hack to say use all samples, DART and non-DART
)

adhoc_fit_LE <- adhoc_fit_LE %>%
  dplyr::rename(
    LE_mutcount = treatment_mutcount,
    LE_trials = treatment_trials,
    LE_se = treatment_se
  )

adhoc_fit_LE <- adhoc_fit_LE %>%
  filter(!gloc %in% c(17, 19, 20))

fwrite(adhoc_fit_LE,
       "example_LE.csv"
)

### High enzyme treatment 2 or more replicates across all samples

adhoc_fit_HE <- id_adhoc_hits(
  lazy_cU,
  treatment = "HE",
  n_hits = 2,
  DART = 2
)

adhoc_fit_HE <- adhoc_fit_HE %>%
  dplyr::rename(
    HE_mutcount = treatment_mutcount,
    HE_trials = treatment_trials,
    HE_se = treatment_se
  )

adhoc_fit_HE <- adhoc_fit_HE %>%
  filter(!gloc %in% c(17, 19, 20))

fwrite(adhoc_fit_HE,
       "example_HE.csv"
)

overlap_table <- inner_join(
  adhoc_fit_LE, adhoc_fit_HE,
  by = c("rname", "gloc")
)

fwrite(overlap_table,
       "example_overlap.csv"
)



