#!/usr/bin/env Rscript

#######  DESCRIPTION   ###########
# Filtering palindromes that are considered as stable through the 'd' score and returns a BED file with
# palindromic regions and features related to these regions for downstream analyzing.

# The 'd' score assess the recombinogenicity of a palindrome given its features. The reference paper:
# "Long inverted repeats in eukaryotic genomes: Recombinogenic motifs determine genomic plasticity"
# Only palindromes having spacer are considered by this score, the other ones can be filtered by their size.

# The following features are taken in account to determine the 'd' score:
# Note that a palindrome is composed of: left_stem + spacer (or not) + rigth_stem with len rigth_stem == left_stem
# 1) Stem length in bp: length of one arm
# 2) Spacer length in bp: length of the sequence in between left and right stems
# 3) Sequence mismatching ratio: ratio mismatch over length between right and left stem

# Then the score is computed as follow:
# d = STEM.LEN/SPACER.LEN

# A palindrome is stated as unstable if this 'd' score is grater or equal to its mismatching ratio. Also in the case of
# a palindrome with a spacer and an exact similarity we consider the miss.rate as 1 at least, so there is a selection
# against those kind of palindromes. By choosing, one, we require that the maximum length for a spacer is the stem
# length of the corresponding palindrome.

#######  Dependencies   ###########
library(dplyr)

#######  FUNCTIONS   ###########
# Check the "d" score introduced above and return TRUE if it is d >= miss.rate
compute_recombinogenic <- function(df){
  # Stem length = (Palindrom length - Spacer length)/2
  # d = ((length.pal - length.spacer)/2)/length.spacer
  # The 'd' score is only computed over spacer present palindromes
  ifelse(df$Sp.len > 0,
         round(((df$Pal.len - df$Sp.len)/2)/df$Sp.len, digits = 2),
         0)
}


# Select the d's that are considered as recombinogenic according to the definition above.
select_recombinogenic <- function(df){
  # Minimum miss.rate is 1 for palindromes with spacer (Maximum length for a spacer is the stem length -> ratio of 1)
  df$Miss.rate <- ifelse(df$Sp.len > 0 & df$Miss.rate == 0,
                         1,
                         df$Miss.rate)
  # Select the lines where the d is equal to 0 (we keep perfect and nearperfect palindromes without sp) and
  # d >= miss.rate
  filtered_df <- df[df$Recomb == 0 | df$Recomb >= df$Miss.rate,]
  return(filtered_df)
}

#######  MAIN INSTANCE FOR FILTERING   ###########
# Data file by chromosome -> BED format
file <- snakemake@input[['chr_bed']]
output <- snakemake@output[[1]]
# Store in a data frame
palindromes <- read.table(file, sep = "\t", header = FALSE)
# Build the header
names(palindromes) <- c("Chr",
                        "Start",
                        "End",
                        "Type",
                        "Pal.len",
                        "Sp.len",
                        "Miss.rate",
                        "AT",
                        "Size")

# Remove Size column because will be newely defined
palindromes.1 <- palindromes[, 1:8]
# Replace percentage signs by '' ans transform value to numeric
palindromes.1$Miss.rate <- as.numeric(sub(pattern = "%",
                                          replacement = "",
                                          x = palindromes$Miss.rate))
# Compute the d
palindromes.1$Recomb <- compute_recombinogenic(palindromes.1)
# Select the ones that recombinogenic
final.palindromes <- select_recombinogenic(palindromes.1)
# Rewrite column miss rate in percentage
final.palindromes$Miss.rate <- paste0(round(final.palindromes$Miss.rate, digits = 2),"%")
# Annotate spacers
final.palindromes$Type <- ifelse(final.palindromes$Sp.len > 0,
                                 "Spacer",
                                 final.palindromes$Type)
# Categorize lengths of palindromes by ranges to easily select subgroups for downstream analyzes
final.palindromes = final.palindromes %>%
  mutate(Size = case_when(
    Pal.len <= 50 ~ "0-50 bp",
    Pal.len %in% 51:99 ~ "51-99 bp",
    Pal.len %in% 100:200 ~ "100-200 bp",
    Pal.len > 200 ~ ">200 bp"
  ))

# Write resulting file -> f = filtered
write.table(final.palindromes,
            output,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)