library(Biostrings)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

#Notes 15Dec23: Issues with alignment noted in df_error file. Unable to align 7043, which aligns in nextclade, need to attempt to align or determine
#cause of error. 8801 aligns in nextclade with frameshift noted. Alignment here produces erroneous sequence.
#Overall need to solve issues with poor alignments.

#upload concatenated sequences
seqs_fasta <- readDNAStringSet("C:/Users/vhogan/Documents/ARIA/SC2_paxlovid_res/fasta/A10-07043-0-Ag.fa")

#User input gene (orf1a)
gene = "orf1a"

#location of fasta reference files
ref_file = case_when(gene == "orf1a" ~ "C:/Users/vhogan/Documents/ARIA/SC2_paxlovid_res/fasta/sc2_ref.fasta")

#read in ref file as DNAStringSet bases on gene
ref_fasta = readDNAStringSet(ref_file)

######### some functions to translate gapped alignments:

## getCodons - a function to split sequences into codons.
# input (myAln) is a DNAStringSet with a gapped alignment
# output is a simple list, one element for each sequence. Each list element is a character vector of each codon
getCodons <- function(myAln) {
  seqs <- as.character(myAln)
  len <- width(myAln)[1]
  starts <- seq(from=1, to=len, by=3)
  ends <- starts + 2
  myViews <- lapply(myAln, function(x) { 
    Views(x, starts, ends)
  })
  myCodons <- lapply(myViews, function(x) {
    as.character(DNAStringSet(x))
  })
  myCodons
}
  
## translateCodons - takes a character vector of codons as input, outputs the corresponding amino acids
translateCodons <- function(myCodons, unknownCodonTranslatesTo="-") {
  ## make new genetic code
  gapCodon <- "-"
  names(gapCodon) <- "---"
  my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)
  
  ## translate the codons
  pep <- my_GENETIC_CODE[myCodons]
  
  ## check for codons that were not possible to translate, e.g. frameshift codons
  if (sum(is.na(pep))>0) {
    cat("\nwarning - there were codons I could not translate. Using this character", unknownCodonTranslatesTo, "\n\n")
    pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
  }
  
  ## prep for output
  pep <- paste(pep, collapse="")
  return(pep)
}

## wrap the getCodons and translateCodons functions together into one:
translateGappedAln <- function(myAln, unknownCodonTranslatesTo="-") {
  myCodons <- getCodons(myAln)
  myAAaln <- AAStringSet(unlist(lapply(myCodons, translateCodons, unknownCodonTranslatesTo=unknownCodonTranslatesTo)))
  return(myAAaln)
}


#Function to trim fastas to specified length based on length of gene segment
trim_to_roi <- function(gene, fasta){
#Select  orf1a portion of gene https://www.ncbi.nlm.nih.gov/nuccore/1798172431
    if (gene == "orf1a"){
      trimmed_fasta <<- subseq(fasta, start = 266, end = 21555)
    }
}

#trim, translate, and put ref fasta into new dataframe
trim_to_roi(gene, ref_fasta)
trimmed_ref <- trimmed_fasta
ref_aa <- translate(trimmed_ref, if.fuzzy.codon = "solve")
ref <- as.character(ref_aa)
ref_commas <- vapply(strsplit(ref, ""), function(x) paste(x, collapse=","), character(1L)) #add commas between aa
ref_sep <- unlist(strsplit(ref_commas, ","))
ref_name = names(ref_fasta) #name of reference sequence

n = width(ref_aa) #number of aa in ref sequence
names = c("name", as.character(1:n)) #names for columns of data frame

df <- data.frame(matrix(NA,    # Create empty data frame
                        nrow = 0,
                        ncol = n+1)) %>%
  rbind(c(ref_name, ref_sep))  #add ref sequence as first row

colnames(df) <- names #add column names to df

df_error <- data.frame(matrix(NA,    # Create empty data frame for alignment error output
                        nrow = 0,
                        ncol = 2))

colnames(df_error) <- c("SampleID", "Error")


#function to align, trim, and translate, and add sequence to df. Input untrimmed fasta and ref_fasta
fasta_to_df <- function(fasta, ref){
  trim_N_fasta <- subseq(fasta, start = 56) #start 56 to removed Ns at beginning of sequence
  alignment <- pairwiseAlignment(ref, trim_N_fasta)
  aligned_seq <- alignedSubject(alignment) #DNAStringSet with aligned sequence
  trim_to_roi(gene, aligned_seq)
  trimmed_seq <- trimmed_fasta
  aligned_seq_aa <- translateGappedAln(trimmed_seq) #translate trimmed, aligned sequence allow for gap characters
  seq_char <- as.character(aligned_seq_aa) #convert translated sequence to character
  seq_commas <- vapply(strsplit(seq_char, ""), function(x) paste(x, collapse=","), character(1L)) #add commas between aa
  seq_sep <- unlist(strsplit(seq_commas, ",")) #unlist to add to df
  seq_name = names(fasta) #name of sequence
  
  #if(substring(seq_char, 1, 1) != "M"){
  #  df_error <<- df_error %>%
  #    rbind(c(seq_name, "Error in alignment"))
  #} else{
      df <<- df %>%
        rbind(c(seq_name, seq_sep)) 
  #}
  
}

#Run through each sequence in concatenated file to add aa to dataframe
for(i in 1:length(seqs_fasta)){
  print(i)
  fasta_to_df(seqs_fasta[i], ref_fasta)
}


#empty data frame for antiviral flags
flag <- data.frame()

#function to find any indicated flag for each gene
flag_func <- function(df){
  if(gene == "orf1a"){
    flag_t <- df %>% 
      mutate(sc2_antiviral = case_when(df$`3313` != "L" & df$`3429` != "E" ~ paste0("L50", df$`3313`, ",E166", df$`3429`),
                                       df$`3313` != "L" ~ paste0("L50", df$`3313`),
                                       df$`3429` != "E"~ paste0("E166", df$`3429`),
                                       TRUE ~ "none" ))
  }
  flag_t <- flag_t %>%
    select(name, sc2_antiviral)
  
  flag <<- flag %>%
    rbind(flag_t)
}


#run through each row of df to find antiviral flags
for(i in 1:nrow(df)){
  flag_func(df[i,])
}

flag <- flag %>%
  rename(SampleID = name)



###SAVE FILE##########
#write.csv(flag, "/Users/vhogan/Documents/ARIA/SC2_paxlovid_res/ARIA_paxlovid_res_07Dec23.csv", row.names = FALSE, na = "")
#write.csv(df_error, "/Users/vhogan/Documents/ARIA/SC2_paxlovid_res/ARIA_paxlovid_res_error_07Dec23.csv", row.names = FALSE, na = "")
