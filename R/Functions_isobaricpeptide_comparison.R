#functions to find unique fragments isobaric peptides
# library(tidyverse)

# based on Spectratable

# fragmentnames like in MQ
# a, b, y, b/y*, b/y(+2), b/y-NH3, b/y-H2O
# * ~ H3PO4_loss
# other neutral loss (HPO3_loss) not in MQ

# watch out: function not optimized for isobaric peptides carrying other modifications then "(Phospho (STY))"

# for a sequence, other isobaric peptides in spectratable
# based on mz (should be identical!!!)

read.delim("/Users/henno/Documents/Skripte/R-Functions/Selbach-Functions/DataAnalysis/20220209_MVB_Rawfiles_library_bio_benchmarking/TransitionTable.txt", 
           header = T, "\t", stringsAsFactors = FALSE, dec = " ", check.names = F) -> SpectraTabletxt

# pairwise comparison function, spits out unique fragments. 
fun_pairwisecomparison <- function(Spectratable_seq1, Spectratable_seq2,SearchMod){
# syntax Spectratable_seq1 _IADPEHDHTGFLTEY(Phospho (STY))VATR_
#read in characteristics sequences
  
if(grepl("Acetyl|Oxidation",Spectratable_seq1)|grepl("Acetyl|Oxidation",Spectratable_seq2)){
  return(NULL)
}

Spectratable_seq1 <- gsub("\\(Phospho \\(STY\\)\\)", 
                          "p", Spectratable_seq1)
seq1 <- gsub("_", "", Spectratable_seq1)

loc_1 <- regexpr("[STY]p", seq1)[[1]][1]
aa1 <- substr(seq1, loc_1, 
              loc_1 +1)

Spectratable_seq2 <- gsub("\\(Phospho \\(STY\\)\\)", 
                          "p", Spectratable_seq2)
seq2 <- gsub("_", "", Spectratable_seq2)

loc_2 <- regexpr("[STY]p", seq2)[[1]][1]
aa2 <- substr(seq1, loc_2, 
              loc_2 +1)

loc_high <- max(loc_1, loc_2)
loc_low <- min(loc_1, loc_2)

seq_len <- nchar(gsub("p", "", seq2))

#find unique y ions: 
ymin <- seq_len - loc_high + 1 
ymax <- seq_len - loc_low
y_unique <- paste0("y", ymin:ymax)

#find unique b ions
bmin <- loc_low
bmax <- loc_high -1
b_unique <- paste0("b", bmin:bmax)

fragmentvector <- c(y_unique, b_unique)

# add other fragment variants
# fragment naming from MQ
fragmentvector <- unique(c(fragmentvector,
                    paste0(fragmentvector, "*"),
                    paste0(fragmentvector, "-NH3"),
                    paste0(fragmentvector, "-H2O"),
                    paste0(fragmentvector, "(2+)"),
                    gsub("b", "a", fragmentvector)))

# pY ion unique identifier? 
if((aa1 == "Yp"|aa2 == "Yp") & aa1 != aa2){
  fragmentvector <- c(fragmentvector, "pY")
  return(fragmentvector)
}

return(fragmentvector) }


# decides if and how many pairwise comparisons need to be ran
# in case of 2 pairwise comparisons, selects intersect as unique fragments. 
function_compareseqs <- function(seqs){
if (length(seqs) > 1){
 unique_frag <- fun_pairwisecomparison(seqs[1],seqs[2])
  if (length(seqs) > 2){
    unique_frag_comp2 <- fun_pairwisecomparison(seqs[1],seqs[3])
  unique_frag <- intersect(unique_frag, unique_frag_comp2)
   print(unique_frag)
  return(unique_frag)
  }
 print(unique_frag)
 return(unique_frag)

}else{print(paste0("No isobaric peptide detected for ", seqs))}
}

stop()

#selects all isobaric peptides for a rownumber in the spectratable. 
#based on identical mz values
#the sequence of interest is first in the 'sequences' vector (important when >2 isobaric peptides)

for (rownumber in 1:nrow(SpectraTabletxt)){
  print(rownumber)
mz_value <- SpectraTabletxt$mz[rownumber]
Modseq <- SpectraTabletxt$Modified.sequence[rownumber]

SpectraTabletxt %>% 
  filter(mz ==	mz_value) %>% # change to selection based on unmodified sequence. Doesn't work now because the LIGHT/HEAVY label isn't working yet. problem with m(ox) peptides as well.  
  pull(Modified.sequence) -> sequences
sequences <- unique(c(Modseq,sequences))

function_compareseqs(sequences)
}



##### testing: ####
mz_value <- "754.67493"  
  SpectraTabletxt %>% View
  filter(mz ==	mz_value) %>% # change to selection based on unmodified sequence. Doesn't work now because the LIGHT/HEAVY label isn't working yet. problem with m(ox) peptides as well.  
  pull(Modified.sequence) -> seqs
length(seqs)
unique_frag_comp1 <- fun_pairwisecomparison("_(Oxidation)IADPEHDHTGFLTEY(Phospho (STY))VATR_","_IADPEHDHTGFLTEYVAT(Phospho (STY))R_")
#pairwise comparison with a third isobaric peptide
unique_frag_comp2 <- fun_pairwisecomparison("_IADPEHDHTGFLTEY(Phospho (STY))VATR_","_IADPEHDHTGFLT(Phospho (STY))EYVATR_")
intersect(unique_frag_comp1,unique_frag_comp2)
unique_frag_comp3 <- fun_pairwisecomparison("_IADPEHDHTGFLT(Phospho (STY))EYVATR_","_IADPEHDHTGFLTEYVAT(Phospho (STY))R_")
#pairwise comparison with a third isobaric peptide
intersect(unique_frag_comp2,unique_frag_comp3)
intersect(unique_frag_comp1,unique_frag_comp2)
# identify primary and secondary identifying fragments. 
function_compareseqs(sequences)


