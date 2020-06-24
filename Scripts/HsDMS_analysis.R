##### NB: This script will work for all HsT5 DMS data
## Currently set for HIV-1 data, with hashed out N-MLV inputs you can switch to as needed



##### LOAD CSV FILE OF RAW COUNTS (compiled by biostatistics core)
# contains each nucleotide (nt) & translated (aa) seq;
# distance (# changes) from WT nt/aa seq;
# the first position (index and aa) deviating from WT;
# and read counts from all libraries for that nt sequence.
HsDMS = read.csv('../Data/Hs-HIV1_summary.csv')
#HsDMS = read.delim('../Data/Hs-NMLV_summary.tsv')
head(HsDMS,3)

## Make easier to work with; get rid of factor class
class(HsDMS[,"seq_nt"])
HsDMS[,"seq_nt"] = as.character(HsDMS[,"seq_nt"]) #convert from factor to character
class(HsDMS[,"seq_aa"])
HsDMS[,"seq_aa"] = as.character(HsDMS[,"seq_aa"]) #convert from factor to character

## Change to 1-based indexing
HsDMS[,"index_aa"] = HsDMS[,"index_aa"] + 1
## Change diff_aa column name to more precise names
colnames(HsDMS)
colnames(HsDMS)[5] = "mut_index"
colnames(HsDMS)[6] = "mut_aa"

## Add column with correct site# of mutation, rather than index
HsDMS[,"mut_site"] = HsDMS[,"mut_index"] + 329

## Load seqinr package; only external package needed below
library(seqinr)



##### GET RID OF NO-GOOD SEQUENCES BASED ON CHARACTERISTICS, NOT ARBITRARY COUNT CUT-OFF
### ANALYZE SEQUENCES (AND DIFFS FROM WT) BASED ON CODONS, NOT SINGLE NT

## Define function to count # of codon differences b/w given seq and WT
CountCodonDiffs2 = function(Testnt) {
  WTnt = "GGGGCACGAGGGACAAGATACCAGACATTTGTGAATTTC" #HsT5 WT v1, NB: 39nt b/c sequenced with Rh; last 6 are after v1 loop
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  CodonDiffs = sum(TestCodons != WTCodons) #finds number of differences from WT for each seq
  return(CodonDiffs)
}

##Apply fxn to whole table, make new column w/ # of codon differences
HsDMS[,"dist_codon"] = sapply(HsDMS[,"seq_nt"], CountCodonDiffs2) 

## Define function to return mutated codon
ReturnMutCodon = function(Testnt) {
  WTnt = "GGGGCACGAGGGACAAGATACCAGACATTTGTGAATTTC" #HsT5 WT v1, NB: 39nt
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  FirstMutCodon = (TestCodons[FirstMutPosition]) #return this codon's sequence
  return(FirstMutCodon)
}

## Apply fxn to whole table, make new column w/ codon sequence at mutated site
HsDMS[,"mut_codon"] = sapply(HsDMS[,"seq_nt"], ReturnMutCodon)

## Make new column with just the wobble position
HsDMS[,"mut_codonWob"] = sapply(strsplit(HsDMS[,"mut_codon"],''), '[', 3)

## Define function to return mutated position (index) based on codon, not amino acid, comparison
ReturnMutIndex = function(Testnt) {
  WTnt = "GGGGCACGAGGGACAAGATACCAGACATTTGTGAATTTC" #HsT5 WT v1 nt seq; NB: 39nt
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  return(FirstMutPosition)
}

## Apply fxn to whole table, make new column
HsDMS[,"mutcod_index"] = sapply(HsDMS[,"seq_nt"], ReturnMutIndex)
HsDMS[,"mutcod_index"] = HsDMS[,"mutcod_index"] + 329
#NB: now have unique identifier for WT-synonymous variants!


### REMOVE ALL ROWS THAT DO NOT MEET THE FOLLOWING 3 CRITERIA: 

## 1) more than 1 codon difference; implicitly removes >3 nt diffs and >1 aa diff
HsDMSfilter1 = HsDMS[which(HsDMS[,"dist_codon"] <= 1), ] #keep all rows that meet this criterion & all columns
# can check that is so with the following double checks:
max(HsDMS[,"dist_nt"]) 
max(HsDMSfilter1[,"dist_nt"])
max(HsDMS[,"dist_aa"])
max(HsDMSfilter1[,"dist_aa"])
dim(HsDMSfilter1) #there are still too many sequences (>32*11=352)
max(HsDMS[which(HsDMS[,"dist_codon"]>1),"InputA"]) # highest read count I just filtered out = 186 (145 for NMLV)
min(HsDMSfilter1[,"InputA"]) #lowest read count left = 0

## 2) mutated codon doesn't end in C or G (non-NNS)
## NB: | indicates OR statements in R
HsDMSfilter2 = HsDMSfilter1[which(HsDMSfilter1[,"mut_codonWob"] == "C" | HsDMSfilter1[,"mut_codonWob"] == "G" | is.na(HsDMSfilter1[,"mut_codonWob"])), ]
# can check that this worked and returned only C/G (or NA)
unique(HsDMSfilter2[,"mut_codonWob"])
dim(HsDMSfilter2) #there are a few too many sequences still (358): extra positions!

## 3) Remove mutations outside of v1 loop (sites 341, 342)
HsDMSfilter3 = HsDMSfilter2[which(HsDMSfilter2[,"mutcod_index"] <= 340 | is.na(HsDMSfilter2[,"mutcod_index"])), ]

## 4) Must be represented at least 1 count in both input libraries (else calculations problematic)
HsDMSfilter4 = HsDMSfilter3[which(HsDMSfilter3[,"InputA"] > 0 & HsDMSfilter3$InputB > 0), ]
dim(HsDMSfilter4) #now 346 (345 for NMLV) sequences, slightly less than expected 352




##### MAKE NORMALIZED COUNTS AND ENRICHMENT SCORES

## Find sum of columns
total_counts = apply(HsDMSfilter4[,c("InputA", "InputB","SortA", "SortB")], 2, sum) #apply sum fxn to all columns in this list
# convert total counts to millions
total_cpm = total_counts/10^6

## Add new columns normalized to total counts per million
for (tempColname in names(total_cpm) ) {
  newColName <- paste(tempColname,"_norm", sep="")
  HsDMSfilter4[,newColName] <- HsDMSfilter4[,tempColname] / total_cpm[tempColname]
}

## Calculate enrichment for each sort relative to input (normalized values)
HsDMSfilter4[,"enrichA"] = HsDMSfilter4[,"SortA_norm"] / HsDMSfilter4[,"InputA_norm"]
HsDMSfilter4[,"enrichB"] = HsDMSfilter4[,"SortB_norm"] / HsDMSfilter4[,"InputB_norm"]
colnames(HsDMSfilter4)
HsDMSfilter4[,"enrichAvg"] = rowMeans(HsDMSfilter4[20:21])

## Calculate the standard deviation of enrichment scores
HsDMSfilter4[,"enrichStdDev"] = apply(HsDMSfilter4[,20:21], 1, sd)


### EXPORT to file so that I can re-analyze / have non-averaged values
write.csv(HsDMSfilter4, file="../Analysis/Hs-HIV_filter3_noAvg.csv",row.names = FALSE)
#write.csv(HsDMSfilter4, file="../Analysis/Hs-NMLV_filter3_noAvg.csv",row.names = FALSE)

### Save point: read in my file again here if I'd like to change the cutoffs below
#HsDMSfilter4 = read.delim("../Analysis/Hs-HIV_filter3_noAvg.csv", sep=",")
#HsDMSfilter4 = read.delim("../Analysis/Hs-NMLV_filter3_noAvg.csv", sep=",")




##### NEED TO FILTER OUT LOW-REPRESENTATION SEQUENCES, WHICH WILL BE NOISY

### A) Decide representation cutoffs
## Plot inputA vs. input B counts ... not that helpful, they're highly correlated
plot(HsDMSfilter4$InputA_norm, HsDMSfilter4$InputB_norm, log="xy")

## Plot StdDev against input counts to choose cutoff
plot(HsDMSfilter4$InputA_norm, HsDMSfilter4$enrichStdDev, log="x")
  # Reasonable cutoff: there's a gap at ~50, where input >50 has low SD across replicates
  # for N-MLV, gap more like ~25 cpm

## Zoom in on the low count range to find cutoff for high std dev
plot(HsDMSfilter4$InputA_norm, HsDMSfilter4$enrichStdDev, xlim= c(1,100))
  # lowest cutoff should be 80? Most high SD are at < 40 input counts, one outlier around 70
  # for N-MLV, cut off of 20cpm would be fine
plot(HsDMSfilter4$InputB_norm, HsDMSfilter4$enrichStdDev, log="x")
  # the 50 cutoff will get rid of the 4s.d. point.
  # for N-MLV, there's one high SD value around 80, but inputA (~20) filter will get rid of it
#50cpm seems like a reasonable cutoff across all these experiments

## Plot enrichment scores against input values
plot(HsDMSfilter4$InputA_norm, HsDMSfilter4$enrichAvg, log="x")
plot(HsDMSfilter4$InputB_norm, HsDMSfilter4$enrichAvg, log="x")
  #no clear correlation here; use the SD-based method above
  #for NMLV, that 1 high SD point is clearly also highly enriched = not reliable data


### B) Filter out sequences based on these cutoffs
HsDMSfilter5 = HsDMSfilter4[which(HsDMSfilter4$InputA_norm > 50 & HsDMSfilter4$InputB_norm > 50), ]
  #now down to 333 seqs from 346




##### AVERAGE ACROSS ALL NT SEQS FOR EACH AA MUTATION

dim(HsDMSfilter5) # check how many unique sequence we have; I expect 32 x 11 = 352; I have 333
table(table(HsDMSfilter5[,"seq_aa"])) # determine how many aa seqs have 1/2/3/etc. unique nt seqs
# NOTE HERE THAT ALL WT SEQS ARE GROUPED INSTEAD OF KEPT SEPARATE, hence 17 WT nt seqs
length(unique(HsDMSfilter5[,"seq_aa"])) # determine # unique aa seqs I have; I expect 20x11 +1 = 221; I have 214 (missing 7 non-WT) (215 for N-MLV)

## Concatenate mutcodon index with aa seq to differentiate WT seqs
HsDMSfilter5[,"mutcodon+seq"] = paste(HsDMSfilter5$mutcod_index, HsDMSfilter5$seq_aa, sep="_")
length(unique(HsDMSfilter5$`mutcodon+seq`)) # determine # unique aa seqs I have; I expect 21 x 11 = 231; I have 223 (224 for N-MLV)
  #also missing 1 WT, b/c there are 2 positions that have no synonymous NNS codons outside of WT

## Split whole dataframe into list of dataframes for eqch unique aa seq
HsDMSfilter5Split <- split(HsDMSfilter5, HsDMSfilter5[,"mutcodon+seq"]) #2nd argument = what to split based upon
head(HsDMSfilter5Split,3)

## Define function to average values; also keeps only useful columns
get_averages = function(datf){  #named function that will act on dataframe
  datf[,"num_nt_seqs"] = dim(datf)[1] #identify # unique nt seqs (rows) for each aa seq (subsetted dataframe)
  datf[,"enrichA_avg"] = mean(datf[,"enrichA"]) #add column w/ mean (across synonymous seqs only)
  datf[,"enrichB_avg"] = mean(datf[,"enrichB"])
  datf[,"enrich_Avg_RepAndAA"] = mean(datf[,"enrichAvg"]) #add column with the mean (across reps & synonymous codons)
  datf[,"enrich_SD_RepAndAA"] = sd(c(datf[,"enrichA"],datf[,"enrichB"]))
  columnsToKeep <- c("dist_codon","seq_aa","mutcod_index","mut_aa", "num_nt_seqs","enrichA_avg", "enrichB_avg", "enrich_Avg_RepAndAA","enrich_SD_RepAndAA")
  datf <- datf[1, columnsToKeep] # overwrite dataframe (w/in fxn only) with the first row & only these columns
  return(datf)
}

## Perform this averaging function on all unique aa seqs
HsDMSfilter5SplitAvgs = lapply (HsDMSfilter5Split, get_averages)
head(HsDMSfilter5SplitAvgs,3)

## Bring all these values back into single table
HsDMSfilter5Avgs = do.call("rbind", HsDMSfilter5SplitAvgs) #row bind each of these rows
rownames(HsDMSfilter5Avgs) = NULL #gets rid of useless row names
head(HsDMSfilter5Avgs,3)

## SAVE OUTPUT (enrichment scores averaged across synonymous codons) AS CSV FILE
write.csv(HsDMSfilter5Avgs, file="../Analysis/Hs-HIV1_avg.csv",row.names = FALSE)
#write.csv(HsDMSfilter5Avgs, file="../Analysis/Hs-NMLV_avg.csv",row.names = FALSE)




##### GENERATE A POSITION X AMINO ACID MATRIX
## make sure nothing is a factor!
class(HsDMSfilter5Avgs$mut_aa) #factor no good
HsDMSfilter5Avgs$mut_aa = as.character(HsDMSfilter5Avgs$mut_aa)
class(HsDMSfilter5Avgs$mutcod_index) #integer OK

##null values are a problem for populating the table
##1) replace null values in mut_aa column (synonymous w/ WT) with "syn"
HsDMSfilter5Avgs$mut_aa[HsDMSfilter5Avgs$mut_aa == ""] = "syn"

##2) replace single null value in index position (WT nt seq)
## I know that HsY336 (and Q337) has only one codon & is thus empty -> put it there
HsDMSfilter5Avgs$mutcod_index[is.na(HsDMSfilter5Avgs$mutcod_index)] = 336

##table populates based on 1/2/3/etc. indexing, won't start at 330 (has 329 empty rows)
##so change it back to indexing starting at 1
HsDMSfilter5Avgs$mutcod_index = HsDMSfilter5Avgs$mutcod_index - 329

## Make empty dataframe
HsDMSmatrix = data.frame()

## Populate dataframe
for (rowindex in 1:dim(HsDMSfilter5Avgs)[1]) { #iterate over the # of rows (first value of dimensions of table)
  RowNum = HsDMSfilter5Avgs[rowindex,"mutcod_index"] # define the row # as the index position
  ColNum = HsDMSfilter5Avgs[rowindex, "mut_aa"] # define the col # as the amino acid
  HsDMSmatrix[RowNum,ColNum] = HsDMSfilter5Avgs[rowindex, "enrich_Avg_RepAndAA"] #populate dataframe at row/col with the enrichment score from that row
}

## Add column with proper position #
HsDMSmatrix[,"Position"] = c("G330","A331","R332","G333","T334","R335","Y336","Q337","T338","F339","V340")

## Reorder rows into my favorite order
HsDMSmatrixOrdered = HsDMSmatrix[,c("Position","G","P","A","I","L","V","F","W","Y","M","C","S","T","N","Q","H","K","R","D","E","*","syn")]

## Save this matrix as csv file
write.csv(HsDMSmatrixOrdered, file="../Analysis/Hs-HIV1_matrix.csv",row.names = FALSE)
#write.csv(HsDMSmatrixOrdered, file="../Analysis/Hs-NMLV_matrix.csv",row.names = FALSE)

