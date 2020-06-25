##### LOAD CSV FILE OF RAW COUNTS (compiled by biostatistics core)
# contains each nucleotide (nt) & translated (aa) seq;
# distance (# changes) from WT nt/aa seq;
# the first position (index and aa) deviating from WT;
# and read counts from all libraries for that nt sequence.
RPDMS = read.csv('../Data/R332P-HIV1_summary.csv')
head(RPDMS,3)

## Make easier to work with; get rid of factor class
RPDMS[,"seq_nt"] = as.character(RPDMS[,"seq_nt"]) #convert from factor to character
RPDMS[,"seq_aa"] = as.character(RPDMS[,"seq_aa"]) #convert from factor to character

## Change to 1-based indexing
RPDMS[,"index_aa"] = RPDMS[,"index_aa"] + 1
## Change diff_aa column name to more precise name
colnames(RPDMS)
colnames(RPDMS)[5] = "mut_index"
colnames(RPDMS)[6] = "mut_aa"

## Load seqinr package; only external package needed below
library(seqinr)




##### GET RID OF NO-GOOD SEQUENCES BASED ON CHARACTERISTICS, NOT COUNT CUT-OFF
### ANALYZE SEQUENCES (AND DIFFS FROM WT) BASED ON CODONS, NOT SINGLE NT

## Define function to count # of codon differences b/w given seq and WT
CountCodonDiffs2 = function(Testnt) {
  WTnt = "GGGGCACCAGGGACAAGATACCAGACATTTGTG" #This is v1 of HsT5-R332P (wt)
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: in R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  CodonDiffs = sum(TestCodons != WTCodons) #finds number of differences from WT for one seq
  return(CodonDiffs)
}

##Apply fxn to whole table, make new column w/ # of codon differences
RPDMS[,"dist_codon"] = sapply(RPDMS[,"seq_nt"], CountCodonDiffs2) 

## Define function to return mutated codon
ReturnMutCodon = function(Testnt) {
  WTnt = "GGGGCACCAGGGACAAGATACCAGACATTTGTG" #This is v1 of HsT5-R332P (wt)
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: in R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  FirstMutCodon = (TestCodons[FirstMutPosition]) #return this codon's sequence
  return(FirstMutCodon)
}

## Apply fxn to whole table, make new column w/ codon sequence at mutated site
RPDMS[,"mut_codon"] = sapply(RPDMS[,"seq_nt"], ReturnMutCodon)

## Make new column with just the wobble position
RPDMS[,"mut_codonWob"] = sapply(strsplit(RPDMS[,"mut_codon"],''), '[', 3)

## Define function to return mutated position (index) based on codon, not amino acid, comparison
ReturnMutIndex = function(Testnt) {
  WTnt = "GGGGCACCAGGGACAAGATACCAGACATTTGTG" #This is v1 of HsT5-R332P (wt)
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: in R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  return(FirstMutPosition)
}

## Apply fxn to whole table, make new column
RPDMS[,"mutcod_index"] = sapply(RPDMS[,"seq_nt"], ReturnMutIndex)
#NB: now have unique identifier for WT-synonymous variants!


### REMOVE ALL ROWS THAT DO NOT MEET THE FOLLOWING 3 CRITERIA: 

## 1) more than 1 codon difference; implicitly removes >3 nt diffs and >1 aa diff
RPDMSfilter1 = RPDMS[which(RPDMS[,"dist_codon"] <= 1), ] #keep all rows that meet this criterion & all columns
## double check this worked:
max(RPDMS[,"dist_nt"])
max(RPDMSfilter1[,"dist_nt"])
max(RPDMS[,"dist_aa"])
max(RPDMSfilter1[,"dist_aa"])
dim(RPDMSfilter1) #there are still too many sequences (669 > 32*10=320)
max(RPDMS[which(RPDMS[,"dist_codon"]>1),"InputA"]) # highest read count I just filtered out = 331
min(RPDMSfilter1[,"InputA"]) #lowest read count left = 0

## 2) mutated codon doesn't end in C or G (non-NNS)
## NB: | indicates OR statements in R
RPDMSfilter2 = RPDMSfilter1[which(RPDMSfilter1[,"mut_codonWob"] == "C" | RPDMSfilter1[,"mut_codonWob"] == "G" | is.na(RPDMSfilter1[,"mut_codonWob"])), ]
# can check that this worked and returned only C/G (or NA)
unique(RPDMSfilter2[,"mut_codonWob"])
max(RPDMSfilter1[which(RPDMSfilter1[,"mut_codonWob"] == "A" | RPDMSfilter1[,"mut_codonWob"] == "T"),"InputA"]) # highest read count I just filtered out = 6962 (yikes!)
min(RPDMSfilter2[,"InputA"]) #lowest read count left = still 0
dim(RPDMSfilter2) #there are a few too many sequences still (336): extra positions!

## 3) Remove mutations at P332 (index 3)
unique(RPDMSfilter2[,"mutcod_index"])
RPDMSfilter3 = RPDMSfilter2[which(RPDMSfilter2[,"mutcod_index"] != 3 | is.na(RPDMSfilter2[,"mutcod_index"])), ]
max(RPDMSfilter2[which(RPDMSfilter2[,"mutcod_index"] ==3), "InputA"]) # highest read count I just filtered out = 62
min(RPDMSfilter3[,"InputA"]) #lowest read count left = 0

## 4) Must be represented at least 1 count in all input libraries (else calculations problematic)
RPDMSfilter4 = RPDMSfilter3[which(RPDMSfilter3$InputA > 0 & 
                                    RPDMSfilter3$InputA > 0), ]




##### MAKE NORMALIZED COUNTS AND ENRICHMENT SCORES
colnames(RPDMSfilter4)

## Find sum of columns
total_counts = apply(RPDMSfilter4[,c("InputA", "InputB", "SortA", "SortB")], 2, sum) #apply sum fxn to all columns in this list
3# convert total counts to millions
total_cpm = total_counts/10^6

## Add new columns normalized to total counts per million
for (tempColname in names(total_cpm) ) {
  newColName <- paste(tempColname,"_norm", sep="")
  RPDMSfilter4[,newColName] <- RPDMSfilter4[,tempColname] / total_cpm[tempColname]
}

## Compare input libraries vs. e/o to determine noisiness b/w replicates
plot(RPDMSfilter4$InputA_norm, RPDMSfilter4$InputB_norm)
plot(RPDMSfilter4$InputA_norm, RPDMSfilter4$InputB_norm, xlim = c(1,20000), ylim= c(1,20000))
reg_in1 = lm(InputB_norm~InputA_norm, data = RPDMSfilter4)
summary(reg_in1) #R2 = 0.991
abline(reg_in1)
# not very noisy b/w replicate inputs

## Calculate enrichment for each sort relative to input (normalized values)
RPDMSfilter4[,"enrichA"] = RPDMSfilter4$SortA_norm / RPDMSfilter4$InputA_norm
RPDMSfilter4[,"enrichB"] = RPDMSfilter4$SortB_norm / RPDMSfilter4$InputB_norm
colnames(RPDMSfilter4)
RPDMSfilter4[,"enrichavg"] = rowMeans(RPDMSfilter4[19:20])

## Compare enrichment across replicates (at nt level)
plot(RPDMSfilter4$enrichA, RPDMSfilter4$enrichB)
#Sadly, this is more noisy b/w reps, even at high enrichment values :(

## Calculate the standard deviation of enrichment scores
RPDMSfilter4[,"enrichSD"] = apply(RPDMSfilter4[,19:20], 1, sd)


### EXPORT to file so that I can re-analyze / have non-averaged values
write.csv(RPDMSfilter4, file="../Analysis/R332P-HIV_filter3_noAvg.csv",row.names = FALSE)
### Save point: read in my file again here if I'd like to change the cutoffs below
#RPDMSfilter4 = read.delim("../Analysis/R332P-HIV_filter3_noAvg.csv", sep=",")




##### NEED TO FILTER OUT LOW-REPRESENTATION SEQUENCES, WHICH WILL BE NOISY

### A) Decide representation cutoffs

## Plot StdDev against input counts to choose cutoff
par(mfrow=c(2,1))
plot(RPDMSfilter4$InputA_norm, RPDMSfilter4$enrichSD, log = "x")
plot(RPDMSfilter4$InputB_norm, RPDMSfilter4$enrichSD, log = "x")
  # Reasonable cutoff: there's a gap at ~150 for normalized average counts? 
## zoom in on the low count range to find cutoff for high std dev
plot(RPDMSfilter4$InputA_norm, RPDMSfilter4$enrichSD, xlim= c(1,300)) #things level out >200ish
plot(RPDMSfilter4$InputB_norm, RPDMSfilter4$enrichSD, xlim= c(1,300)) #things level out >200ish
  # Lowest cutoff should be 150-200 in both inputs

## Plot enrichment scores against input values
plot(RPDMSfilter4$InputA_norm, RPDMSfilter4$enrichavg, log="x")
plot(RPDMSfilter4$InputB_norm, RPDMSfilter4$enrichavg, log="x")
  # no clear pattern, not very informative


### B) Filter out sequences into separate experiments based on these cutoffs
RPDMSfilter5 = RPDMSfilter4[which(RPDMSfilter4$InputA_norm > 200 & RPDMSfilter4$InputB_norm > 200), ]

#Compare replicates (at nt level, post-filter)
plot(RPDMSfilter5$enrichA, RPDMSfilter5$enrichB) # Much better but only 297 seqs




##### AVERAGE ACROSS ALL NT SEQS FOR EACH AA MUTATION

dim(RPDMSfilter5) # check how many unique sequence we have; I expect 32 x 10 = 320; I have 297
table(table(RPDMSfilter5[,"seq_aa"])) # determine how many aa seqs have 1/2/3/etc. unique nt seqs
# NOTE HERE THAT ALL WT SEQS ARE GROUPED INSTEAD OF KEPT SEPARATE, hence 14 wt nt seqs
length(unique(RPDMSfilter5[,"seq_aa"])) # determine # unique aa seqs I have; I expect 20x10 +1 = 201; I have 192 (missing 9 non-WT)

## Concatenate mutcodon index with aa seq to differentiate WT seqs
RPDMSfilter5[,"mutcodon+seq"] = paste(RPDMSfilter5$mutcod_index, RPDMSfilter5$seq_aa, sep="_")
length(unique(RPDMSfilter5$`mutcodon+seq`)) # determine # unique aa seqs I have; I expect 21 x 10 = 210; I have 200 (missing 9 non-WT, as expected)

## Split whole dataframe into list of dataframes for eqch unique aa seq
RPDMSfilter5Split <- split(RPDMSfilter5, RPDMSfilter5[,"mutcodon+seq"]) #2nd argument = what to split based upon
head(RPDMSfilter5Split,3)

## Define function to average values; also keeps only useful columns
get_averages = function(datf){  #named function that will act on dataframe
  datf[,"num_nt_seqs"] = dim(datf)[1] #identify # unique nt seqs (rows) for each aa seq (subsetted dataframe)
  datf[,"enrich_Avg_RepAndAA"] = mean(datf[,"enrichavg"]) # add column with the mean (across reps & synonymous codons)
  datf[,"enrichA_Avg"] = mean(datf[,"enrichA"]) #add column w/ mean (across synonymous seqs only)
  datf[,"enrichB_Avg"] = mean(datf[,"enrichB"])
  datf[,"enrich_SD_RepAndAA"] = sd(c(datf[,"enrichA"],datf[,"enrichB"]))
  columnsToKeep <- c("dist_codon","seq_aa","mutcod_index","mut_aa", "num_nt_seqs","enrichA_Avg","enrichB_Avg","enrich_Avg_RepAndAA","enrich_SD_RepAndAA")
  datf <- datf[1, columnsToKeep] # overwrite dataframe (w/in fxn only) with the first row & only these columns
  return(datf)
}

## Perform this averaging function on all unique aa seqs
RPDMSfilter5SplitAvgs = lapply (RPDMSfilter5Split, get_averages)
head(RPDMSfilter5SplitAvgs,3)

## Bring all these values back into single table
RPDMSfilter5Avgs = do.call("rbind", RPDMSfilter5SplitAvgs) #row bind each of these rows
rownames(RPDMSfilter5Avgs) = NULL #gets rid of useless row names
head(RPDMSfilter5Avgs,3)

## SAVE OUTPUT (enrichment scores averaged across synonymous codons) AS CSV FILE
write.csv(RPDMSfilter5Avgs, file="../Analysis/R332P-HIV1_avg.csv",row.names = FALSE)




##### GENERATE A POSITION X AMINO ACID MATRIX
#NB: make sure nothing is a factor!
class(RPDMSfilter5Avgs$mut_aa) #factor no good
RPDMSfilter5Avgs$mut_aa = as.character(RPDMSfilter5Avgs$mut_aa)
class(RPDMSfilter5Avgs$mutcod_index) #integer OK

##null values are a problem for populating the table
##1) replace null values in mut_aa column (synonymous w/ WT) with "syn"
RPDMSfilter5Avgs$mut_aa[RPDMSfilter5Avgs$mut_aa == ""] = "syn"
table(RPDMSfilter5Avgs$mut_aa) #9 synonymous codons (as expected)
RPDMSfilter5Avgs[is.na(RPDMSfilter5Avgs$mutcod_index),] #1 WT nt seq

##2) replace single null value in index position (WT nt seq)
## I know that Y336 (and Q337) has only one codon & is thus empty -> put it there (position 7)
unique(RPDMSfilter5Avgs$mutcod_index) #no more position 3!
RPDMSfilter5Avgs$mutcod_index[is.na(RPDMSfilter5Avgs$mutcod_index)] = 7

## Make empty dataframe
RPDMSmatrix = data.frame()

## Populate dataframe
for (rowindex in 1:dim(RPDMSfilter5Avgs)[1]) { #iterate over the # of rows (first value of dimensions of table)
  RowNum = RPDMSfilter5Avgs[rowindex,"mutcod_index"] # define the row # as the index position
  ColNum = RPDMSfilter5Avgs[rowindex, "mut_aa"] # define the col # as the amino acid
  RPDMSmatrix[RowNum,ColNum] = RPDMSfilter5Avgs[rowindex, "enrich_Avg_RepAndAA"] #populate dataframe at row/col with the Sort value from that row
}

## Add column with proper position #
## NB: keep R332 in b/c this is a categorical/factor-based filling
RPDMSmatrix[,"Position"] = c("G330","A331","R332","G333","T334","R335","Y336","Q337","T338","F339","V340")

## Reorder rows into my favorite order
RPfilter5MatrixOrdered = RPDMSmatrix[,c("Position","G","P","A","I","L","V","F","W","Y","M","C","S","T","N","Q","H","K","R","D","E","*","syn")]

## Save this matrix as csv file
write.csv(RPfilter5MatrixOrdered, file="../Analysis/R332P-HIV1_matrix.csv",row.names = FALSE)
