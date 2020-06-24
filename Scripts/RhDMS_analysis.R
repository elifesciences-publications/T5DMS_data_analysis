##### LOAD CSV FILE OF RAW COUNTS (compiled by biostatistics core)
# contains each nucleotide (nt) & translated (aa) seq;
# distance (# changes) from WT nt/aa seq;
# the first position (index and aa) deviating from WT;
# and read counts from all libraries for that nt sequence.
RhDMS = read.csv('../Data/RhDMS-HIV1-NMLV.csv')
head(RhDMS,3)

## Make easier to work with; get rid of factor class
class(RhDMS[,"seq_nt"])
RhDMS[,"seq_nt"] = as.character(RhDMS[,"seq_nt"]) #convert from factor to character
RhDMS[,"seq_aa"] = as.character(RhDMS[,"seq_aa"]) #convert from factor to character

## Change to 1-based indexing
RhDMS[,"index_aa"] = RhDMS[,"index_aa"] + 1
## Change diff_aa column name to more precise names
colnames(RhDMS)
colnames(RhDMS)[5] = "mut_index"
colnames(RhDMS)[6] = "mut_aa"

## Load seqinr package; only external package needed below
library(seqinr)




##### GET RID OF NO-GOOD SEQUENCES BASED ON CHARACTERISTICS, NOT ARBITRARY COUNT CUT-OFF

### (0) Because this is a very large file, get rid of some obviously bad sequences up-front (save processing time).
## Remove >3 nt diffs and >1 aa diff
RhDMSfilter0 = RhDMS[which(RhDMS[,"dist_nt"] <= 3 & RhDMS[,"dist_aa"] <=1), ]

## Double check this worked:
max(RhDMS[,"dist_nt"])
max(RhDMSfilter0[,"dist_nt"]) #should be 3
max(RhDMS[,"dist_aa"])
max(RhDMSfilter0[,"dist_aa"]) #should be 1

## Get rough estimate of reads for sequences I just filtered = 2413 counts!?!?!
max(RhDMS[which(RhDMS[,"dist_nt"] > 3 & RhDMS[,"dist_aa"] > 1),"uninfected"])
head(sort(RhDMS[which(RhDMS[,"dist_nt"] > 3 & RhDMS[,"dist_aa"] > 1),"uninfected"], decreasing = TRUE))
#after that first value, it goes to 670, 544, 422, 181, 176 = mostly very low counts
#So what's that first one?
RhDMS[which(RhDMS[,"uninfected"]==2413),]
#Some random variant (qapVtlftRpslt). Likely it really was part of the library 
#Exclude anyway as not single mutant = can't assign effect to one amino acid change


### ANALYZE SEQUENCES (AND DIFFS FROM WT) BASED ON CODONS, NOT SINGLE NT

## Define function to count # of codon differences b/w given seq and WT
CountCodonDiffs2 = function(Testnt) {
  WTnt = "CAGGCACCAGGGACATTATTTACGTTTCCGTCACTCACG" #This is v1 of WT Rhesus TRIM5a.
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  CodonDiffs = sum(TestCodons != WTCodons) #finds number of codon differences from WT for each seq
  return(CodonDiffs)
}

##Apply fxn to whole table, make new column w/ # of codon differences
RhDMSfilter0[,"dist_codon"] = sapply(RhDMSfilter0[,"seq_nt"], CountCodonDiffs2) 

### Define function to return mutated codon
ReturnMutCodon = function(Testnt) {
  WTnt = "CAGGCACCAGGGACATTATTTACGTTTCCGTCACTCACG" #This is v1 of WT Rhesus TRIM5a.
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  FirstMutCodon = (TestCodons[FirstMutPosition]) #return this codon's sequence
  return(FirstMutCodon)
}

## Apply fxn to whole table, make new column w/ codon sequence at mutated site
RhDMSfilter0[,"mut_codon"] = sapply(RhDMSfilter0[,"seq_nt"], ReturnMutCodon)

## Make new column with just the wobble position
RhDMSfilter0[,"mut_codonWob"] = sapply(strsplit(RhDMSfilter0[,"mut_codon"],''), '[', 3)

## Define function to return mutated position (index) based on codon, not amino acid, comparison
ReturnMutIndex = function(Testnt) {
  WTnt = "CAGGCACCAGGGACATTATTTACGTTTCCGTCACTCACG" #This is v1 of WT Rhesus TRIM5a.
  WTntSplit = strsplit(WTnt,"")[[1]] #split into vector of characters. NB: In R, need to pull first element from the list of 1 to get a character, not a list
  TestntSplit = strsplit(Testnt,"")[[1]]
  WTCodons = splitseq(WTntSplit, frame = 0, word = 3) #requires seqinr package; splits into codons
  TestCodons = splitseq(TestntSplit, frame = 0, word = 3) 
  FirstMutPosition = which(TestCodons != WTCodons)[1] #find 1st codon where test seq differs from WT
  return(FirstMutPosition)
}

## Apply fxn to whole table, make new column
RhDMSfilter0[,"mutcod_index"] = sapply(RhDMSfilter0[,"seq_nt"], ReturnMutIndex)
#NB: now have unique identifier for WT-synonymous variants!


### REMOVE ALL ROWS THAT DO NOT MEET THE FOLLOWING 3 CRITERIA: 

## 1) more than 1 codon difference; implicitly removes >3 nt diffs and >1 aa diff
RhDMSfilter1 = RhDMSfilter0[which(RhDMSfilter0[,"dist_codon"] <= 1), ] #keep all rows that meet this criterion & all columns
dim(RhDMSfilter1) #there are still too many sequences (817 > 32*13=416)
max(RhDMSfilter0[which(RhDMSfilter0[,"dist_codon"]>1),"uninfected"]) # highest read count I just filtered out = 327
min(RhDMSfilter1[,"uninfected"]) #lowest read count left = 0

## 2) mutated codon doesn't end in C or G (non-NNS)
RhDMSfilter2 = RhDMSfilter1[which(RhDMSfilter1[,"mut_codonWob"] == "C" | RhDMSfilter1[,"mut_codonWob"] == "G" | is.na(RhDMSfilter1[,"mut_codonWob"])), ]
# can check that this worked and returned only C/G (or NA = WT)
unique(RhDMSfilter2[,"mut_codonWob"])
max(RhDMSfilter1[which(RhDMSfilter1[,"mut_codonWob"] == "A" | RhDMSfilter1[,"mut_codonWob"] == "T"),"uninfected"]) # highest read count I just filtered out = 316
min(RhDMSfilter2[,"uninfected"]) #lowest read count left = 8
dim(RhDMSfilter2) #now have 411 < 416 unique sequences

## 3) Must be represented at least 1 count in all input libraries (else calculations problematic)
RhDMSfilter3 = RhDMSfilter2[which(RhDMSfilter2$uninfected > 0 & 
                                    RhDMSfilter2$HIV1_inA > 0 & 
                                    RhDMSfilter2$HIV1_inB > 0 & 
                                    RhDMSfilter2$NMLV_inA > 0 & 
                                    RhDMSfilter2$NMLV_inB > 0), ]




##### MAKE NORMALIZED COUNTS AND ENRICHMENT SCORES
colnames(RhDMSfilter3)

## Find sum of columns
total_counts = apply(RhDMSfilter3[,c("uninfected", "HIV1_inA", "HIV1_inB", "HIV1_sortA", "HIV1_sortB",
                                     "NMLV_inA", "NMLV_inB", "NMLV_sortA", "NMLV_sortB")], 2, sum) #apply sum fxn to all columns in this list
## convert total counts to millions
total_cpm = total_counts/10^6

## Add new columns for read counts normalized to total counts per million
for (tempColname in names(total_cpm) ) {
  newColName <- paste(tempColname,"_norm", sep="")
  RhDMSfilter3[,newColName] <- RhDMSfilter3[,tempColname] / total_cpm[tempColname]
}

## Calculate enrichment for each sort relative to input (normalized values)
RhDMSfilter3[,"enrich1A"] = RhDMSfilter3$HIV1_sortA_norm / RhDMSfilter3$HIV1_inA_norm
RhDMSfilter3[,"enrich1B"] = RhDMSfilter3$HIV1_sortB_norm / RhDMSfilter3$HIV1_inB_norm
RhDMSfilter3[,"enrich2A"] = RhDMSfilter3$NMLV_sortA_norm / RhDMSfilter3$NMLV_inA_norm
RhDMSfilter3[,"enrich2B"] = RhDMSfilter3$NMLV_sortB_norm / RhDMSfilter3$NMLV_inB_norm

## Calculate enrichment averaged across replicates for each nt seq
colnames(RhDMSfilter3)
RhDMSfilter3[,"enrich1avg"] = rowMeans(RhDMSfilter3[29:30])
RhDMSfilter3[,"enrich2avg"] = rowMeans(RhDMSfilter3[31:32])

## Compare replicates (at nt level)
par(mfrow=c(2,1)) #view 2 graphs in same panel
plot(RhDMSfilter3$enrich1A, RhDMSfilter3$enrich1B)
plot(RhDMSfilter3$enrich2A, RhDMSfilter3$enrich2B)
  #high correlation b/w reps at nt level

## Calculate the standard deviation of enrichment scores at nt level
colnames(RhDMSfilter3)
RhDMSfilter3[,"enrich1SD"] = apply(RhDMSfilter3[,29:30], 1, sd)
RhDMSfilter3[,"enrich2SD"] = apply(RhDMSfilter3[,31:32], 1, sd)

## Compare input libraries vs. e/o & vs. uninfected
## to determine noisiness b/w replicates and if infection made a difference
plot(RhDMSfilter3$HIV1_inA_norm, RhDMSfilter3$HIV1_inB_norm)
plot(RhDMSfilter3$HIV1_inA_norm, RhDMSfilter3$HIV1_inB_norm, xlim = c(1,5000), ylim= c(1,5000))
  #NB: Expt sorts were old samples resequenced; uninfected was fresh thaw of library
  #Sorts had very high counts for WT nt; not so for uninfected
reg_in1 = lm(HIV1_inB_norm~HIV1_inA_norm, data = RhDMSfilter3)
summary(reg_in1) #R2 = 0.9875
abline(reg_in1)
  # very good correlation b/w replicate inputs @ nt level

## Compare uninfected to infected (NB: different thaws of library)
plot(RhDMSfilter3$uninfected_norm, RhDMSfilter3$HIV1_inA_norm)
plot(RhDMSfilter3$uninfected_norm, RhDMSfilter3$HIV1_inA_norm, ylim = c(1,5000))
reg_in1Avs0 = lm(HIV1_inA_norm~uninfected_norm, data = RhDMSfilter3)
summary(reg_in1Avs0) #R2 = 0.1828
abline(reg_in1Avs0)
  #Poor correlation; however in unpublished data (removed here), compared infected input library from same thaw
  #In that case, got R2 = 0.9893
#Conclusion: noise between thawed vials, but infection doesn't introduce noise




##### EXPORT to a csv file so I can graph and/or reopen; Also gives non-averaged values
write.csv(RhDMSfilter3, file="../Analysis/Rh-HIV1-NMLV_filter3_noAvg.csv",row.names = FALSE)

### Save point: read in my file again here if I'd like to change the cutoffs below
#RhDMSfilter3 = read.delim('../Analysis/Rh-HIV1-NMLV_filter3_noAvg.csv', sep=",")





##### NEED TO FILTER OUT LOW-REPRESENTATION SEQUENCES, WHICH WILL BE NOISY

### A) Decide representation cutoffs
## Calculate average input counts
RhDMSfilter3[,"in1_avg"] = rowMeans(RhDMSfilter3[21:22])
RhDMSfilter3[,"in2_avg"] = rowMeans(RhDMSfilter3[25:26])

## Plot StdDev against input counts to choose cutoff
plot(RhDMSfilter3$in1_avg, RhDMSfilter3$enrich1SD, log = "x")
plot(RhDMSfilter3$in2_avg, RhDMSfilter3$enrich2SD, log = "x")
  # Lower counts generally noisier, but not a clear cutoff/pattern

## Try against individual, not averaged, input counts
plot(RhDMSfilter3$HIV1_inA_norm, RhDMSfilter3$enrich1SD, log = "x")
plot(RhDMSfilter3$HIV1_inB_norm, RhDMSfilter3$enrich1SD, log = "x")
plot(RhDMSfilter3$HIV1_inA_norm, RhDMSfilter3$enrich1SD, xlim = c(0,1000))
plot(RhDMSfilter3$HIV1_inB_norm, RhDMSfilter3$enrich1SD, xlim = c(0,1000))
plot(RhDMSfilter3$NMLV_inA_norm, RhDMSfilter3$enrich1SD, log = "x")
plot(RhDMSfilter3$NMLV_inB_norm, RhDMSfilter3$enrich2SD, log = "x")
plot(RhDMSfilter3$NMLV_inA_norm, RhDMSfilter3$enrich1SD, xlim = c(0,1000))
plot(RhDMSfilter3$NMLV_inB_norm, RhDMSfilter3$enrich1SD, xlim = c(0,1000))

## Plot enrichment scores against input values
plot(RhDMSfilter3$in1_avg, RhDMSfilter3$enrich1avg, log="x")
plot(RhDMSfilter3$in2_avg, RhDMSfilter3$enrich2avg, log="x")
plot(RhDMSfilter3$in1_avg, RhDMSfilter3$enrich1avg, xlim=c(0,1000))
plot(RhDMSfilter3$in2_avg, RhDMSfilter3$enrich2avg, xlim=c(0,1000))
# still no clear pattern; the lowest representation are not those with the largest scores

## Let's do >50 in each input (normalized) library, to keep consistent w/ Hs cutoffs
  # That won't actually filter any data, fine.


### B) Filter out sequences into separate experiments based on these cutoffs
RhDMS_HIV1 = RhDMSfilter3[which(RhDMSfilter3$HIV1_inA_norm > 50 & RhDMSfilter3$HIV1_inB_norm > 50), ]
RhDMS_NMLV = RhDMSfilter3[which(RhDMSfilter3$NMLV_inA_norm > 50 & RhDMSfilter3$NMLV_inB_norm > 50), ]

## Rename columns for the relevant columns to a generic name in each dataset
colnames(RhDMS_HIV1)
colnames(RhDMS_HIV1)[29:30] = c("enrichA", "enrichB")
colnames(RhDMS_NMLV)[31:32] = c("enrichA", "enrichB")
colnames(RhDMS_HIV1)[33] = "enrichAvg"
colnames(RhDMS_NMLV)[34] = "enrichAvg"
colnames(RhDMS_HIV1)[35] = "enrichSD"
colnames(RhDMS_NMLV)[36] = "enrichSD"




##### AVERAGE ACROSS ALL NT SEQS FOR EACH AA MUTATION

dim(RhDMS_HIV1) # check how many unique sequence we have; I expect 32 x 13 = 416; I have 411
table(table(RhDMS_HIV1[,"seq_aa"])) # determine how many aa seqs have 1/2/3/etc. unique nt seqs
length(unique(RhDMS_HIV1[,"seq_aa"])) # determine # unique aa seqs I have; I expect 20x13 +1 = 261; I have all 261!
## NB: Here, all WT seqs grouped, not kept separate, hence 21 WT nt seqs
length(unique(RhDMS_NMLV[,"seq_aa"])) #261

## Concatenate mutcodon index with aa seq to differentiate WT seqs
RhDMS_HIV1[,"mutcodon+seq"] = paste(RhDMS_HIV1$mutcod_index, RhDMS_HIV1$seq_aa, sep="_")
RhDMS_NMLV[,"mutcodon+seq"] = paste(RhDMS_NMLV$mutcod_index, RhDMS_NMLV$seq_aa, sep="_")
length(unique(RhDMS_HIV1$`mutcodon+seq`)) # determine # unique aa seqs I have; I expect 21 x 13 = 273; I have all 273

## Split whole dataframe into list of dataframes for eqch unique aa seq
RhDMS_HIV1Split <- split(RhDMS_HIV1, RhDMS_HIV1[,"mutcodon+seq"]) #2nd argument = what to split based upon
head(RhDMS_HIV1Split,3)
RhDMS_NMLVSplit <- split(RhDMS_NMLV, RhDMS_NMLV[,"mutcodon+seq"])

## Define function to average values; also keeps only useful columns
get_averages = function(datf){  #named function that will act on dataframe
  datf[,"num_nt_seqs"] = dim(datf)[1] #identify # unique nt seqs (rows) for each aa seq (subsetted dataframe)
  datf[,"enrich_Avg_RepAndAA"] = mean(datf[,"enrichAvg"]) # add column with the mean (across reps & synonymous seqs)
  datf[,"enrichA_Avg"] = mean(datf[,"enrichA"]) # add column w/ mean (across synonymous seqs only)
  datf[,"enrichB_Avg"] = mean(datf[,"enrichB"])
  datf[,"enrich_SD_RepAndAA"] = sd(c(datf[,"enrichA"],datf[,"enrichB"]))
  columnsToKeep <- c("dist_codon","seq_aa","mutcod_index","mut_aa", "num_nt_seqs","enrichA_Avg","enrichB_Avg","enrich_Avg_RepAndAA","enrich_SD_RepAndAA")
  datf <- datf[1, columnsToKeep] # overwrite dataframe (w/in fxn only) with the first row & only these columns
  return(datf)
}

## Perform this averaging function on all unique aa seqs
RhDMS_HIV1SplitAvgs = lapply (RhDMS_HIV1Split, get_averages)
head(RhDMS_HIV1SplitAvgs,3)
RhDMS_NMLVSplitAvgs = lapply (RhDMS_NMLVSplit, get_averages)

## Bring all these values back into single table
RhDMS_HIV1Avgs = do.call("rbind", RhDMS_HIV1SplitAvgs) #row bind each of these rows
rownames(RhDMS_HIV1Avgs) = NULL #gets rid of useless row names
head(RhDMS_HIV1Avgs,3)
RhDMS_NMLVAvgs = do.call("rbind", RhDMS_NMLVSplitAvgs)
rownames(RhDMS_NMLVAvgs) = NULL 

## Save output (enrichment scores averaged across synonymous codons) as CSV file
write.csv(RhDMS_HIV1Avgs, file="../Analysis/Rh-HIV1_avg.csv",row.names = FALSE)
write.csv(RhDMS_NMLVAvgs, file="../Analysis/Rh-NMLV_avg.csv",row.names = FALSE)




##### GENERATE A POSITION X AMINO ACID MATRIX
## make sure nothing is a factor!
class(RhDMS_HIV1Avgs$mut_aa) #factor no good
RhDMS_HIV1Avgs$mut_aa = as.character(RhDMS_HIV1Avgs$mut_aa)
RhDMS_NMLVAvgs$mut_aa = as.character(RhDMS_NMLVAvgs$mut_aa)
class(RhDMS_HIV1Avgs$mutcod_index) #integer OK

##null values are a problem for populating the table
## 1) replace null values in mut_aa column (synonymous w/ WT) with "syn"
RhDMS_HIV1Avgs$mut_aa[RhDMS_HIV1Avgs$mut_aa == ""] = "syn"
RhDMS_NMLVAvgs$mut_aa[RhDMS_NMLVAvgs$mut_aa == ""] = "syn"

## 2) replace single null value in index position (truly WT nt seq)
## Only Q332 (index1) has no synonymous codons -> put it there
RhDMS_HIV1Avgs$mutcod_index[is.na(RhDMS_HIV1Avgs$mutcod_index)] = 1
RhDMS_NMLVAvgs$mutcod_index[is.na(RhDMS_NMLVAvgs$mutcod_index)] = 1

## Make empty dataframe
RhDMSmatrix_1 = data.frame()
RhDMSmatrix_2 = data.frame()

## Populate dataframe
for (rowindex in 1:dim(RhDMS_HIV1Avgs)[1]) { #iterate over the # of rows (first value of dimensions of table)
  RowNum = RhDMS_HIV1Avgs[rowindex,"mutcod_index"] #define the row # as the index position
  ColNum = RhDMS_HIV1Avgs[rowindex, "mut_aa"] #define the col # as the amino acid
  RhDMSmatrix_1[RowNum,ColNum] = RhDMS_HIV1Avgs[rowindex, "enrich_Avg_RepAndAA"] #populate dataframe at row/col with the enrichment score from that row
}

for (rowindex in 1:dim(RhDMS_NMLVAvgs)[1]) { 
  RowNum = RhDMS_NMLVAvgs[rowindex,"mutcod_index"] 
  ColNum = RhDMS_NMLVAvgs[rowindex, "mut_aa"] 
  RhDMSmatrix_2[RowNum,ColNum] = RhDMS_NMLVAvgs[rowindex, "enrich_Avg_RepAndAA"] 
}

## Add column with proper position #
RhDMSmatrix_1[,"Position"] = c("Q332","A333","P334","G335","T336","L337","F338", "T339","F340","P341","S342","L343","T344")
RhDMSmatrix_2[,"Position"] = c("Q332","A333","P334","G335","T336","L337","F338", "T339","F340","P341","S342","L343","T344")

## Reorder rows into my favorite order
Rh_Expt1Matrix = RhDMSmatrix_1[,c("Position","G","P","A","I","L","V","F","W","Y","M","C","S","T","N","Q","H","K","R","D","E","*","syn")]
Rh_Expt2Matrix = RhDMSmatrix_2[,c("Position","G","P","A","I","L","V","F","W","Y","M","C","S","T","N","Q","H","K","R","D","E","*","syn")]

## Save this matrix as csv file
write.csv(Rh_Expt1Matrix, file="../Analysis/Rh-HIV1_matrix.csv",row.names = FALSE)
write.csv(Rh_Expt2Matrix, file="../Analysis/Rh-NMLV_matrix.csv",row.names = FALSE)
