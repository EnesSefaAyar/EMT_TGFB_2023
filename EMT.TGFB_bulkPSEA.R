### Script to run Protein Set enrichment Analysis on Pre-formatted DIANN outputs
## Also included code at the bottom of this script to generate output for this format using DIANN files

### Load libraries
library(ggplot2)
library(gplots)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(UniProt.ws)
library(progress)


### Data Import

# Your data should be formatted in terms of a melted list with the following characteristics:
# 1. Condition column
# 2. Relative Protein Abundance Column
# 3. Column of Uniprot IDs
# 4. (optional) Column of Genes
# Make sure to title the relevant columns: [Protein, Intensity, Condition, Gene]


# Hi friend, I'm sorry, this is an old script and is super inefficient, across the board. 
# In it's current unfortunate state, you'll have to reassign each BioRep to the "data" variable below, re-run and assign to output etc. I sorry. Much shames. 

## Assign data paths
DataPath <- ".../EMT.Bulk1_ProcessedInputForPSEA.txt"
DataPath <- ".../EMT.Bulk2_ProcessedInputForPSEA.txt"

GO_DB_Path <- ".../GOA_db.txt"


## Import and process data

dataOne <- read.table(DataPath, header = TRUE, stringsAsFactors = F, sep = "\t")
dataTwo <- read.table(DataPath, header = TRUE, stringsAsFactors = F, sep = "\t")


dataOne$Condition <- str_replace_all(dataOne$Condition, c("d0" = "d0_brep1", "d3" = "d3_brep1","d9" = "d9_brep1"))

dataTwo$Condition <- str_replace_all(dataTwo$Condition, c("d0" = "d0_brep2", "d3" = "d3_brep2","d9" = "d9_brep2"))

dataOne <- dataOne[,c(1:4)]
dataTwo <- dataTwo[,c(1:4)]

data <- dataOne
data <- dataTwo


data <- rbind(dataOne,dataTwo)



## GO db
GO_db <- read.delim(GO_DB_Path, sep = " ")

GO_db <- GO_db %>% dplyr::select(Gene,GO_term_name)
GO_db <- unique(GO_db)
GO_db$Gene <- toupper(GO_db$Gene)

## Prot to Gene Map
Prot_to_Gene_map <- read.delim('.../uniprot-proteome_UP000005640.txt')
colnames(Prot_to_Gene_map) <- c("Protein","Gene")
Prot_to_Gene_map <- Prot_to_Gene_map[,c("Protein", "Gene")]

###### ffs, fix this
## again, I'm so sorry. 

forIdiotMapping <- data
data$ForJoins <-gsub("-.*","",forIdiotMapping$Protein)
Prot_to_Gene_map$ForJoins <- Prot_to_Gene_map$Protein
data <- data %>% left_join(Prot_to_Gene_map, by = "ForJoins")
data <- data[ ,-which(names(data) %in% c("ForJoins","Protein.y"))]
data$Gene <-gsub("_.*","",data$Gene)

colnames(data) <- c("Protein","Intensity","Condition", "Gene")
data <- data[,1:4]



### PSEA 
unique_GO_terms <- unique(GO_db$GO_term_name)
KW_output <- data.frame(GO_term = character(), pVal = numeric(),numberOfMatches= numeric(),fractionOfDB_Observed = numeric(),stringsAsFactors = FALSE)
for (i in 1:length(unique_GO_terms)){
  GO_term <- unique_GO_terms[i]
  GO_db_lim <- GO_db %>% dplyr::filter(GO_term_name == GO_term)
  
  data_matches <- data %>% dplyr::filter(Gene %in% GO_db_lim$Gene)
  matches_number <- length(unique(data_matches$Protein))
  total_Proteins_in_GOterm <- length(unique(GO_db_lim$Gene))
  fractionObserved <- matches_number/total_Proteins_in_GOterm
  dataProt_matches_medInt <- data_matches %>% dplyr::group_by(Condition)  %>% dplyr::summarise(medianInt = median(Intensity, na.rm = T))
  dataProt_matches_medInt$GO_term <- GO_term
  
  
  ## Kruskall Wallis
  # due to poor designs, if this fails out, uncomment matches_number > 0, run stop, re-run. Super easy fix to code, but...
  #if(matches_number > 0){
  if(fractionObserved > 0.3 & total_Proteins_in_GOterm < 55 & !is.na(total_Proteins_in_GOterm) &   !is.na(fractionObserved)){
    KW_out <- kruskal.test(Intensity ~ Condition, data = data_matches)
    pVal <- KW_out$p.value
  }
  if(matches_number == 0){
    pVal = NA
  }
  KW_output[i,1] <- GO_term
  KW_output[i,2] <- pVal 
  KW_output[i,3] <- matches_number
  KW_output[i,4] <- fractionObserved
  
  if(i == 1){
    KW_Int_output <- dataProt_matches_medInt[0,]
    
  }
  if(i > 1){
    KW_Int_output <- rbind(KW_Int_output, dataProt_matches_medInt)
    
  }
  
}  

### To check distributions of proteins and GO terms
hist(KW_output$numberOfMatches/KW_output$fractionOfDB_Observed, breaks = 50, xlab = 'Proteins per GO term', main = 'Distribution of proteins per GO term')
abline(v = 25, col = 'red')

hist(KW_output$fractionOfDB_Observed, breaks = 50, xlab = 'Fraction of GO term matching Sample', main = 'Distribution of fraction matching GO term')
abline(v = 0.3, col = 'red')



# Multiple hypothesis testing correction
KW_output$qVal <- p.adjust(KW_output$pVal, method = "BH")
# Join Effect size
KW_output_comb <- KW_Int_output %>% left_join(KW_output, by = "GO_term")

KW_output_comb_effectSize <- KW_output_comb %>% group_by(GO_term) %>% dplyr::summarise(EffectSize = max(abs(medianInt)))
KW_output_comb <- KW_output_comb %>% left_join(KW_output_comb_effectSize, by = "GO_term")

# filter below qVal
KW_output_filt_comb <- KW_output_comb %>% filter(qVal <= .05)
# bounding effect size into quantiles
KW_Effect_quantile <- quantile(KW_output_comb_effectSize$EffectSize,c(.95,.75,.50,.20), na.rm = T)

### check distribution of p, q vals
par(mfrow = c(1,2))
hist(log10(KW_output_filt_comb$pVal), main = "Distribution of p values - KW", xlab = "log10(pvalue")
hist(log10(KW_output_filt_comb$qVal), main = "Distribution of q values - BH, 5%", xlab = "log10(qvalue")



## Heatmap
# here we are, the nastiness. 
filtered_brep1_pre <- KW_output_filt_comb
filtered_brep1_effect <- KW_output_comb_effectSize

#filtered_brep2_pre <- KW_output_filt_comb
#filtered_brep2_effect <- KW_output_comb_effectSize

filtered_brep1_pre$Condition <- str_replace_all(filtered_brep1_pre$Condition, c("d0" = "d0_brep1", "d3" = "d3_brep1","d9" = "d9_brep1"))
filtered_brep2_pre$Condition <- str_replace_all(filtered_brep2_pre$Condition, c("d0" = "d0_brep2", "d3" = "d3_brep2","d9" = "d9_brep2"))

KW_Effect_quantile_brep1 <- quantile(filtered_brep1_effect$EffectSize,c(.80,.50,.40,.20))
KW_Effect_quantile_brep2 <- quantile(filtered_brep2_effect$EffectSize,c(.80,.50,.40,.20))


filtered_brep1 <- filtered_brep1_pre %>% filter(EffectSize > KW_Effect_quantile_brep1[[1]])
filtered_brep2 <- filtered_brep2_pre %>% filter(EffectSize > KW_Effect_quantile_brep2[[1]])

filtered <- rbind(filtered_brep1,filtered_brep2)

plottingFormat <- na.omit(data.matrix(acast(filtered, GO_term ~ Condition, value.var = 'medianInt')))

## These were the GO terms we selected based on an intermediary between effect size and biology 
bulk_GOTerms <- c("podosome", "cortical cytoskeleton","filamentous actin", "superoxide metabolic process", "stress fiber", "extracellular matrix disassembly", "substrate adhesion-dependent cell spreading","fibril organization", "spindle organization", "ligand-gated ion channel activity", "hemidesmosome assembly", "L-ascorbic acid binding","DNA-dependent DNA replication initiation","deoxyribonucleoside monophosphate biosynthetic process", "desmosome", "platelet dense granule organization", "negative regulation of fatty acid biosynthetic process", "cell-cell junction organization", "mitotic cell cycle spindle assembly checkpoint", "heat shock protein binding")


reordered_plottingFormat <- plottingFormat[match(bulk_GOTerms, rownames(plottingFormat)), ] %>% na.omit

redBlue <- c("navy","white", "red") 
pal <- colorRampPalette(redBlue)(50)


heatmap.2(reordered_plottingFormat, col = pal, trace = "none", density.info = "none", margins = c(8, 22), key.xlab = NA, key.title = "", Colv = F, dendrogram = "row", key.par=list(mar=c(8,1,1,2)), breaks=seq(-1,1,length.out=51), cex.row = 1.5)


### END PSEA Script



### code to generate files to input in this pipeline
# load packages
library(plyr)
library(dplyr)
library(stringr)
library(diann)

# assign functions
cr_norm_log<-function(dat){
  
  for(k in 1:ncol(dat)){
    
    dat[,k]<-dat[,k]-median(dat[,k], na.rm = T)
    
    
  }
  
  
  for(k in 1:nrow(dat)){
    
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
    
  }
  
  return(dat)
}

## Again, repetitivly for each Brep *sigh*
## file paths
data_path <- ".../EMT.Bulk1_DiannReport.tsv" 
#data_path <- ".../EMT.Bulk2_DiannReport.tsv.tsv" 

## preprocess
evi <- read.delim(data_path)
evi <- evi[which(evi$Lib.PG.Q.Value < 0.01),]
#evi <- evi[which(evi$CScore > 0.95),]
#evi <- evi[which(evi$Ms1.Profile.Corr > 0),]
evi <- evi[grepl("HUMAN", evi$Protein.Names),]

evi$run <- paste0(evi$Run)
evi[,grepl("Run",colnames(evi))] <- evi %>% dplyr::select(run)

evi <- evi[which(evi$Proteotypic == T),]
evi$Protein.Group <- gsub(";.*", "", evi$Protein.Group)

evi$seqcharge <- paste0(evi$Precursor.Id)

evi <- evi[which(evi$Precursor.Translated > 0),]

## collapse to protein
evi_proteinLevel <- diann::diann_maxlfq(evi, sample.header = "Run", group.header = "Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")

## norm at protein level
evi_proteinLevel <- log2(evi_proteinLevel)
evi_proteinLevel <- cr_norm_log(evi_proteinLevel)
evi_proteinLevel <- na.omit(evi_proteinLevel)
evi_proteinLevel <- reshape2::melt(evi_proteinLevel)

colnames(evi_proteinLevel) <-c("Protein", "Condition", "Intensity")

## Reassign file names to timepoint
# brep 1
evi_proteinLevel$Condition <- str_replace_all(evi_proteinLevel$Condition, c("eSK154-001" = "d0", "eSK155-002" = "d3","eSK156-003" = "d9"))
# brep 2
#evi_proteinLevel$Condition <- str_replace_all(evi_proteinLevel$Condition, c("eSK327_notpass" = "d0", "eSK328" = "d3","eSK329" = "d9"))

## write out tables
write.table(evi_proteinLevel, ".../EMT.Bulk1_ProcessedInputForPSEA.txt", sep = "\t", row.names = FALSE)
#write.table(evi_proteinLevel, ".../EMT.Bulk2_ProcessedInputForPSEA.txt", sep = "\t", row.names = FALSE)




