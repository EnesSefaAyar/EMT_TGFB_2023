source("functions_parameters.R")

# User specific

# Reference channel number (1-11, or 1-16)
ref_channel<-2

# Add your cell type labels, must match those used in experimental design
your_labels<-c("d0","d3", "d9", "neg")
your_control_label<-"neg"

# Import ------------------------------------------------------------------

# Load raw data 
# NOTE: have commented all preDART commands out where there appear to be more options

  ev<-read.delim("ev_updated_WallETMTDIA.txt")

# Parse sp|Q00000|HUMAN_XXX into just Uniprot accession: Q00000  

parse_row<-grep("|",ev$Leading.razor.protein, fixed=T)

if(length(parse_row)>0){
  split_prot<-str_split(ev$Leading.razor.protein[parse_row], pattern = fixed("|"))
  split_prot2<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  ev$Leading.razor.protein[parse_row]<-split_prot2
}
# Load experimental design and batches

design<-read.csv("annotation.csv")

batch<-read.csv("batch.csv")

# Attach batch data to protein data
ev[,colnames(batch)[-1]]<-NA
for(X in batch$set){
  
  ev$lcbatch[ev$Raw.file==X] <- as.character(batch$lcbatch[batch$set%in%X])
  ev$sortday[ev$Raw.file==X] <- as.character(batch$sortday[batch$set%in%X])
  ev$digest[ev$Raw.file==X] <- as.character(batch$digest[batch$set%in%X])
  
}


# Create unique peptide+charge column:
ev$modseq<-paste0(ev$Modified.sequence,ev$Charge)

# Add X in front of experiment names because R doesn't like column names starting with numbers
ev$Raw.file<-paste0("X",ev$Raw.file)
design$Set<-paste0("X",design$Set)


# Which columns hold the TMT Reporter ion (RI) data
ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.",1:16))
#ri.index<-which(colnames(ev)%in%paste0("Reporter.intensity.corrected.",1:16))

# Make sure all runs are described in design, if not, print and remove them:
not.described<-unique(ev$Raw.file)[ !unique(ev$Raw.file) %in% paste0(design$Set) ]
ev<-ev[!ev$Raw.file%in%not.described,]

# Filter out reverse hits, contaminants, and contaminated spectra...
ev<-ev[-which(ev$Reverse=="+"),]
if(length(grep("REV", ev$Leading.razor.protein))>0){ ev<-ev[-grep("REV", ev$Leading.razor.protein),] }
if(length(grep("CON", ev$Leading.razor.protein))>0){ ev<-ev[-grep("CON", ev$Leading.razor.protein),] }
if(length(which(ev$Potential.contaminant=="+"))>0){ ev<-ev[-which(ev$Potential.contaminant=="+"),] }

#ev<-ev[!is.na(ev$Raw.file),]

ev<-ev[!is.na(ev$PIF),]
ev<-ev[ev$PIF>0.8,]

# Remove peptides that are more the 10% the intensity of the carrier in the single cell runs (only)
ev<-as.data.frame(ev)
ev$mrri<-0
ev$mrri <- rowMeans(ev[, ri.index[4:length(ri.index)]] / ev[, ri.index[1]], na.rm = T)
ev<-ev[ev$mrri < 0.1, ]


# Filter by PEP or FDR: CHOOSE ONE
## Regular
# ev<-ev[ev$PEP<0.02, ]
# ev<-ev[calc_fdr(ev$PEP)<0.01, ]
# ev<-Empirical_FDR_PEP(ev) 
# ev<- ev[ev$qVal < 0.01,] 

## DART
ev<-ev[calc_fdr(ev$dart_PEP)<0.01, ]
# ev<-ev[ev$dart_PEP<0.02, ]
# ev<-Empirical_FDR_dPEP(ev) 
# ev<- ev[ev$qVal < 0.01,] 


# Normalize single cell runs to normalization channel
ev<-as.data.frame(ev)
ev[, ri.index] <- ev[, ri.index] / ev[, ri.index[ref_channel]]


# Organize data into a more convenient data structure:
# Create empty data frame
ev.melt<-melt(ev[0, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev)[ri.index]) ],
              id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"))

colnames(ev.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")


# Record mapping of cell type to Channel:
ct.v<-c()
qt.v<-c()
rt.v <-c()

# Create a unique ID string
unique.id.numeric<-1:length(ri.index)
unique.id<-paste0("i",unique.id.numeric)

RI_keep<-ri.index

# Give each sample a unique identifier
for(X in unique(ev$Raw.file)){
  
  # Subset data by X'th experiment
  ev.t<-ev[ev$Raw.file%in%X, ]
  
  # Name the RI columns by what sample type they are: carrier, single cell, unused, etc...
  colnames(ev.t)[ri.index]<-paste0(as.character(unlist(design[design$Set==X,-1])),"-", unique.id)
  
  # Melt it! and combine with other experimental sets
  ev.t.melt<-melt(ev.t[, c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest", colnames(ev.t)[RI_keep]) ],
                  id.vars = c("Raw.file","modseq","Leading.razor.protein","lcbatch","sortday","digest"));
  
  # Record mapping of cell type to Channel:
  ct.v<-c(ct.v, unique.id[which(ri.index%in%RI_keep)] )
  qt.v<-c(qt.v, colnames(ev)[RI_keep] )
  rt.v <- c(rt.v, rep(X, times = length(ri.index)))
  
  colnames(ev.t.melt)<-c("Raw.file","sequence","protein","lcbatch","sortday","digest","celltype","quantitation")
  
  ev.melt<-rbind(ev.melt, ev.t.melt)
  
  # Update unique ID string
  unique.id.numeric<-unique.id.numeric + length(ri.index)
  unique.id<-paste0("i", unique.id.numeric)
  
}

# Grab the unique number associate to each and every cell, carrier channel, and empty channel
ev.melt$id<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(2,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$celltype<-unlist(strsplit(as.character(ev.melt$celltype),"-"))[seq(1,length(unlist(strsplit(as.character(ev.melt$celltype),"-"))),2)]
ev.melt$id<-as.factor(ev.melt$id)

# Remove duplicate observations of peptides from a single experiment
ev.melt<-remove.duplicates(ev.melt,c("sequence","id") )
ev.melt<-ev.melt[!is.na(ev.melt$protein), ]

# Create additional meta data matrices
ev.melt.uniqueID<-remove.duplicates(ev.melt,"id")
ev.melt.pep<-remove.duplicates(ev.melt, c("sequence","protein") )

# Create data frame of peptides x cells, populated by quantitation
ev.unmelt<-dcast(ev.melt, sequence ~ id, value.var = "quantitation", fill=NA)

# Also create matrix of same shape
ev.matrix<-as.matrix(ev.unmelt[,-1]); row.names(ev.matrix)<-ev.unmelt$sequence

# Replace all 0s with NA
ev.matrix[ev.matrix==0]<-NA
ev.matrix[ev.matrix==Inf]<-NA
ev.matrix[ev.matrix==-Inf]<-NA

# Divide matrix into single cells (including intentional blanks) and carriers
sc_cols<-unique(ev.melt$id[(ev.melt$celltype%in%c(your_labels))])
ev.matrix.sc<-ev.matrix[, colnames(ev.matrix)%in%sc_cols]

#write.table(ev.matrix.sc, "EpiToMesen.TGFB.nPoP_trial1_filteredButNotNormalised.txt", sep = "\t", row.names = TRUE)


# Filter single cells ----------------------------------------------------------------------


sc.melt<-ev.melt

xd<-as_tibble( sc.melt )

xd <- xd %>% group_by(id) %>% mutate(med_per_c = median(quantitation, na.rm=T)); length(unique(xd$id))

length(unique(xd$id))

xd$quantitation[(xd$quantitation)==Inf]<-NA
xd$quantitation[(xd$quantitation)==0]<-NA

xd <- xd %>% mutate_if(is.factor, as.character)

xd1 <- xd %>%
  group_by(id) %>%
  mutate(norm_q1 = quantitation / median(quantitation, na.rm=T))

xd2 <- xd1 %>%
  group_by(sequence, Raw.file) %>%
  mutate(norm_q = quantitation / mean(norm_q1, na.rm=T))

xd3<- xd2 %>%
  filter(celltype%in%c(your_labels))

xd4<- xd3 %>%
  group_by(protein, id) %>%
  mutate(cvq = cv(norm_q))

xd5<- xd4 %>%
  group_by(protein, id) %>%
  mutate(cvn = cvna(norm_q))

xd6<- xd5 %>%
  filter(cvn > 5)

xd7<-xd6 %>% group_by(id) %>% mutate(cvm=median(cvq, na.rm=T))

xdf<-xd7

cvPar <- 0.40

print("Number of unique proteins used in calculation:", length(unique(xdf$protein)))

hist(unique(xdf$cvm[xdf$celltype!=your_control_label]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(xdf$cvm[xdf$celltype==your_control_label]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=40)


# USER TUNED

hist(unique(xdf$cvm[xdf$celltype!=your_control_label]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(xdf$cvm[xdf$celltype==your_control_label]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=40)
abline (v = cvPar)
# Filter out variable wells and controls
sc_kept<-unique( xdf$id[xdf$celltype!=your_control_label & xdf$cvm < cvPar])
sc0_kept<-unique( xdf$id[xdf$celltype==your_control_label & xdf$cvm > cvPar])

# Which wells to keep
keep_these<-unique( xdf$id)

sc_total<-unique( xdf$id[xdf$celltype!=your_control_label])
sc0_total<-unique( xdf$id[xdf$celltype==your_control_label])
scrate<-round(length(sc_kept) / length(sc_total),2)*100

ev.matrix.sc.f<-ev.matrix.sc[,colnames(ev.matrix.sc)%in%sc_kept]; dim(ev.matrix.sc.f)
ev.matrix.sc.f[ev.matrix.sc.f==Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==-Inf]<-NA
ev.matrix.sc.f[ev.matrix.sc.f==0]<-NA

xdf$control<-"sc"
xdf$control[xdf$celltype==your_control_label]<-"ctl"


my_col3<-c( "black", "purple2")

# Plot!
ggplot(data=xdf, aes(x=cvm)) + geom_density(aes(fill=control, alpha=0.5), adjust=4) + theme_pubr() +
  scale_fill_manual(values=my_col3[c(1,2)]) +
  xlab("Quantification variability") + ylab("Density") + rremove("y.ticks") + rremove("y.text") +
  font("xylab", size=35) +
  font("x.text", size=30) +
  coord_cartesian(xlim=c(0,1.5))+
  #xlim(c(-0.15, 0.35)) +
  # annotate("text", x=0.27, y= 14, label=paste0(scrate,"% single cells passed"), size=8, color=my_col3[c(2)])+
  # annotate("text", x=0.27, y= 12.5, label=paste0(sc0rate,"% control wells passed"), size=8, color=my_col3[c(1)])+
  annotate("text", x=0.172, y= 14, label=paste0(length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.6, y= 12, label=paste0(length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  annotate("text", x=0.6, y= 14, label=paste0(length(sc_total) -length(sc_kept)," single cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.165, y= 12, label=paste0(length(sc0_total) - length(sc0_kept)," control wells"), size=10, color=my_col3[c(1)])+
  #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) +
  rremove("legend") + geom_vline(xintercept=0.4, lty=2, size=2, color="gray50")




# Data transformations ----------------------------------------------------



# Perform normalizations / transformations in multiple steps with visual sanity checks:
b.t<-"FD"
xlim.t<-c(-2,2)
par(mfrow=c(3,3))

# Original data, normalized to reference channel, filtered for failed wells:
t0<-ev.matrix.sc.f

hist(c(t0), breaks=b.t, xlim=xlim.t)

# Column then row normalize by median or mean (see source functions):
#t1<-cr_norm(t0)
t1<-medInt_cr_norm(t0)

hist(c(t1), breaks=b.t, xlim=xlim.t)


# Filter for missing data:
t2<-filt.mat.rc(t1, na.row, na.col)
hist(c(t2), breaks=b.t, xlim=xlim.t)

# Log2 transform:
t3<-log2(t2)
t3[t3==Inf]<-NA
t3[t3==-Inf]<-NA
t3[t3==0]<-NA
hist(c(t3), breaks=b.t, xlim=xlim.t)


# # Collapse to protein level by median:
t3m<-data.frame(t3)
t3m$pep<-rownames(t3)
t3m$prot <- ev.melt.pep$protein[match(t3m$pep, ev.melt.pep$sequence)]

t3m<-melt(t3m, variable.names = c("pep", "prot"))
colnames(t3m) <-c("pep","prot","id","quantitation")
t3m2<- t3m %>% group_by(prot, id) %>% dplyr::summarize(qp = median(quantitation, na.rm=T))
t4m<-dcast(t3m2, prot ~ id, value.var = "qp", fill=NA)
t4<-as.matrix(t4m[,-1]); row.names(t4)<-t4m[,1]
hist(c(t4), breaks=b.t, xlim=xlim.t)

# Re-column and row normalize:
t4b<-medInt_cr_norm_log(t4)
hist(c(t4b), breaks=b.t, xlim=xlim.t)

# Assign to a final variable name:
ev.matrix.sc.f.n<-t4b

## Impute single celldata
imp.input<-ev.matrix.sc.f.n
sc.imp <- hknn(imp.input, k.t)
t5<-sc.imp
sum(is.na(sc.imp))
dim(sc.imp)

#write.csv(t4b, "EpiToMesen.TGFB.nPoP_trial1_1PercDartFDRTMTBulkDIA.WallE_unimputed.txt", row.names=T)

matrix.sc.batch <- t5

# visual sanity checks post-imputation:
hist(c(t5), breaks=b.t, xlim=xlim.t)


# PC's to display:
PCx<-"PC1"
PCy<-"PC2"

mat.sc.imp<-medInt_cr_norm_log(matrix.sc.batch)

# Dot product of each protein correlation vector with itself
r1<-rcorr(t(matrix.sc.batch))$r
rsum<-rowSums(r1^2)

# Calculate the weighted data matrix:
X.m <- mat.sc.imp
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- rcorr(X.m)$r

# PCA
sc.pca<-eigen(pca.imp.cor)
scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

# Percent of variance explained by each principle component
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
plot(1:length(percent_var), percent_var, xlab="PC", ylab="% of variance explained")

# Map meta data
pca.melt <- melt(scx); colnames(pca.melt)<-c("id","pc","value")
#add.cols<-colnames(ev.melt)[c(4:7,10,11)]
add.cols<-colnames(ev.melt)[c(4:7)]

pca.melt[,add.cols]<-NA

for(X in unique(pca.melt$id)){
  
  pca.melt[pca.melt$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]
  
}


# Re map ...
pca.display <- dcast(pca.melt, id ~ pc, value.var = "value", fill=NA)

pca.display[,add.cols]<-NA

for(X in unique(pca.display$id)){
  
  pca.display[pca.display$id==X, add.cols]<-ev.melt.uniqueID[ev.melt.uniqueID$id==X, add.cols]
  
}


### Make PCA for publication
# basically, get timepoint to cell ID map, join and plot
pca.display_forPub <- pca.display
timepointToID<-read.table("cellIDToTimepoint_WallETMTDIADART.txt", header = TRUE)

pca.display_forPub <- pca.display_forPub %>% left_join(timepointToID, by = "id")
pca.display_forPub$celltype <- str_replace_all(pca.display_forPub$celltype, c("d0" = "Day 0", "d3" = "Day 3","d9" = "Day 9"))

# Display 
ggscatter(pca.display_forPub,color = 'celltype' ,x =PCx, y = PCy, size = 5, alpha=0.5, palette = c("#F8766D", "grey", "#619CFF")) +
  xlab(paste0(PCx,"  (", round(percent_var[1],0),"%)")) +
  ylab(paste0(PCy,"  (", round(percent_var[2],0),"%)")) +
  font("ylab",size=30) +
  font("xlab",size=30) +
  font("xy.text", size=20)+
  font("legend.text", size = 14) +
  rremove("legend.title") +
  #scale_color_manual(values = my_colors[2:3]) +
  #annotate("text", x=-0.025, y=0.21,label="Macrophage", color=my_colors[2], size=10)  +
  #annotate("text", x=0.03, y=0.21, label="Monocyte", color=my_colors[3], size=10) +
  annotate("text", x=-0.042, y=-0.11, label=paste0(dim(mat.sc.imp)[1], " proteins"), size=8) +
  annotate("text", x=-0.05, y=-0.14, label=paste0(dim(mat.sc.imp)[2], " cells"), size=8) 
