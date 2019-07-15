# set working directory
setwd("/path_to_directory/")
source("source.R")
# set the number of cpus 
BPPARAM=MulticoreParam(4)
# Load VCF file
myData<-read.table("file.vcf", sep="\t")
myInfo<-read.delim("file.vcf", sep="\t", nrows=1, skip=9, header=F, stringsAsFactors = FALSE)
myInfo[,1]<-"CHROM"
names(myData)<-as.character(myInfo[1,])
#variable exploration
table(myData$CHROM)
# check duplicate SNPs
length(unique(paste(myData$CHR, myData$POS, sep=":")))
length(unique(paste(myData$CHR, myData$POS, myData$ALT,sep=":")))
# Exist more than one SNP in the same position. Indeed, there are SNPs identified
#with different ID and different REF because the REF is defined using the criteria
# REF allele is which have the highest frequency.
length(which(duplicated(paste(myData$CHR, myData$POS,sep=":"))))
#check alternative and reference allele frequencies
summary(myData$REF)
summary(myData$ALT)
# check quality
table(myData$QUAL)
# check filter flag
table(myData$FILTER)
# check SNP info
length(unique(myData$INFO))
## Filter redundant columns
myData<-myData[,c(-6, -7,-9)]
# separate the info column
infoExpand<-do.call(rbind, lapply(1:nrow(myData), function(x){
  aux<-strsplit(as.character(myData$INFO[x]), split=";")[[1]]
  NS<-strsplit(aux[1], split="=")[[1]][2]
  aux2<-strsplit(aux[2], split="=")[[1]][2]
  AF_REF<-strsplit(aux2, split=",")[[1]][1]
  AF_ALT<-strsplit(aux2, split=",")[[1]][2]
  return(cbind(NS=NS, AF_REF=AF_REF, AF_ALT=AF_ALT))
}))
infoExpand<-as.data.frame(infoExpand)
infoExpand$AF_REF<-as.numeric(as.character(infoExpand$AF_REF))
infoExpand$AF_ALT<-as.numeric(as.character(infoExpand$AF_ALT))
infoExpand$NS<-as.numeric(as.character(infoExpand$NS))
# Check number of genotiped samples
table(infoExpand$NS)
# check frequencies consistency
all((infoExpand$AF_REF + infoExpand$AF_ALT) ==1)
# build the expanded dataset
myDataFormat<-cbind(myData[, 1:5], infoExpand, myData[, 7:ncol(myData)])
## Filter older samples
old<-c("Plate2_PMA126.sorted", "Plate2_PMA132.sorted","Plate1_PMA77.sorted",
       "Plate1_PMA64.sorted","Plate2_PMA123.sorted","Plate2_PMA125.sorted")
myDataFormat<-myDataFormat[,!names(myDataFormat) %in% old]
# recompute allele frequencies and NS field
AF_ALT_NS<-do.call(rbind, bplapply(1:nrow(myDataFormat), function(x){
  AF_ALT<-calculateAFALT(myDataFormat[x,9:ncol(myDataFormat)])
  NS<-getNS(myDataFormat[x,9:ncol(myDataFormat)])
  return(cbind(AF_ALT=AF_ALT,NS=NS))
}, BPPARAM = BPPARAM))

myDataFormat$AF_ALT<-AF_ALT_NS[,"AF_ALT"]
myDataFormat$AF_REF<-1-AF_ALT_NS[,"AF_ALT"]
myDataFormat$NS<-AF_ALT_NS[,"NS"]
# remove SNPS without samples
myDataFormat<-myDataFormat[which(!is.na(myDataFormat$AF_ALT)),]
##Duplicated SNPs
aux<-paste(myDataFormat$CHR, myDataFormat$POS,sep=":")
dups<-which(duplicated(aux))
# combine duplicated SNPs
myDataFormatNoDup<-mergeSNPs(dups, myDataFormat,ncolInfo=8,patNA="./.",
                             BPPARAM=MulticoreParam(4))
# invert frequencies
idxInv<-which(myDataFormatNoDup$AF_REF<0.5) 
invertedGTs<-exchangeGT(myDataFormatNoDup,idxInv, BPPARAM = MulticoreParam(4))

for(i in 9:ncol(myDataFormatNoDup)){
    myDataFormatNoDup[,i]<-as.character(myDataFormatNoDup[,i])
    
}
myDataFormatNoDup[idxInv,]<-invertedGTs
# filter low frequency SNPs
MAF<-0.05
myDataFormatNoDup<-myDataFormatNoDup[myDataFormatNoDup$AF_ALT >=MAF & 
                                         myDataFormatNoDup$AF_REF >=MAF,]
summary(myDataFormatNoDup$AF_REF)
summary(myDataFormatNoDup$AF_ALT)

####### Filtering SNPs with high percentage of missing data
# total samples
N<-135
myDataFormatNoDup$NS<-as.numeric(as.character(myDataFormatNoDup$NS))
maxPerc<-80
table(100-myDataFormatNoDup$NS/N*100 > maxPerc)
myDataFormatNoDup<-myDataFormatNoDup[100-myDataFormatNoDup$NS/N*100 <= maxPerc,]
# Filter SNPs coming from IDs with more than 5 SNPs
myDataFormatNoDup<-myDataFormatNoDup[myDataFormatNoDup$ID %in% rownames(IDS_SNP[IDS_SNP$total <5,]),]

##### split sample columns GT:DP:AD:GL

sampleInfo<-(do.call(cbind, lapply(9:ncol(myDataFormatNoDup), function(x){
  sample<-strsplit(colnames(myDataFormatNoDup)[x], split="_")[[1]][2]
  sample<-strsplit(sample, split="[.]")[[1]][1]
  aux<-strsplit(as.character(myDataFormatNoDup[,x]), split=":")
  infoExp<-do.call(rbind, aux)
  colnames(infoExp)<-paste(sample, c("GT", "DP", "AD", "GL"), sep="_")
  return(infoExp)
})))

rownames(myDataFormatNoDup)<-1:nrow(myDataFormatNoDup)
myDataFormat<-myDataFormatNoDup[,1:8]
# construct genotypes table
charNA<-"./."
nofCol<-4 # number of columns in the info field
GT_table<-as.data.frame(do.call(rbind, lapply(1:nrow(sampleInfo), function(x){
  aux<-as.numeric(table(sampleInfo[x,seq(1, ncol(sampleInfo), by=nofCol)]))
  genot<-names(table(sampleInfo[x,seq(1, ncol(sampleInfo), by=nofCol)]))
  if(!(charNA %in% genot)){
    aux<-c(0, aux)
    genot<-c(charNA, genot)
  }
  return(cbind(SNP=x, genotypes=genot, freq=aux, perc=aux/sum(aux)*100))
  })))
levels(GT_table$genotypes)
# [1] "./." "0/0" "0/1" "1/0" "1/1"
levels(GT_table$genotypes)[levels(GT_table$genotypes) == charNA]<-"NA"
GT_table$perc<-round(as.numeric(as.character(GT_table$perc)),3)
# explore NA percentage distribution
summary(GT_table$perc[GT_table$genotypes == "NA"])
# number of complete SNPs
table(GT_table$freq[GT_table$genotypes=="NA"] == 0 )  
# genotype table of complete SNPs
GT_table_complete<-GT_table[GT_table$SNP %in% GT_table$SNP[GT_table$freq== 0 & GT_table$genotypes=="NA"],]
GT_table_complete<-GT_table_complete[GT_table_complete$genotypes !="NA",]
GT_table_complete$genotypes<-factor(as.character(GT_table_complete$genotypes))
# info matrix  of complete SNPs
sampleInfoComplete<-sampleInfo[as.numeric(as.character(unique(GT_table_complete$SNP))),]
# data matrix of complete SNPs
myDataFormatComplete<-myDataFormat[as.numeric(as.character(unique(GT_table_complete$SNP))),]

save.image("step1.RData", compress="xz")

# Set working directory
setwd("/path_to_directory/sim1")
# load the source's files
load("../step1.RData")
# define the SNPs with missing data
withNA<-as.character(GT_table$SNP[GT_table$genotypes=="NA" & !(GT_table$SNP %in% GT_table_complete$SNP)])
# compute NA percentages
percNAs<-100-myDataFormat$NS[myDataFormat$NS!= 105] *100/105
# extract genotypes
wholeGT<-as.data.frame(do.call(rbind, lapply(1:nrow(sampleInfo), function(x){
    aux<-as.character(sampleInfo[x,seq(1, ncol(sampleInfo), by=3)])
    return(aux)
})))
colnames(wholeGT)<-do.call(rbind,strsplit(colnames(sampleInfo)[seq(1, ncol(sampleInfo), by=3)], split="_"))[,1]
#genotypes of complete SNPs
wholeGT_complete<-wholeGT[unique(as.character(GT_table_complete$SNP)),]
wholeGT_complete<-as.data.frame(do.call(cbind, lapply(1:ncol(wholeGT_complete), function(x){
    column<-as.character(wholeGT_complete[,x])
    column[column=="0/0"]<-"0"
    column[column=="0/1"]<-"1"
    column[column =="1/0"]<-"1"
    column[column=="1/1"]<-"2"
    return(as.numeric(column))
    
})))
# redefine genotypes
wholeGT<-as.data.frame(do.call(cbind, lapply(1:ncol(wholeGT), function(x){
    column<-as.character(wholeGT[,x])
    column[column=="0/0"]<-"0"
    column[column=="0/1"]<-"1"
    column[column =="1/0"]<-"1"
    column[column=="1/1"]<-"2"
    column[column=="./."]<-"NA"
    return(as.numeric(column))
    
})))

# RF prediction
predDataAll<-predictSNP(SNPsCor=NULL, wholeGT, withNA, RF=TRUE, cor=FALSE, num.threads=4)
# Extract predicted genotypes
predRFAll<-as.data.frame(do.call(rbind, lapply(predDataAll, function(x){
    return(as.character(x[[1]]))
})))
# Extract associated OOB
oobAll<-do.call(c, lapply(predDataAll, function(x){
    return(as.numeric(as.character(x[[2]])))
}))
# RF imputation
RFImputed<-imputeSNP(sampleInfo, predictedData=predRFAll, myDataFormat, myDataFormatNoDup,OOB=FALSE, estimOOB = NULL)
# RFOOB imputation
RFOOBImputed<-imputeSNP(sampleInfo, predictedData=predRFAll, myDataFormat, myDataFormatNoDup,OOB=TRUE, estimOOB = oobAll)

# Correlated SNPs identification
SNPsCor<-corSNPs(withNA, wholeGT,SNPsData=NULL,LD=FALSE )
# RFCor prediction
predRFCorData<-predictSNP(SNPsCor, wholeGT,withNA, RF=TRUE, cor=TRUE)
# Extract predicted genotypes
predRFcor<-as.data.frame(do.call(rbind, lapply(predRFCorData, function(x){
    return(as.character(x[[1]]))
})))
# Extract associated OOB
oob<-do.call(c, lapply(predRFCorData, function(x){
    return(as.numeric(as.character(x[[2]])))
}))
# RFCor imputation
RFCorImputed<-imputeSNP(sampleInfo, predictedData=predRFcor, myDataFormatNoDup,ncolInfo=9, OOB=FALSE, estimOOB = NULL)
# RFCorOOB imputation
RFCorOOBImputed<-imputeSNP(sampleInfo, predictedData=predRFcor, myDataFormatNoDup,ncolInfo=9,OOB=TRUE, estimOOB = oob)

# Correlated SNPs with LD identification
SNPsCorLD<-corSNPs(withNA, wholeGT, SNPsData=myDataFormat,LD=TRUE )
# RFCorLD prediction
predDataLD<-predictSNP(SNPsCor=SNPsCorLD, wholeGT,withNA, RF=TRUE, cor=TRUE)
# Extract predicted genotypes
predRFcorLD<-as.data.frame(do.call(rbind, lapply(predDataLD, function(x){
    return(as.character(x[[1]]))
})))
# Extract associated OOB
oobLD<-do.call(c, lapply(predDataLD, function(x){
    return(as.numeric(as.character(x[[2]])))
}))
# RFCorLD imputation
RFCorLDImputed<-imputeSNP(sampleInfo, predictedData=predRFcorLD, myDataFormatNoDup,ncolInfo=9,OOB=FALSE, estimOOB = NULL)
# RFCorLDOOB imputation
RFCorLDOOBImputed<-imputeSNP(sampleInfo, predictedData=predRFcorLD, myDataFormatNoDup,ncolInfo=9,OOB=TRUE, estimOOB = oobLD)


# Mode imputation
predDataMed<-predictSNP(SNPsCor, wholeGT,withNA, RF=FALSE, cor=FALSE)
predMode<-do.call(rbind, predDataMed)
# RFCorLDOOB imputation
ModeImputed<-imputeSNP(sampleInfo, predictedData=predMode, myDataFormatNoDup,ncolInfo=9,OOB=FALSE, estimOOB = NULL)

save.image("step2.RData", compress="xz")

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
