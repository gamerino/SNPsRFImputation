library("BiocParallel")
library("ranger")
library("Rcpp")
# calculateAFALT
#Function to compute the frequency of the alternative allele.
#@param sampleData: data.frame containing sample information
#@param ngt: numeric indicating which field correspond to the genotype
#@param patNA: character indicating the missing genotype pattern
#@param phased: logical indicating if the genotype is phased. If it is TRUE, the
#allele separator is "|", else, the allele separator is "/"
calculateAFALT<-function(sampleData, ngt=1, phased=FALSE,patNA="./."){
    if (phased){allSep<-"|"}else{allSep<-"/"}
    GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(sampleData)), 
                               split="[:]"))[,1]
    GTTrue<-paste(do.call(rbind,strsplit(GT[GT!=patNA], allSep)))
    AAF<-round(length(which(GTTrue ==1))/length(GTTrue),3)
    return(AAF)
}
# getNS
# Function to obtain the number of samples with the SNP genotyped
# @param sampleData: data.frame containing samples' information
# @param patNA: character indicating the missing genotype pattern
getNS<-function(sampleData,patNA="./."){
    GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(sampleData)), 
                               split="[:]"))[,1]
    NS<-length(which(GT != patNA))
    return(NS)
}
# mergeSNPs
# Function to merge the information of SNPs identified in more than one tag
# @param dups: numeric vector indicating the position of duplicated CHROM:POS:ALT
# information in sampleData
# @param sampleData: data.frame containing the VCF indormation. The first 8 
# columns are CHROM, POS,ID,REF,ALT, NS AF_REF and AF_ALT. The next columns are
# the sample's data
#@param ncolInfo: number of columns of sampleData containing information about
#the SNPs. The next columnt should contain sample data
#@param patNA: character indicating the pattern of the missing genotype
mergeSNPs<-function(dups, sampleData,ncolInfo,patNA="./.",BPPARAM=bpparam()){
    if (any(!(c("CHROM", "POS", "ALT", "ID","NS","REF") %in% colnames(
        sampleData)))){
        stop("sampleData should have at least 'CHROM', 'POS', 'ALT', 'ID',
             'NS' and 'REF' columns")
    } 
    dupsSNP<-unique(sampleData[dups,c("CHROM","POS")])
    ncolSamples<-ncolInfo+1
    withoutRep<-as.data.frame(do.call(rbind, bplapply(1:nrow(dupsSNP), 
                                                      function(x){
                                                          repet<-sampleData[sampleData[,"CHROM"] == dupsSNP[x,"CHROM"] & 
                                                                                sampleData[,"POS"] == dupsSNP[x,"POS"] ,]
                                                          if(length(unique(repet[,"ALT"]))>1){
                                                              if (length(unique(repet[,"REF"])) > 1 & all(unique(
                                                                  repet[,"REF"]) %in% unique(repet[,"ALT"]))){
                                                                  IDMAX<-repet[which.max(repet$NS)[1],"ID"]
                                                                  REF<-repet[repet[,"ID"] ==IDMAX,"REF"]    
                                                                  ALT<-repet[repet[,"ID"] ==IDMAX,"ALT"]    
                                                                  repet[repet$REF != REF,"ALT"]<-ALT
                                                                  auxAFREF<-repet[repet$REF != REF,"AF_REF"]
                                                                  repet[repet$REF != REF,"AF_REF"]<-repet[
                                                                      repet$REF != REF,"AF_ALT"]
                                                                  repet[repet$REF != REF,"AF_ALT"]<-auxAFREF
                                                                  GT<-do.call(rbind,strsplit(as.character(
                                                                      as.matrix.data.frame(repet[repet$REF != REF,
                                                                                                 ncolSamples:ncol(repet)])), split="[:]"))[,1]
                                                                  GTcorr<-GT
                                                                  GTcorr[GT =="0/0"]<-"1/1"
                                                                  GTcorr[GT =="1/1"]<-"0/0"
                                                                  for(k in ncolSamples:ncol(repet)){
                                                                      repet[,k]<-as.character(repet[,k ])
                                                                  } 
                                                                  repet[repet$REF != REF,ncolSamples:ncol(repet) ]<-paste(
                                                                      GTcorr, do.call(rbind,strsplit(as.character(
                                                                          as.matrix.data.frame(repet[repet$REF != REF,
                                                                                                     ncolSamples:ncol(repet)])), split="[:]"))[,2],
                                                                      do.call(rbind,strsplit(as.character(
                                                                          as.matrix.data.frame(repet[repet$REF != REF,
                                                                                                     ncolSamples:ncol(repet)])), split="[:]"))[,3],
                                                                      do.call(rbind,strsplit(as.character(
                                                                          as.matrix.data.frame(repet[repet$REF != REF,
                                                                                                     ncolSamples:ncol(repet)])), split="[:]"))[,4], 
                                                                      sep=":")
                                                                  
                                                                  repet[repet$REF != REF,"REF"]<-REF
                                                                  
                                                              }
                                                              
                                                          }
                                                          IDMAX<-repet[which.max(as.numeric(as.character(repet[,"NS"])))[1], "ID"]
                                                          ref<-repet[repet[,"ID"]==IDMAX,]
                                                          alt<-repet[repet[,"ID"] !=IDMAX,]
                                                          ref[,"ID"]<-paste(sort(repet[,"ID"]), collapse="+")
                                                          for( k in 1:nrow(alt)){
                                                              concord<-NULL
                                                              for(j in ncolSamples:ncol(repet)){
                                                                  concord<-c(concord,length(unique(repet[,j]))==1)
                                                              }
                                                              if(any(!concord)){
                                                                  gts<-do.call(rbind, strsplit(as.matrix.data.frame(ref[,ncolInfo+
                                                                                                                            which(!concord),drop=FALSE]), ":"))[,1]
                                                                  gtsAlt<-do.call(rbind, strsplit(as.matrix.data.frame(alt[k,ncolInfo+
                                                                                                                               which(!concord),drop=FALSE]), ":"))[,1]
                                                                  if (any(gts==patNA)){
                                                                      NAref<-which(gts==patNA)
                                                                      altInNAret<-alt[k,ncolInfo+which(!concord)]
                                                                      toReplace<-which(gtsAlt[NAref] !=patNA)
                                                                      ref[,ncolInfo+which(!concord)][, NAref[toReplace]]<-alt[k,ncolInfo+
                                                                                                                                  which(!concord),drop=FALSE][,NAref[toReplace]]
                                                                  }
                                                                  if(any(gts!=patNA)){
                                                                      idx<-which(gts!=patNA)
                                                                      genotConRef<-gts[gts!=patNA]
                                                                      genotConAlt<-gtsAlt[gts!=patNA]
                                                                      if(any(genotConRef != genotConAlt & genotConAlt!="./." & 
                                                                             genotConAlt!="1/0" & genotConAlt!="0/1")){
                                                                          idxHetero<-which(genotConRef != genotConAlt & 
                                                                                               genotConAlt!="./." & genotConAlt!="1/0" & 
                                                                                               genotConAlt!="0/1")
                                                                          ref[,ncolInfo+which(!concord)][, which(gts!=patNA)][,
                                                                                                                              idxHetero]<-alt[k,ncolInfo+which(!concord),drop=FALSE][, 
                                                                                                                                                                                     which(gts!=patNA),
                                                                                                                                                                                     drop=FALSE][,idxHetero,drop=FALSE]
                                                                      }
                                                                  }
                                                              }
                                                          }
                                                          AAF<-calculateAFALT(ref[1,ncolSamples:ncol(ref)])
                                                          ref[,"AF_ALT"]<-AAF
                                                          ref[,"AF_REF"]<-1-AAF
                                                          GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(ref[1,
                                                                                                                           ncolSamples:ncol(ref)])), split="[:]"))[,1]
                                                          
                                                          ref[,"NS"]<-length(which(GT != "./."))
                                                          
                                                          return(ref)
                                                      }, BPPARAM=BPPARAM)))
    toRemove<-which(paste(sampleData[,"CHROM"] , sampleData[,"POS"] ,sep=":") 
                    %in% paste(withoutRep[,"CHROM"], withoutRep[,"POS"],
                               sep=":") )
    sampleDataNew<-sampleData[-toRemove,]
    sampleDataNew<-rbind(sampleDataNew,withoutRep)
    sampleDataNew[,"CHROM"]<-factor(as.character(sampleDataNew[,"CHROM"]), 
                                    levels=unique(sampleData[,"CHROM"]))
    sampleDataNew<-sampleDataNew[order(sampleDataNew[,"CHROM"], 
                                       sampleDataNew[,"POS"]),]
    return(sampleDataNew)
    }
# Mode
# Function to identify the SNP mode and impute missing genotypes with it
# @param x: character vector of detected SNPs
# @param method: character indicating the way to mode calculation
# @param n2imp: amount of SNPs to impute
Mode <- function(x, method = "one", n2imp=NULL) {
    x <- unlist(x)
    
    # Get unique values
    ux <- unique(x)
    n <- length(ux)
    
    # Get frequencies of all unique values
    frequencies <- tabulate(match(x, ux))
    modes <- frequencies == max(frequencies)
    
    # Determine number of modes
    nmodes <- sum(modes)
    nmodes <- ifelse(nmodes==n, 0L, nmodes)
    
    if (method %in% c("one", "mode", "") | is.na(method)) {
        # Return NA if not exactly one mode, else return the mode
        if (nmodes != 1) {
            return(sample(ux[which(modes)],n2imp, replace=T))
        } else {
            return(ux[which(modes)])
        }
    } else if (method %in% c("n", "nmodes")) {
        # Return the number of modes
        return(nmodes)
    } else if (method %in% c("all", "modes")) {
        # Return NA if no modes exist, else return all modes
        if (nmodes > 0) {
            return(ux[which(modes)])
        } else {
            return(NA)
        }
    }
    warning("Warning: method not recognised.  Valid methods are 'one'/'mode' [default], 'n'/'nmodes' and 'all'/'modes'")
}
#invertGT
# Function to reverse the genotipe orders when the reference allele is changed
# @param myDataFormatNoDup: data.frame containing the expanded VCF file without
# duplicated SNPs
# @param idxInv: numeric with the indexes of SNPs (rows of myDataFormatNoDup) 
# to invert
# @param ncolInfo: numeric indicating the last column of SNPs information 
# @param patNA: character indicating the pattern of the missing genotype
# @param BPPARAM: An optional BiocParallelParam instance defining the parallel
# back-end to be used during evaluation
invertGT<-function(myDataFormatNoDup,idxInv, ncolInfo=8, patNA="./.",
                   BPPARAM=bpparam()){
    infoCol<-myDataFormatNoDup[idxInv,1:ncolInfo]
    GTFinal<-as.data.frame(do.call(rbind,bplapply(1:length(idxInv),function(i){
        samplesGT<-myDataFormatNoDup[idxInv[i], (ncolInfo+1):ncol(
            myDataFormatNoDup)]
        GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(samplesGT)),
                                   split="[:]"))[,1]
        GTcorr<-GT
        GTcorr[GT =="0/0"]<-"1/1"
        GTcorr[GT =="1/1"]<-"0/0"
        for(k in 1:ncol(samplesGT)){
            samplesGT[,k]<-as.character(samplesGT[,k ])
        } 
        samplesGTInv<-paste(GTcorr, do.call(rbind,strsplit(as.character(
            as.matrix.data.frame(samplesGT)), split="[:]"))[,2], do.call(rbind,
                                                                         strsplit(as.character(as.matrix.data.frame(samplesGT)), 
                                                                                  split="[:]"))[,3],  sep=":")
        return(samplesGTInv)
        
        
    },BPPARAM=BPPARAM)))
    colnames(GTFinal)<-colnames(myDataFormatNoDup[, (ncolInfo+1):ncol(
        myDataFormatNoDup)])
    for(i in 1:ncol(GTFinal)){
        GTFinal[,i]<-as.character(GTFinal[,i])
    }
    auxREF<-infoCol[,"REF"]
    infoCol[,"REF"]<-infoCol[,"ALT"]
    infoCol[,"ALT"]<-auxREF
    auxAFREF<-infoCol[,"AF_REF"]
    infoCol[,"AF_REF"]<-infoCol[,"AF_ALT"]
    infoCol[,"AF_ALT"]<-auxAFREF
    
    return(cbind(infoCol, GTFinal))
}
# simulateNAs
# Function to simulate missing genotypes
# @param myDataFormatNoDup: data.frame containing the expanded VCF file without
# duplicated SNPs
# @param NS: numeric indicating the number of samples to simulate with missing
# genotypes
# @param ncolInfo: numeric indicating the last column of SNPs information 
# @param patNA: character indicating the pattern of the missing genotype
# @param BPPARAM: An optional BiocParallelParam instance defining the parallel
# back-end to be used during evaluation
# @param seed: numeric indicating the seed of the random selection
simulateNAs<-function(myDataFormatNoDup, NS,  ncolInfo=8, patNA="./.",
                      BPPARAM=bpparam(),seed=1234){
    set.seed(seed)
    GTFinal<-as.data.frame(do.call(rbind,bplapply(1:length(NS),function(i){
        samplesGT<-myDataFormatNoDup[i, (ncolInfo+1):ncol(myDataFormatNoDup)]
        samples2NA<-sample(1:length(samplesGT), NS[i])
        GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(samplesGT)
        ), split="[:]"))[,1]
        GT[samples2NA]<-patNA
        samplesGTNA<-paste(GT, do.call(rbind,strsplit(as.character(
            as.matrix.data.frame(samplesGT)), split="[:]"))[,2], do.call(rbind,
                                                                         strsplit(as.character(as.matrix.data.frame(samplesGT)), 
                                                                                  split="[:]"))[,3],  sep=":")
        return(samplesGTNA)
    }, BPPARAM=BPPARAM)))
    colnames(GTFinal)<-colnames(myDataFormatNoDup)[(ncolInfo+1):ncol(
        myDataFormatNoDup)]
    GTFinal<-cbind(myDataFormatNoDup[1:ncolInfo], GTFinal)
    return(GTFinal)
}
# corSNPs
# Function to determine which of complete SNPs are correlated with each of the 
# incomplete SNPs
# @param withNA: character vector containing the name of incomplete SNPs
# @param wholeGT: genotype matrix of incomplete and complete SNPs
# @param wholeGT_complete: genotype matrix of complete SNPs
# @param myDataFormat: information matrix of incomplete and complete SNPs
# @param myDataFormatComplete: information matrix of complete SNPs
# @param LD: logical indicating if linkage disequilibrium should be considered
# for SNPs correlation
# @param pval: numeric indicating the p-value used in the independence test
corSNPs<-function(withNA, wholeGT, wholeGT_complete,myDataFormat, 
                  myDataFormatComplete, LD=TRUE, pval=0.05){
    highCorr<-bplapply(1:length(withNA),function(x){
        snpToImpute<-withNA[x]
        snpRead<-wholeGT[snpToImpute,]
        detected<-snpRead[!is.na(snpRead)]
        miss<-snpRead[is.na(snpRead)]
        if(LD==FALSE){
            dfw<-as.data.frame(t(wholeGT_complete[, !is.na(snpRead)]))
            colnames(dfw)<-paste("SNP", colnames(dfw), sep="")
            pvalIndep<-callmychi(detected,as.matrix(dfw))
            names(pvalIndep)<-colnames(dfw)
            if( length(which(pvalIndep< pval)) >=1){
                SNPsCorr<-names(which(pvalIndep< pval))
            }else{
                SNPsCorr<-NA
            }
        }else{
            chromOK<-myDataFormat[snpToImpute, "CHROM"]
            if(any(myDataFormatComplete$CHROM==chromOK)){
                dfw<-as.data.frame(t(wholeGT_complete[, !is.na(snpRead)]))
                colnames(dfw)<-paste("SNP", colnames(dfw), sep="")
                dfw<-dfw[,myDataFormatComplete$CHROM==chromOK, drop=FALSE]
                pvalIndep<-callmychi(detected,as.matrix(dfw))
                names(pvalIndep)<-colnames(dfw)
                if( length(which(pvalIndep< pval)) >=1){
                    SNPsCorr<-names(which(pvalIndep< pval))
                }else{
                    SNPsCorr<-NA
                }
            }else{
                SNPsCorr<-NA
            }
        }
        return(SNPsCorr)
    }, BPPARAM=MulticoreParam(4))
    return(highCorr)
}
# predictSNP
# Function to predict the missing genotyps of incomplete SNPs
# @param SNPsCor: list containing the correlated 'complete SNPs' for each
# 'incomplete SNP'
# @param withNA: character vector containing the name of incomplete SNPs
# @param wholeGT: genotype matrix of incomplete and complete SNPs
# @param wholeGT_complete: genotype matrix of complete SNPs
# @param traspWholeGT_complete: transpose wholeGT_complete matrix
# @param RF: logical indicating if predictios should be based on random forest.
# If it was set to FALSE, mode imputation will be used.
# @param cor: logical indicating if SNPs correlation should be considered
# @param num.threads: number of threads used to random forest construction
predictSNP<-function(SNPsCor, withNA, wholeGT, wholeGT_complete, 
                     traspWholeGT_complete, RF=TRUE, cor=TRUE, num.threads=4){
    imputedData<-lapply(1:length(withNA), function(x){
        snpToImpute<-withNA[x]
        snpRead<-wholeGT[snpToImpute,]
        detected<-snpRead[!is.na(snpRead)]
        if(RF){
            miss<-snpRead[is.na(snpRead)]
            df<-data.frame(detected=factor(as.character(detected), 
                                           levels=c("0","1","2")))
            dfw<-cbind(df, traspWholeGT_complete[!is.na(snpRead), ])
            colnames(dfw)[-1]<-paste("SNP", colnames(dfw)[-1], sep="")
            
            if(length(unique(detected))>1){
                if(cor){
                    if(!any(is.na(SNPsCor[[x]]))){
                        dfCor<-dfw[,c("detected", SNPsCor[[x]])]
                    }else{dfCor<-dfw
                    }
                    rf<-ranger(formula(dfCor), dfCor,num.trees=500,
                               num.threads = num.threads)
                }else{
                    rf<-ranger(formula(dfw), dfw,num.trees=500, 
                               num.threads = num.threads)
                }
                OOBER<-rf$prediction.error
                dfM<-data.frame(missing=factor(as.character(miss), levels=c(
                    "NA","0","1","2")))
                dfwM<-cbind(dfM, traspWholeGT_complete[is.na(snpRead),])
                colnames(dfwM)[-1]<-paste("SNP", colnames(dfwM)[-1], sep="")
                if(cor){
                    if(!any(is.na(SNPsCor[[x]]))){
                        dfCorM<-dfwM[,SNPsCor[[x]], drop=FALSE] 
                    }else{dfCorM<-dfwM
                    }
                    imputed<-predictions(predict(rf, dfCorM))
                    
                }else{
                    
                    imputed<-predictions(predict(rf, dfwM))
                }
         }else{
                imputed<-as.character(unique(dfw$detected))
                OOBER<-NA
            }
            snpReadImputed<-do.call(c,lapply(snpRead[1,,drop=T], as.character))
            
            snpReadImputed[is.na(snpReadImputed)]<-as.character(imputed)
            
            snpReadImputed<-factor(as.character(snpReadImputed), levels=c("0","1","2"))
            
            return(list(imp=as.character(snpReadImputed), error=OOBER))
            
        }else{
            imputed<-Mode(detected, n2imp=length(snpRead[is.na(snpRead)]))
            snpReadImputed<-snpRead
            snpReadImputed[is.na(snpRead)]<-imputed
            return(snpReadImputed)
            
        }
    })    
    return(imputedData)
    
}
# imputeSNP
# Method to impute the missing genotypes with predicted values obtained using the
# predictSNP method.
# @param sampleInfo: data.frame containing the sample information columns from 
# the VCF file
# @param myDataFormat: information matrix of incomplete and complete SNPs
# @param myDataFormatNoDup: data.frame containing the expanded VCF file without
# duplicated SNPs
# @param OOB: logical indicating if OOB estimation should be used
# @param estimOOB: numeric vector containing OOB estimations
# @param OOBThres: numeric indicating the OOB Threshold to be used in order to 
# decide if the incomplete SNP should be imputed
imputeSNP<-function(sampleInfo, predictedData,myDataFormat,myDataFormatNoDup,
    OOB=TRUE,estimOOB, OOBThres=0.2){
    sampleInfoImp<-sampleInfo
    imputedMedGT<-as.data.frame(do.call(cbind, lapply(1:ncol(predictedData), 
        function(x){
        aux<-predictedData[,x]
        aux[aux==0]<-"0/0"
        aux[aux==1]<-"1/0"
        aux[aux==2]<-"1/1"
        return(aux)
    })))
    if(OOB){
        idxOOB<-which(is.na(estimOOB) | estimOOB < 0.2)
        
        sampleInfoImp[as.numeric(withNA),seq(1, ncol(sampleInfoImp), by=3)][
            idxOOB,]<-as.matrix(imputedMedGT[idxOOB,])
    }else{    
        sampleInfoImp[as.numeric(withNA),seq(1, ncol(sampleInfoImp), by=3)]<-
            as.matrix(imputedMedGT)
    }
    myDataFormatImputedM<-myDataFormat
    
    mdfIM<-(do.call(cbind, lapply(seq(1,ncol(sampleInfoImp),by=3), function(x){
        infoCompres<-paste0(sampleInfoImp[,x],":", sampleInfoImp[,x+1],":", 
                            sampleInfoImp[,x+2], sep="")
        return(infoCompres)
    })))
    colnames(mdfIM)<-colnames(myDataFormatNoDup[,9:ncol(myDataFormatNoDup)])
    myDataFormatImputedM<-cbind(myDataFormatImputedM, mdfIM)
    # recompute NS, AF_REF y AF_ALT
    AF_ALT_NS<-do.call(rbind, bplapply(1:nrow(myDataFormatImputedM), function(x){
        AF_ALT<-calculateAFALT(myDataFormatImputedM[x,9:ncol(myDataFormatImputedM)])
        NS<-getNS(myDataFormatImputedM[x,9:ncol(myDataFormatImputedM)])
        return(cbind(AF_ALT=AF_ALT,NS=NS))
    }, BPPARAM = MulticoreParam(4)))
    
    myDataFormatImputedM$AF_ALT<-AF_ALT_NS[,"AF_ALT"]
    myDataFormatImputedM$AF_REF<-1-AF_ALT_NS[,"AF_ALT"]
    myDataFormatImputedM$NS<-AF_ALT_NS[,"NS"]
    
    myDataImputed<-myDataFormatImputedM[,1:5]
    
    mdIM<-paste0("NS=", myDataFormatImputedM[,"NS"],";AF=", 
                 myDataFormatImputedM[,"AF_REF"],",", 
                 myDataFormatImputedM[,"AF_ALT"], sep="")
    
    myDataImputed<-cbind(myDataImputed,QUAL=".", FILTER="PASS", INFO=mdIM, 
                         FORMAT="GT:DP:AD:GL", myDataFormatImputedM[,9:ncol(
                             myDataFormatImputedM)])
    return(myDataImputed)
}
