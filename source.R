library("BiocParallel")
library("ranger")
library("Rcpp")
sourceCpp("source.cpp")
#######################################################
## Computing the frequency of the alternative allele ##
#######################################################
# Parameters:
# sampleData, df containing sample information
# ngt, numeric indicating which field correspond to the genotype
# patNA, character indicating the missing genotype patterns
# phased, logical indicating if the genotype is phased. If it is TRUE, the 
# allele separator is "|", else, the allele separator is "/"
calculateAFALT<-function(sampleData, ngt=1, phased=FALSE,patNA="./."){
    if (phased){allSep<-"|"}else{allSep<-"/"}
    GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(sampleData)), 
    split="[:]"))[,1]
  GTTrue<-paste(do.call(rbind,strsplit(GT[GT!=patNA], allSep)))
  
  AAF<-round(length(which(GTTrue ==1))/length(GTTrue),3)
  return(AAF)
}

###############################################
## Obtaining the number of genotyped samples ##
###############################################

# Parameters
# sampleData, data.frame  containing samples' information
# patNA, character indicating the pattern used to refer to missing genotype pattern

getNS<-function(sampleData,patNA="./."){
  GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(sampleData)), 
                             split="[:]"))[,1]
  NS<-length(which(GT != patNA))

  return(NS)
}  

#############################
## Merging duplicated SNPs ##
#############################
# Parameters:
# dups, numeric vector indicating the position of duplicated CHROM:POS:ALT 
# in sampleData
# sampleData, data.frame containing the VCF information. The first 8 columns are CHROM,
# POS,ID,REF,ALT, NS AF_REF and AF_ALT. The next columns are the sample's data
# ncolInfo: number of columns of sampleData containing information about 
# the SNPs. The next columnt should contain sample data
# patNA: character indicating the pattern used to refer to missing genotype pattern

mergeSNPs<-function(dups, sampleData,ncolInfo,patNA="./.",BPPARAM=bpparam()){
    if (any(!(c("CHROM", "POS", "ALT", "ID","NS","REF") %in% colnames(
        sampleData)))){stop("sampleData should have at least 'CHROM', 'POS', 
                            'ALT', 'ID','NS' and 'REF' columns")} 
    # dupsSNP<-unique(sampleData[dups,c("CHROM","POS", "ALT")])
    dupsSNP<-unique(sampleData[dups,c("CHROM","POS")])
    #aca tengo que evaluar que pasa si son iguales los REF y que pasa sin son distintos
    ncolSamples<-ncolInfo+1
    withoutRep<-as.data.frame(do.call(rbind, bplapply(1:nrow(dupsSNP), function(x){
        repet<-sampleData[sampleData[,"CHROM"] == dupsSNP[x,"CHROM"] & 
                              sampleData[,"POS"] == dupsSNP[x,"POS"] ,]
        
        if(length(unique(repet[,"ALT"]))>1){
            
            if (length(unique(repet[,"REF"])) > 1 & all(unique(repet[,"REF"]) %in% 
                                                        unique(repet[,"ALT"]))){
                #Tengo SNPs que tienen cruzada la info de ref y alt, tomo el más
                #frecuente como de referencia
                IDMAX<-repet[which.max(repet$NS)[1],"ID"]
                REF<-repet[repet[,"ID"] ==IDMAX,"REF"]    
                ALT<-repet[repet[,"ID"] ==IDMAX,"ALT"]    
                repet[repet$REF != REF,"ALT"]<-ALT
                auxAFREF<-repet[repet$REF != REF,"AF_REF"]
                repet[repet$REF != REF,"AF_REF"]<-repet[repet$REF != REF,"AF_ALT"]
                repet[repet$REF != REF,"AF_ALT"]<-auxAFREF
                GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(
                    repet[repet$REF != REF,ncolSamples:ncol(repet)])), split="[:]"))[,1]
                GTcorr<-GT
                GTcorr[GT =="0/0"]<-"1/1"
                GTcorr[GT =="1/1"]<-"0/0"
                for(k in ncolSamples:ncol(repet)){
                    repet[,k]<-as.character(repet[,k ])
                } 
                readAll<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(
                    repet[repet$REF != REF,ncolSamples:ncol(repet)])), split="[:]"))[,3]
                readAllInv<-do.call(c, lapply(readAll, function(l){
                    return(paste(strsplit(l, split="[,]")[[1]][2],strsplit(l, split="[,]")[[1]][1], sep=","))
                }))
                if(length(strsplit(as.character(as.matrix.data.frame(repet[
                    repet$REF != REF,ncolSamples:ncol(repet)])), split="[:]")[[1]])==4){
                    repet[repet$REF != REF,ncolSamples:ncol(repet) ]<-paste(GTcorr, do.call(rbind,strsplit(as.character(as.matrix.data.frame(
                    repet[repet$REF != REF,ncolSamples:ncol(repet)])), split="[:]"))[,2],
                    readAllInv,
                    do.call(rbind,strsplit(as.character(as.matrix.data.frame(
                        repet[repet$REF != REF,ncolSamples:ncol(repet)])), split="[:]"))[,4], sep=":")
                }else{
                    repet[repet$REF != REF,ncolSamples:ncol(repet) ]<-paste(GTcorr, do.call(rbind,strsplit(as.character(as.matrix.data.frame(
                        repet[repet$REF != REF,ncolSamples:ncol(repet)])), split="[:]"))[,2],
                        readAllInv,
                        do.call(rbind,strsplit(as.character(as.matrix.data.frame(
                            repet[repet$REF != REF,ncolSamples:ncol(repet)])), split="[:]"))[,4], sep=":")
                    
                }
                repet[repet$REF != REF,"REF"]<-REF
                
            }
            
        }
        IDMAX<-repet[which.max(as.numeric(as.character(repet[,"NS"])))[1], "ID"]
        ref<-repet[repet[,"ID"]==IDMAX,]
        alt<-repet[repet[,"ID"] !=IDMAX,]
        ref[,"ID"]<-paste(sort(repet[,"ID"]), collapse="+")
        for( k in 1:nrow(alt)){
            concord<-NULL
            for(j in ncolSamples:ncol(repet)){concord<-c(concord,length(unique(
                repet[,j]))==1)}
            if(any(!concord)){
                #debo ver de los que no tengo datos iguales
                # primero me fijo si alguno esNA en rel
                gts<-do.call(rbind, strsplit(as.matrix.data.frame(ref[,ncolInfo+
                    which(!concord),drop=FALSE]), ":"))[,1]
                gtsAlt<-do.call(rbind, strsplit(as.matrix.data.frame(
                    alt[k,ncolInfo+ which(!concord),drop=FALSE]), ":"))[,1]
                
                if (any(gts==patNA)){
                    NAref<-which(gts==patNA)
                    altInNAret<-alt[k,ncolInfo+which(!concord)]
                    toReplace<-which(gtsAlt[NAref] !=patNA)
                    ref[,ncolInfo+which(!concord)][, NAref[toReplace]]<-alt[k,ncolInfo+
                                                                                which(!concord),drop=FALSE][,NAref[toReplace]]
                }
                if(any(gts!=patNA)){
                    # primer caso posible, mismo genotipo (o ref homocigoto y alt heteroc)
                    # y diferente valor de los restantes campos, me quedo con elref
                    idx<-which(gts!=patNA)
                    genotConRef<-gts[gts!=patNA]
                    genotConAlt<-gtsAlt[gts!=patNA]
                    rdRef<-as.numeric(do.call(rbind, strsplit(as.matrix.data.frame(ref[,ncolInfo+
                            which(!concord),drop=FALSE]), ":"))[,2])[gts!=patNA]
                    rdAlt<-as.numeric(do.call(rbind, strsplit(as.matrix.data.frame(
                        alt[k,ncolInfo+ which(!concord),drop=FALSE]), ":"))[,2])[gts!=patNA]
                    
                    if (any((genotConRef != genotConAlt | rdRef<rdAlt) & genotConAlt!="./." )){
                        
                        idxHetero<-which(genotConRef != genotConAlt & 
                                             genotConAlt!="./." & rdRef<rdAlt)
                        ref[,ncolInfo+which(!concord)][, which(gts!=patNA)][,
                         idxHetero]<-alt[k,ncolInfo+which(!concord),drop=FALSE][, 
                            which(gts!=patNA),drop=FALSE][,idxHetero,drop=FALSE]
                        idxHomo<-which(genotConAlt!="./." & genotConRef == genotConAlt  & rdRef<rdAlt)
                        ref[,ncolInfo+which(!concord)][, which(gts!=patNA)][,idxHomo]<-alt[k,
                            ncolInfo+which(!concord),drop=FALSE][, 
                            which(gts!=patNA),drop=FALSE][,idxHomo,drop=FALSE]
                        
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

Mode <- function(x, method = "one", na.rm = FALSE, n2imp=NULL) {
  x <- unlist(x)
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  
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
##################################
## exchanging genotypes of SNPs ##
##################################
# Parameters:
# idxInv, numeric vector indicating the position of the SNPs which reference and
# alternative alleles information must be exchanged 
# myDataFormatNoDup, data.frame containing the VCF information. The first 8 columns are CHROM,
# POS,ID,REF,ALT, NS AF_REF and AF_ALT. The next columns are the sample's data
# ncolInfo: number of columns of sampleData containing information about 
# the SNPs. The next columnt should contain sample data
# patNA: character indicating the pattern used to refer to missing genotype pattern

exchangeGT<-function(myDataFormatNoDup,idxInv, ncolInfo=8, patNA="./.",BPPARAM=bpparam()){
    infoCol<-myDataFormatNoDup[idxInv,1:ncolInfo]
    GTFinal<-as.data.frame(do.call(rbind,bplapply(1:length(idxInv),function(i){
        samplesGT<-myDataFormatNoDup[idxInv[i], (ncolInfo+1):ncol(myDataFormatNoDup)]
        GT<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(samplesGT)), split="[:]"))[,1]
        GTcorr<-GT
        GTcorr[GT =="0/0"]<-"1/1"
        GTcorr[GT =="1/1"]<-"0/0"
        for(k in 1:ncol(samplesGT)){
            samplesGT[,k]<-as.character(samplesGT[,k ])
        } 
        readAll<-do.call(rbind,strsplit(as.character(as.matrix.data.frame(
            samplesGT[,])), split="[:]"))[,3]
        readAll[which(GT!= patNA)]<-do.call(c, lapply(readAll[which(GT!= patNA)], function(l){
            return(paste(strsplit(l, split="[,]")[[1]][2],strsplit(l, split="[,]")[[1]][1], sep=","))
        }))
        
        samplesGTInv<-paste(GTcorr, do.call(rbind,strsplit(as.character(as.matrix.data.frame(
            samplesGT)), split="[:]"))[,2], 
            readAll, sep=":")
        if(length(strsplit(as.character(samplesGT), split=":")[[1]])==4){
            samplesGTInv<-paste(samplesGTInv, do.call(rbind,strsplit(as.character(as.matrix.data.frame(
                samplesGT)), split="[:]"))[,4],  sep=":")
        }
        return(samplesGTInv)
        
        
    },BPPARAM=BPPARAM)))
    colnames(GTFinal)<-colnames(myDataFormatNoDup[, (ncolInfo+1):ncol(myDataFormatNoDup)])
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

#################################
## Identifying correleted SNPs ##
#################################
# Function to determine which of the complete SNPs are correlated with each of 
# the incompletely genotyped SNPs
# Parameters:
# withNA, character vector containing the name of incomplete SNPs
# wholeGT, genotype matrix of incomplete and complete SNPs
# wholeGT_complete, genotype matrix of complete SNPs
# SNPsData, information matrix of incomplete and complete SNPs
# myDataFormatComplete, information matrix of complete SNPs
# LD, logical indicating if linkage disequilibrium should be considered
# for SNPs correlation
# pval, numeric indicating the significance threshold used in the independence 
# test
corSNPs<-function(withNA, wholeGT,SNPsData, LD=TRUE, pval=0.05, BPPARAM=bpparam()){
    wholeGT_complete<-wholeGT[-as.numeric(withNA),]
    completeSNPsData<-SNPsData[-as.numeric(withNA),]
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
            chromOK<-SNPsData[snpToImpute, "CHROM"]
            if(any(completeSNPsData$CHROM==chromOK)){
                dfw<-as.data.frame(t(wholeGT_complete[, !is.na(snpRead)]))
                colnames(dfw)<-paste("SNP", colnames(dfw), sep="")
                dfw<-dfw[,completeSNPsData$CHROM==chromOK, drop=FALSE]
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
    }, BPPARAM=BPPARAM)
    return(highCorr)
}
# predictSNP
# Function to predict the missing genotyps of incomplete SNPs
# @param wholeGT: genotype matrix of incomplete and complete SNPs
# @param wholeGT_complete: genotype matrix of complete SNPs
# @param withNA: character vector containing the name of incomplete SNPs
# @param RF: logical indicating if predictios should be based on random forest.
# If it was set to FALSE, mode imputation will be used.
# @param cor: logical indicating if SNPs correlation should be considered
# @param SNPsCor: list containing the correlated 'complete SNPs' for each
# 'incomplete SNP'
# @param num.threads: number of threads used to random forest construction
predictSNP<-function(SNPsCor, wholeGT, withNA, 
                     RF=TRUE, cor=TRUE, num.threads=4){
    # Traspose the complete SNPs' genotypes
    wholeGT_complete<-wholeGT[-as.numeric(withNA),]
    traspWholeGT_complete<-as.data.frame(t(wholeGT_complete))
    for(j in 1:ncol(traspWholeGT_complete)){
        traspWholeGT_complete[,j]<-factor(as.character(traspWholeGT_complete[,j]), levels=c("0", "1", "2"))
    }
    imputedData<-lapply(1:length(withNA), function(x){
        snpToImpute<-withNA[x]
        snpRead<-wholeGT[snpToImpute,]
        detected<-snpRead[!is.na(snpRead)]
        if(RF){
            miss<-snpRead[is.na(snpRead)]
            df<-data.frame(detected=factor(as.character(detected), 
                                           levels=c("0","1","2")))
            dfw<-cbind(df, traspWholeGT_complete[!is.na(snpRead), , drop=FALSE])
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
                dfwM<-cbind(dfM, traspWholeGT_complete[is.na(snpRead),,drop=FALSE])
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
# @param SNPsData: information matrix of incomplete and complete SNPs
# @param OOB: logical indicating if OOB estimation should be used
# @param estimOOB: numeric vector containing OOB estimations
# @param OOBThres: numeric indicating the OOB Threshold to be used in order to 
# decide if the incomplete SNP should be imputed
imputeSNP<-function(sampleInfo, predictedData,myDataFormatNoDup, ncolInfo, 
                    OOB=TRUE,estimOOB, OOBThres=0.2, INFOcols=c("NS", "AF_REF", "AF_ALT"), 
                    nfields=4, BPPARAM=bpparam()){
    sampleInfoImp<-sampleInfo
    SNPsData<-myDataFormatNoDup[,1:(ncolInfo-1)]
    imputedMedGT<-as.data.frame(do.call(cbind, lapply(1:ncol(predictedData), 
            function(x){
                aux<-as.character(predictedData[,x])
                aux[aux==0]<-"0/0"
                aux[aux==1]<-"1/0"
                aux[aux==2]<-"1/1"
                return(aux)
                })))
    if(OOB){
        idxOOB<-which(is.na(estimOOB) | estimOOB < OOBThres)
        
        sampleInfoImp[as.numeric(withNA),seq(1, ncol(sampleInfoImp), by=nfields)][
            idxOOB,]<-as.matrix(imputedMedGT[idxOOB,])
    }else{    
        sampleInfoImp[as.numeric(withNA),seq(1, ncol(sampleInfoImp), by=nfields)]<-
            as.matrix(imputedMedGT)
    }
    myDataFormatImputedM<-SNPsData
    len<-seq(1,ncol(sampleInfoImp),by=nfields)
    mdfIM<-(do.call(cbind, lapply(1:length(len), function(x){
        infoCompres<-paste0(sampleInfoImp[,len[x]], sep="")
        for(j in 1:(nfields-1)){
            infoCompres<-paste0(infoCompres,":",sampleInfoImp[,len[x]+j], sep="")
        }
        return(infoCompres)
    })))

    colnames(mdfIM)<-colnames(myDataFormatNoDup[,ncolInfo:ncol(myDataFormatNoDup)])
    myDataFormatImputedM<-cbind(myDataFormatImputedM, mdfIM)
    # recompute NS, AF_REF y AF_ALT
    AF_ALT_NS<-do.call(rbind, bplapply(1:nrow(myDataFormatImputedM), function(x){
        AF_ALT<-calculateAFALT(myDataFormatImputedM[x,ncolInfo:ncol(myDataFormatImputedM)])
        NS<-getNS(myDataFormatImputedM[x,ncolInfo:ncol(myDataFormatImputedM)])
        return(cbind(AF_ALT=AF_ALT,NS=NS))
    }, BPPARAM = BPPARAM))
    
    myDataFormatImputedM$AF_ALT<-AF_ALT_NS[,"AF_ALT"]
    myDataFormatImputedM$AF_REF<-1-AF_ALT_NS[,"AF_ALT"]
    myDataFormatImputedM$NS<-AF_ALT_NS[,"NS"]
    
    myDataImputed<-myDataFormatImputedM[,1:5]
    if("AF_REF" %in% INFOcols){
        AF=paste0("AF=",myDataFormatImputedM[,"AF_REF"],",",
                 myDataFormatImputedM[,"AF_ALT"],sep="")
    }else{
        AF=paste0("AF=",myDataFormatImputedM[,"AF_ALT"],sep="")
    }
    if("LOCORI" %in% INFOcols){
        LOCORI=paste0(";locori=",myDataFormatImputedM[,"LOCORI"],sep="")
        }else{LOCORI=NULL}
    mdIM<-paste0("NS=", myDataFormatImputedM[,"NS"],";",AF,LOCORI, sep="")
    if(nfields==4){
        FORMAT="GT:DP:AD:GL"
    }else{
        FORMAT="GT:DP:AD"
    }
    myDataImputed<-cbind(myDataImputed,QUAL=".", FILTER="PASS", INFO=mdIM, 
                         FORMAT, myDataFormatImputedM[,ncolInfo:ncol(
                             myDataFormatImputedM)])
    return(myDataImputed)
}
