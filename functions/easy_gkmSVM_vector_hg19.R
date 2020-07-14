easy_gkmSVM_vector_hg19 <- function(TS,Trail,xfold){
  #TS=training set in bed format
  #Trial=times of trial finding negative control
  #xfold=number of fold
  POS_fa<-paste(TS,".pos.fa",sep="")
  NEG_bed<-paste(TS,".neg.bed",sep="")
  NEG_fa<-paste(TS,".neg.fa",sep="")
  genNullSeqs(TS,nMaxTrials=20,xfold=1,
              genome=BSgenome.Hsapiens.UCSC.hg19.masked,   
              outputPosFastaFN=POS_fa, 
              outputBedFN=NEG_bed, 
              outputNegFastaFN=NEG_fa)
  kernel_out<-paste(TS,".",xfold,"x.out",sep="")
  gkmsvm_kernel(POS_fa,
                NEG_fa, 
                kernel_out,
                K=6,
                L=10)
  
  cvres=gkmsvm_trainCV(kernel_out,
                       POS_fa,
                       NEG_fa,
                       svmfnprfx=TS,
                       nCV=5,
                       outputCVpredfn=paste(TS,".",xfold,"x.cvpred.out",sep=""), 
                       outputROCfn=paste(TS,".",xfold,"x.ROC.out",sep=""),
                       showPlots = TRUE,
                       outputPDFfn = paste(TS,".",xfold,"x.CURVE.pdf",sep=""))
  

}


auPRC <- function(perf) {
  rec <- perf@x.values
  prec <- perf@y.values
  result <- list()
  for (i in 1:length(perf@x.values)) {
    result[i] <- list(sum((rec[[i]][2:length(rec[[i]])] - 
                             rec[[i]][2:length(rec[[i]]) - 1]) * prec[[i]][-1]))
  }
  return(result)
}
