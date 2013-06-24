isopattern <-
function(
  isotopes,
  chemforms,
  threshold=.001,
  charge=FALSE,
  emass=0.00054858,
  plotit=FALSE,
  algo=2
  ){

    ############################################################################
    # (1) issue warnings #######################################################
    if(length(isotopes)!=4){stop("WARNING: invalid isotope list\n")}
    if(threshold>100 || threshold<0){stop("WARNING: invalid threshold; 0<=threshold<100.\n")}  
    if(plotit!="TRUE"&plotit!="FALSE"){stop("WARNING: plotit invalid. TRUE, FALSE.\n")}
    if(emass!=0.00054858){cat("NOTE: You are sure that is the mass of an electrone?")}
    if(length(charge)!=length(chemforms) & length(charge)>1){stop("length of charge does not match number of chemforms!\n")}
    if(charge==0 & charge!=FALSE){stop("WARNING: charge=0?")}
    if(is.numeric(charge)==FALSE & charge!=FALSE){stop("WARNING: charge either numeric or FALSE!")}
    if(algo!=1 & algo!=2){stop("invalid algo argument!")}
    options(digits=10);
    ############################################################################
    # (2) parse isolist, set charge ############################################
    cat("\n Calculate isotope patterns ...");
    isolist<-""
    for(i in 1:length(isotopes[,1])){
      getit<-as.character(i)
      for(j in 1:4){
        getit<-paste(getit,"  ",as.character(isotopes[i,j]),sep="")
      }
      isolist<-paste(isolist,getit,"$",sep="");
    } 
    substr(isolist,nchar(isolist),nchar(isolist))<-"@"
    if(length(charge)!=length(chemforms)){
      charge<-rep(charge,length(chemforms))
    }
    # (3) run isotope pattern generator ######################################## 
    pattern<-list(0)
    for(i in 1:length(chemforms)){
      # sum:   chemical formula
      # threshold: peaks with abundance < threshold are neglegted. 
      #              value is relative to the maximum abundance
      if(algo==1){
        out <- .Call( "iso_pattern_Call",
      s1 = as.character(chemforms[i]),   # chemical formula
      pl = as.integer(1E6),             # number of peaks to be reserved for
      t1 = as.double(threshold),        # relative intensity cutoff
      iso_list1 = as.character(isolist) # parsed isotope list
    )
      }else{
        out <- .Call( "iso_pattern_Call_2",
      s1 = as.character(chemforms[i]),   # chemical formula
      pl = as.integer(1E6),             # number of peaks to be reserved for
      t1 = as.double(threshold),        # relative intensity cutoff
      iso_list1 = as.character(isolist) # parsed isotope list
    )
  }
  # parse output ###########################################################
  if(length(out[[1]])==0){
         pattern[[i]]<-"error"
         names(pattern)[i]<-chemforms[i]
      }else{
        out2<-as.data.frame(out[[1]])
        for(j in 2:(length(out)-1)){
      out2<-cbind(out2,out[[j]])
    }
    out2<-out2[order(out2[,1],decreasing=FALSE),]
        names(out2)<-out[[length(out)]][-length(out[[length(out)]])]
        names(out2)[1]<-"m/z"
        if(charge[i]!=FALSE){
          out2[,1]<-c(out2[,1]-(charge[i]*emass));  # electrone mass
          out2[,1]<-c(out2[,1]/abs(charge[i]));  # /charge=z
        }
        pattern[[i]]<-out2
    names(pattern)[i]<-as.character(chemforms[i]);
    if(plotit==TRUE){
        plot(out2[,1],out2[,2],type="h",
              xlab="m/z",ylab="Relative abundance",main=names(pattern)[i])
    }
  }
}
cat(" done.");
    ############################################################################
    # (3) output ###############################################################
    return(pattern) 
    ############################################################################

}
