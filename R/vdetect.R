vdetect <-
function(profiles,detect="centroid",plotit=TRUE){

    ############################################################################
    # (1) issue warnings #######################################################
    if(length(profiles[[1]])<2 || names(profiles[[1]])[1:2]!=c("m/z","abundance")){stop("WARNING: profile has invalid entries\n")}
    if(any(detect==c("centroid","intensoid","valley"))==FALSE){stop("WARNING: invalid detect argument")}
    if(plotit!="TRUE"&plotit!="FALSE"){stop("WARNING: plotit invalid. TRUE, FALSE.\n")}
    options(digits=10);
    ############################################################################
    # (2) detect ###############################################################
    cat(paste("\n Detect ",detect,"s... ",sep=""));
    getlist<-list(0);
    for(i in 1:length(profiles)){
        p_m<-profiles[[i]][[1]]   
        p_a<-profiles[[i]][[2]]
        if(detect=="centroid"){return_type<-0}  
        if(detect=="intensoid"){return_type<-1}
        if(detect=="valley"){return_type<-2}  
        # m:     double array of profile masses
        # a:     double array of profile abundances
        # return type: 0=centroid / 1=intensoid / 2=valley
      out <- .Call("iso_centroid_Call",
      f1 = as.double(p_m),
      a1 = as.double(p_a),
      t1 = as.integer(return_type),
	  PACKAGE="enviPat"
     )
       if(length(out[[1]])==0){
           getlist[[i]]<-"error"
           names(getlist)[i]<-names(profiles)[i]
        }else{
          out2<-data.frame(out[[1]],out[[2]])
          names(out2)<-c("m/z","abundance")
          getlist[[i]]<-out2
          names(getlist)[i]<-names(profiles)[i]
      if(plotit==TRUE){
        plot(p_m,p_a,type="h",xlab="m/z",ylab="Relative abundance",main=names(profiles)[i])
        if(detect=="valley"){
          points(out[[1]],out[[2]],col="red",type="h")
        }else{
          points(out[[1]],out[[2]],col="red",type="h",lwd=2);
        }
          }
    }
    }
    cat(" done.");
    ############################################################################
    # (3) output ###############################################################
    #if(detect=="centroid"){assign(detect,getlist);return(centroid);}
    #if(detect=="intensoid"){assign(detect,getlist);return(centroid);}
    #if(detect=="valley"){assign(detect,getlist);return(centroid);}
    return(getlist);
    ############################################################################

}
