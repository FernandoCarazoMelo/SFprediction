sacartranscritos <- function(edgetr,events){
  
  
  p1 <- events$Path1
  p2 <- events$Path2
  pref <- events$PathR
  
  
  p1<-sapply(strsplit(p1,","),function(X){
    n<-length(X)
    trans <- ""
    for (i in 1:n){
      ss <- strsplit(X[i],"-")[[1]][1]
      ee <- strsplit(X[i],"-")[[1]][2]
      if (ss == ee){
        ss <- paste0(ss,".a")
        ee <- paste0(ee,".b")
      }else{
        ss <- paste0(ss,".b")
        ee <- paste0(ee,".a")
      }
      if (i==1){
        trans<-paste0(trans,edgetr$transcripts[ which(edgetr$From == ss & edgetr$To == ee)])
      }else{
        trans<-paste(trans,edgetr$transcripts[ which(edgetr$From == ss & edgetr$To == ee)],sep = "|")
      }
      
    }
    trans <- sapply(strsplit(trans,"\\|"),function(X){
      tt<-unique(X)
      if(length(tt)==1){
        return(tt)
      }else{
        tran <- tt[1]
        for (j in 2:length(tt)){
          tran <- paste(tran,tt[j],sep="|")
        }
        return(tran)
      }
    })
    
    return(trans)
  })
  
  p2<-sapply(strsplit(p2,","),function(X){
    n<-length(X)
    trans <- ""
    for (i in 1:n){
      ss <- strsplit(X[i],"-")[[1]][1]
      ee <- strsplit(X[i],"-")[[1]][2]
      if (ss == ee){
        ss <- paste0(ss,".a")
        ee <- paste0(ee,".b")
      }else{
        ss <- paste0(ss,".b")
        ee <- paste0(ee,".a")
      }
      if (i==1){
        trans<-paste0(trans,edgetr$transcripts[ which(edgetr$From == ss & edgetr$To == ee)])
      }else{
        trans<-paste(trans,edgetr$transcripts[ which(edgetr$From == ss & edgetr$To == ee)],sep = "|")
      }
      
    }
    
    trans <- sapply(strsplit(trans,"\\|"),function(X){
      tt<-unique(X)
      if(length(tt)==1){
        return(tt)
      }else{
        tran <- tt[1]
        for (j in 2:length(tt)){
          tran <- paste(tran,tt[j],sep="|")
        }
        return(tran)
      }
    })
    
    return(trans)
  })
  
  pref<-sapply(strsplit(pref,","),function(X){
    n<-length(X)
    trans <- ""
    for (i in 1:n){
      ss <- strsplit(X[i],"-")[[1]][1]
      ee <- strsplit(X[i],"-")[[1]][2]
      if (ss == ee){
        ss <- paste0(ss,".a")
        ee <- paste0(ee,".b")
      }else{
        ss <- paste0(ss,".b")
        ee <- paste0(ee,".a")
      }
      if (i==1){
        trans<-paste0(trans,edgetr$transcripts[ which(edgetr$From == ss & edgetr$To == ee)])
      }else{
        trans<-paste(trans,edgetr$transcripts[ which(edgetr$From == ss & edgetr$To == ee)],sep = "|")
      }
      
    }
    
    trans <- sapply(strsplit(trans,"\\|"),function(X){
      tt<-unique(X)
      if(length(tt)==1){
        return(tt)
      }else{
        tran <- tt[1]
        for (j in 2:length(tt)){
          tran <- paste(tran,tt[j],sep="|")
        }
        return(tran)
      }
    })
    
    
    return(trans)
  })
  

  return(data.frame(p1=p1,p2=p2,ref=pref))

}