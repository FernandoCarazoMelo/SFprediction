transfromedge <- function(SG,SG_Gene){
  # La salida es un data.frame que tiene el numero del edge, su start y end y a los transcritos a los que pertenece
  
  
  A <- SG$Edges # lo que viene de SG_Creation
  
  SG_Gene_SoloE <- SG_Gene[type(SG_Gene)=="E"]
  B <- SG_Gene_SoloE # El correspondiente a un Gen
  
  nn<- length(txName(B))
  aa <- c()
  for (i in 1:nn){
    aa <- c(aa,txName(B)[[i]])
  }
  aa <- unique(aa)
  
  edge <- 1:dim(A)[1]
  transcripts <- vector(mode="character",length = length(edge))
  
  iixe <-which(A$Type=="E")
  
  matrixexons <- matrix(0,nrow=length(aa),ncol=length(iixe))
  colnames(matrixexons)<-iixe
  rownames(matrixexons)<-aa
  
  for (ii in iixe){
    s <- A$Start[ii]
    e <- A$End[ii]
    ix <- which(start(B)==s & end(B)==e)
    trans <- txName(B)[[ix]]
    n <- length(trans)
    if (n == 0){
      transcripts[ii]<-"no mapeado en la referencia"
    }else if (n==1){
      matrixexons[trans,as.character(ii)]<-1
      transcripts[ii]<-trans
    }else {
      matrixexons[trans,as.character(ii)]<-1
      tt <-trans
      trans <- trans[1]
      for (kk in 2:n){
        trans <- paste(trans,tt[kk],sep="|")
      }
      transcripts[ii] <- trans
    }
    
  }
  
  
  iix <-which(A$Type=="V")
  transcripts[iix] <- ""
  
  
  
  #matrixexons
  iix <-which(A$Type=="J")
  for(ii in iix){
    s <- as.character(A$From[ii])
    e <- as.character(A$To[ii])
    
    # ex1 <- as.character(iixe[which(A$To[iixe]==s)])
    # ex2 <- as.character(iixe[which(A$From[iixe]==e)])
    
    ex1 <- which(A$To[iixe]==s)
    ex2 <- which(A$From[iixe]==e)
    
    
    trans <- rownames(matrixexons)[which(rowSums(matrixexons[,ex1:ex2])==2 & matrixexons[,ex1]==1 & matrixexons[,ex2]==1)]
    
    n <- length(trans)
    if (n == 0){
      transcripts[ii]<-"no mapeado en la referencia"
    }else if (n==1){
      
      transcripts[ii]<-trans
    }else {
      
      tt <-trans
      trans <- trans[1]
      for (kk in 2:n){
        trans <- paste(trans,tt[kk],sep="|")
      }
      transcripts[ii] <- trans
    }
    
    
    
    
  }
  
  From <- A$From
  To <- A$To
  Start <- A$Start
  End <- A$End
  D <- data.frame(edge=edge,From=From,To=To,Start=Start, End = End,transcripts = transcripts)
  return(D)

}











