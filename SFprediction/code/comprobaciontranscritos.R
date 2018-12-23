comprobaciontranscritos <- function(Result){
  
  
  Events <- Result[,c("tran_P1","tran_P2","tran_Ref")]
  
  comprobacion <- apply(Events,1,function(X){
    
    p1p2 <- paste(X[1],X[2],sep="|")
    p1p2 <- sapply(strsplit(p1p2,"\\|"),function(XX) return(XX))
    
    p3<-as.character(X[3])
    p3 <- sapply(strsplit(p3,"\\|"),function(XX) return(XX))
    
    
    r1 <- (!(any(p1p2%in%p3==FALSE) | any(p3%in%p1p2==FALSE)))
    
  })
  
}

comprobaciontranscritos2 <- function(Result){
  Events <- Result[,c("tran_P1","tran_P2","tran_Ref")]
  
  comprobacion <- apply(Events,1,function(X){
    #condicion 0: los transcritos de un path no pueden ser unicamente "no mapeado en la referencia" (hay q eliminarlos)
    p1 <- unlist(strsplit(as.character(X[1]),"\\|"))
    p2 <- unlist(strsplit(as.character(X[2]),"\\|"))
    p3 <- unlist(strsplit(as.character(X[3]),"\\|"))
    
    iix1 <- which(p1=="no mapeado en la referencia")
    if(length(iix1)>0){
      p1<-p1[-iix1]
    }
    
    iix2 <- which(p2=="no mapeado en la referencia")
    if(length(iix2)>0){
      p2<-p2[-iix2]
    }
    
    iix3 <- which(p3=="no mapeado en la referencia")
    if(length(iix3)>0){
      p3<-p3[-iix3]
    }
    
    if(length(p1)>0 & length(p2)>0 & length(p3)>0){
      
      
      #condicion 1: los que estan en p1 no pueden estar en p2
      
      cond1 <- (any(p1%in%p2==T) | any(p2%in%p1==T)) #Falso si se cumple la condicion (True si no se cumple)
      
      #condicion 2: la suma de los que estan en p1 y p2 tienen q ser los mismos que los q estan en p3 
      
      #p1p2 <- paste(X[1],X[2],sep="|")
      #p1p2 <- sapply(strsplit(p1p2,"\\|"),function(XX) return(XX))
      
      #p3<-as.character(X[3])
      #p3 <- sapply(strsplit(p3,"\\|"),function(XX) return(XX))
      
      p1p2 <- c(p1,p2)
      
      cond2 <- (any(p1p2%in%p3==FALSE) | any(p3%in%p1p2==FALSE)) #Falso si se cumple la condicon (True si no se cumple)
      
      return(!(cond1 | cond2))
      
      
    }else{
      return(FALSE)
    }
    
    
    
  })
  
  
}








