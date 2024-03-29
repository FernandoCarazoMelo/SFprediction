AnnotateEvents_KLL <- function (Events,Gxx,GenI){
  {
    #
    #Gxx <- GeneName
    #
    GeneName<-Gxx
    GeneID<-GenI
    Chrom<-gsub("chr","",as.vector(Events[[1]]$P1[1,"Chr"]))
    Result<-vector("list")
    #Flat<-vector("list")
    mm<-0
    
    for(ii in seq_along(Events))
    {
      if (!any(c(identical(unique(Events[[ii]]$P1$Type),"V"),identical(unique(Events[[ii]]$P2$Type),"V"),identical(unique(Events[[ii]]$Ref$Type),"V"))==T))
      {
        mm<-mm+1
        
        EventNumber<-ii
        
        EventType<-Events[[ii]]$Type
        
        Positions<-rbind(Events[[ii]]$P1,Events[[ii]]$P2)[,4:5]
        Start<-as.numeric(Positions[,1])
        End<-as.numeric(Positions[,2])
        Start<-Start[which(Start!=0)]
        End<-End[which(End!=0)]
        
        # browser()
        minGPos<-min(Start)
        maxGPos<-max(End)
        GPos<-paste(Chrom,":",minGPos,"-",maxGPos,sep="")
        
        CP1s<-which(Events[[ii]]$P1[,1]=="S")
        CP1e<-which(Events[[ii]]$P1[,2]=="E")
        
        if(length(CP1s)>0|length(CP1e)>0)
        {
          CC<-c(CP1s,CP1e)
          Events[[ii]]$P1<-Events[[ii]]$P1[-CC,]
        }
        
        PS1<-as.numeric(gsub(".[ab]","",Events[[ii]]$P1[,1]))
        PE1<-as.numeric(gsub(".[ab]","",Events[[ii]]$P1[,2]))
        Path1<-as.matrix(cbind(PS1,PE1))
        Path1<-Path1[order(Path1[,1],Path1[,2]),,drop=FALSE]
        
        CP2s<-which(Events[[ii]]$P2[,1]=="S")
        CP2e<-which(Events[[ii]]$P2[,2]=="E")
        
        if(length(CP2s)>0|length(CP2e)>0)
        {
          CC<-c(CP2s,CP2e)
          Events[[ii]]$P2<-Events[[ii]]$P2[-CC,]
        }
        
        PS2<-as.numeric(gsub(".[ab]","",Events[[ii]]$P2[,1]))
        PE2<-as.numeric(gsub(".[ab]","",Events[[ii]]$P2[,2]))
        Path2<-as.matrix(cbind(PS2,PE2))
        Path2<-Path2[order(Path2[,1],Path2[,2]),,drop=FALSE]
        
        CPRs<-which(Events[[ii]]$Ref[,1]=="S")
        CPRe<-which(Events[[ii]]$Ref[,2]=="E")
        
        if(length(CPRs)>0|length(CPRe)>0)
        {
          CC<-c(CPRs,CPRe)
          Events[[ii]]$Ref<-Events[[ii]]$Ref[-CC,]
        }
        
        PSR<-as.numeric(gsub(".[ab]","",Events[[ii]]$Ref[,1]))
        PER<-as.numeric(gsub(".[ab]","",Events[[ii]]$Ref[,2]))
        PathR<-as.matrix(cbind(PSR,PER))
        PathR<-PathR[order(PathR[,1],PathR[,2]),,drop=FALSE]
        
        
        Path1<-paste(Path1[,1],"-",Path1[,2],sep="",collapse=",")
        Path2<-paste(Path2[,1],"-",Path2[,2],sep="",collapse=",")
        PathR<-paste(PathR[,1],"-",PathR[,2],sep="",collapse=",")
        
        
        NEv<-data.frame(GeneName,GeneID,EventNumber,EventType,GPos,Path1,Path2,PathR,stringsAsFactors = FALSE)
        Result[[mm]]<-NEv
        
        
      }
      
    }
    
    Result<-do.call(rbind,Result)
    
    return(Result)
  }
  
}