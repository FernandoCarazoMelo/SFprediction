ggmatplot <- function(dat, xnames=TRUE, highlight = NULL){
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(tibble)
  
  if(!is.matrix(dat)){
    dat <- as.data.frame(dat)
    df <- dat %>% 
      rownames_to_column() %>% 
      gather(reading, value, -rowname) %>% 
      group_by(rowname) %>% 
      mutate(x=1:n()) 
    colnames(df)[1] <- "reading"
    colnames(df)[2] <- "Legend"
    
    if(xnames){
      a <- ggplot(data = df, aes(x=1:nrow(df), y=value)) +
        geom_line(aes(color = Legend), size = 1.5) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(limits = df$reading)
    }else{
      a <- ggplot(data = df, aes(x=1:nrow(df), y=value)) +
        geom_line(aes(color = Legend), size = 1.5)
    }
    
  }else{
    dat <- as.data.frame(t(dat))
    df <- dat %>% 
    rownames_to_column() %>% 
    gather(reading, value, -rowname) %>% 
    group_by(rowname) %>% 
    mutate(x=1:n()) 
    colnames(df)[1] <- "Legend"
    
  
    if(xnames){
      df$reading <- factor(df$reading,levels = unique(df$reading))
      a <- ggplot(data = df, aes(x=reading, y=value, group=Legend)) +
        geom_line(aes(color=Legend)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }else{
      a <- ggplot(df, aes(x=x, y=value, group=Legend)) +
        geom_line(aes(color=Legend))
    }
    
    if(length(highlight > 0)){a <- a + geom_line(data = df[which(df$Legend %in% highlight), ],aes(color=Legend), size = 1.5)}
}
    return(a)
}
