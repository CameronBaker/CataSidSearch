
#Creates a list of Amino acid susbstitutions by row in the dataframe
#For example ASP-ASP would not show up in this list but ASP-ARG would
#Load into memory for actual use, not called normally
aasubs <- function(df){
  listsubs <- c()
  for(i in 1:nrow(df)){
    for(j in 4:ncol(df)){
      tmp = TRUE
      if(!is.na(df[i,j]) && tmp){
        y <- unlist(strsplit(df[i,j],"-"))
        if(y[1] != y[2]){
          listsubs <- c(listsubs,rownames(df)[i])
          tmp = FALSE
        }
      }
    }
  }
  return(listsubs)
}

readmotifs <- function(){
  files = list.files(pattern="cdb*")
  
  df <- data.frame(matrix(ncol=15,nrow=1))
  for(i in files){
    nr <- unlist(read.csv(i,header=FALSE)[5,])
    if(!is.na(nr[1])){
      if(nr[1]=="Target"){
        nr <- unlist(read.csv(i,header=FALSE)[6,])
      }
      newrow <- c()
      for(j in nr){
        newrow <- c(newrow,j)
      }
      for(j in length(newrow):14){
        newrow <- c(newrow,NA)
      }
      df <- rbind(df,newrow)
    }
  }
  df <- subset(df,select = -c(X5,X6,X15))
  df <- df[rowSums(is.na(df)) != ncol(df),]
  df[is.na(df)] <- ""
  rownames(df) <- df[,1]
  colnames(df) <- c("Res","Match Score","Best Match","EC","AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8")
  df <- subset(df[,2:ncol(df)])
  return(df)
}

df <- readmotifs()