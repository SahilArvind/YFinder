#' Depth Finder
#'
#' Find depth of a subclade on the YFull json tree (not actual depth, but depth
#'  in the json which is always more)
#' @param hg The subclade eg R-Z2124
#' @param jsonnest the loaded JSON yfull file (use jsonlite::fromJSON() in package 
#' 'jsonlite')
#' @return The depth in the YFull v11.02 json nested list
#' @examples 
#' depth <- depthfinder("R-Z2124",loadedjsonobject);


depthfinder <- function(hg,jsonnest){
  abcd=rrapply::rrapply(
    jsonnest, 
    condition = \(x, .xname) (hg %in% x) && .xname == "id",
    f = \(x, .xpos) length(.xpos), 
    how = "unlist"
  ) %>% unname(.) %>% { ifelse(is.null(.), 0, . ) }
  
  return(abcd)
}

#' Call Children
#'
#' Find the immediate children of a subclade as per provided yfull json tree
#' @param hg The subclade eg R-Z2124
#' @param jsonnest the loaded JSON yfull file (use fromJSON() in package 
#' 'jsonlite')
#' @return ancestor of the subclade
#' @examples 
#' anc <- call_children("R-Z2124",loadedjsonobject);

call_children <- function(hg,jsonnest){
  
  abcd=rrapply::rrapply(
    jsonnest, 
    condition = \(x, .xname) (hg %in% x) && .xname == "id",
    f = \(x, .xparents) .xparents, 
    how = "flatten"
  )
  #index1=as.numeric(tail(abcd[[1]][1:(max(1,length(abcd[[1]])-1))],n=1))
  #abcd[[1]] <- abcd[[1]][1:(max(1,length(abcd[[1]])-3))]
  abcd1=sapply(abcd[[1]],function(x){ifelse(nchar(x)>3 || x == "id",paste('$',x,sep=""),paste('[[',x,']]',sep=""))})
  abcd1=stringr::str_flatten(abcd1)
  abcd1 <- paste(deparse(substitute(jsonnest)),abcd1,sep="")
  lista=eval(parse(text=abcd1))
  index1=which(lista==hg)
  abcd[[1]] <- abcd[[1]][1:(max(1,length(abcd[[1]])-1))]
  abcd=sapply(abcd[[1]],function(x){ifelse(nchar(x)>3 || x == "id",paste('$',x,sep=""),paste('[[',x,']]',sep=""))})
  abcd=stringr::str_flatten(abcd)
  abcd <- paste(deparse(substitute(jsonnest)),abcd,sep="")
  
  return(eval(parse(text=paste(abcd,"$children[[",index1,"]]$id",sep=""))))
}




#' Call Ancestor
#'
#' Find the immediate ancestor of a subclade as per provided yfull json tree
#' @param hg The subclade eg R-Z2124
#' @param jsonnest the loaded JSON yfull file (use fromJSON() in package 
#' 'jsonlite')
#' @param yfullsnp additional yfull snp list file provided in/data
#' @return ancestor of the subclade
#' @examples 
#' anc <- call_ancestor("R-Z2124",loadedjsonobject,loadedyfullsnpfile);


call_ancestor <- function(hg,jsonnest,yfullsnp){
  
if(yfullsnp$depth[which(yfullsnp$HG == hg)][[1]]<=2){
  return(hg)
}
  
  abcd=rrapply::rrapply(
    jsonnest, 
    condition = \(x, .xname) (hg %in% x) && .xname == "id",
    f = \(x, .xparents) .xparents, 
    how = "flatten"
  )
if (length(abcd)==0){return(NA)}
index1=as.numeric(tail(abcd[[1]][1:(max(1,length(abcd[[1]])-1))],n=1))
abcd[[1]] <- abcd[[1]][1:(max(1,length(abcd[[1]])-3))]
abcd=sapply(abcd[[1]],function(x){ifelse(nchar(x)>3 || x == "id",paste('$',x,sep=""),paste('[[',x,']]',sep=""))})
abcd=stringr::str_flatten(abcd)
abcd <- paste(deparse(substitute(jsonnest)),abcd,sep="")

return(eval(parse(text=paste(abcd,"$id[[",index1,"]]",sep=""))))
}


#' Check if subclade is +ve and return path
#'
#' Check if subclade is +ve and return path
#' @param hg The subclade eg R-Z2124
#' @param jsonnest the loaded JSON yfull file (use jsonlite::fromJSON() in package 
#' 'jsonlite')
#' @param calltable object provided by calling function
#' @param callquality value provided by calling function, default 0.9
#' @param count value provided by calling function, count of ancestors who 
#' are +ve as well
#' @param depth value provided by calling function, default 2
#' @param path list provided by calling function
#' @param yfullsnp table provided by calling function
#' @param numsnptotal value provided by calling function. default 2. 
#' Min number of SNPs needed to give +ve call. 1 gives false calls.
#' @return List object with Path and derived status TRUE or FALSE


checkpath <- function(hg,jsonnest,calltable,callquality,count=1,depth=2,path=NULL,yfullsnp,numsnptotal){

if (!hg %in% yfullsnp$HG){
    newList <- list("Path" = c(), "Derived" = FALSE)
    return(newList)
} else if(yfullsnp$depth[which(yfullsnp$HG == hg)][[1]] == 2){
  newList <- list("Path" = path, "Derived" = TRUE)
  return(newList)
} else if (yfullsnp$depth[which(yfullsnp$HG == hg)][[1]] == 0){
  newList <- list("Path" = c(), "Derived" = FALSE)
  return(newList)
} 

anc_hg <- call_ancestor(hg,jsonnest,yfullsnp)
if (is.na(anc_hg)){
  newList <- list("Path" = c(), "Derived" = FALSE)
  return(newList)
}
newtable <- calltable %>% filter((HG %in% anc_hg) & (!sample %in% NA))

call1 <- nrow(newtable[newtable$sample  == "D",])/nrow(newtable)
if (is.null(call1) || is.nan(call1) || is.na(call1)) {
  path=paste(path,"<",anc_hg,"? ",sep="")
  return(checkpath(anc_hg,jsonnest,calltable,callquality=callquality,
                   count=count,depth=depth,path=path,yfullsnp = yfullsnp,
                   numsnptotal))
} else if (call1<callquality){
  newList <- list("Path" = c(), "Derived" = FALSE)
  return(newList)
} else if (call1 >= callquality){
  count=count+1
  path=paste(path,"<",anc_hg,"+ ",sep="")
  if (count>depth || nrow(newtable) >= numsnptotal ){
    newList <- list("Path" = path, "Derived" = TRUE)
    return(newList)
  }
  return(checkpath(anc_hg,jsonnest,calltable,callquality=callquality,count=count
                   ,depth=depth,path=path,yfullsnp = yfullsnp,numsnptotal))
} 
}

#' YFinder - Get haplogroups from VCF File
#'
#' Get list of all Y Hg calls based on YFull from VCF file of Y chr 
#' with multiple samples
#' @param vcffile The location of the vcf file
#' @param yfullpositionsfile The location of the yfullpositions file (
#' provided in /data folder)
#' @param indivlistfile the location of the file with individual sample ids 
#' to be tested. 1 sampleid per line. Last line should be empty. Ensure that
#' each id matches exactly with the sampleid column names in the vcf file.
#' @param yfulljsonfile location of the yfulljsonfile (provided in /data)
#' @param transversionsonly TRUE or FALSE. TRUE ignore C/T and G/A positions
#' Use TRUE for ancient samples to avoid false calls
#' @param depth number of ancestors to check in case of a single +ve SNP. 
#' default 2
#' @param callquality default 0.9. Percentage of Derived calls per subclade for 
#' assigning derived status to subclade.
#' @param numsnptotal Min number of SNPs to decide +ve subclade 
#' (without checking for ancestors). Default 2. Avoid 1.
#' @param viewrealtime Updates and shows the output table after every individual 
#' sample is processed for realtime view. Deafult is FALSE.
#' @return Data frame object with list of assigned haplogroups and 
#' a call table to crossverify manually
#' @export
#' @import tidyverse
#' @import dplyr
#' @import rrapply
#' @import jsonlite
#' @import purrr
#' @import stringr





YFinder <- function(vcffile,yfullpositionsfile,indivlistfile,
                          yfulljsonfile,transversionsonly=F,depth=1,
                          callquality=0.9,numsnptotal=2,viewrealtime=FALSE){

library(dplyr)
vcf=readLines(vcffile)
jsonnest <- jsonlite::fromJSON(yfulljsonfile,simplifyVector = T,simplifyDataFrame = T)

for (i in 1:length(vcf)){
  if(substr(vcf[[i]][1],1,6) == "#CHROM"){
    break
  }
}

vcf <- vcf[i:length(vcf)]
write.table(vcf, file = "tempvcf.txt", sep = "",row.names = FALSE,
            col.names = FALSE)
vcf <- gsub("#CHROM","CHROM",vcf)
vcf <- gsub("0/0","0",vcf)
vcf <- gsub("1/1","1",vcf)
vcf <- gsub("./.",".",vcf)
vcf <- gsub("0|0","0",vcf)
vcf <- gsub("1|1","1",vcf)

write.table(vcf, file = "tempvcf.txt", sep = "",row.names = FALSE,
            col.names = FALSE)
a=read.table("tempvcf.txt",sep = "\t",quote ="",colClasses = "character",
             header = TRUE )
rm(vcf)
file.remove("tempvcf.txt")
yfullsnp=read.table(yfullpositionsfile,sep=" ",header = TRUE)

newf <- merge(yfullsnp[,2:8],a[2:ncol(a)],by='POS')
colnames(newf)[ncol(newf)]=substr(colnames(newf)[ncol(newf)],1,-1+nchar(colnames(newf)[ncol(newf)]))
test <- colnames(newf[15:ncol(newf)])
#test <- sapply(test,function(x){gsub("\\.","",as.character(x))})
#colnames(newf)[15:ncol(newf)] <- test
newf[,ncol(newf)] <- sapply(newf[,ncol(newf)], 
                          function(x){gsub('"',"",as.character(x))})
#newf[,ncol(newf)] <- sapply(newf[,ncol(newf)], function(x){as.numeric(x)})

indivlist <- as.matrix(read.table(indivlistfile))

if (!purrr::is_empty(setdiff(indivlist,test))){
print(setdiff(indivlist,test))
print("These samples are not present in VCF, ignoring them")
}

temp <- indivlist[!indivlist %in% setdiff(indivlist,test)]
newf3 <- newf[,1:14]
newf2 = newf %>% select(all_of(temp))
newf <- cbind(newf3,newf2)
indivlist <- colnames(newf)[15:ncol(newf)]
#colnames(newf)[15:ncol(newf)] <- temp


tempdf <- data.frame()
j=0
for (j in 15:ncol(newf)){
  m=j-14
  tempdf <- as.data.frame(newf[,j])
  colnames(tempdf) <- c("V1")
  tempdf$V1 <- ifelse(tempdf$V1 == 0,newf$REF,ifelse(tempdf$V1 == 1,newf$ALT,NA))
  tempdf$V1 <- ifelse(tempdf$V1 == newf$ANC,"A",ifelse(tempdf$V1 == newf$DER,"D",NA))
  newf[,j] <- tempdf$V1
  colnames(newf)[j] <- indivlist[[m]]
}

finaloutput <- data.frame(SampleID=rep(NA,length(indivlist)),
                        HG=rep(NA,length(indivlist)),
                        Path=rep(NA,length(indivlist)))


#Progress bar initialization 
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(indivlist),   # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar


for (k in 1:length(indivlist)){
  
  if(viewrealtime==TRUE){
   view(finaloutput)
  }
  
  temp = newf %>% select(all_of(c("POS","SNP","HG","depth","MUTATION",c(indivlist[k]))))
  temp$hgdepth <- paste(substr(temp$HG,1,1),temp$depth,sep="")
  colnames(temp)[6] <- "sample"

  if (transversionsonly == TRUE){
    temp1 <- temp %>% filter(sample == "D" & !MUTATION %in% c("C->T","G->A") )
    temp2 <- temp %>% filter(!MUTATION %in% c("C->T","G->A") )
  } else{
    temp1 <- temp %>% filter(sample == "D")
    temp2 <- temp
  }
  if (dim(temp1)[1] == 0){
    finaloutput$SampleID[[k]] <- indivlist[[k]]
    finaloutput$HG[[k]] <- "NoSNP/lowquality"
    finaloutput$Path[[k]] <- c("")
    next 
    
  }
    
  temp1 <- temp1[order(temp1$hgdepth,decreasing = TRUE),]
  
  
  for (j in 1:nrow(temp1)){
    
    hg <- temp1$HG[[j]]
    
    newtable <- temp2 %>% filter((HG %in% hg) & (!sample %in% NA))
    call1 <- nrow(newtable[newtable$sample  == "D",])/nrow(newtable)

    if(call1 >= callquality && nrow(newtable) >= numsnptotal){
      finaloutput$SampleID[[k]] <- indivlist[[k]]
      finaloutput$HG[[k]] <- hg
      finaloutput$Path[[k]] <- c("")
      break
    } else if(call1 < callquality){
      next
      
    }
    
    out_hg <- checkpath(hg=hg,jsonnest=jsonnest,calltable=temp2,
                      callquality=callquality,depth=depth,yfullsnp=yfullsnp,numsnptotal)
    if (out_hg$Derived == TRUE){
      finaloutput$SampleID[[k]] <- indivlist[[k]]
      finaloutput$HG[[k]] <- hg
      finaloutput$Path[[k]] <- stringr::str_flatten(c(paste(hg,"+ ",sep=""),out_hg$Path))
      break
    } else {
      finaloutput$SampleID[[k]] <- indivlist[[k]]
      finaloutput$HG[[k]] <- "NoSNP/lowquality"
      finaloutput$Path[[k]] <- c("")
      next
    }
  }
  setTxtProgressBar(pb, k)
  
}
close(pb)
return(list(Final_Assignment=finaloutput,CallsTable=newf))
}

