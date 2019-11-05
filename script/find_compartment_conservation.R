#!usr/bin/env/Rscript
#find conservation location using P/M PC1 value#
#mean correlation#
setwd("F:/HiC_garen/HiC_result/development/compartment/PM")
maternal = read.csv("dev_m.csv",header = TRUE)
paternal = read.csv("dev_p.csv",header = TRUE)
maternal[is.na(maternal)] = 0
paternal[is.na(paternal)] = 0
#maternal
data = maternal
data[is.na(data)] = 0
out = data.frame(chr=character(0),start=numeric(0),end=numeric(0),correlation=numeric(0),stringsAsFactors=FALSE)
for (k in 1:19){
  name = paste("chr",k,sep="")
  part = subset(data,data$chr == name)
  n = floor(nrow(part)/5)
  res = data.frame(chr=character(0),start=numeric(0),end=numeric(0),correlation=numeric(0),stringsAsFactors=FALSE)
  for (i in 1:n){
    n_row = i*5
    start = (i-1)*5+1
    part1 = part[start:n_row,]
    r = cor(part1[,-c(1:3)])
    M = as.matrix(r)
    correlation = mean(abs(M))
    res[i,1] = name
    res[i,2] = part1[1,2]
    res[i,3] = part1[5,3]
    res[i,4] = correlation
  }
  if (nrow(part) != n*50){
    part1 = part[(n*5+1):nrow(part),]
    r = cor(part1[,-c(1:3)])
    M = as.matrix(r)
    correlation = mean(abs(M))
    res[(n+1),1] = name
    res[(n+1),2] = part1[1,2]
    res[(n+1),3] = part1[nrow(part1),3]
    res[(n+1),4] = correlation		
  }
  res = res[which(res[,4]>0.6),]
  out = rbind(out,res)
}
write.csv(out,"maternal_PC1_correlation.csv",quote=FALSE,row.names = FALSE)
#paternal
data = paternal
data[is.na(data)] = 0
out = data.frame(chr=character(0),start=numeric(0),end=numeric(0),correlation=numeric(0),stringsAsFactors=FALSE)
for (k in 1:19){
  name = paste("chr",k,sep="")
  part = subset(data,data$chr == name)
  n = floor(nrow(part)/5)
  res = data.frame(chr=character(0),start=numeric(0),end=numeric(0),correlation=numeric(0),stringsAsFactors=FALSE)
  for (i in 1:n){
    n_row = i*5
    start = (i-1)*5+1
    part1 = part[start:n_row,]
    r = cor(part1[,-c(1:3)])
    M = as.matrix(r)
    correlation = mean(abs(M))
    res[i,1] = name
    res[i,2] = part1[1,2]
    res[i,3] = part1[5,3]
    res[i,4] = correlation
  }
  if (nrow(part) != n*50){
    part1 = part[(n*5+1):nrow(part),]
    r = cor(part1[,-c(1:3)])
    M = as.matrix(r)
    correlation = mean(abs(M))
    res[(n+1),1] = name
    res[(n+1),2] = part1[1,2]
    res[(n+1),3] = part1[nrow(part1),3]
    res[(n+1),4] = correlation		
  }
  res = res[which(res[,4]>0.6),]
  out = rbind(out,res)
}
write.csv(out,"paternal_PC1_correlation.csv",quote=FALSE,row.names = FALSE)
