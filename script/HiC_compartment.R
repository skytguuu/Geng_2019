#compartment A/B#
#Rscript_mean#
library("HiTC")
name = list.files(pattern=".matrix")##批量读取符合某种匹配的数据
n = length(name)
pas = function(x){
  a = paste(x[,1],x[,2],sep = ":")
  b = paste(a,x[,3],sep = "-")
  return(b)
}
M = c()
#chr1-19为1:19,chrX为20
for (j in 1:19){
  chr.name = paste("chr",j,sep = "")
  chr_matrix = c()
  for (i in 1:n){
    Exp = importC(name[i],xgi.bed = "genome_400kb.bed",allPairwise = F,rm.trans = F,lazyload = F)		
    if (i == 1){
      chr.name <- pca.hic(Exp[isIntraChrom(Exp)][[j]],normPerExpected=TRUE,method = "mean",npc=1)
      chr.name = as.matrix(as.data.frame(chr.name))
      chr.name = chr.name[,c(3:5,8)]
      colnames(chr.name)[4] = name[i]
      chr_matrix = chr.name
    }
    else {
      chr.name <- pca.hic(Exp[isIntraChrom(Exp)][[j]],normPerExpected=TRUE,method = "mean",npc=1)
      chr.name = as.matrix(as.data.frame(chr.name))
      chr.name = chr.name[,c(3:5,8)]
      colnames(chr.name)[4] = name[i]
      chr_matrix = merge(chr_matrix, chr.name,by=c("seqnames","start","end"))
    }
  }
  M = rbind(M,chr_matrix)
}
write.csv(M,"chr1-19_mean_400kb.csv",row.names = TRUE)