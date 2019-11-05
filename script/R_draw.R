#!usr/bin/env/Rscript
#1.Compartment A/B plot#
setwd("F:/HiC_result/compartment/400kb")
pca = read.csv("R_compartment.csv",header = TRUE)
cop = pca[pca$Chr=="chr8" & pca$Start>=2800000 & pca$Start<=19700000,]
plot(x=cop$Start,y=cop$iPS,type="h",lwd=3,col=ifelse(cop$iPS>0,"darkblue","orange"))
plot(x=cop$Start,y=cop$SSC,type="h",lwd=3,col=ifelse(cop$SSC>0,"darkblue","orange"))
plot(x=cop$Start,y=cop$FGSC,type="h",lwd=3,col=ifelse(cop$FGSC>0,"darkblue","orange"))
plot(x=cop$Start,y=cop$NPC,type="h",lwd=3,col=ifelse(cop$NPC>0,"darkblue","orange"))
plot(x=cop$Start,y=cop$STO,type="h",lwd=3,col=ifelse(cop$STO>0,"darkblue","orange"))
cop = pca[pca$Chr=="chr6" & pca$Start>=112400001 & pca$Start<=137600000,]
#2.Compartment PCA Heatmap#
setwd("F:/HiC_result/compartment/400kb/result")
T = read.table("R_pca_input.txt",header = TRUE)
data = T[,-1]
rownames(data) = T[,1]
data[is.na(data)]=0
library("circlize")
library("ComplexHeatmap")
set.seed(11)
hm=colorRampPalette(c('yellow','darkblue'))
Heatmap(as.matrix(data),cluster_columns = FALSE,km=2,show_row_names = FALSE,col=hm(20),name="PC1 Value")
#3.chrX TAD signal#
setwd("F:/HiC_result/TAD/chrX/normalization")
T = read.csv("chrX_normalize.csv",header = TRUE)
data = subset(T,T$Start>=100000000 & T$End<=102000000)
plot(x=data$Start,y=data$FGSC,type="l",lwd=3,col="#E69F00",bty="n",ylim=c(-0.8,0.8))
lines(x=data$Start,y=data$SSC,type="l",lwd=3,col="#0072B2")
lines(x=data$Start,y=data$FGSC.SSC,type="l",lwd=3,col="red")
abline(h=0,col="black",lty=2,lwd=2)
#4.chrX Compartment signal#
setwd("F:/HiC_result/compartment/400kb/chrX/merge")
E = read.table("R_chrx_input.txt",header = TRUE)
plot(x=E[,2],y=E$FGSC,type="l",col="#E69F00",lwd=3,ylab="PC1 Score",xlab = "Genome Region",bty="l")
abline(h=0,col="black",lty=2,lwd=2)
lines(x=E[,2],y=E$SSC,type="l",col="#0072B2",lwd=3)
lines(x=E[,2],y=E$fESC,type="l",col="#009E73",lwd=3)
legend("top",legend=paste(c("FGSC","SSC","female_ESC")),col=c("#E69F00","#0072B2","#009E73"),lty=1)
#autosome#
setwd("F:/HiC_result/compartment/400kb/pca_loess")
T = read.csv("pca.csv",header = TRUE)
tiff("autosome1.tiff",height=1000,width = 1000)
layout(matrix(1:20,10,2))
par(oma=c(2,4,4,2),mar=c(1,4,1,4),cex=1,las=2)
for (i in 1:19){
  ch = paste("chr",i,sep="")
  chr = paste(ch,":",sep="")
  data = T[grep(pattern = chr,x=T$location),c(3,4)]
  data = na.omit(data)
  corr = cor(data$FGSC,data$SSC)
  co = paste("r",corr,sep="=")
  tit = paste(ch,corr,sep=":")
  plot(x=row.names(data),y=data$FGSC,type="l",col="#E69F00",lwd=3,ylab="PC1 value",cex.lab=1.2,font.lab=2,xlab = "Genome Region",axes = FALSE,ylim=c(-0.1,0.1),main=tit)
  axis(2,lwd=2,cex.axis=1,font.axis=2)
  abline(h=0,col="black",lty=2,lwd=2)
  lines(x=row.names(data),y=data$SSC,type="l",col="#0072B2",lwd=3)
}
dev.off()
#5. TAD signal 
library("ggplot2")
setwd("F:/HiC_garen/HiC_result/TAD")
tad = read.table("FGSC_20kb_DI_TAD.txt")
filelist <- list.files(pattern=".bedGraph")
datalist <- lapply(filelist, function(x) read.table(x,header=TRUE,stringsAsFactors=F))
m = length(filelist)
for (k in 1:m){
  data = as.data.frame(datalist[k])
  data[is.na(data)]=0
  f_loc = c()
  f_value = c()
  for (i in 1:nrow(tad)){
    length = tad[i,3]-tad[i,2]
    up = 20000*floor((tad[i,2]-0.5*length)/20000)
    down = 20000*ceiling((tad[i,3]+0.5*length)/20000)
    part = subset(data,data$Chr == as.character(tad[i,1]) & up<=data$Start & down>data$End)
    fgsc = c()
    fgsc_loc = c()
    for (j in 1:nrow(part)){
      if (part[j,4]>0.25){
        fgsc_loc = c(fgsc_loc,j*20000/length)
        fgsc = c(fgsc,part[j,4])
      }
    }
    f_loc = c(f_loc,fgsc_loc)
    f_value = c(f_value,fgsc)
  }
  name = strsplit(filelist[k],split="_")
  cell = rep(name[[1]][1],length(f_loc))
  if (k == 1){
    mat = data.frame(location=f_loc,type=cell,value=f_value)
  }
  else {
    tmp = data.frame(location=f_loc,type=cell,value=f_value)
    mat = rbind(mat,tmp)
  }
}
mat = mat[mat$type=="FGSC" | mat$type=="SSC",]
ggplot(mat,aes(x=location,y=value))+geom_smooth(aes(fill=type),method="loess",span=0.05,se=FALSE)+theme_bw()+theme(panel.grid=element_blank())+scale_fill_manual(values = c("#CC79A7", "#D55E00", "#0072B2","#009E73","#E69F00","darkred"))
#final figure#
setwd("F:/HiC_garen/HiC_result/development/TAD_signal")
filelist <- list.files(pattern=".bedGraph")
datalist <- lapply(filelist, function(x) read.table(x,header=TRUE,stringsAsFactors=F))
m = length(filelist)
for (k in 1:m){
  data = as.data.frame(datalist[k])
  data[is.na(data)]=0
  f_loc = c()
  f_value = c()
  for (i in 1:nrow(tad)){
    length = tad[i,3]-tad[i,2]
    up = 20000*floor((tad[i,2]-0.5*length)/20000)
    down = 20000*ceiling((tad[i,3]+0.5*length)/20000)
    part = subset(data,data$Chr == as.character(tad[i,1]) & up<=data$Start & down>data$End)
    fgsc = c()
    fgsc_loc = c()
    for (j in 1:nrow(part)){
      if (part[j,4]>0.25){
        fgsc_loc = c(fgsc_loc,j*20000/length)
        fgsc = c(fgsc,part[j,4])
      }
    }
    f_loc = c(f_loc,fgsc_loc)
    f_value = c(f_value,fgsc)
  }
  name = strsplit(filelist[k],split="_")
  cell = rep(name[[1]][1],length(f_loc))
  if (k == 1){
    mat = data.frame(location=f_loc,type=cell,value=f_value)
  }
  else {
    tmp = data.frame(location=f_loc,type=cell,value=f_value)
    mat = rbind(mat,tmp)
  }
}
r1 = ggplot(mat,aes(x=location,y=value))+theme_bw()+theme(panel.grid=element_blank())
r2 = r1 + geom_smooth(data = mat[mat$type=="FGSC",],aes(color=type),size=2,method="loess",color="black",span=0.05,se=FALSE,show.legend = TRUE)#scale_color_manual(values = "black")
r3 = r2 + geom_smooth(data = mat[mat$type=="GV",],aes(color=type),size=2,method="loess",span=0.05,se=FALSE)
r4 = r3 + geom_smooth(data = mat[mat$type=="MII",],aes(color=type),size=2,method="loess",span=0.05,se=FALSE)
r5 = r4 + geom_smooth(data = mat[mat$type=="zygote",],aes(color=type),size=2,method="loess",span=0.05,se=FALSE)
r6 = r5 + geom_smooth(data = mat[mat$type=="early2",],aes(color=type),size=2,method="loess",span=0.05,se=FALSE)
r7 = r6 + geom_smooth(data = mat[mat$type=="late2",],aes(color=type),size=2,method="loess",span=0.05,se=FALSE)
r8 = r7 + geom_smooth(data = mat[mat$type=="8cell",],aes(color=type),size=2,method="loess",span=0.05,se=FALSE)
r8

#6.genome contact distance of development#
setwd("F:/HiC_garen/HiC_result/development/contact_distance")
T = read.table("r_input_dis.txt",header = TRUE,stringsAsFactors=FALSE)
Exp = log10(T)
model1 = loess(FGSC~distance,data=Exp,span=0.1)
model2 = loess(GV~distance,data=Exp,span=0.1)
model3 = loess(MII~distance,data=Exp,span=0.1)
model4 = loess(Zygote~distance,data=Exp,span=0.1)
model5 = loess(Early2~distance,data=Exp,span=0.1)
model6 = loess(Late2~distance,data=Exp,span=0.1)
model7 = loess(cell8~distance,data=Exp,span=0.1)
plot(Exp$distance,model1$fitted,lty=1,lwd=4,type="l",col=2,ylim=c(-9,-5),yaxs="i",xlim=c(5,8),xaxs="i",xlab = "log10(Genomic distance(bp))",ylab = "log10(P(s))",bty="n")
lines(Exp$distance,model2$fitted,lty=1,lwd=4,col=3)
lines(Exp$distance,model3$fitted,lty=1,lwd=4,col=4)
lines(Exp$distance,model4$fitted,lty=1,lwd=4,col=5)
lines(Exp$distance,model5$fitted,lty=1,lwd=4,col=6)
lines(Exp$distance,model6$fitted,lty=1,lwd=4,col=7)
lines(Exp$distance,model7$fitted,lty=1,lwd=4,col=8)
#formula: y = ks^(L)#
L = -1
s = Exp$distance #log10(T$distance)
fit = lm(I(Exp$FGSC-L*s)~1)
b = fit$coefficients
y = L*s+b
lines(Exp$distance,y,lty=4,lwd=3,col=1)
L = -0.5
s = Exp$distance #log10(T$distance)
fit = lm(I(Exp$FGSC-L*s)~1)
b = fit$coefficients
y = L*s+b
lines(Exp$distance,y,lty=4,lwd=3,col=1)
legend("topright",legend=c("FGSC","GV","MII","Zygote","Early 2-Cell","Late 2-Cell","8-Cell","S^(-1)","S^(-0.5)"),col=c(2:8,1,1),lty=c(rep(1,7),4,4),bg="gray90",bty="n",merge=TRUE)
box(lwd=3)
