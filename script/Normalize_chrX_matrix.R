#!usr/bin/env/Rscript
#Normalize chrx matrix#
#!/usr/bin/env/Rescript
library("limma")
ref = read.table("genome_20kb.bed")
SSC = read.table("SSC_merge_20000_iced.matrix")
FGSC = read.table("FGSC_Li_merge_20000_iced.matrix")
ref = ref[ref$V1 =="chrX",]
start = ref[1,4]
end = ref[nrow(ref),4]
fgsc = FGSC[FGSC$V1>=start & FGSC$V2<=end,]
ssc = SSC[SSC$V1>=start & SSC$V2<=end,]
Exp = merge(fgsc,ssc,by=c("V1","V2"),all=TRUE)
colnames(Exp) = c("rowbin","colbin","FGSC","SSC")
data = Exp[,c(3,4)]
norm = normalizeQuantiles(data)
ssc_xa = norm$SSC/2
fgsc_xi = norm$FGSC-ssc_xa
Xi_norm = cbind(Exp[,1:2],fgsc_xi)
Xa_norm = cbind(Exp[,1:2],ssc_xa)
write.table(Xi_norm,"Xi_NA_norm.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
#HiTC draw pic#
library("HiTC")
FGSC = importC("FGSC_NA_chrx_norm.txt",xgi.bed = "genome_20kb.bed",allPairwise = F,rm.trans = F,lazyload = F)
SSC = importC("SSC_NA_chrx_norm.txt",xgi.bed = "genome_20kb.bed",allPairwise = F,rm.trans = F,lazyload = F)
Xi = importC("Xi_NA_norm.txt",xgi.bed = "genome_20kb.bed",allPairwise = F,rm.trans = F,lazyload = F)
fgsc = binningC(FGSC$chrXchrX,binsize=500000,step=1)
ssc = binningC(SSC$chrXchrX,binsize=500000,step=1)
xi = binningC(Xi$chrXchrX,binsize=500000,step=1)
mapC(forceSymmetric(HTClist(fgsc$chrXchrX)),col.pos=c("white","red"),trim.range=0.9)
mapC(forceSymmetric(HTClist(ssc$chrXchrX)),col.pos=c("white","red"),trim.range=0.9)
mapC(forceSymmetric(HTClist(xi$chrXchrX)),col.pos=c("white","red"),col.neg=c("blue","white"),trim.range=0.9)
f_dxz4 = extractRegion(FGSC$chrXchrX,chr="chrX",from=70000000,to=75000000)
s_dxz4 = extractRegion(SSC$chrXchrX,chr="chrX",from=70000000,to=75000000)
xi_dxz4 = extractRegion(Xi$chrXchrX,chr="chrX",from=70000000,to=75000000)
f_xist = extractRegion(FGSC$chrXchrX,chr="chrX",from=100000000,to=102000000)
s_xist = extractRegion(SSC$chrXchrX,chr="chrX",from=100000000,to=102000000)
xi_xist = extractRegion(Xi$chrXchrX,chr="chrX",from=100000000,to=102000000)
mapC(forceSymmetric(HTClist(xi_dxz4$chrXchrX)),col.pos=c("white","red"),col.neg=c("blue","white"),trim.range=0.9)
mapC(forceSymmetric(HTClist(xi_xist$chrXchrX)),col.pos=c("white","red"),col.neg=c("blue","white"),trim.range=0.9)