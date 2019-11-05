#!usr/bin/env/Rscript
dir = getwd()
files_name = list.files(path=dir,pattern="chr*")
data = lapply(files_name,function(x){read.table(x)})
Exp = data.frame(dis=numeric(0),ESC=numeric(0))
for (k in 1:1533) {
	count = 0
	total = 0
	bins = 0
	for (m in 1:19){
		T = as.matrix(data[[m]])
		for (i in 1:(nrow(T)-k)){
       			count = count + T[i,i+k]
  		}
  		bins = bins + i
  		total = total + sum(T)/2
  	}
  	Exp[k,1] = k*40000
  	Exp[k,2] = count/bins
	print(k)
}
proba = Exp[,2]/total
Exp$proba = proba
pro = proba/sum(proba)
Exp$pro = pro
write.table(Exp,"contact_probability.txt",quote=FALSE,row.name=FALSE)
