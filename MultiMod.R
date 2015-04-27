
options <- commandArgs(trailingOnly = TRUE)
print(options)
File=options[1] # file name
Zero=options[2] # Whether exclude 0 when looking for modes
LOD=as.numeric(options[3]) # lower limit of detection
Norm=options[4] # whether perform normalization; if "T" is specified, median-by-ratio normalization will be performed.
if(length(options)<2)Zero="T"
if(length(options)<3)LOD=0
if(length(options)<4)Norm="F"

#Text=T
#n=5
#File="PCA_example.csv"

# csv or txt
tmp=strsplit(File, split="\\.")[[1]]
FileType=tmp[length(tmp)]

if(FileType=="csv"){
	cat("\n Read in csv file \n")
	prefix=strsplit(File,split="\\.csv")[[1]][1]
	In=read.csv(File,stringsAsFactors=F,row.names=1)
}
if(FileType!="csv"){
	cat("\n Read in tab delimited file \n")
	prefix=strsplit(File,split=paste0("\\.",FileType))[[1]][1]
	In=read.table(File,stringsAsFactors=F,row.names=1)
}



Matraw=data.matrix(In)

Max=apply(Matraw,1,max)
NumNo0=apply(Matraw,1,function(i)length(which(i>0)))
WhichRM=union(which(Max<LOD),which(NumNo0<5))
print(paste(length(WhichRM),"genes with less than 5 expressed cells and genes with max expression < ", LOD, "are removed"))

Mat=Matraw
if(length(WhichRM)>0)Mat=Matraw[-WhichRM,]
print(str(Mat))

if(Norm=="T"){
cat("\n ==== Performing normalization ==== \n")
library(EBSeq)
Sizes=MedianNorm(Mat)
if(is.na(Sizes[1]))cat("\n Warning: all genes have 0(s), normalization is not performed \n")
else Mat=GetNormalizedMat(Mat, MedianNorm(Mat))
}

#Log2(data+1)
MatLog=log2(Mat+1)
library(mclust)
if(Zero=="T")Gs=apply(MatLog,1,function(i)Mclust(i)$G)
if(Zero=="F")Gs=apply(MatLog,1,function(i){
											j=i[i>1]
											Mclust(j)$G})

names(Gs)=rownames(Mat)
str(Gs)

More1=which(Gs>1)
GMore=Gs[More1]
MoreNames=rownames(Mat)[More1]
names(GMore)=MoreNames
GMore=sort(GMore,decreasing=T)

GOut=matrix(Gs,ncol=1);rownames(GOut)=rownames(Gs);colnames(GOut)="Num_clusters"
GOutM=matrix(GMore,ncol=1);rownames(GOutM)=rownames(GMore);colnames(GOutM)="Num_clusters"
write.csv(GOut,file=paste0(prefix,"_Zero_",Zero,"_Num_clusters.csv"))
write.csv(GOutM,file=paste0(prefix,"_Zero_",Zero,"_Multi_mod.csv"))




