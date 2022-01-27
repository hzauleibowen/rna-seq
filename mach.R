rm(list = ls())
a=read.table("jobs.txt",header=F,na.strings = c("NA"),sep="\t")
colnames(a)="Sample Id "
library(readxl)
b=read_xlsx("Duroc(RNA-seq)(1).xlsx")
A<-as.vector(a$`Sample Id`)
B<-as.vector(b$`Sample Id`)
setdiff(A,B)

