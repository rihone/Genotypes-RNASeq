#####################Extracted SNPs from the paper#######################
#~/riho/Processed_data/Genotypes
genotypes<-read.table("SNP-inds_open3_res.txt",header=T,sep="\t",row.names=1)
maf<-as.data.frame(matrix(rep(NA,length(genotypes[,1])*7),ncol=7))
rownames(maf)<-rownames(genotypes)
colnames(maf)<-c(0,1,2,"#REF","#ALT","MA","MAF(%)")
al<-c("REF","ALT")
sam<-length(genotypes[1,])*2
for (i in 1:length(genotypes[,1])){
   maf[i,1]<-length(which(genotypes[i,]==0))
   maf[i,2]<-length(which(genotypes[i,]==1))
   maf[i,3]<-length(which(genotypes[i,]==2))
   maf[i,4]<-2*maf[i,1]+maf[i,2]
   maf[i,5]<-2*maf[i,3]+maf[i,2]
   if (maf[i,4]==maf[i,5]){
   maf[i,6]<-"REF/ALT"
   maf[i,7]<-(maf[i,4]/sam) *100
   }else{
   maf[i,6]<-al[which(maf[i,4:5] == min(maf[i,4],maf[i,5]))]
   maf[i,7]<-(min(maf[i,4],maf[i,5])/sam) *100
   }
}

######################maf###########################
                 0  1  2 #REF #ALT      MA    MAF(%)
1_109634341_C_T 54 15  0  123   15     ALT 10.869565
1_109782037_A_G 38 24  7  100   38     ALT 27.536232
1_109880721_T_C 40 25  4  105   33     ALT 23.913043
1_110119732_C_T 49 18  2  116   22     ALT 15.942029
1_110733794_A_G  9 33 27   51   87     REF 36.956522
1_110950295_C_A 63  6  0  132    6     ALT  4.347826
1_111812709_G_A 19 31 19   69   69 REF/ALT 50.000000
1_117529458_G_A 27 29 13   83   55     ALT 39.855072
1_117603156_G_T 63  6  0  132    6     ALT  4.347826
1_11852300_C_T  33 32  4   98   40     ALT 28.985507
#####################################################

write.table(maf,"SNP-inds_open3_MAF.txt",sep="\t",quote=F)

maf<-read.table("SNP-inds_open3_MAF.txt",header=T,sep="\t",row.names=1)
genotypes<-read.table("SNP-inds_open3_res.txt",header=T,sep="\t",row.names=1)
cnames<-colnames(genotypes)
writeLines(cnames,"SNP-inds_open3_noRare.txt",sep="\t")
i<-1
for (i in 1:length(genotypes[,1])){
   if (maf[i,7]>5){
      print(maf[i,7])
      print(length(genotypes[i,]))
#      write.table(genotypes[i,],file="SNP-inds_open3_noRare.txt",sep="\t",append=T,col.names=F, row.names=T, quote=F)
   }
}

gene1<-read.table("SNP-inds_open3_noRare.txt",header=T,sep="\t")

#####################All SNPs#######################
#~/riho/Processed_data/Genotypes
genotypes<-read.table("SNP-inds_open4_res.txt",header=T,sep="\t",row.names=1)
maf<-as.data.frame(matrix(rep(NA,length(genotypes[,1])*7),ncol=7))
rownames(maf)<-rownames(genotypes)
colnames(maf)<-c(0,1,2,"#REF","#ALT","MA","MAF(%)")
al<-c("REF","ALT")
sam<-length(genotypes[1,])*2
for (i in 1:length(genotypes[,1])){
   maf[i,1]<-length(which(genotypes[i,]==0))
   maf[i,2]<-length(which(genotypes[i,]==1))
   maf[i,3]<-length(which(genotypes[i,]==2))
   maf[i,4]<-2*maf[i,1]+maf[i,2]
   maf[i,5]<-2*maf[i,3]+maf[i,2]
   if (maf[i,4]==maf[i,5]){
   maf[i,6]<-"REF/ALT"
   maf[i,7]<-(maf[i,4]/sam) *100
  }else{
   maf[i,6]<-al[which(maf[i,4:5] == min(maf[i,4],maf[i,5]))]
   maf[i,7]<-(min(maf[i,4],maf[i,5])/sam) *100
   }
}

write.table(maf,"SNP-inds_open4_MAF.txt",sep="\t",quote=F)

maf<-read.table("SNP-inds_open4_MAF.txt",header=T,sep="\t",row.names=1)
genotypes<-read.table("SNP-inds_open4_res.txt",header=T,sep="\t",row.names=1)
cnames<-colnames(genotypes)
writeLines(cnames,"SNP-inds_open4_noRare.txt",sep="\t")
i<-1
for (i in 31:length(genotypes[,1])){
   if (maf[i,7]>5){
      write.table(genotypes[i,],file="SNP-inds_open4_noRare.txt",sep="\t",append=T,col.names=F, row.names=T, quote=F)
   }
}


gene<-read.table("SNP-inds_open4_noRare.txt",header=T,sep="\t",row.names=1)


