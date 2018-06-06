#####Find_Matching_Samples#############
setwd("~/riho/Processed_data")
rid<-read.table("RNASeq_open_list.txt")
gid<-read.table("Genotype_open_list.txt")
for (i in 1:length(rid[,1])){
   index<-charmatch(rid[i,1],gid[,1])
   if(is.na(gid[index,1])==F){
   write.table(gid[index,1], file="match_open.txt",append=T,col.names=F, row.names=F, quote=F)
   }
}

#####################RNASeq#################
#Rscript
#Find unique target_id
#setwd("~/riho/Processed_data/RNASeq/open")
name<-list.files("open")
infile<-read.table(name[1],header=T)
for( i in 2:length(name)){
  infile2<-read.table(name[i],header=T)
  infile <- rbind(infile, infile2)
  infile <- infile[!duplicated(infile$target_id),]
}
#Choose extracted sample files
mat<-scan("match_open.txt",what = character(), sep = "\n")
files<-list.files("open") #~/RNASeq
for (i in 1:length(mat)){
   line<-charmatch(mat[i],files)
   write.table(files[line], file="match_open_files.txt",append=T,col.names=F, row.names=F, quote=F)
}
#Make a table
mat<-sub("-",".",mat)
ID<-infile[,1]
tab<-matrix(rep(NA,length(ID)*length(mat)),ncol=length(mat))
colnames(tab)<-mat
rownames(tab)<-ID
files<-scan("match_open_files.txt",what=character(),sep="\n")
setwd("open")
for (i in 1:length(files)){
  print(i)
  line<-read.table(files[i],header=T)
  for(j in 1:length(line[,1])){
    index<-match(line[j,1],ID)
    if(is.na(index)== F){
      tab[index,i] <- line[j,4]
    }
  } 
}
setwd("../")
write.table(tab,file="RNA-celline_extracted(open).txt",quote=F)

###############Remove missing values#############
rna<-read.table("RNA-celline_extracted(open).txt",header=T,sep="\t")
dl<-0
for (i in 1:length(rna[,1])){
  if (sum(rna[i,])== 0) {
     dl<-append(-i,dl,after=dl[1])
  }
}

rna2<-rna[dl[2:length(dl)],]
write.table(rna2,"RNA-celline_extracted(open)_no0.txt",sep="\t",quote=F)

#############Extract Individuals#################
#~riho/Processed_data/RNASeq/eQTL_open"
rna<-read.table("RNA-celline_extracted(open)_no0.txt",header=T,sep="\t",row.names=1)
cnames<-colnames(rna)
np<-strsplit(cnames,"_")
sep<-rep("A",length(np))
for (i in 1:length(np)){
   sep[i]<-np[[i]][1]
}

dup<-sep[duplicated(sep)]
ndup<-sep[!duplicated(sep)]
for (i in 1:length(dup)){
   index<-grep(dup[i],cnames)
   col1<-rna[,index[1]]
   col2<-rna[,index[2]]
   print(dup[i])
   print(cor(col1,col2))
}

al<-length(unique(sep))
dat<-matrix(rep(NA,length(rna[,1])*al),ncol=al)
rownames(dat)<-rownames(rna)
colnames(dat)<-unique(sep)
cnames<-colnames(rna)
rnames<-rownames(rna)
dup<-sep[duplicated(sep)]
ndup<-sep[!duplicated(sep)]
sdup<-c(3,4,5,6,9,10,11,12,15,16,19,20,22,23,27,28,31,32,34,35,37,38,43,44,58,59)
id<-1:82
id<-id[-sdup] # columns of no duplicated columns
idx<-c(3,5,9,12,16,20,22,27,31,34,38,44,58) #extracted samples based on the "count" table
#"HPSI0114i.lexy_1","HPSI0115i.qoog_4","HPSI0115i.vazt_2","HPSI0214i.feec_3","HPSI0214i.pelm_3","HPSI0414i.oikd_2","HPSI0414i.seru_1","HPSI0514i.lako_2","HPSI0614i.lepk_4","HPSI0914i.laey_4

idx2<-sort(append(id,idx)) #all single columns
dnames<-cnames[idx2]
ex<-strsplit(dnames,"_")
exsep<-rep("A",length(ex))
for (i in 1:length(ex)){
   exsep[i]<-ex[[i]][1]
}

#unique individuals
for (i in 1:length(idx2)){
  dat[rnames,exsep[i]]<-rna[rnames,idx2[i]]
}

write.table(dat,"RNASeq-inds_open_res.txt",sep="\t",quote=F)

#############Genotypes########################
#####Choose_Matching_Sample_files#############
mat<-scan("match_open.txt",what = character(), sep = "\n")
files<-list.files("open_genlight") #~/Genotypes
for (i in 1:length(mat)){
   line<-charmatch(mat[i],files,value=T)
   write.table(line, file="match_open_files3.txt",append=T,col.names=F, row.names=F, quote=F)
}

###########Concatenate SNP Attributes#############
filenames <- list.files(path = ".", pattern = ".snpids$", recursive = TRUE)
library(dplyr)
currentfile <- filenames[1]
snps <- read.table(currentfile, header = TRUE)
snps <-snps %>% as.data.frame() %>% mutate(UNI = paste(!!!rlang::syms(c("CHROM", "POS", "REF", "ALT")), sep="_"))
for (i in 2:length(filenames)) {
  currentfile <- filenames[i]
  newfile<-read.table(currentfile, header = TRUE)
  newfile <- newfile %>% as.data.frame() %>% mutate(UNI = paste(!!!rlang::syms(c("CHROM", "POS", "REF", "ALT")), sep="_"))
  snps <- rbind(snps, newfile)
  snps <- snps[!duplicated(snps$UNI),]
}
snpmatrix_out <- as.matrix(snps[,c(6,1,2,3,4,5)])
write.table(snpmatrix_out, file = "all_snps_attribute_open2.txt", quote = FALSE, row.names = FALSE)
#unix
sort -k2,2n -k3,3 all_snps_attribute_open2.txt > all_snps_attribute_open3.txt 

#################Extract_SNPS_from_Paper#################
library(dplyr)
hip<-read.table("all_snps_attribute_open3.txt",header=T)
hipdat<-as.data.frame(hip)
tb<-read.csv("nature22403-s5.csv")
for (i in 1:length(tb[,1])){
  chr<-tb[i,1]
  pos<-tb[i,2]
  nchr<-filter(hipdat,CHROM==chr,POS==pos)
  write.table(nchr, file="extractedSNPs3.txt",append=T,col.names=F, row.names=F, quote=F)
}

#####################Make_the_Table#####################
name<-read.table("match_open.txt")
name<-sort(name[,1])
extra<-read.table("extractedSNPs3.txt",header=T)
ex<-extra[!duplicated(extra$UNI),]
tab<-matrix(rep(NA,length(name)*length(ex[,1])),ncol=length(name))
file<-scan("match_open_files3.txt", what = character(), sep = "\n", blank.lines.skip = F)
colnames(tab)<-name
rownames(tab)<-ex[,1]
setwd("open")
for (i in 1:length(file)){
  print(i)
  infile<-read.table(file[i],header=T)
  for(j in 1:length(infile[,1])){
    index<-match(infile[j,1],ex[,1])
    if(is.na(index)== F){
      tab[index,i] <- infile[j,2]
    }
  } 
}

setwd("../")
write.table(tab,"SNP-cellines_open3.txt",quote=F)

#all SNPs
tab<-matrix(rep(NA,length(name)*length(hip[,1])),ncol=length(name))
file<-scan("match_open_files_non.txt",what=character(),sep="\n")
colnames(tab)<-name
rownames(tab)<-hip[,1]
setwd("open_genlight")
for (i in 1:length(file)){
  print(i)
  infile<-read.table(file[i],header=T)
  for(j in 1:length(infile[,1])){
    index<-match(infile[j,1],hip[,1])
    if(is.na(index)== F){
      tab[index,i] <- infile[j,2]
    }
  } 
}

#all SNPs for extracted samples
mat<-scan("match_open.txt", what = character(), sep = "\n", blank.lines.skip = F)
files<-list.files("snpids") #~/Genotypes
for (i in 1:length(mat)){
   line<-charmatch(mat[i],files)
   write.table(files[line], file="match_open_files4.txt",append=T,col.names=F, row.names=F, quote=F)
}
filenames <- scan("match_open_files4.txt", what = character(), sep = "\n")

# Concatenate SNP Attributes
library(dplyr)
currentfile <- filenames[1]
snps <- read.table(currentfile, header = TRUE)
snps <-snps %>% as.data.frame() %>% mutate(UNI = paste(!!!rlang::syms(c("CHROM", "POS", "REF", "ALT")), sep="_"))
for (i in 2:length(filenames)) {
  currentfile <- filenames[i]
  newfile<-read.table(currentfile, header = TRUE)
  newfile <- newfile %>% as.data.frame() %>% mutate(UNI = paste(!!!rlang::syms(c("CHROM", "POS", "REF", "ALT")), sep="_"))
  snps <- rbind(snps, newfile)
  snps <- snps[!duplicated(snps$UNI),]
}
snpmatrix_out <- as.matrix(snps[,c(6,1,2,3,4,5)])
write.table(snpmatrix_out, file = "match_snps_attribute_open.txt", quote = FALSE, row.names = FALSE)

#unix
sort -k2,2n -k3,3n match_snps_attribute_open.txt > match_snps_attribute_open2.txt 

name<-scan("match_open.txt",what=character(),sep="\n")
hip<-read.table("match_snps_attribute_open2.txt",header=T)
tab<-matrix(rep(NA,length(name)*length(hip[,1])),ncol=length(name))
file<-scan("match_open_files3.txt",what=character(),sep="\n")
colnames(tab)<-name
rownames(tab)<-hip[,1]
i<-1
setwd("open_genlight")
for (i in 1:length(file)){
  print(i)
  infile<-as.matrix(read.table(file[i],header=T,row.names=1))
  ID<-rownames(infile)
  tab[ID,i] <- infile[ID,1] 
}
setwd("../")
write.table(tab,"SNP-cellines_open4.txt",quote=F)

###############Remove missing values#############
gene<-read.table("SNP-cellines_open3.txt",header=T,sep=" ")
gene2<-read.table("SNP-cellines_open4.txt",header=T,sep=" ")

library(dplyr)
n<-gene2 %>% as.data.frame() %>% filter(complete.cases(gene2)==F) %>% select(id)
cha<-as.character(n[[1]])
dl<-0
for (i in 1:length(cha)){
  num<-as.numeric(rownames(gene2[gene2$id==cha[i],]))
  dl<-append(-num,dl,after=dl[1])
}

res2<-gene2[dl[2:length(dl)],]
write.table(res,"SNP-cellines_open3_res.txt",sep="\t",quote=F,row.names=F)
write.table(res2,"SNP-cellines_open4_res.txt",sep="\t",quote=F,row.names=F)

#############Extract Individuals#################
gene<-read.table("SNP-cellines_open3_res.txt",header=T,sep="\t",row.names=1)
cnames<-colnames(gene)
np<-strsplit(cnames,"_")
sep<-rep("A",length(np))
for (i in 1:length(np)){
   sep[i]<-np[[i]][1]
}

dup<-sep[duplicated(sep)]
ndup<-sep[!duplicated(sep)]
for (i in 1:length(dup)){
   index<-grep(dup[i],cnames)
   print(length(index))
   col1<-gene[,index[1]]
   col2<-gene[,index[2]]
   print(dup[i])
   print(sum(col1-col2)==0)
   
}

#open3 
"HPSI0115i.vazt"
11	114263860	rs3887314	A	G	.	.	.	GT:LRR:BAF:IA:IB:GC 1/1:-0.4592:0.9224:0.079:0.552:0.8953
11	114263860	rs3887314	A	G	.	.	.	GT:LRR:BAF:IA:IB:GC 0/1:-0.098:0.5045:0.401:0.472:0.9407
http://www.hipsci.org/lines/#/lines/HPSI0115i-vazt_1
http://www.hipsci.org/lines/#/lines/HPSI0115i-vazt_2

dn <- rownames(gene[gene[,11]!=gene[,12],11:12])
idx <- grep(dn,rownames(gene))

snps<-gene[,11]
col1<-11
col2<-12
ridx <- which(gene[,col1]!=gene[,col2]) #rownum of different SNP
count<-matrix(rep(0,5*length(ridx)),ncol=5) #1="0",2="1",3="2"
rownames(count)<-ridx
colnames(count)<-c("al1","al2","num1","num2","max")
for (i in 1:length(ridx)){
  count[i,"al1"]<-gene[ridx[i],col1]
  count[i,"al2"]<-gene[ridx[i],col2]
  count[i,"num1"]<-length(which(gene[ridx[i],]==gene[ridx[i],col1]))
  count[i,"num2"]<-length(which(gene[ridx[i],]==gene[ridx[i],col2]))
  count[i,"max"]<-count[i,which.max(count[i,3:4])]
}

#########count###########
   al1 al2 num1 num2 max
69   2   1   11   41   1
#########################

snps[ridx]<-count[1,"max"]

al<-length(unique(sep))
dat<-matrix(rep(NA,length(gene[,1])*al),ncol=al)
rownames(dat)<-rownames(gene)
colnames(dat)<-unique(sep)
cnames<-colnames(gene)
rnames<-rownames(gene)
dup<-sep[duplicated(sep)]
ndup<-sep[!duplicated(sep)]
sdup<-c(3,4,5,6,9,10,11,12,15,16,19,20,22,23,27,28,31,32,34,35,37,38,43,44,58,59)
id<-1:82
id<-id[-sdup] # columns of no duplicated columns
dcol<- c(3,4,5,6,9,10,15,16,19,20,22,23,27,28,31,32,34,35,37,38,43,44,58,59) #duplicated columns of no different genotypes
idx<-c(3,5,9,11,15,19,22,27,31,34,37,43,58)
idx2<-sort(append(id,idx)) #all single columns
idx3<- c(3,5,9,15,19,22,27,31,34,37,43,58) #columns of no different genotypes
idx4<-sort(append(id,idx3))#without the column of diferrent genotype

dnames<-cnames[idx4]
ex<-strsplit(dnames,"_")
exsep<-rep("A",length(ex))
for (i in 1:length(ex)){
   exsep[i]<-ex[[i]][1]
}

#unique individuals
for (i in 1:length(idx4)){
  dat[rnames,exsep[i]]<-gene[rnames,idx4[i]]
}

dat2<-dat
#dat2[,8] = HPSI0115i.vazt
dat2[,8]<-snps

write.table(dat2,"SNP-inds_open3_res.txt",sep="\t",quote=F)

#open4
gene<-read.table("SNP-cellines_open4_res.txt",header=T,sep="\t",row.names=1)
cnames<-colnames(gene)
np<-strsplit(cnames,"_")
sep<-rep("A",length(np))
for (i in 1:length(np)){
   sep[i]<-np[[i]][1]
}

dup<-sep[duplicated(sep)]
ndup<-sep[!duplicated(sep)]
lst<-rep("A",10)
num<-1
for (i in 1:length(dup)){
   index<-grep(dup[i],n)
   col1<-gene[,index[1]]
   col2<-gene[,index[2]]
   if(sum(col1-col2)!=0){
    lst[num]<-dup[i]
    num<-num+1
   }
}

               
# "HPSI0114i.lexy"
gene[gene[,5]!=gene[,6],5:6]
#unix (~michoel_grp/HipSci/Genotype/Open/ftp.sra.ebi.ac.uk/vol1/)
find . -type f | grep "HPSI0114i.lexy"
zcat HPSI0114i_lexy_1.wec.gtarray.HumanCoreExome-12_v1_0.20141111.genotypes.vcf.gz | grep "115227215"
zcat HPSI0114i-lexy_2.wec.gtarray.HumanCoreExome-12_v1_0.20141111.genotypes.vcf.gz | grep "115227215"
4	115227215	rs4597902	G	T	.	.	AC=0;AN=2	GT:LRR:BAF:IA:IB:GC	0/0:-0.2869:0:0.623:0.019:0.9306
4	115227215	rs4597902	G	T	.	.	AC=1;AN=2	GT:LRR:BAF:IA:IB:GC	1/0:-0.1091:0.484:0.44:0.41:0.9306

16	84613405	rs2967846	T	C	.	.	AC=1;AN=2	GT:LRR:BAF:IA:IB:GC	0/1:0.0489:0.6619:0.54:1.169:0.7591
16	84613405	rs2967846	T	C	.	.	AC=2;AN=2	GT:LRR:BAF:IA:IB:GC	1/1:-0.0032:0.996:0.027:1.375:0.7809

21	46313371	exm1578770	G	A	.	.	AC=1;AN=2	GT:LRR:BAF:IA:IB:GC	1/0:0.2143:0.352:0.659:0.427:0.2825
21	46313371	exm1578770	G	A	.	.	AC=0;AN=2	GT:LRR:BAF:IA:IB:GC	0/0:-0.0113:0.0013:0.882:0.018:0.9072

###############count###################
                 al1 al2 num1 num2 max
4_115227215_G_T   1   0   46   16   1
16_84613405_T_C   1   2   27   10   1
21_46313371_G_A   1   0    1   81   0
########################################

snps<-gene[,c(5,9,12,15,19,31,34,37,43,58)]
spnames<-colnames(snps)
spn<-strsplit(spnames,"_")
spcol<-rep("A",length(spn))
for (i in 1:length(spn)){
   spcol[i]<-spn[[i]][1]
}

colnames(snps)<-spcol

col1<-58
col2<-59
nm<-rownames(gene[gene[,col1]!=gene[,col2],])
count<-as.data.frame(matrix(rep(0,6*length(nm)),ncol=6)) #1="0",2="1",3="2"
rownames(count)<-nm
colnames(count)<-c("al1","al2","num1","num2","max","ale")
al<-c("al1","al2")
for (i in 1:length(nm)){
  count[i,"al1"]<-gene[nm[i],col1]
  count[i,"al2"]<-gene[nm[i],col2]
  count[i,"num1"]<-length(which(gene[nm[i],]==gene[nm[i],col1]))
  count[i,"num2"]<-length(which(gene[nm[i],]==gene[nm[i],col2]))
  count[i,"max"]<-count[i,which.max(count[i,3:4])]
  count[i,"ale"]<-al[which(count[i,1:2] == count[i,"max"])]
}

i<-10
snps[nm,i]<-count[nm,"max"]
al<-length(unique(sep))
dat<-matrix(rep(NA,length(gene[,1])*al),ncol=al)
rownames(dat)<-rownames(gene)
colnames(dat)<-unique(sep)
cnames<-colnames(gene)
rnames<-rownames(gene)
dup<-sep[duplicated(sep)]
ndup<-sep[!duplicated(sep)]
sdup<-c(3,4,5,6,9,10,11,12,15,16,19,20,22,23,27,28,31,32,34,35,37,38,43,44,58,59)
id<-1:82
id<-id[-sdup] # columns of no duplicated columns
idx<-c(3,5,9,12,16,20,22,27,31,34,38,44,58) #extracted columnes based on the allele freq (all)
idx2<-sort(append(id,idx)) #all single columns
idx3<-c(3,22,27) #duplicated columns of no different genotypes
idx4<-sort(append(id,idx3)) #all columns without the column of diferrent genotype
idx5<-c(5,9,12,16,20,31,34,38,44,58) #extracted columnes based on the allele freq (dup)

dnames<-cnames[idx4]
ex<-strsplit(dnames,"_")
exsep<-rep("A",length(ex))
for (i in 1:length(ex)){
   exsep[i]<-ex[[i]][1]
}

#unique individuals
for (i in 1:length(idx4)){
  dat[rnames,exsep[i]]<-gene[rnames,idx4[i]]
}

dat2<-dat
for (i in 1:length(spcol)){
   dat2[rnames,spcol[i]]<-gene[rnames,idx5[i]]
}

write.table(dat2,"SNP-inds_open4_res.txt",sep="\t",quote=F)

#"HPSI0115i.qoog"
gene[gene[,9]!=gene[,10],9:10]
                al1 al2 num1 num2 max ale
4_182226071_A_G   1   2   25    5   1 al1
5_160047528_C_A   0   1   81    1   0 al1
7_34607206_A_C    1   0   22   59   0 al2
9_8504321_G_A     0   1   81    1   0 al1
10_3738103_G_A    1   0   24   53   0 al2
10_75783340_T_C   1   0   38   32   1 al1
10_75868114_C_T   1   0   36    9   1 al1
10_75881347_G_T   1   0   39   23   1 al1
10_75900462_C_A   1   0   36    9   1 al1
10_75902177_G_A   1   0   36   23   1 al1
10_75946564_T_C   1   0   34   22   1 al1
10_75948420_C_A   1   0   34    9   1 al1
10_76071203_G_A   1   0   36    9   1 al1
10_76092691_A_G   1   0   35   22   1 al1
10_76131244_C_T   1   0   36    9   1 al1
10_76180335_T_C   1   0   36   23   1 al1
10_76214709_C_A   1   0   36    9   1 al1
10_76221239_T_C   1   0   36    9   1 al1
10_76225597_G_A   1   0   36    9   1 al1
10_76273579_G_A   1   0   34    9   1 al1
10_76295789_A_G   1   0   34   36   0 al2
10_76323150_C_A   1   2   26    6   1 al1
13_99049254_A_G   1   2   34   43   2 al2
13_99051172_C_A   1   2   38   28   1 al1
13_99059793_C_T   1   2   39   33   1 al1

#"HPSI0115i.vazt"
gene[gene[,11]!=gene[,12],11:12]
                 al1 al2 num1 num2 max ale
3_66287056_G_A     2   1    9   34   1 al2
3_66329205_C_T     2   1    2   35   1 al2
10_54799931_C_T    1   0   27   53   0 al2
10_96229166_C_T    2   1   32   42   1 al2
11_113570385_C_T   0   1   58   22   0 al1
11_113578677_G_A   2   1    1   18   1 al2
11_113684809_A_G   0   1   42   31   0 al1
11_113710905_T_C   0   1   26   51   1 al2
11_113773946_C_T   2   1   11   51   1 al2
11_113803028_A_C   2   1    3   44   1 al2
11_113812733_G_A   2   1   25   48   1 al2
11_113815354_T_G   2   1    1    1   2 al1
11_113833152_A_G   0   1   16   50   1 al2
11_113875853_T_C   0   1    6   45   1 al2
11_113890125_T_C   2   1    7   39   1 al2
11_113957197_G_A   0   1   21   46   1 al2
11_113969561_A_G   0   1    7   34   1 al2
11_113986162_G_A   2   1   14   32   1 al2
11_113997675_C_T   2   1   15   34   1 al2
11_114000727_C_T   2   1    9   38   1 al2
11_114005213_C_T   2   1    8   36   1 al2
11_114053768_C_T   2   1   21   28   1 al2
11_114063362_G_T   2   1   30   28   2 al1
11_114072875_C_T   2   1    6   27   1 al2
11_114074681_G_A   2   1    9   34   1 al2
11_114074756_C_T   2   1   12   36   1 al2
11_114077940_C_T   0   1    7   29   1 al2
11_114079713_T_C   0   1   28   32   1 al2
11_114101709_T_C   0   1    9   43   1 al2
11_114105627_C_T   0   1    8   35   1 al2
11_114107173_G_T   0   1   37   35   0 al1
11_114115298_A_G   2   1    2   20   1 al2
11_114120774_A_G   0   1   33   34   1 al2
11_114167939_G_C   2   1   37   30   2 al1
11_114179701_T_C   0   1    2   32   1 al2
11_114186645_T_C   2   1    2   22   1 al2
11_114189367_C_T   2   1   31   37   1 al2
11_114191016_A_G   2   1   31   37   1 al2
11_114199213_A_G   0   1   24   44   1 al2
11_114214237_T_C   0   1   14   52   1 al2
11_114233785_A_G   2   1   30   36   1 al2
11_114235304_A_G   2   1   16   49   1 al2
11_114241636_G_T   2   1   14   46   1 al2
11_114263860_A_G   2   1   11   41   1 al2
11_114308684_C_T   2   1   20   42   1 al2
11_114319760_A_G   2   1   10   39   1 al2
11_114324026_A_G   2   1   17   49   1 al2
11_114333457_A_G   2   1    9   28   1 al2
11_114393652_C_T   2   1   11   36   1 al2
11_114397499_T_C   2   1    9   24   1 al2
11_114413446_G_T   2   1   11   36   1 al2
11_114436944_A_G   2   1   16   43   1 al2
11_114442103_A_G   2   1   16   43   1 al2
11_114471076_C_T   0   1   53   27   0 al1
11_114480735_T_C   2   1    2   19   1 al2
11_114512627_A_C   2   1    1    8   1 al2
11_114535196_C_T   0   1    1    8   1 al2
11_114566047_T_C   0   1    5   25   1 al2
11_114604149_T_G   2   1    1   12   1 al2
11_114622352_C_T   0   1   31   42   1 al2
11_114644778_T_C   2   1   12   33   1 al2
11_114649770_C_T   0   1   34   30   0 al1
11_114702663_A_G   0   1   47   26   0 al1
11_114727155_T_C   0   1   44   24   0 al1
11_114738157_G_A   2   1   24   35   1 al2
11_114781237_G_A   0   1   20   34   1 al2
11_114803157_A_G   0   1   14   31   1 al2
11_114811893_T_C   2   1   14   32   1 al2
11_114839858_T_C   2   1   27   44   1 al2
11_114945485_T_G   0   1    4   27   1 al2
11_114946119_C_T   0   1   27   42   1 al2
11_114958994_C_T   0   1   21   45   1 al2
11_114961537_T_C   0   1   21   47   1 al2
11_115022404_A_G   2   1   18   47   1 al2
11_115032240_A_G   2   1   44   30   2 al1
11_115043574_T_G   0   1   47   29   0 al1
11_115043888_C_A   2   1   23   39   1 al2
11_115045237_T_C   2   1   25   36   1 al2
11_115081563_T_C   0   1   22   40   1 al2
11_115083875_G_A   0   1   22   44   1 al2
11_115111711_G_A   0   1   51   28   0 al1
11_115189228_T_C   2   1    7   27   1 al2
11_115212725_C_T   2   1   50   31   2 al1
11_115215743_A_C   2   1   32   42   1 al2
11_115220254_T_G   2   1    6   34   1 al2
11_115268121_A_G   2   1   17   30   1 al2
11_115283535_C_T   2   1   55   24   2 al1
11_115285376_G_A   2   1   17   33   1 al2
11_115301871_C_T   2   1   52   27   2 al1
11_115304363_A_G   2   1   17   33   1 al2
11_115308919_C_T   0   1   54   28   0 al1
11_115397017_A_G   2   1    7   37   1 al2
11_115412993_T_C   2   1   36   39   1 al2
11_115463786_G_A   2   1   39   36   2 al1
11_115469253_C_T   2   1    3   37   1 al2
11_115498833_C_T   0   1    5   19   1 al2
11_115501343_C_T   0   1   13   35   1 al2
11_115515326_C_T   0   1   11   35   1 al2
11_115523910_A_G   0   1    4   14   1 al2
11_115542933_T_C   0   1    2   21   1 al2
11_115570540_G_A   0   1   14   41   1 al2
11_115598377_T_G   0   1   60   21   0 al1
11_115618087_G_T   0   1   57   25   0 al1
11_115685475_A_G   2   1   35   32   2 al1
11_115743487_C_T   0   1   19   41   1 al2
11_115753077_A_C   0   1   24   40   1 al2
11_115757874_C_T   0   1   50   25   0 al1
11_115764045_C_T   0   1   34   32   0 al1
11_115768256_C_A   0   1   32   33   1 al2
11_115774774_G_A   0   1   49   29   0 al1
11_115777037_G_A   0   1   31   37   1 al2
11_115964809_A_C   2   1    1   31   1 al2
11_115982197_G_A   2   1    8   29   1 al2
11_116002530_G_A   0   1   50   26   0 al1
11_116034732_G_A   0   1   57   25   0 al1
11_116042048_G_A   2   1   22   42   1 al2
11_116060076_C_T   0   1    4   27   1 al2
11_116076882_C_T   0   1    8   34   1 al2
11_116077726_C_A   0   1    8   34   1 al2
11_116108902_G_A   0   1   61   20   0 al1
11_116118912_T_G   2   1   21   46   1 al2
11_116130782_G_A   0   1   34   42   1 al2
11_116139119_G_A   2   1   31   43   1 al2
11_116184920_A_G   0   1   53   29   0 al1
11_116233857_T_C   0   1   34   45   1 al2
11_116236775_C_T   0   1   23   45   1 al2
11_116245350_A_G   2   1   49   30   2 al1
12_32189006_A_G    0   1   21   37   1 al2
12_32193344_A_G    2   1    8   33   1 al2
12_32194375_C_T    0   1   57   22   0 al1
12_32197877_G_T    0   1   54   21   0 al1
12_32204749_C_T    2   1    6   33   1 al2
12_32206560_C_T    0   1   18   40   1 al2
15_53441177_A_G    2   1   46   32   2 al1
16_81204172_A_C    1   2   42   10   1 al1

#"HPSI0214i.feec"
gene[gene[,15]!=gene[,16],15:16]
                al1 al2 num1 num2 max ale
6_17676470_A_G    2   1   18   40   1 al2
6_17744751_C_A    0   1   46   33   0 al1
6_18106076_T_G    0   1   30   44   1 al2
10_67996275_A_G   0   1   27   41   1 al2
16_78713526_G_A   2   1   29   43   1 al2
16_78715586_T_G   2   1   54   23   2 al1
16_78803121_T_C   0   1    9   53   1 al2
16_78808512_G_A   2   1    9   33   1 al2
16_78814629_G_A   2   1    8   32   1 al2
16_78817229_G_A   0   1   14   51   1 al2
16_78818322_T_C   0   1    8   34   1 al2
16_78822393_C_T   0   1   51   30   0 al1
16_78827057_T_C   0   1   44   34   0 al1
16_78840745_G_A   2   1   41   33   2 al1
16_78854487_A_G   0   1    1   18   1 al2
22_36231357_G_T   2   1   10   37   1 al2

#"HPSI0214i.pelm"
gene[gene[,19]!=gene[,20],19:20]
               al1 al2 num1 num2 max ale
2_99209314_G_A   2   1    7   36   1 al2

#"HPSI0414i.oikd"
gene[gene[,31]!=gene[,32],31:32]
               al1 al2 num1 num2 max ale
1_56458937_A_G   1   0   38   30   1 al1

#"HPSI0414i.seru"
gene[gene[,34]!=gene[,35],34:35]
                al1 al2 num1 num2 max ale
2_154278240_C_T   1   2   30    7   1 al1
2_154731055_T_C   1   2   28   23   1 al1
2_154748983_G_A   1   2   28   42   2 al2
2_154767281_T_C   1   2   28   45   2 al2
7_110661337_C_A   1   0   34    8   1 al1
7_110681482_A_C   1   0   34    8   1 al1
7_110843795_A_G   1   0   40    8   1 al1
7_110867535_A_G   1   0   34    8   1 al1
7_110879776_A_G   1   0   47   10   1 al1
7_110921994_G_T   1   2   47   10   1 al1
7_110949566_C_A   1   2   32    6   1 al1
9_109912752_A_C   1   0   31    5   1 al1
16_78722387_T_C   1   0   30   40   0 al2
21_15447212_T_G   2   1   59   20   2 al1

#"HPSI0514i.lako"
gene[gene[,37]!=gene[,38],37:38]
                 al1 al2 num1 num2 max ale
2_226243056_T_C    2   1   25   39   1 al2
5_29996418_A_C     1   2   45    7   1 al1
6_161082695_T_C    0   1    5   34   1 al2
6_161084642_C_T    0   1   11   53   1 al2
6_161111113_A_G    2   1    8   28   1 al2
6_161115245_A_T    0   1   21   46   1 al2
6_161136058_T_A    2   1   15   43   1 al2
8_432146_A_G       1   0   36   37   0 al2
12_133168145_C_T   0   1   19   35   1 al2

#"HPSI0614i.lepk"
gene[gene[,43]!=gene[,44],43:44]
                al1 al2 num1 num2 max ale
MT_1811_A_G       1   0    1   71   0 al2
X_34148918_G_A    0   1   81    1   0 al1
1_246347153_G_A   1   2   30    3   1 al1
2_71891480_G_A    0   1   81    1   0 al1
3_162266475_G_A   2   1   14   35   1 al2
4_96546107_T_C    1   2   35   28   1 al1
5_135391462_A_G   0   1   79    3   0 al1
6_10270377_G_A    0   1   16   43   1 al2
6_68838783_C_A    0   1   48   31   0 al1
7_123516944_C_T   1   0    1   81   0 al2
20_14096380_T_G   0   1   13   40   1 al2
20_14106881_T_G   2   1    3   26   1 al2
20_14111226_G_A   2   1    3   26   1 al2
20_14205048_C_T   0   1    7   37   1 al2
20_14244971_A_G   0   1    5   41   1 al2
20_14293522_T_C   0   1   32   39   1 al2
20_14295645_C_A   0   1    5   40   1 al2
20_14306953_G_T   0   1   32   39   1 al2
20_14317390_T_C   0   1    5   40   1 al2
20_14321880_A_C   2   1    2   29   1 al2
20_14327899_G_A   0   1   23   44   1 al2
20_14338493_A_G   0   1   12   39   1 al2
20_14353386_C_T   2   1    9   39   1 al2
20_14358117_T_C   2   1   27   42   1 al2
20_14363827_T_C   0   1   43   33   0 al1
20_14388215_C_T   2   1    9   43   1 al2
20_14400165_G_A   2   1   13   45   1 al2
20_14425337_A_G   0   1   47   29   0 al1
20_14466618_C_T   2   1    2   29   1 al2
20_14466997_T_G   0   1   51   26   0 al1
20_14514820_G_A   2   1   15   34   1 al2
20_14521106_C_T   2   1   10   37   1 al2
20_14524548_A_G   2   1    7   33   1 al2
20_14530220_C_A   2   1   25   41   1 al2
20_14562027_G_A   2   1   14   39   1 al2
20_14568758_G_A   2   1   35   33   2 al1
20_14576448_C_T   2   1   30   32   1 al2


#"HPSI0914i.laey"
gene[gene[,58]!=gene[,59],58:59]

#RNA expression correration 
rna<-read.table("RNA-celline_extracted(open)_no0.txt",header=T,sep="\t")
cor(rna[,5],rna[,6])
#0.9800941
cor(rna[,9],rna[,10])
#0.9818828
cor(rna[,11],rna[,12])
#0.9603363
cor(rna[,15],rna[,16])
#0.9204957
cor(rna[,19],rna[,20])
#0.9921373
cor(rna[,31],rna[,32])
#0.9700321
cor(rna[,34],rna[,35])
#0.9943587
cor(rna[,37],rna[,38])
#0.9725006
cor(rna[,43],rna[,44])
#0.9937367
cor(rna[,58],rna[,59])
#0.9502427



########Python###########
import pandas as pd
import os
#RNASeq
files = os.listdir("open")
name = os.listdir("open")
os.chdir("open")

for i,line in enumerate(files):
    stn=line.split("._")
    if stn[0] == "":
        name[i]="".join(stn)

name=list(set(name))

num = 0
for i, f in enumerate(name):
    print i
    infile = pd.read_table(f)
    for j in range(len(infile)):
            if ID[j] != infile.iloc[j,0]:
                num += 1
    if num != 0: 
            print f
            num = 0


for i, f in enumerate(name):
    print ("file=",i)
    infile = pd.read_table(f)
    stn=f.split(".")
    if i == 0:
        dat=infile.iloc[:,[0,4]]
        dat.columns = ["target_id",stn[0]]
        print ("dat =",len(dat.iloc[0,:])) 
    else:
            dat[stn[0]]=infile.iloc[:,4] 
            print ("dat=",len(dat.iloc[0,:]))       

dat.to_csv("table.tsv",sep="\t",index=False)

files = os.listdir("open")
name = os.listdir("open")
os.chdir("open")

for i,line in enumerate(files):
    stn=line.split("._")
    if stn[0] == "":
        name[i]="".join(stn)

name=list(set(name))
ID=infile.iloc[:,0]

aray=[]
for i, f in enumerate(name):
    print i
    infile = pd.read_table(f)
    aray.append(len(infile))
#name[2]=max
infile = pd.read_table(name[2],sep=" ")
stn=f.split(".")
ID=infile.iloc[:,0]
dat = infile
dat.columns = ["id",stn[0]]
ar=[]
for i, f in enumerate(name):
    print i
    infile = pd.read_table(f,sep=" ")
    stn=f.split(".")
    infile.columns = ["id",stn[0]]
    ar.append(infile.iloc[:,0])
    
ar=list(set(ar))          

f = open("SNP.txt","a")
for i in range(len(ar)):
    f.write(ar[i]+"\n")
f.close()









