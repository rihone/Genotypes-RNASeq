##############Open_Samples###############
library(dplyr)
#non-labeled
tab<-read.table("hipsci.proteomics.maxquant.uniprot.MQ1_MQ2.20151023.proteinGroups.txt",header=T,sep="\t")
mtx<-tab[,c(7,165:174,176:179,181:184)]
samples<-colnames(mtx[2:19])
samples<-strsplit(samples,".",fixed=T)
tmt<-c("HPSI0114i-eipl_2","HPSI0215i-oilg_1","HPSI0215i-oilg_3","HPSI0914i-suop_2","HPSI0914i-suop_5","HPSI0114i-eipl_1") #from labeled data
#labeled
tab<-read.table("hipsci.proteomics.maxquant.uniprot.TMT_batch_1_1_QE+.20161201.proteinGroups.txt",header=T,sep="\t",fill=T)
samples2<-colnames(tab)
tab2<-read.table("hipsci.proteomics.maxquant.uniprot.TMT_batch_1_2_fusion.20161201.proteinGroups.txt",header=T,sep="\t",fill=T)
samples3<-colnames(tab2)
lst<-rep("A",6)
tmt<-strsplit(tmt,"-")
for (i in 1:length(tmt)){ 
  lst[i]<-tmt[[i]][2]
}

names<-scan("Genotype_open_list.txt",what=character(),sep="\n")
all<-strsplit(names,"-")

sep<-rep("A",length(all))
for (i in 1:length(all)){
   sep[i]<-all[[i]][2]
}

sam<-rep("A",3)
n<-1
for (i in 1:length(lst)){
   line<-match(lst[i],sep)
   if (is.na(line) == F){
   sam[n]<-sep[line]
   n<-n+1
   }
}
#"eipl_2" "oilg_1" "suop_2" from labeled
smb://cmvm.datastore.ed.ac.uk/cmvm/eb/groups/michoel_grp/riho/Proteomics/HipSci_PXD003903_PXD005506/hipsci.proteomics.maxquant.uniprot.TMT_batch_1_3_fusion_DDA.20161201.proteinGroups.txt
tab<-read.table("hipsci.proteomics.maxquant.uniprot.TMT_batch_1_3_fusion_DDA.20161201.proteinGroups.txt",header=T,sep="\t",fill=T)
tbl<-tab[,c(7,68:73)]
rn<-as.character(tab[,1])
rownames(tbl) <-rn
tbl
#############################
gnames<-tab[,7]
x<-0
y<-0
for (i in 1:length(gnames)){
 if (gnames[i] == ""){
    y<-append(i,y,after=y[1])
   }
}

#include 0 genes
mtx<-matrix(rep(NA,(length(x)-1)*6),ncol=6)
gnames<-as.character(gnames)
rownames(mtx)<-gnames[x[2:8800]]
gn<-gnames[x[2:8800]]
tn<-colnames(tab)
colnames(mtx)<-tn[68:73]
cn<-tn[68:73]

for (i in 1:6){
  for (j in 1:length(mtx[,1])){ 
    mtx[j,i]<-tbl[x[j+1],i]
  }
}
#################################
write.table(tbl,"gnames_corrected_values_DDA.txt",quote=F)
tbl2<-read.table("gnames_corrected_values_DDA.txt",header=T,sep=" ",fill=T)

#remove 0 genes
dl<-0
for (i in 1:length(tbl2[,1])){
  if (sum(tbl2[i,3:8])== 0) {
     dl<-append(-i,dl,after=dl[1])
  }
}

tbl3<-tbl2[dl[2:192],]

write.table(tbl3,"gnames_corrected_values_no0_DDA.txt",quote=F,row.names=F)
write.table(tbl3[,c(1,3:8)],"ID_values_no0_DDA.txt",quote=F,row.names=F)
######################sort names#################
(head -n 1 gnames_corrected_values_no0.txt && tail -n +1 gnames_corrected_va
lues_no0.txt | sort -k1) >  gnames_sorted_values.txt
sed -n 6132,6132p gnames_sorted_values.txt
#################################################


