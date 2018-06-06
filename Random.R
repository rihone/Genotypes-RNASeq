geno<-read.table('RNASeq-inds_open_res.txt', sep='\t', stringsAsFactors=F,header=T,row.names=1)
shID<-sample(1:length(geno[1,]),replace=F)
geno2<-geno[,shID]
write.table(geno2,"RNASeq-inds_open_res_ram.txt",sep="\t",quote=F)
