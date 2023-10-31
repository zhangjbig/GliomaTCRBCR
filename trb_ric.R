library(stringr)
library(vegan)

CGGA <- read.table("/Users/wanglu/Documents/2021_TC/CGGA_total_TRUST.cdr3.tsv", stringsAsFactors = F)
dim(CGGA)

tmp <- CGGA
names(tmp) <- c("CGGA_ID", "Clones", "Proportion", "CDR3.nt", "CDR3.aa", "V.name", "D.name", "J.name", "C.name")

tmp$V.end <- NA
tmp$D.start <- NA
tmp$D.end <- NA
tmp$J.start <- NA
tmp$VJ.ins <- NA
tmp$VD.ins <- NA
tmp$DJ.ins <- NA
tmp$Sequence <- NA

clinical1 <- read.csv("/Users/wanglu/Documents/2021_TC/clinical_info/IDH-MUT_X1p19q-Codel.csv")
clinical2 <- read.csv("/Users/wanglu/Documents/2021_TC/clinical_info/IDH-MUT_X1p19q-Noncodel.csv")
clinical3 <- read.csv("/Users/wanglu/Documents/2021_TC/clinical_info/IDH-WT.csv")

C_clinical <- rbind(clinical1, clinical2, clinical3)

intersect_data <- tmp[which(tmp$CGGA_ID %in% C_clinical$CGGA_ID),] 

Sample<-c()
for(ID in intersect_data$CGGA_ID) {
  print(ID)
  col_sample = C_clinical[which(C_clinical$CGGA_ID == ID),c("Sample")]
  print(col_sample)
  Sample<-c(Sample,col_sample)
}

#save.image("/Users/wanglu/Documents/2021_TC/tb.RData") 
load("/Users/wanglu/Documents/2021_TC/tb.RData")

#Clinical information update
c_info <- read.csv("/Users/wanglu/Documents/2021_TC/XZY_database/CGGA_clinical_TCRBCR.csv") 

tmps <- cbind(Sample,intersect_data[,1:ncol(tmp)])  

tmps <- tmps[which(tmps$CDR3.aa != "out_of_frame"),]
tmps <- tmps[which(tmps$CDR3.aa != "partial"),]
tmps <- tmps[which(tmps$CDR3.aa != "?"),]
tmps <- tmps [!grepl("_",tmps$CDR3.aa),] 
tmps <- tmps [!grepl("\\?",tmps$CDR3.aa),] 

#-----------------------------#
#----------  wt ----------#
wt <- c_info[which(c_info$Sample =="IDH-WT"),]

category <-  tmps[which(tmps$CGGA_ID %in% wt$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRB")),] 
Chain <- Chain[which(Chain$D.name!="*"),] 
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","D.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            D.name = Chain$D.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[7] <- "Clones"

Chain <- Chain[,-1]
Chain <- split(Chain[,-1], Chain$CGGA_ID ) 
Chain <- list(Chain)
names(Chain) <- c("data")

Chain$data[[1]]$Clones <- round(Chain$data[[1]]$Clones)
est <- estimateR(Chain$data[[1]]$Clones)
est <-as.data.frame(est)
colnames(est) <- names(Chain$data)[1]
est <- t(est)
for (i in 2:length(names(Chain$data))){
  Chain$data[[i]]$Clones <- round(Chain$data[[i]]$Clones)
  est1 <- estimateR(Chain$data[[i]]$Clones)
  est1 <-as.data.frame(est1)
  colnames(est1) <- names(Chain$data)[i]
  est1 <- t(est1)
  est <- rbind(est,est1)
}
wt_ric <- est
wt_ric <- as.data.frame(wt_ric)
wt_ric$log2obs <- log2(wt_ric$S.obs)
wt_ric$log2chao1 <- log2(wt_ric$S.chao1)
wt_ric$log2ace <- log2(wt_ric$S.ACE)
write.csv(wt_ric,"/Users/wanglu/Documents/2021_TC/final_version2/evaluating_indicator/trb_wt_ric.csv",row.names = TRUE ) 


#----------  mut codel ----------#
mut_codel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Codel"),]

category <-  tmps[which(tmps$CGGA_ID %in% mut_codel$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRB")),]  
Chain <- Chain[which(Chain$D.name!="*"),] 
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","D.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            D.name = Chain$D.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[7] <- "Clones"

Chain <- Chain[,-1]
Chain <- split(Chain[,-1], Chain$CGGA_ID ) 
Chain <- list(Chain)
names(Chain) <- c("data")

Chain$data[[1]]$Clones <- round(Chain$data[[1]]$Clones)
est <- estimateR(Chain$data[[1]]$Clones)
est <-as.data.frame(est)
colnames(est) <- names(Chain$data)[1]
est <- t(est)
for (i in 2:length(names(Chain$data))){
  Chain$data[[i]]$Clones <- round(Chain$data[[i]]$Clones)
  est1 <- estimateR(Chain$data[[i]]$Clones)
  est1 <-as.data.frame(est1)
  colnames(est1) <- names(Chain$data)[i]
  est1 <- t(est1)
  est <- rbind(est,est1)
}
mutcodel_ric <- est
mutcodel_ric <- as.data.frame(mutcodel_ric)
mutcodel_ric$log2obs <- log2(mutcodel_ric$S.obs)
mutcodel_ric$log2chao1 <- log2(mutcodel_ric$S.chao1)
mutcodel_ric$log2ace <- log2(mutcodel_ric$S.ACE)
write.csv(mutcodel_ric,"/Users/wanglu/Documents/2021_TC/final_version2/evaluating_indicator/trb_mutcodel_ric.csv",row.names = TRUE ) 


#----------  mut noncodel----------#
mut_noncodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Noncodel"),]

category <-  tmps[which(tmps$CGGA_ID %in% mut_noncodel$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRB")),]  
Chain <- Chain[which(Chain$D.name!="*"),] 
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","D.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            D.name = Chain$D.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[7] <- "Clones"

Chain <- Chain[,-1]
Chain <- split(Chain[,-1], Chain$CGGA_ID ) 
Chain <- list(Chain)
names(Chain) <- c("data")

Chain$data[[1]]$Clones <- round(Chain$data[[1]]$Clones)
est <- estimateR(Chain$data[[1]]$Clones)
est <-as.data.frame(est)
colnames(est) <- names(Chain$data)[1]
est <- t(est)
for (i in 2:length(names(Chain$data))){
  Chain$data[[i]]$Clones <- round(Chain$data[[i]]$Clones)
  est1 <- estimateR(Chain$data[[i]]$Clones)
  est1 <-as.data.frame(est1)
  colnames(est1) <- names(Chain$data)[i]
  est1 <- t(est1)
  est <- rbind(est,est1)
}
mutnoncodel_ric <- est
mutnoncodel_ric <- as.data.frame(mutnoncodel_ric)
mutnoncodel_ric$log2obs <- log2(mutnoncodel_ric$S.obs)
mutnoncodel_ric$log2chao1 <- log2(mutnoncodel_ric$S.chao1)
mutnoncodel_ric$log2ace <- log2(mutnoncodel_ric$S.ACE)
write.csv(mutnoncodel_ric,"/Users/wanglu/Documents/2021_TC/final_version2/evaluating_indicator/trb_mutnoncodel_ric.csv",row.names = TRUE ) 
