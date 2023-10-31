library(stringr)
library(tidyr)

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
Chain <- category[which(str_detect(category$V.name, "IGL")),]  
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[6] <- "Clones"


wt_Chain <- Chain[,c(3,4,5)]
#以下是所有重复的行，只计算一次
repeat_Chain <- unique(wt_Chain[duplicated(wt_Chain),])
repeat_Chain1 <- unite(repeat_Chain,"newcolname","CDR3.aa","V.name","J.name",sep = ",",remove = FALSE)

wt_Chain1 <- unite(wt_Chain,"newcolname","CDR3.aa","V.name","J.name",sep = ",",remove = FALSE)

public <- wt_Chain1[which(wt_Chain1$newcolname %in% repeat_Chain1$newcolname),]
public <- public[,-1]
public <- unique(public) #2940种克隆类型

private <- wt_Chain1[-which(wt_Chain1$newcolname %in% repeat_Chain1$newcolname),]
private <- private[,-1] 
private <- unique(private) #41030种克隆类型

#计算样本数目
Chain1 <- Chain
Chain1 <- unite(Chain1,"newcolname","CDR3.aa","V.name","J.name",sep = ",",remove = FALSE)

public <- wt_Chain1[which(wt_Chain1$newcolname %in% repeat_Chain1$newcolname),]
public_sample <- Chain1[which(Chain1$newcolname %in% public$newcolname),]

private <- wt_Chain1[-which(wt_Chain1$newcolname %in% repeat_Chain1$newcolname),]
private_sample <- Chain1[which(Chain1$newcolname %in% private$newcolname),]

#计算一共有多少个样本
unique_sample <- c(public_sample$CGGA_ID,private_sample$CGGA_ID)
unique_sample <- unique(unique_sample) #399个样本
