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
Chain <- category[which(str_detect(category$V.name, "IGK")),] 
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
public <- unique(public) 

public$len <- nchar(public$CDR3.aa)
public <- public[,c(1,4)]
public <- unique(public)
public$log2len <- log2(public$len)
write.csv(public,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igk_public.csv",row.names = FALSE )

private <- wt_Chain1[-which(wt_Chain1$newcolname %in% repeat_Chain1$newcolname),]
private <- private[,-1] 

private$len <- nchar(private$CDR3.aa)
private <- private[,c(1,4)]
private <- unique(private)
private$log2len <- log2(private$len)
write.csv(private,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igk_private.csv",row.names = FALSE )
