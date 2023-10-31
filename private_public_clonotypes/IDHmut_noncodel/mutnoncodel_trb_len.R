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
#----------  mutnoncodel ----------#
mutnoncodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Noncodel"),]

category <-  tmps[which(tmps$CGGA_ID %in% mutnoncodel$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRB")),]
Chain <- Chain[which(Chain$D.name!="*"),]
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","D.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            D.name = Chain$D.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[7] <- "Clones"


mutnoncodel_Chain <- Chain[,c(3,4,5,6)]
#以下是所有重复的行，只计算一次
repeat_Chain <- unique(mutnoncodel_Chain[duplicated(mutnoncodel_Chain),])

repeat_line <- data.frame(CDR3.aa = 0,V.name = 0,D.name = 0,J.name = 0)
repeat_line <- repeat_line[-1,]
row = data.frame(row = 0)
row = row[-1,]
for (i in c(1:nrow(repeat_Chain))){
  for (j in c(1:nrow(mutnoncodel_Chain))){
    if (sum(repeat_Chain[i,] == mutnoncodel_Chain[j,]) == ncol(mutnoncodel_Chain)){
      repeat_line <- rbind(repeat_line,mutnoncodel_Chain[j,])
      row = rbind(row,j)
    }
  }
}
row <- as.data.frame(row)
public <- Chain[row$V1,]
private <- Chain[-row$V1,]

public <- public[,c(3,4,5,6)]
public$len <- nchar(public$CDR3.aa)
public <- public[,c(1,5)]
public <- unique(public)
public$log2len <- log2(public$len)
write.csv(public,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_trb_public.csv",row.names = FALSE )

private <- private[,c(3,4,5,6)]
private$len <- nchar(private$CDR3.aa)
private <- private[,c(1,5)]
private <- unique(private)
private$log2len <- log2(private$len)
write.csv(private,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_trb_private.csv",row.names = FALSE )
