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
#----------  mutcodel ----------#
mutcodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Codel"),]

category <-  tmps[which(tmps$CGGA_ID %in% mutcodel$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRA")),] 
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[6] <- "Clones"


mutcodel_Chain <- Chain[,c(3,4,5)]
#以下是所有重复的行，只计算一次
repeat_Chain <- unique(mutcodel_Chain[duplicated(mutcodel_Chain),])

repeat_line <- data.frame(CDR3.aa = 0,V.name = 0,D.name = 0)
repeat_line <- repeat_line[-1,]
row = data.frame(row = 0)
row = row[-1,]
for (i in c(1:nrow(repeat_Chain))){
  for (j in c(1:nrow(mutcodel_Chain))){
    if (sum(repeat_Chain[i,] == mutcodel_Chain[j,]) == ncol(mutcodel_Chain)){
      repeat_line <- rbind(repeat_line,mutcodel_Chain[j,])
      row = rbind(row,j)
    }
  }
}
row <- as.data.frame(row)
public <- Chain[row$V1,]
private <- Chain[-row$V1,]

public <- public[,c(3,4,5)]
public$len <- nchar(public$CDR3.aa)
public <- public[,c(1,4)]
public <- unique(public)
public$log2len <- log2(public$len)
write.csv(public,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_tra_public.csv",row.names = FALSE )

private <- private[,c(3,4,5)]
private$len <- nchar(private$CDR3.aa)
private <- private[,c(1,4)]
private <- unique(private)
private$log2len <- log2(private$len)
write.csv(private,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_tra_private.csv",row.names = FALSE )
