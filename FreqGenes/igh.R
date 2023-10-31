library(immunarch)
library(stringr)

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
Chain <- category[which(str_detect(category$V.name, "IGH")),] 
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

wt_ighv <- geneUsage(Chain$data, "hs.ighv",.quant = "count",.type = "segment",.norm = TRUE)
wt_ighv <- wt_ighv[apply(wt_ighv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_ighv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighv_wt.csv",row.names = FALSE ) 

wt_ighd <- geneUsage(Chain$data, "hs.ighd",.quant = "count",.type = "segment",.norm = TRUE)
wt_ighd <- wt_ighd[apply(wt_ighd, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_ighd,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighd_wt.csv",row.names = FALSE ) 

wt_ighj <- geneUsage(Chain$data, "hs.ighj",.quant = "count",.type = "segment",.norm = TRUE)
wt_ighj <- wt_ighj[apply(wt_ighj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_ighj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighj_wt.csv",row.names = FALSE ) 

####----mut_codel-------####
mut_codel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Codel"),]
category <-  tmps[which(tmps$CGGA_ID %in% mut_codel$CGGA_ID),]

Chain <- category[which(str_detect(category$V.name, "IGH")),] 
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

mutcodel_ighv <- geneUsage(Chain$data, "hs.ighv",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_ighv <- mutcodel_ighv[apply(mutcodel_ighv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_ighv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighv_mutcodel.csv",row.names = FALSE ) 

mutcodel_ighd <- geneUsage(Chain$data, "hs.ighd",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_ighd <- mutcodel_ighd[apply(mutcodel_ighd, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_ighd,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighd_mutcodel.csv",row.names = FALSE ) 

mutcodel_ighj <- geneUsage(Chain$data, "hs.ighj",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_ighj <- mutcodel_ighj[apply(mutcodel_ighj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_ighj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighj_mutcodel.csv",row.names = FALSE ) 

####----mut_noncodel-------####
mut_noncodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Noncodel"),]
category <-  tmps[which(tmps$CGGA_ID %in% mut_noncodel$CGGA_ID),] 

Chain <- category[which(str_detect(category$V.name, "IGH")),] 
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

mutnoncodel_ighv <- geneUsage(Chain$data, "hs.ighv",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_ighv <- mutnoncodel_ighv[apply(mutnoncodel_ighv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_ighv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighv_mutnoncodel.csv",row.names = FALSE ) 

mutnoncodel_ighd <- geneUsage(Chain$data, "hs.ighd",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_ighd <- mutnoncodel_ighd[apply(mutnoncodel_ighd, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_ighd,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighd_mutnoncodel.csv",row.names = FALSE ) 

mutnoncodel_ighj <- geneUsage(Chain$data, "hs.ighj",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_ighj <- mutnoncodel_ighj[apply(mutnoncodel_ighj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_ighj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/ighj_mutnoncodel.csv",row.names = FALSE ) 










