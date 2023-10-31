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
Chain <- category[which(str_detect(category$V.name, "IGK")),] 
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[6] <- "Clones"

Chain <- Chain[,-1]
Chain <- split(Chain[,-1], Chain$CGGA_ID ) 
Chain <- list(Chain)
names(Chain) <- c("data")

wt_igkv <- geneUsage(Chain$data, "hs.igkv",.quant = "count",.type = "segment",.norm = TRUE)
wt_igkv <- wt_igkv[apply(wt_igkv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_igkv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/igkv_wt.csv",row.names = FALSE ) 

wt_igkj <- geneUsage(Chain$data, "hs.igkj",.quant = "count",.type = "segment",.norm = TRUE)
wt_igkj <- wt_igkj[apply(wt_igkj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_igkj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/igkj_wt.csv",row.names = FALSE ) 

####----mut_codel-------####
mut_codel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Codel"),]
category <-  tmps[which(tmps$CGGA_ID %in% mut_codel$CGGA_ID),]

Chain <- category[which(str_detect(category$V.name, "IGK")),] 
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[6] <- "Clones"

Chain <- Chain[,-1]
Chain <- split(Chain[,-1], Chain$CGGA_ID )
Chain <- list(Chain)
names(Chain) <- c("data")

mutcodel_igkv <- geneUsage(Chain$data, "hs.igkv",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_igkv <- mutcodel_igkv[apply(mutcodel_igkv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_igkv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/igkv_mutcodel.csv",row.names = FALSE ) 

mutcodel_igkj <- geneUsage(Chain$data, "hs.igkj",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_igkj <- mutcodel_igkj[apply(mutcodel_igkj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_igkj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/igkj_mutcodel.csv",row.names = FALSE ) 

####----mut_noncodel-------####
mut_noncodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Noncodel"),]
category <-  tmps[which(tmps$CGGA_ID %in% mut_noncodel$CGGA_ID),] 

Chain <- category[which(str_detect(category$V.name, "IGK")),] 
Chain <- Chain[which(Chain$J.name!="*"),]  

Chain <- Chain[,c("Sample","CGGA_ID","Clones","CDR3.aa","V.name","J.name")]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            CDR3.aa = Chain$CDR3.aa, V.name = Chain$V.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[6] <- "Clones"

Chain <- Chain[,-1]
Chain <- split(Chain[,-1], Chain$CGGA_ID )
Chain <- list(Chain)
names(Chain) <- c("data")

mutnoncodel_igkv <- geneUsage(Chain$data, "hs.igkv",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_igkv <- mutnoncodel_igkv[apply(mutnoncodel_igkv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_igkv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/igkv_mutnoncodel.csv",row.names = FALSE ) 

mutnoncodel_igkj <- geneUsage(Chain$data, "hs.igkj",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_igkj <- mutnoncodel_igkj[apply(mutnoncodel_igkj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_igkj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/igkj_mutnoncodel.csv",row.names = FALSE ) 
