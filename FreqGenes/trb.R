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
tmps <- tmps [!grepl("\\?",tmps$CDR3.aa),] #369529，18

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

wt_trbv <- geneUsage(Chain$data, "hs.trbv",.quant = "count",.type = "segment",.norm = TRUE)
wt_trbv <- wt_trbv[apply(wt_trbv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_trbv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbv_wt.csv",row.names = FALSE ) 

wt_trbd <- geneUsage(Chain$data, "hs.trbd",.quant = "count",.type = "segment",.norm = TRUE)
wt_trbd <- wt_trbd[apply(wt_trbd, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_trbd,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbd_wt.csv",row.names = FALSE ) 

wt_trbj <- geneUsage(Chain$data, "hs.trbj",.quant = "count",.type = "segment",.norm = TRUE)
wt_trbj <- wt_trbj[apply(wt_trbj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_trbj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbj_wt.csv",row.names = FALSE ) 

####----mut_codel-------####
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

mutcodel_trbv <- geneUsage(Chain$data, "hs.trbv",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_trbv <- mutcodel_trbv[apply(mutcodel_trbv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_trbv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbv_mutcodel.csv",row.names = FALSE ) 

mutcodel_trbd <- geneUsage(Chain$data, "hs.trbd",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_trbd <- mutcodel_trbd[apply(mutcodel_trbd, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_trbd,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbd_mutcodel.csv",row.names = FALSE ) 

mutcodel_trbj <- geneUsage(Chain$data, "hs.trbj",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_trbj <- mutcodel_trbj[apply(mutcodel_trbj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_trbj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbj_mutcodel.csv",row.names = FALSE ) 

####----mut_noncodel-------####
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

mutnoncodel_trbv <- geneUsage(Chain$data, "hs.trbv",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_trbv <- mutnoncodel_trbv[apply(mutnoncodel_trbv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_trbv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbv_mutnoncodel.csv",row.names = FALSE ) 

mutnoncodel_trbd <- geneUsage(Chain$data, "hs.trbd",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_trbd <- mutnoncodel_trbd[apply(mutnoncodel_trbd, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_trbd,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbd_mutnoncodel.csv",row.names = FALSE ) 

mutnoncodel_trbj <- geneUsage(Chain$data, "hs.trbj",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_trbj <- mutnoncodel_trbj[apply(mutnoncodel_trbj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_trbj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/trbj_mutnoncodel.csv",row.names = FALSE ) 
