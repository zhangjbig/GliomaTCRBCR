library(immunarch)
library(stringr)


load("/Users/wanglu/Documents/2021_TC/tb.RData")

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
Chain <- category[which(str_detect(category$V.name, "IGL")),] 
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

wt_iglv <- geneUsage(Chain$data, "hs.iglv",.quant = "count",.type = "segment",.norm = TRUE)
wt_iglv <- wt_iglv[apply(wt_iglv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_iglv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/iglv_wt.csv",row.names = FALSE ) 

wt_iglj <- geneUsage(Chain$data, "hs.iglj",.quant = "count",.type = "segment",.norm = TRUE)
wt_iglj <- wt_iglj[apply(wt_iglj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(wt_iglj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/iglj_wt.csv",row.names = FALSE ) 

####----mut_codel-------####
mut_codel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Codel"),]
category <-  tmps[which(tmps$CGGA_ID %in% mut_codel$CGGA_ID),]

Chain <- category[which(str_detect(category$V.name, "IGL")),] 
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

mutcodel_iglv <- geneUsage(Chain$data, "hs.iglv",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_iglv <- mutcodel_iglv[apply(mutcodel_iglv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_iglv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/iglv_mutcodel.csv",row.names = FALSE ) 

mutcodel_iglj <- geneUsage(Chain$data, "hs.iglj",.quant = "count",.type = "segment",.norm = TRUE)
mutcodel_iglj <- mutcodel_iglj[apply(mutcodel_iglj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutcodel_iglj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/iglj_mutcodel.csv",row.names = FALSE ) 

####----mut_noncodel-------####
mut_noncodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Noncodel"),]
category <-  tmps[which(tmps$CGGA_ID %in% mut_noncodel$CGGA_ID),] 

Chain <- category[which(str_detect(category$V.name, "IGL")),]  
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

mutnoncodel_iglv <- geneUsage(Chain$data, "hs.iglv",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_iglv <- mutnoncodel_iglv[apply(mutnoncodel_iglv, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_iglv,"/Users/wanglu/Documents/2021_TC/final_version2/gene/iglv_mutnoncodel.csv",row.names = FALSE ) 

mutnoncodel_iglj <- geneUsage(Chain$data, "hs.iglj",.quant = "count",.type = "segment",.norm = TRUE)
mutnoncodel_iglj <- mutnoncodel_iglj[apply(mutnoncodel_iglj, 1, function(row) sum(!is.na(row)) > 3),] #如果非NA的数目大于3个，就留下来，这样才有统计学意义
write.csv(mutnoncodel_iglj,"/Users/wanglu/Documents/2021_TC/final_version2/gene/iglj_mutnoncodel.csv",row.names = FALSE ) 
