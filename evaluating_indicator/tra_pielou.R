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
tmps <- tmps [!grepl("\\?",tmps$CDR3.aa),] #369529，18

#-----------------------------#
#----------  wt ----------#
wt <- c_info[which(c_info$Sample =="IDH-WT"),]

category <-  tmps[which(tmps$CGGA_ID %in% wt$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRA")),] 
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

wt_shdi <- unlist(lapply(1:length(names(Chain$data)), FUN = function(x){diversity(Chain$data[[x]]$Clones,"shannon",base = 2)}))
N <- unlist(lapply(1:length(names(Chain$data)), FUN = function(x){specnumber(Chain$data[[x]]$Clones)}))
wt_pielou <- wt_shdi/log(N,2)
wt_pielou <- as.data.frame(wt_pielou)
wt_pielou$CGGA_ID <- names(Chain$data)
wt_pielou$log2 <- log2(wt_pielou$wt_pielou)
write.csv(wt_pielou,"/Users/wanglu/Documents/2021_TC/final_version2/evaluating_indicator/tra_wt_pielou.csv",row.names = FALSE ) 

#----------  mut codel ----------#
mut_codel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Codel"),]

category <-  tmps[which(tmps$CGGA_ID %in% mut_codel$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRA")),] 
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

mutcodel_shdi <- unlist(lapply(1:length(names(Chain$data)), FUN = function(x){diversity(Chain$data[[x]]$Clones,"shannon",base = 2)}))
N <- unlist(lapply(1:length(names(Chain$data)), FUN = function(x){specnumber(Chain$data[[x]]$Clones)}))
mutcodel_pielou <- mutcodel_shdi/log(N,2)
mutcodel_pielou <- as.data.frame(mutcodel_pielou)
mutcodel_pielou$CGGA_ID <- names(Chain$data)
mutcodel_pielou$log2 <- log2(mutcodel_pielou$mutcodel_pielou)
write.csv(mutcodel_pielou,"/Users/wanglu/Documents/2021_TC/final_version2/evaluating_indicator/tra_mutcodel_pielou.csv",row.names = FALSE ) 

#----------  mut noncodel----------#
mut_noncodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Noncodel"),]

category <-  tmps[which(tmps$CGGA_ID %in% mut_noncodel$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRA")),] 
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

mutnoncodel_shdi <- unlist(lapply(1:length(names(Chain$data)), FUN = function(x){diversity(Chain$data[[x]]$Clones,"shannon",base = 2)}))
N <- unlist(lapply(1:length(names(Chain$data)), FUN = function(x){specnumber(Chain$data[[x]]$Clones)}))
mutnoncodel_pielou <- mutnoncodel_shdi/log(N,2)
mutnoncodel_pielou <- as.data.frame(mutnoncodel_pielou)
mutnoncodel_pielou$CGGA_ID <- names(Chain$data)
mutnoncodel_pielou$log2 <- log2(mutnoncodel_pielou$mutnoncodel_pielou)
write.csv(mutnoncodel_pielou,"/Users/wanglu/Documents/2021_TC/final_version2/evaluating_indicator/tra_mutnoncodel_pielou.csv",row.names = FALSE ) 


