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
#####---------mut codel----------#####
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

mutcodel <- data.frame(len = 0,prop = 0,CGGA_ID = 0)
mutcodel <- mutcodel[-1,]
for (i in 1:length(names(Chain$data))){
  myvec <- Chain$data[[i]]
  myvec$len <- nchar(myvec$CDR3.aa)
  myvec$prop <- myvec$Clones/sum(myvec$Clones)
  myvec <- myvec[,c("len","prop")]
  myvec <- aggregate(x =myvec$prop, by=list(len = myvec$len),FUN=sum)
  colnames(myvec)[2] <- "prop"
  myvec$CGGA_ID <- names(Chain$data)[i]
  mutcodel <-rbind(mutcodel,myvec)
}
mutcodel <- mutcodel[order(mutcodel$len),]
write.csv(mutcodel,"/Users/wanglu/Documents/2021_TC/final_version2/length/tra_mutcodel_len.csv",row.names = FALSE ) 
