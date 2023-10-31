library(stringr)
library(tidyr)
library(voronoiTreemap)
library(randomcoloR)

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
#----------mutcodel----------#
mutcodel <- c_info[which(c_info$Sample =="IDH-MUT_X1p19q-Codel"),]

category <-  tmps[which(tmps$CGGA_ID %in% mutcodel$CGGA_ID),] 
Chain <- category[which(str_detect(category$V.name, "TRA")),] 
Chain <- Chain[which(Chain$J.name!="*"),]  


Chain <- Chain[,c("Sample","CGGA_ID","Clones","V.name","J.name")]
Chain <- separate(Chain,V.name, into= c("V.name","delete1"),sep= "\\*")
Chain <- separate(Chain,J.name, into= c("J.name","delete2"),sep= "\\*")
Chain <- Chain[,-c(5,7)]
Chain <- aggregate(x =Chain$Clones, by=list(Sample = Chain$Sample, CGGA_ID = Chain$CGGA_ID,
                                            V.name = Chain$V.name,
                                            J.name = Chain$J.name),FUN=sum)
colnames(Chain)[5] <- "Clones"

treedf <- Chain

treedf <- unite(treedf,"h2",c("V.name","J.name"), sep=",", remove = F)
treedf <- treedf[,c(2,3,6)]

#克隆类型的种类
clones <- as.data.frame(unique(treedf$h2))
colnames(clones) <- "clones"

DF <- data.frame(clones = 0,weight = 0)
DF <- DF[-1,]
for (i in c(1:nrow(clones))){
  df1 <- treedf[which(treedf$h2 == clones[i,1]),]
  weight <- sum(df1$Clones)/nrow(df1)
  df <- cbind(clones[i,1],weight)
  df <- as.data.frame(df)
  colnames(df)[1] <- "clones"
  DF <- rbind(DF,df)
}

treemap_data <- DF
treemap_data$weight <- as.numeric(treemap_data$weight)
treemap_data$weight <- treemap_data$weight/sum(treemap_data$weight)*100

palette <- distinctColorPalette(337) 
palette <- as.data.frame(palette)
palette1 <- palette$palette[duplicated(palette$palette)] 
treemap_data$color <- palette$palette

treemap_data$h1 <- "mutcodel"
treemap_data$h2 <- treemap_data$clones
treemap_data$h3 <- treemap_data$clones
treemap_data$codes <- treemap_data$clones

treemap_data <- treemap_data[,c("h1","h2","h3","color","weight","codes")]
treemap <- vt_export_json(vt_input_from_df(treemap_data,scaleToPerc = FALSE))
vt_d3(treemap, label=FALSE,width = 400,height = 400)
