library(tidyr)
# mutnoncodel igl private
mutnoncodel_igl_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_igl_private.csv")
mutnoncodel_igl_private <- mutnoncodel_igl_private[,c(1,2)]
mutnoncodel_igl_private <- unique(mutnoncodel_igl_private) 

mutnoncodel_igl_private$first <- floor(mutnoncodel_igl_private$len/2)
mutnoncodel_igl_private$position1 <- substr(mutnoncodel_igl_private$CDR3.aa,mutnoncodel_igl_private$first,mutnoncodel_igl_private$first)
mutnoncodel_igl_private$position2 <- substr(mutnoncodel_igl_private$CDR3.aa,mutnoncodel_igl_private$first+1,mutnoncodel_igl_private$first+1)
mutnoncodel_igl_private$position3 <- substr(mutnoncodel_igl_private$CDR3.aa,mutnoncodel_igl_private$first+2,mutnoncodel_igl_private$first+2)

p1_private <- sum(mutnoncodel_igl_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igl_private$position1)
p2_private <- sum(mutnoncodel_igl_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igl_private$position2)
p3_private <- sum(mutnoncodel_igl_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igl_private$position3)

#mutnoncodel igl public
mutnoncodel_igl_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_igl_public.csv")
mutnoncodel_igl_public <- mutnoncodel_igl_public[,c(1,2)]
mutnoncodel_igl_public <- unique(mutnoncodel_igl_public) 

mutnoncodel_igl_public$first <- floor(mutnoncodel_igl_public$len/2)
mutnoncodel_igl_public$position1 <- substr(mutnoncodel_igl_public$CDR3.aa,mutnoncodel_igl_public$first,mutnoncodel_igl_public$first)
mutnoncodel_igl_public$position2 <- substr(mutnoncodel_igl_public$CDR3.aa,mutnoncodel_igl_public$first+1,mutnoncodel_igl_public$first+1)
mutnoncodel_igl_public$position3 <- substr(mutnoncodel_igl_public$CDR3.aa,mutnoncodel_igl_public$first+2,mutnoncodel_igl_public$first+2)

p1_public <- sum(mutnoncodel_igl_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igl_public$position1)
p2_public <- sum(mutnoncodel_igl_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igl_public$position2)
p3_public <- sum(mutnoncodel_igl_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igl_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/proportion_hydro/mutnoncodel_igl_line_data.csv",row.names = FALSE ) 

#p1
binom.test(5230,17740,p=0.2522936) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(4364,17740,p=0.1399083) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(7579,17740,p=0.565367) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("2.2e-16","2.2e-16","2.2e-16"), "BH")
fdr
