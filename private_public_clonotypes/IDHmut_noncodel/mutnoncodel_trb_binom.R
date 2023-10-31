library(tidyr)
# mutnoncodel trb private
mutnoncodel_trb_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_trb_private.csv")
mutnoncodel_trb_private <- mutnoncodel_trb_private[,c(1,2)]
mutnoncodel_trb_private <- unique(mutnoncodel_trb_private) 

mutnoncodel_trb_private$first <- floor(mutnoncodel_trb_private$len/2)
mutnoncodel_trb_private$position1 <- substr(mutnoncodel_trb_private$CDR3.aa,mutnoncodel_trb_private$first,mutnoncodel_trb_private$first)
mutnoncodel_trb_private$position2 <- substr(mutnoncodel_trb_private$CDR3.aa,mutnoncodel_trb_private$first+1,mutnoncodel_trb_private$first+1)
mutnoncodel_trb_private$position3 <- substr(mutnoncodel_trb_private$CDR3.aa,mutnoncodel_trb_private$first+2,mutnoncodel_trb_private$first+2)

p1_private <- sum(mutnoncodel_trb_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_trb_private$position1)
p2_private <- sum(mutnoncodel_trb_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_trb_private$position2)
p3_private <- sum(mutnoncodel_trb_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_trb_private$position3)

#mutnoncodel trb public
mutnoncodel_trb_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_trb_public.csv")
mutnoncodel_trb_public <- mutnoncodel_trb_public[,c(1,2)]
mutnoncodel_trb_public <- unique(mutnoncodel_trb_public) 

mutnoncodel_trb_public$first <- floor(mutnoncodel_trb_public$len/2)
mutnoncodel_trb_public$position1 <- substr(mutnoncodel_trb_public$CDR3.aa,mutnoncodel_trb_public$first,mutnoncodel_trb_public$first)
mutnoncodel_trb_public$position2 <- substr(mutnoncodel_trb_public$CDR3.aa,mutnoncodel_trb_public$first+1,mutnoncodel_trb_public$first+1)
mutnoncodel_trb_public$position3 <- substr(mutnoncodel_trb_public$CDR3.aa,mutnoncodel_trb_public$first+2,mutnoncodel_trb_public$first+2)

p1_public <- sum(mutnoncodel_trb_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_trb_public$position1)
p2_public <- sum(mutnoncodel_trb_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_trb_public$position2)
p3_public <- sum(mutnoncodel_trb_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_trb_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/proportion_hydro/mutnoncodel_trb_line_data.csv",row.names = FALSE ) 

#p1
binom.test(1688,3006,p=0.375) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(1813,3006,p=0.625) #以public为期望的概率 #p-value = 0.01359
#p3
binom.test(1477,3006,p=0.5) #以public为期望的概率 #p-value = 0.3523

fdr <- p.adjust(c("2.2e-16","0.01359","0.3523"), "BH")
fdr
#6.6000e-16 2.0385e-02 3.5230e-01