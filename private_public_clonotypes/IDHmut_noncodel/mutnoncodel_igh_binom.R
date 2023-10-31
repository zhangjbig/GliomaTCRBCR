library(tidyr)
# mutnoncodel igh private
mutnoncodel_igh_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_igh_private.csv")
mutnoncodel_igh_private <- mutnoncodel_igh_private[,c(1,2)]
mutnoncodel_igh_private <- unique(mutnoncodel_igh_private) 

mutnoncodel_igh_private$first <- floor(mutnoncodel_igh_private$len/2)
mutnoncodel_igh_private$position1 <- substr(mutnoncodel_igh_private$CDR3.aa,mutnoncodel_igh_private$first,mutnoncodel_igh_private$first)
mutnoncodel_igh_private$position2 <- substr(mutnoncodel_igh_private$CDR3.aa,mutnoncodel_igh_private$first+1,mutnoncodel_igh_private$first+1)
mutnoncodel_igh_private$position3 <- substr(mutnoncodel_igh_private$CDR3.aa,mutnoncodel_igh_private$first+2,mutnoncodel_igh_private$first+2)

p1_private <- sum(mutnoncodel_igh_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igh_private$position1)
p2_private <- sum(mutnoncodel_igh_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igh_private$position2)
p3_private <- sum(mutnoncodel_igh_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igh_private$position3)

#mutnoncodel igh public
mutnoncodel_igh_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_igh_public.csv")
mutnoncodel_igh_public <- mutnoncodel_igh_public[,c(1,2)]
mutnoncodel_igh_public <- unique(mutnoncodel_igh_public) 

mutnoncodel_igh_public$first <- floor(mutnoncodel_igh_public$len/2)
mutnoncodel_igh_public$position1 <- substr(mutnoncodel_igh_public$CDR3.aa,mutnoncodel_igh_public$first,mutnoncodel_igh_public$first)
mutnoncodel_igh_public$position2 <- substr(mutnoncodel_igh_public$CDR3.aa,mutnoncodel_igh_public$first+1,mutnoncodel_igh_public$first+1)
mutnoncodel_igh_public$position3 <- substr(mutnoncodel_igh_public$CDR3.aa,mutnoncodel_igh_public$first+2,mutnoncodel_igh_public$first+2)

p1_public <- sum(mutnoncodel_igh_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igh_public$position1)
p2_public <- sum(mutnoncodel_igh_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igh_public$position2)
p3_public <- sum(mutnoncodel_igh_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igh_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/proportion_hydro/mutnoncodel_igh_line_data.csv",row.names = FALSE ) 

#p1
binom.test(16701,34299,p=0.5694444) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(16889,34299,p=0.5138889) #以public为期望的概率 #p-value = 1.752e-15
#p3
binom.test(15924,34299,p=0.4027778) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("2.2e-16","1.752e-15","2.2e-16"), "BH")
fdr
#3.300e-16 1.752e-15 3.300e-16