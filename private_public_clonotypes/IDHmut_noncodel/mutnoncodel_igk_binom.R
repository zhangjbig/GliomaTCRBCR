library(tidyr)
# mutnoncodel igk private
mutnoncodel_igk_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_igk_private.csv")
mutnoncodel_igk_private <- mutnoncodel_igk_private[,c(1,2)]
mutnoncodel_igk_private <- unique(mutnoncodel_igk_private) 

mutnoncodel_igk_private$first <- floor(mutnoncodel_igk_private$len/2)
mutnoncodel_igk_private$position1 <- substr(mutnoncodel_igk_private$CDR3.aa,mutnoncodel_igk_private$first,mutnoncodel_igk_private$first)
mutnoncodel_igk_private$position2 <- substr(mutnoncodel_igk_private$CDR3.aa,mutnoncodel_igk_private$first+1,mutnoncodel_igk_private$first+1)
mutnoncodel_igk_private$position3 <- substr(mutnoncodel_igk_private$CDR3.aa,mutnoncodel_igk_private$first+2,mutnoncodel_igk_private$first+2)

p1_private <- sum(mutnoncodel_igk_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igk_private$position1)
p2_private <- sum(mutnoncodel_igk_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igk_private$position2)
p3_private <- sum(mutnoncodel_igk_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igk_private$position3)

#mutnoncodel igk public
mutnoncodel_igk_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_igk_public.csv")
mutnoncodel_igk_public <- mutnoncodel_igk_public[,c(1,2)]
mutnoncodel_igk_public <- unique(mutnoncodel_igk_public) 

mutnoncodel_igk_public$first <- floor(mutnoncodel_igk_public$len/2)
mutnoncodel_igk_public$position1 <- substr(mutnoncodel_igk_public$CDR3.aa,mutnoncodel_igk_public$first,mutnoncodel_igk_public$first)
mutnoncodel_igk_public$position2 <- substr(mutnoncodel_igk_public$CDR3.aa,mutnoncodel_igk_public$first+1,mutnoncodel_igk_public$first+1)
mutnoncodel_igk_public$position3 <- substr(mutnoncodel_igk_public$CDR3.aa,mutnoncodel_igk_public$first+2,mutnoncodel_igk_public$first+2)

p1_public <- sum(mutnoncodel_igk_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igk_public$position1)
p2_public <- sum(mutnoncodel_igk_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igk_public$position2)
p3_public <- sum(mutnoncodel_igk_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_igk_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/proportion_hydro/mutnoncodel_igk_line_data.csv",row.names = FALSE ) 

#p1
binom.test(4925,19001,p=0.2401392) #以public为期望的概率 #p-value = 1.074e-09
#p2
binom.test(4273,19001,p=0.149652) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(10608,19001,p=0.5156613) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("1.074e-09","2.2e-16","2.2e-16"), "BH")
fdr
#1.074e-09 3.300e-16 3.300e-16