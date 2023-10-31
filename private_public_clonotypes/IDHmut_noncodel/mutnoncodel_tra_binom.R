library(tidyr)
# mutnoncodel tra private
mutnoncodel_tra_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_tra_private.csv")
mutnoncodel_tra_private <- mutnoncodel_tra_private[,c(1,2)]
mutnoncodel_tra_private <- unique(mutnoncodel_tra_private) 

mutnoncodel_tra_private$first <- floor(mutnoncodel_tra_private$len/2)
mutnoncodel_tra_private$position1 <- substr(mutnoncodel_tra_private$CDR3.aa,mutnoncodel_tra_private$first,mutnoncodel_tra_private$first)
mutnoncodel_tra_private$position2 <- substr(mutnoncodel_tra_private$CDR3.aa,mutnoncodel_tra_private$first+1,mutnoncodel_tra_private$first+1)
mutnoncodel_tra_private$position3 <- substr(mutnoncodel_tra_private$CDR3.aa,mutnoncodel_tra_private$first+2,mutnoncodel_tra_private$first+2)

p1_private <- sum(mutnoncodel_tra_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_tra_private$position1)
p2_private <- sum(mutnoncodel_tra_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_tra_private$position2)
p3_private <- sum(mutnoncodel_tra_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_tra_private$position3)

#mutnoncodel tra public
mutnoncodel_tra_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/mutnoncodel_tra_public.csv")
mutnoncodel_tra_public <- mutnoncodel_tra_public[,c(1,2)]
mutnoncodel_tra_public <- unique(mutnoncodel_tra_public) 

mutnoncodel_tra_public$first <- floor(mutnoncodel_tra_public$len/2)
mutnoncodel_tra_public$position1 <- substr(mutnoncodel_tra_public$CDR3.aa,mutnoncodel_tra_public$first,mutnoncodel_tra_public$first)
mutnoncodel_tra_public$position2 <- substr(mutnoncodel_tra_public$CDR3.aa,mutnoncodel_tra_public$first+1,mutnoncodel_tra_public$first+1)
mutnoncodel_tra_public$position3 <- substr(mutnoncodel_tra_public$CDR3.aa,mutnoncodel_tra_public$first+2,mutnoncodel_tra_public$first+2)

p1_public <- sum(mutnoncodel_tra_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_tra_public$position1)
p2_public <- sum(mutnoncodel_tra_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_tra_public$position2)
p3_public <- sum(mutnoncodel_tra_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutnoncodel_tra_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutnoncodel/proportion_hydro/mutnoncodel_tra_line_data.csv",row.names = FALSE ) 

#p1
binom.test(640,1422,p=0.3571429) #以public为期望的概率 #p-value = 6.498e-13
#p2
binom.test(711,1422,p=0.3571429) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(635,1422,p=0.4285714) #以public为期望的概率 #p-value = 0.1719

fdr <- p.adjust(c("6.498e-13","2.2e-16","0.1719"), "BH")
fdr
#9.747e-13 6.600e-16 1.719e-01