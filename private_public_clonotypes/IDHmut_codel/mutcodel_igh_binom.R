library(tidyr)
# mutcodel igh private
mutcodel_igh_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_igh_private.csv")
mutcodel_igh_private <- mutcodel_igh_private[,c(1,2)]
mutcodel_igh_private <- unique(mutcodel_igh_private) 

mutcodel_igh_private$first <- floor(mutcodel_igh_private$len/2)
mutcodel_igh_private$position1 <- substr(mutcodel_igh_private$CDR3.aa,mutcodel_igh_private$first,mutcodel_igh_private$first)
mutcodel_igh_private$position2 <- substr(mutcodel_igh_private$CDR3.aa,mutcodel_igh_private$first+1,mutcodel_igh_private$first+1)
mutcodel_igh_private$position3 <- substr(mutcodel_igh_private$CDR3.aa,mutcodel_igh_private$first+2,mutcodel_igh_private$first+2)

p1_private <- sum(mutcodel_igh_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igh_private$position1)
p2_private <- sum(mutcodel_igh_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igh_private$position2)
p3_private <- sum(mutcodel_igh_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igh_private$position3)

#mutcodel igh public
mutcodel_igh_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_igh_public.csv")
mutcodel_igh_public <- mutcodel_igh_public[,c(1,2)]
mutcodel_igh_public <- unique(mutcodel_igh_public) 

mutcodel_igh_public$first <- floor(mutcodel_igh_public$len/2)
mutcodel_igh_public$position1 <- substr(mutcodel_igh_public$CDR3.aa,mutcodel_igh_public$first,mutcodel_igh_public$first)
mutcodel_igh_public$position2 <- substr(mutcodel_igh_public$CDR3.aa,mutcodel_igh_public$first+1,mutcodel_igh_public$first+1)
mutcodel_igh_public$position3 <- substr(mutcodel_igh_public$CDR3.aa,mutcodel_igh_public$first+2,mutcodel_igh_public$first+2)

p1_public <- sum(mutcodel_igh_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igh_public$position1)
p2_public <- sum(mutcodel_igh_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igh_public$position2)
p3_public <- sum(mutcodel_igh_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igh_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/proportion_hydro/mutcodel_igh_line_data.csv",row.names = FALSE ) 

#p1
binom.test(2922,5870,p=0.5882353) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(2790,5870,p=0.4117647) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(2821,5870,p=0.5294118) #以public为期望的概率 #p-value = 7.195e-14

fdr <- p.adjust(c("2.2e-16","2.2e-16","7.195e-14"), "BH")
fdr

#3.300e-16 3.300e-16 7.195e-14
