library(tidyr)
# mutcodel trb private
mutcodel_trb_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_trb_private.csv")
mutcodel_trb_private <- mutcodel_trb_private[,c(1,2)]
mutcodel_trb_private <- unique(mutcodel_trb_private) 

mutcodel_trb_private$first <- floor(mutcodel_trb_private$len/2)
mutcodel_trb_private$position1 <- substr(mutcodel_trb_private$CDR3.aa,mutcodel_trb_private$first,mutcodel_trb_private$first)
mutcodel_trb_private$position2 <- substr(mutcodel_trb_private$CDR3.aa,mutcodel_trb_private$first+1,mutcodel_trb_private$first+1)
mutcodel_trb_private$position3 <- substr(mutcodel_trb_private$CDR3.aa,mutcodel_trb_private$first+2,mutcodel_trb_private$first+2)

p1_private <- sum(mutcodel_trb_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_trb_private$position1)
p2_private <- sum(mutcodel_trb_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_trb_private$position2)
p3_private <- sum(mutcodel_trb_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_trb_private$position3)

#mutcodel trb public
mutcodel_trb_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_trb_public.csv")
mutcodel_trb_public <- mutcodel_trb_public[,c(1,2)]
mutcodel_trb_public <- unique(mutcodel_trb_public) 

mutcodel_trb_public$first <- floor(mutcodel_trb_public$len/2)
mutcodel_trb_public$position1 <- substr(mutcodel_trb_public$CDR3.aa,mutcodel_trb_public$first,mutcodel_trb_public$first)
mutcodel_trb_public$position2 <- substr(mutcodel_trb_public$CDR3.aa,mutcodel_trb_public$first+1,mutcodel_trb_public$first+1)
mutcodel_trb_public$position3 <- substr(mutcodel_trb_public$CDR3.aa,mutcodel_trb_public$first+2,mutcodel_trb_public$first+2)

p1_public <- sum(mutcodel_trb_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_trb_public$position1)
p2_public <- sum(mutcodel_trb_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_trb_public$position2)
p3_public <- sum(mutcodel_trb_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_trb_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/proportion_hydro/mutcodel_trb_line_data.csv",row.names = FALSE ) 

#p1
binom.test(453,796,p=0) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(442,796,p=0.5) #以public为期望的概率 #p-value = 0.002026
#p3
binom.test(399,796,p=0.5) #以public为期望的概率 #p-value = 0.9717

fdr <- p.adjust(c("2.2e-16","0.002026","0.9717"), "BH")
fdr
#6.600e-16 3.039e-03 9.717e-01