library(tidyr)
# mutcodel igk private
mutcodel_igk_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_igk_private.csv")
mutcodel_igk_private <- mutcodel_igk_private[,c(1,2)]
mutcodel_igk_private <- unique(mutcodel_igk_private) 

mutcodel_igk_private$first <- floor(mutcodel_igk_private$len/2)
mutcodel_igk_private$position1 <- substr(mutcodel_igk_private$CDR3.aa,mutcodel_igk_private$first,mutcodel_igk_private$first)
mutcodel_igk_private$position2 <- substr(mutcodel_igk_private$CDR3.aa,mutcodel_igk_private$first+1,mutcodel_igk_private$first+1)
mutcodel_igk_private$position3 <- substr(mutcodel_igk_private$CDR3.aa,mutcodel_igk_private$first+2,mutcodel_igk_private$first+2)

p1_private <- sum(mutcodel_igk_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igk_private$position1)
p2_private <- sum(mutcodel_igk_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igk_private$position2)
p3_private <- sum(mutcodel_igk_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igk_private$position3)

#mutcodel igk public
mutcodel_igk_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_igk_public.csv")
mutcodel_igk_public <- mutcodel_igk_public[,c(1,2)]
mutcodel_igk_public <- unique(mutcodel_igk_public) 

mutcodel_igk_public$first <- floor(mutcodel_igk_public$len/2)
mutcodel_igk_public$position1 <- substr(mutcodel_igk_public$CDR3.aa,mutcodel_igk_public$first,mutcodel_igk_public$first)
mutcodel_igk_public$position2 <- substr(mutcodel_igk_public$CDR3.aa,mutcodel_igk_public$first+1,mutcodel_igk_public$first+1)
mutcodel_igk_public$position3 <- substr(mutcodel_igk_public$CDR3.aa,mutcodel_igk_public$first+2,mutcodel_igk_public$first+2)

p1_public <- sum(mutcodel_igk_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igk_public$position1)
p2_public <- sum(mutcodel_igk_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igk_public$position2)
p3_public <- sum(mutcodel_igk_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igk_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/proportion_hydro/mutcodel_igk_line_data.csv",row.names = FALSE ) 

#p1
binom.test(1080,4385,p=0.2703704) #以public为期望的概率 #p-value = 0.0003121
#p2
binom.test(942,4385,p=0.1) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(2402,4385,p=0.437037) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("0.0003121","2.2e-16","2.2e-16"), "BH")
fdr
#3.121e-04 3.300e-16 3.300e-16