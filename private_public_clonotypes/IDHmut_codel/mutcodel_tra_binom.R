library(tidyr)
# mutcodel tra private
mutcodel_tra_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_tra_private.csv")
mutcodel_tra_private <- mutcodel_tra_private[,c(1,2)]
mutcodel_tra_private <- unique(mutcodel_tra_private) 

mutcodel_tra_private$first <- floor(mutcodel_tra_private$len/2)
mutcodel_tra_private$position1 <- substr(mutcodel_tra_private$CDR3.aa,mutcodel_tra_private$first,mutcodel_tra_private$first)
mutcodel_tra_private$position2 <- substr(mutcodel_tra_private$CDR3.aa,mutcodel_tra_private$first+1,mutcodel_tra_private$first+1)
mutcodel_tra_private$position3 <- substr(mutcodel_tra_private$CDR3.aa,mutcodel_tra_private$first+2,mutcodel_tra_private$first+2)

p1_private <- sum(mutcodel_tra_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_tra_private$position1)
p2_private <- sum(mutcodel_tra_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_tra_private$position2)
p3_private <- sum(mutcodel_tra_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_tra_private$position3)

#mutcodel tra public
mutcodel_tra_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_tra_public.csv")
mutcodel_tra_public <- mutcodel_tra_public[,c(1,2)]
mutcodel_tra_public <- unique(mutcodel_tra_public) 

mutcodel_tra_public$first <- floor(mutcodel_tra_public$len/2)
mutcodel_tra_public$position1 <- substr(mutcodel_tra_public$CDR3.aa,mutcodel_tra_public$first,mutcodel_tra_public$first)
mutcodel_tra_public$position2 <- substr(mutcodel_tra_public$CDR3.aa,mutcodel_tra_public$first+1,mutcodel_tra_public$first+1)
mutcodel_tra_public$position3 <- substr(mutcodel_tra_public$CDR3.aa,mutcodel_tra_public$first+2,mutcodel_tra_public$first+2)

p1_public <- sum(mutcodel_tra_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_tra_public$position1)
p2_public <- sum(mutcodel_tra_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_tra_public$position2)
p3_public <- sum(mutcodel_tra_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_tra_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/proportion_hydro/mutcodel_tra_line_data.csv",row.names = FALSE ) 

#p1
binom.test(174,391,p=0.75) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(201,391,p=0.25) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(143,391,p=0.25) #以public为期望的概率 #p-value = 4.285e-07

fdr <- p.adjust(c("2.2e-16","2.2e-16","4.285e-07"), "BH")
fdr
#3.300e-16 3.300e-16 4.285e-07
