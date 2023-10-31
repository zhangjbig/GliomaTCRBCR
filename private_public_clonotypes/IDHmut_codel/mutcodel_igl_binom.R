library(tidyr)
# mutcodel igl private
mutcodel_igl_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_igl_private.csv")
mutcodel_igl_private <- mutcodel_igl_private[,c(1,2)]
mutcodel_igl_private <- unique(mutcodel_igl_private) 

mutcodel_igl_private$first <- floor(mutcodel_igl_private$len/2)
mutcodel_igl_private$position1 <- substr(mutcodel_igl_private$CDR3.aa,mutcodel_igl_private$first,mutcodel_igl_private$first)
mutcodel_igl_private$position2 <- substr(mutcodel_igl_private$CDR3.aa,mutcodel_igl_private$first+1,mutcodel_igl_private$first+1)
mutcodel_igl_private$position3 <- substr(mutcodel_igl_private$CDR3.aa,mutcodel_igl_private$first+2,mutcodel_igl_private$first+2)

p1_private <- sum(mutcodel_igl_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igl_private$position1)
p2_private <- sum(mutcodel_igl_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igl_private$position2)
p3_private <- sum(mutcodel_igl_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igl_private$position3)

#mutcodel igl public
mutcodel_igl_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/mutcodel_igl_public.csv")
mutcodel_igl_public <- mutcodel_igl_public[,c(1,2)]
mutcodel_igl_public <- unique(mutcodel_igl_public) 

mutcodel_igl_public$first <- floor(mutcodel_igl_public$len/2)
mutcodel_igl_public$position1 <- substr(mutcodel_igl_public$CDR3.aa,mutcodel_igl_public$first,mutcodel_igl_public$first)
mutcodel_igl_public$position2 <- substr(mutcodel_igl_public$CDR3.aa,mutcodel_igl_public$first+1,mutcodel_igl_public$first+1)
mutcodel_igl_public$position3 <- substr(mutcodel_igl_public$CDR3.aa,mutcodel_igl_public$first+2,mutcodel_igl_public$first+2)

p1_public <- sum(mutcodel_igl_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igl_public$position1)
p2_public <- sum(mutcodel_igl_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igl_public$position2)
p3_public <- sum(mutcodel_igl_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(mutcodel_igl_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_mutcodel/proportion_hydro/mutcodel_igl_line_data.csv",row.names = FALSE ) 

#p1
binom.test(886,3257,p=0.1584158) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(787,3257,p=0.1188119) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(1359,3257,p=0.5643564) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("2.2e-16","2.2e-16","2.2e-16"), "BH")
fdr

