library(tidyr)
# wt igk private
wt_igk_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igk_private.csv")
wt_igk_private <- wt_igk_private[,c(1,2)]
wt_igk_private <- unique(wt_igk_private) 

wt_igk_private$first <- floor(wt_igk_private$len/2)
wt_igk_private$position1 <- substr(wt_igk_private$CDR3.aa,wt_igk_private$first,wt_igk_private$first)
wt_igk_private$position2 <- substr(wt_igk_private$CDR3.aa,wt_igk_private$first+1,wt_igk_private$first+1)
wt_igk_private$position3 <- substr(wt_igk_private$CDR3.aa,wt_igk_private$first+2,wt_igk_private$first+2)

p1_private <- sum(wt_igk_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igk_private$position1)
p2_private <- sum(wt_igk_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igk_private$position2)
p3_private <- sum(wt_igk_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igk_private$position3)

#wt igk public
wt_igk_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igk_public.csv")
wt_igk_public <- wt_igk_public[,c(1,2)]
wt_igk_public <- unique(wt_igk_public) 

wt_igk_public$first <- floor(wt_igk_public$len/2)
wt_igk_public$position1 <- substr(wt_igk_public$CDR3.aa,wt_igk_public$first,wt_igk_public$first)
wt_igk_public$position2 <- substr(wt_igk_public$CDR3.aa,wt_igk_public$first+1,wt_igk_public$first+1)
wt_igk_public$position3 <- substr(wt_igk_public$CDR3.aa,wt_igk_public$first+2,wt_igk_public$first+2)

p1_public <- sum(wt_igk_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igk_public$position1)
p2_public <- sum(wt_igk_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igk_public$position2)
p3_public <- sum(wt_igk_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igk_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/proportion_hydro/wt_igk_line_data.csv",row.names = FALSE ) 

#p1
binom.test(10287,38977,p=0.2520968) #以public为期望的概率 #p-value = 8.851e-08
#p2
binom.test(9130,38977,p=0.1548047) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(21647,38977,p=0.5291157) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("8.851e-08","2.2e-16","2.2e-16"), "BH")
fdr
