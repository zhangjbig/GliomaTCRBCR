library(tidyr)
# wt igh private
wt_igh_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igh_private.csv")
wt_igh_private <- wt_igh_private[,c(1,2)]
wt_igh_private <- unique(wt_igh_private) 

wt_igh_private$first <- floor(wt_igh_private$len/2)
wt_igh_private$position1 <- substr(wt_igh_private$CDR3.aa,wt_igh_private$first,wt_igh_private$first)
wt_igh_private$position2 <- substr(wt_igh_private$CDR3.aa,wt_igh_private$first+1,wt_igh_private$first+1)
wt_igh_private$position3 <- substr(wt_igh_private$CDR3.aa,wt_igh_private$first+2,wt_igh_private$first+2)

p1_private <- sum(wt_igh_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igh_private$position1)
p2_private <- sum(wt_igh_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igh_private$position2)
p3_private <- sum(wt_igh_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igh_private$position3)

#wt igh public
wt_igh_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igh_public.csv")
wt_igh_public <- wt_igh_public[,c(1,2)]
wt_igh_public <- unique(wt_igh_public) 

wt_igh_public$first <- floor(wt_igh_public$len/2)
wt_igh_public$position1 <- substr(wt_igh_public$CDR3.aa,wt_igh_public$first,wt_igh_public$first)
wt_igh_public$position2 <- substr(wt_igh_public$CDR3.aa,wt_igh_public$first+1,wt_igh_public$first+1)
wt_igh_public$position3 <- substr(wt_igh_public$CDR3.aa,wt_igh_public$first+2,wt_igh_public$first+2)

p1_public <- sum(wt_igh_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igh_public$position1)
p2_public <- sum(wt_igh_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igh_public$position2)
p3_public <- sum(wt_igh_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igh_public$position3)
#binom.test(445,500,p=0.85)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/proportion_hydro/wt_igh_line_data.csv",row.names = FALSE ) 

#p1
binom.test(40487,83129,p=0.5014493) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(41038,83129,p=0.4666667) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(38633,83129,p=0.4869565) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("2.2e-16","2.2e-16","2.2e-16"), "BH")
fdr
