library(tidyr)
# wt trb private
wt_trb_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_trb_private.csv")
wt_trb_private <- wt_trb_private[,c(1,2)]
wt_trb_private <- unique(wt_trb_private) 

wt_trb_private$first <- floor(wt_trb_private$len/2)
wt_trb_private$position1 <- substr(wt_trb_private$CDR3.aa,wt_trb_private$first,wt_trb_private$first)
wt_trb_private$position2 <- substr(wt_trb_private$CDR3.aa,wt_trb_private$first+1,wt_trb_private$first+1)
wt_trb_private$position3 <- substr(wt_trb_private$CDR3.aa,wt_trb_private$first+2,wt_trb_private$first+2)

p1_private <- sum(wt_trb_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_trb_private$position1)
p2_private <- sum(wt_trb_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_trb_private$position2)
p3_private <- sum(wt_trb_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_trb_private$position3)

#wt trb public
wt_trb_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_trb_public.csv")
wt_trb_public <- wt_trb_public[,c(1,2)]
wt_trb_public <- unique(wt_trb_public) 

wt_trb_public$first <- floor(wt_trb_public$len/2)
wt_trb_public$position1 <- substr(wt_trb_public$CDR3.aa,wt_trb_public$first,wt_trb_public$first)
wt_trb_public$position2 <- substr(wt_trb_public$CDR3.aa,wt_trb_public$first+1,wt_trb_public$first+1)
wt_trb_public$position3 <- substr(wt_trb_public$CDR3.aa,wt_trb_public$first+2,wt_trb_public$first+2)

p1_public <- sum(wt_trb_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_trb_public$position1)
p2_public <- sum(wt_trb_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_trb_public$position2)
p3_public <- sum(wt_trb_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_trb_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/proportion_hydro/wt_trb_line_data.csv",row.names = FALSE ) 

#p1
binom.test(3329,5957,p=0.3571429) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(3564,5957,p=0.7142857) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(2944,5957,p=0.2857143) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("2.2e-16","2.2e-16","2.2e-16"), "BH")
