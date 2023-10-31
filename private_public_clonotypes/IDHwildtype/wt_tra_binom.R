library(tidyr)
# wt tra private
wt_tra_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_tra_private.csv")
wt_tra_private <- wt_tra_private[,c(1,2)]
wt_tra_private <- unique(wt_tra_private) 

wt_tra_private$first <- floor(wt_tra_private$len/2)
wt_tra_private$position1 <- substr(wt_tra_private$CDR3.aa,wt_tra_private$first,wt_tra_private$first)
wt_tra_private$position2 <- substr(wt_tra_private$CDR3.aa,wt_tra_private$first+1,wt_tra_private$first+1)
wt_tra_private$position3 <- substr(wt_tra_private$CDR3.aa,wt_tra_private$first+2,wt_tra_private$first+2)

p1_private <- sum(wt_tra_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_tra_private$position1)
p2_private <- sum(wt_tra_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_tra_private$position2)
p3_private <- sum(wt_tra_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_tra_private$position3)

#wt tra public
wt_tra_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_tra_public.csv")
wt_tra_public <- wt_tra_public[,c(1,2)]
wt_tra_public <- unique(wt_tra_public) 

wt_tra_public$first <- floor(wt_tra_public$len/2)
wt_tra_public$position1 <- substr(wt_tra_public$CDR3.aa,wt_tra_public$first,wt_tra_public$first)
wt_tra_public$position2 <- substr(wt_tra_public$CDR3.aa,wt_tra_public$first+1,wt_tra_public$first+1)
wt_tra_public$position3 <- substr(wt_tra_public$CDR3.aa,wt_tra_public$first+2,wt_tra_public$first+2)

p1_public <- sum(wt_tra_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_tra_public$position1)
p2_public <- sum(wt_tra_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_tra_public$position2)
p3_public <- sum(wt_tra_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_tra_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/proportion_hydro/wt_tra_line_data.csv",row.names = FALSE ) 

#p1
binom.test(1254,2757,p=0.4242424) #以public为期望的概率 #p-value = 0.001205
#p2
binom.test(1401,2757,p=0.3333333) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(1143,2757,p=0.4242424) #以public为期望的概率 #p-value = 0.3072

fdr <- p.adjust(c("0.001205","2.2e-16","0.3072"), "BH")
fdr
#1.8075e-03 6.6000e-16 3.0720e-01
