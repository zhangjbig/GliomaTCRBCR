library(tidyr)
# wt igl private
wt_igl_private <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igl_private.csv")
wt_igl_private <- wt_igl_private[,c(1,2)]
wt_igl_private <- unique(wt_igl_private) 

wt_igl_private$first <- floor(wt_igl_private$len/2)
wt_igl_private$position1 <- substr(wt_igl_private$CDR3.aa,wt_igl_private$first,wt_igl_private$first)
wt_igl_private$position2 <- substr(wt_igl_private$CDR3.aa,wt_igl_private$first+1,wt_igl_private$first+1)
wt_igl_private$position3 <- substr(wt_igl_private$CDR3.aa,wt_igl_private$first+2,wt_igl_private$first+2)

p1_private <- sum(wt_igl_private$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igl_private$position1)
p2_private <- sum(wt_igl_private$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igl_private$position2)
p3_private <- sum(wt_igl_private$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igl_private$position3)

#wt igl public
wt_igl_public <- read.csv("/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/wt_igl_public.csv")
wt_igl_public <- wt_igl_public[,c(1,2)]
wt_igl_public <- unique(wt_igl_public) 

wt_igl_public$first <- floor(wt_igl_public$len/2)
wt_igl_public$position1 <- substr(wt_igl_public$CDR3.aa,wt_igl_public$first,wt_igl_public$first)
wt_igl_public$position2 <- substr(wt_igl_public$CDR3.aa,wt_igl_public$first+1,wt_igl_public$first+1)
wt_igl_public$position3 <- substr(wt_igl_public$CDR3.aa,wt_igl_public$first+2,wt_igl_public$first+2)

p1_public <- sum(wt_igl_public$position1 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igl_public$position1)
p2_public <- sum(wt_igl_public$position2 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igl_public$position2)
p3_public <- sum(wt_igl_public$position3 %in% c("M","F","W","V","L","I","A","P","G"))/length(wt_igl_public$position3)

line_data <- cbind(p1_private,p2_private,p3_private,p1_public,p2_public,p3_public)

write.csv(line_data,"/Users/wanglu/Documents/2021_TC/final_version2/private_public/private_public_wt/proportion_hydro/wt_igl_line_data.csv",row.names = FALSE ) 

#p1
binom.test(10999,36909,p=0.25) #以public为期望的概率 #p-value < 2.2e-16
#p2
binom.test(9310,36909,p=0.1824913) #以public为期望的概率 #p-value < 2.2e-16
#p3
binom.test(15756,36909,p=0.5313589) #以public为期望的概率 #p-value < 2.2e-16

fdr <- p.adjust(c("2.2e-16","2.2e-16","2.2e-16"), "BH")
fdr
