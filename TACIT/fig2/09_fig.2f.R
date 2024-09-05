


################################ state expression for Figure 2f ########################################

setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/")

c2_1_exp <- read.csv("./02_chromHMM_2cell/2cell_1/2cell_1_all_state_expression.txt" ,sep = "\t")
c2_2_exp <- read.csv("./02_chromHMM_2cell/2cell_2/2cell_2_all_state_expression.txt" ,sep = "\t")
c4_1_exp <- read.csv("./03_chromHMM_4cell/4cell_1/4cell_1_all_state_expression.txt" ,sep = "\t")
c4_2_exp <- read.csv("./03_chromHMM_4cell/4cell_2/4cell_2_all_state_expression.txt" ,sep = "\t")
c4_3_exp <- read.csv("./03_chromHMM_4cell/4cell_3/4cell_3_all_state_expression.txt" ,sep = "\t")
c4_4_exp <- read.csv("./03_chromHMM_4cell/4cell_4/4cell_4_all_state_expression.txt" ,sep = "\t")

####### multivalent ############
multi <- c2_1_exp[which(c2_1_exp$state == "E10"),]
multi$group <- "Multi"

###### promoter weak ##########
b <- c4_3_exp[which(c4_3_exp$state == "E1"),]
c <- c4_4_exp[which(c4_4_exp$state == "E5"),]
pr_w <- rbind(b,c)
pr_w$group <- "Pr-W"

###### promoter strong ##########
a <- c2_1_exp[which(c2_1_exp$state == "E9" |c2_1_exp$state == "E1"),]
b <- c2_2_exp[which(c2_2_exp$state == "E1" |c2_1_exp$state == "E11"),]
c <- c4_1_exp[which(c4_1_exp$state == "E8"),]
d <- c4_2_exp[which(c4_2_exp$state == "E12"),]
e <- c4_4_exp[which(c4_4_exp$state == "E4"),]
pr_s <- rbind(a,b,c,d,e)
pr_s$group <- "Pr-S"

###### enhancer weak ##########
a <- c4_1_exp[which(c4_1_exp$state == "E12" ),]
b <- c4_2_exp[which(c4_2_exp$state == "E2" | c4_2_exp$state == "E3" |c4_2_exp$state == "E4"),]
en_w <- rbind(a,b)
en_w$group <- "En-W"

###### enhancer strong ##########
a <- c2_1_exp[which(c2_1_exp$state == "E12" |c2_1_exp$state == "E11" ),]
b <- c2_2_exp[which(c2_2_exp$state == "E12"|c2_2_exp$state == "E9"),]
c <- c4_1_exp[which(c4_1_exp$state == "E9"),]
d <- c4_2_exp[which(c4_2_exp$state == "E1"),]
e <- c4_3_exp[which(c4_3_exp$state == "E11" | c4_3_exp$state == "E12"),]
f <- c4_4_exp[which(c4_4_exp$state == "E1" | c4_4_exp$state == "E4" |c4_4_exp$state == "E2"),]
en_s <- rbind(a,b,c,d,e,f)
en_s$group <- "En-S"

##### genebody weak ############
a <- c2_1_exp[which(c2_1_exp$state == "E2" |c2_1_exp$state == "E8" ),]
b <- c4_1_exp[which(c4_1_exp$state == "E11"),]
c <- c4_3_exp[which(c4_3_exp$state == "E10"),]
ge_w <- rbind(a,b,c)
ge_w$group <- "Tr-W"

##### genebody strong ############
a <- c2_2_exp[which(c2_2_exp$state == "E10"|c2_2_exp$state == "E8"),]
b <- c4_1_exp[which(c4_1_exp$state == "E10"),]
c <- c4_2_exp[which(c4_2_exp$state == "E11"),]
d <- c4_3_exp[which(c4_3_exp$state == "E9"),]
e <- c4_4_exp[which(c4_4_exp$state == "E10"),]
ge_s <- rbind(a,b,c,d,e)
ge_s$group <- "Tr-S"

##### hetero k27me3 ###########
a <- c2_1_exp[which(c2_1_exp$state == "E5"|c2_1_exp$state == "E6"),]
b <- c2_2_exp[which(c2_2_exp$state == "E5"|c2_1_exp$state == "E6"),]
c <- c4_1_exp[which(c4_1_exp$state == "E1"),]
d <- c4_2_exp[which(c4_2_exp$state == "E10"),]
e <- c4_3_exp[which(c4_3_exp$state == "E4"),]
he_p <- rbind(a,b,c,d,e)
he_p$group <- "He-P"

##### hetero k9me3 ###########
c <- c4_1_exp[which(c4_1_exp$state == "E4"|c4_1_exp$state == "E6"),]
d <- c4_2_exp[which(c4_2_exp$state == "E8"),]
e <- c4_3_exp[which(c4_3_exp$state == "E6"),]
he_k9 <- rbind(c,d,e)
he_k9$group <- "He-K9"

##### hetero both ###########
a <- c2_1_exp[which(c2_1_exp$state == "E7"),]
b <- c2_2_exp[which(c2_2_exp$state == "E7"),]
c <- c4_1_exp[which(c4_1_exp$state == "E5"|c4_1_exp$state == "E2"),]
d <- c4_2_exp[which(c4_2_exp$state == "E9"|c4_2_exp$state == "E7"),]
e <- c4_3_exp[which(c4_3_exp$state == "E3"|c4_3_exp$state == "E2"),]
f <- c4_4_exp[which(c4_4_exp$state == "E12"|c4_4_exp$state == "E6"|c4_4_exp$state == "E7"),]
he_b <- rbind(a,b,c,d,e,f)
he_b$group <- "He"


################################ plot together #####################################
all <- rbind(multi, pr_w, pr_s, en_w, en_s, ge_w, ge_s, he_p, he_k9, he_b)
all$group <- factor(all$group, levels = c("Multi", "Pr-W", "Pr-S", "En-W", "En-S", "Tr-W", "Tr-S", "He-P", "He-K9", "He"))

pdf("Figure_2d_all_state_exp.pdf")
ggplot(data=all,mapping = aes(x = group, y = exp))+
  #geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(aes(fill = group),width = 0.6)+
  #ylim(0,14)+
  theme_bw()
dev.off()



########################################## only 4 cell ####################################
####### multivalent ############
multi <- c2_1_exp[which(c2_1_exp$state == "E10"),]
multi$group <- "Multi"

###### promoter weak ##########
b <- c4_3_exp[which(c4_3_exp$state == "E1"),]
c <- c4_4_exp[which(c4_4_exp$state == "E5"),]
pr_w <- rbind(b,c)
pr_w$group <- "Pr-W"

###### promoter strong ##########
b <- c4_1_exp[which(c4_1_exp$state == "E8"),]
c <- c4_2_exp[which(c4_2_exp$state == "E12"),]
d <- c4_4_exp[which(c4_4_exp$state == "E4"),]
pr_s <- rbind(b,c,d)
pr_s$group <- "Pr-S"

###### enhancer weak ##########
a <- c4_1_exp[which(c4_1_exp$state == "E12" ),]
b <- c4_2_exp[which(c4_2_exp$state == "E2" | c4_2_exp$state == "E3" |c4_2_exp$state == "E4"),]
en_w <- rbind(a,b)
en_w$group <- "En-W"

###### enhancer strong ##########
c <- c4_1_exp[which(c4_1_exp$state == "E9"),]
d <- c4_2_exp[which(c4_2_exp$state == "E1"),]
e <- c4_3_exp[which(c4_3_exp$state == "E11" | c4_3_exp$state == "E12"),]
f <- c4_4_exp[which(c4_4_exp$state == "E1" | c4_4_exp$state == "E4" |c4_4_exp$state == "E2"),]
en_s <- rbind(c,d,e,f)
en_s$group <- "En-S"

##### genebody weak ############
c <- c4_1_exp[which(c4_1_exp$state == "E11" ),]
d <- c4_3_exp[which(c4_3_exp$state == "E10"),]
e <- c4_4_exp[which(c4_4_exp$state == "E10"),]
ge_w <- rbind(d,e)
ge_w$group <- "Tr-W"

##### genebody strong ############
c <- c4_2_exp[which(c4_2_exp$state == "E11"),]
d <- c4_1_exp[which(c4_1_exp$state == "E10" ),]
e <- c4_3_exp[which(c4_3_exp$state == "E9"),]
ge_s <- rbind(c,d,e)
ge_s$group <- "Tr-S"

##### hetero k27me3 ###########
a <- c4_1_exp[which(c4_1_exp$state == "E1"),]
b <- c4_2_exp[which(c4_2_exp$state == "E10"),]
c <- c4_3_exp[which(c4_3_exp$state == "E4"),]
he_p <- rbind(a,b,c)
he_p$group <- "He-P"

##### hetero k9me3 ###########
c <- c4_1_exp[which(c4_1_exp$state == "E4"|c4_1_exp$state == "E6"),]
d <- c4_2_exp[which(c4_2_exp$state == "E8"),]
e <- c4_3_exp[which(c4_3_exp$state == "E6"),]
he_k9 <- rbind(c,d,e)
he_k9$group <- "He-K9"

##### hetero both ###########
b <- c4_1_exp[which(c4_1_exp$state == "E5"|c4_1_exp$state == "E2"),]
c <- c4_2_exp[which(c4_2_exp$state == "E9"|c4_2_exp$state == "E7"),]
d <- c4_3_exp[which(c4_3_exp$state == "E3"|c4_3_exp$state == "E2"),]
e <- c4_4_exp[which(c4_4_exp$state == "E12"|c4_4_exp$state == "E6"|c4_4_exp$state == "E7"),]
he_b <- rbind(b,c,d,e)
he_b$group <- "He"


################################ plot together #####################################
he_p <- he_p[which(he_p$exp <10),]
he_k9 <- he_k9[which(he_k9$exp <10),]
he_b <- he_b[which(he_b$exp <10),]

all <- rbind(multi, pr_w, pr_s, en_w, en_s, ge_w, ge_s, he_p, he_k9, he_b)
all$group <- factor(all$group, levels = c("Multi", "Pr-W", "Pr-S", "En-W", "En-S", "Tr-W", "Tr-S", "He-P", "He-K9", "He"))

pdf("Figure_2d_all_state_exp_no2cell.pdf")
ggplot(data=all,mapping = aes(x = group, y = exp))+
  #geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(aes(fill = group),width = 0.6)+
  ylim(0,15)+
  theme_bw()
dev.off()

