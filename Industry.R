options(download.file.method = "wininet")
devtools::install_github("stagedtrees/stagedtrees")
setwd("C:/Users/manuele.leonelli/OneDrive - IE University/Documentos/GitHub/causal_stagedtrees")
library(stagedtrees)
library(bnlearn)
library(qgraph)
library(igraph)
library(dplyr)
library(caret)
library(ggplot2)

palette("Tableau")


albn <- function(tree,...){
  pl <- as_parentslist(tree)
  G <- graph_from_adjacency_matrix(as_adj_matrix(pl))
  el <- as_edgelist(G)
  col <- c()
  for (i in 1:nrow(el)){
    col[i] <- "black"
    if (el[i,1] %in% pl[[el[i, 2]]]$context){
      if (el[i,1] %in% pl[[el[i, 2]]]$partial){
        col[i] <- "purple"
      }else{
        col[i] <- "red"
      }
    }else{
      if (el[i,1] %in% pl[[el[i, 2]]]$partial){
        col[i] <- "blue"
      }
      if (el[i,1] %in% pl[[el[i, 2]]]$local){
        col[i]<- "green"
      }
    }
  }
  qgraph(el, edge.color = col,...)
}


data <-  read.delim("C:/Users/manuele.leonelli/Dropbox/Desktop/Dati/CIS_Microdati_Anno_2012.txt")

## Data Preparation

data$GP <- factor(data$GP,levels = c(0,1), labels = c("No","Yes"))
data$LARMAR <- factor(data$LARMAR,levels = c("A","B","C","D"),labels = c("Regional","National","International","International")) 
data$INPDGD <- factor(data$INPDGD, levels = c(0,1),labels = c("No","Yes"))
data$INPDSV <- factor(data$INPDSV, levels = c(0,1),labels= c("No","Yes"))
data$INPD <- rep(0,nrow(data))
for(i in 1:nrow(data)) data$INPD[i] <- any(data[i,32:34]) == 1
data$INPD <- factor(data$INPD, levels = c(0,1),labels = c("No","Yes"))
data$INABA <- factor(data$INABA, levels = c(0,1), labels = c("No","Yes"))
data$INONG <- factor(data$INONG, levels = c(0,1), labels = c("No","Yes"))
data$CO <- factor(data$CO, levels = c(0,1), labels = c("No","Yes"))
data$ORG <- rep(0,nrow(data))
for(i in 1:nrow(data)) data$ORG[i] <- any(data[i,123:125]) == 1
data$ORG <- factor(data$ORG, levels = c(0,1),labels = c("No","Yes"))
data$MKT <- rep(0,nrow(data))
for(i in 1:nrow(data)) data$MKT[i] <- any(data[i,126:129]) == 1
data$MKT <- factor(data$MKT, levels = c(0,1),labels = c("No","Yes"))
data$PUB <- rep(0,nrow(data))
for(i in 1:nrow(data)) data$PUB[i] <- any(data[i,130:131]) == 1
data$PUB <- factor(data$PUB, levels = c(0,1),labels = c("No","Yes"))
data$emp12 <- factor(data$emp12, levels = c(10,50,250), labels = c("10-49","50-249",">250"))
data$EMPUD <- factor(data$EMPUD, levels = c(1,2,3,4,5,6,7), labels = c("0%","<10%","<10%",">10%",">10%",">10%",">10%"))
data$GROWTH <- factor(ifelse(data$Gturn12 >1,"Yes","No" ))

data$RR<- rep(0,nrow(data))
for(i in 1:nrow(data)) data$RR[i] <- any(data[i,c(42,44:50)]) == 1
data$RR <- factor(data$RR, levels = c(0,1),labels = c("No","Yes"))
data$EMP <- data$emp12
data$RD <- data$RR
dad <- data %>% select(c('GP',  'LARMAR','INPDGD','INPDSV','INPD','INABA','INONG','CO','ORG','MKT','PUB','EMP','EMPUD','GROWTH','RD'))
dad <- na.omit(dad)

## Learning BN
grafo <- hc(dad)
dag <- bn.fit(hc(dad),dad)
-2*BIC(hc(dad),dad)

## Learning ALDAG BHC
tree3 <- stages_bhc(sevt_fit(as_sevt(dag),data = dad,lambda = 0))
albn(tree3)

tree4 <- stages_csbhc(sevt_fit(as_sevt(dag),data = dad,lambda = 0))
albn(tree4)


png("test1.png", units="in", width=6, height=2.6, res=1200)
albn(tree4)
dev.off()


try <- dad
set.seed(2023)
rows <- createFolds(try[,1],k=5,list = T,returnTrain = F)

loglik_bn <- loglik_hc <- loglik_csbhc <- rep(0,5)
for(i in 1:5){
  bn <- bn.fit(hc(try[-rows[[i]],]), try[-rows[[i]],])
  tree <- stages_bhc(sevt_fit(as_sevt(bn),data = try,lambda = 0))
  tree2 <- stages_csbhc(sevt_fit(as_sevt(bn),data = try,lambda = 0))
  loglik_hc[i] <- logLik(sevt_fit(tree,data = try[rows[[i]],], lambda = 0))
  loglik_csbhc[i] <- logLik(sevt_fit(tree2,data = try[rows[[i]],], lambda = 0))
  loglik_bn[i] <- logLik(bn.fit(hc(try[-rows[[i]],]), try[rows[[i]],]), try[rows[[i]],])
  print(i)
}

value <- c(loglik_hc,loglik_csbhc,loglik_bn,pred_hc,pred_csbhc,pred_bn,pred_hc_ov,pred_csbhc_ov,pred_bn_ov)
model <- rep(rep(c("ST_HC","ST_CSBHC","BN"),each =5),3)
score <- rep(c("Pred_LogLik","Acc_Life","Overall_Acc"), each = 15)
results <- data.frame(value,model,score)

ggplot(results, aes(model,value, color=model)) + geom_boxplot()+  facet_wrap( ~ score,scales = "free_y")
