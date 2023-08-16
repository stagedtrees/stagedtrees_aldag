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


## Data Preparation
data <- read.csv("D_SN147.dat")

try <- data %>% select("FIDUCIA","SESSO","SPOCON","AMICI2","VOTOVI") 
try <- na.omit(try)
try[,1] <- factor(try[,1],levels = c(1,2), labels = c("Yes","No"))
try[,2] <- factor(try[,2],levels = c(1,2), labels = c("Male","Female"))
try[,3] <- factor(try[,3],levels = c(1,2),labels = c("No","Yes"))
try[,4] <- factor(try[,4],levels = c(1,2,3),labels = c("No","Yes","Unsure"))
try[,5] <- factor(try[,5],levels = c(0,1,2,3,4,5,6,7,8,9,10),labels = c("Little","Little","Little","Little","Little","Little","Average","Average","Average","A lot","A lot"))

# Learning DAG
wl = matrix(c("VOTOVI", "FIDUCIA", "VOTOVI", "AMICI2","VOTOVI","SPOCON","VOTOVI","SESSO"), ncol = 2, byrow = TRUE,      dimnames = list(NULL, c("from", "to")))
-2*BIC(hc(try,blacklist = wl),try)
plot(hc(try,blacklist=wl))

# Learning generic staged tree (Figure 7 left)
dat <- try[,1:4]
best <- search_best(data = dat, alg = stages_hc, lambda = 0, join_unobserved = FALSE)
tree <- stages_hc(full(try,order = c(names(best$tree),"VOTOVI"), join_unobserved = FALSE))
BIC(tree)

# Learning generic context-specific tree (Figure 7 right)
best2 <- search_best(data = dat, alg = stages_csbhc, lambda = 0, join_unobserved = FALSE)
tree2 <- stages_csbhc(full(try,order = c(names(best2$tree),"VOTOVI"), join_unobserved = FALSE))
BIC(tree2)

## ALDAGS
albn(tree)
albn(tree2)

## Cross-Validation
set.seed(2023)
rows <- createFolds(try[,1],k=5,list = T,returnTrain = F)
pred_bn <- pred_hc <- pred_csbhc <-  rep(0,5)
pred_bn_ov <- pred_hc_ov <- pred_csbhc_ov <- rep(0,5)
loglik_bn <- loglik_hc <- loglik_csbhc <- rep(0,5)
for(i in 1:5){
  bn <- bn.fit(hc(try[-rows[[i]],],blacklist = wl), try[-rows[[i]],])
  best <- search_best(data = dat[-rows[[i]],], alg = stages_hc, lambda = 0, join_unobserved = FALSE)
  tree <- stages_hc(full(try[-rows[[i]],],order = c(names(best$tree),"VOTOVI"), join_unobserved = FALSE))
  best2 <- search_best(data = dat[-rows[[i]],], alg = stages_csbhc, lambda = 0, join_unobserved = FALSE)
  tree2 <- stages_csbhc(full(try[-rows[[i]],],order = c(names(best2$tree),"VOTOVI"), join_unobserved = FALSE))
  loglik_hc[i] <- logLik(sevt_fit(tree,data = try[rows[[i]],], lambda = 0))
  loglik_csbhc[i] <- logLik(sevt_fit(tree2,data = try[rows[[i]],], lambda = 0))
  loglik_bn[i] <- logLik(bn.fit(hc(try[-rows[[i]],],blacklist = wl), try[rows[[i]],]), try[rows[[i]],])
  tree_bn <- sevt_fit(as_sevt(bn), try[-rows[[i]],], lambda = 0)
  pred_bn[i] <- mean(predict(tree_bn, newdata = try[rows[[i]],], class= "VOTOVI") == try[rows[[i]],"VOTOVI"])
  pred_hc[i] <- mean(predict(tree, newdata = try[rows[[i]],], class= "VOTOVI") == try[rows[[i]],"VOTOVI"])
  pred_csbhc[i] <- mean(predict(tree, newdata = try[rows[[i]],], class= "VOTOVI") == try[rows[[i]],"VOTOVI"])
  for(j in colnames(try)){
    pred_bn_ov[i] <- pred_bn_ov[i] + mean(predict(tree_bn, newdata = try[rows[[i]],], class= j) == try[rows[[i]],j])
    pred_hc_ov[i] <- pred_hc_ov[i] + mean(predict(tree, newdata = try[rows[[i]],], class= j) == try[rows[[i]],j])
    pred_csbhc_ov[i] <- pred_csbhc_ov[i] + mean(predict(tree2, newdata = try[rows[[i]],], class= j) == try[rows[[i]],j])
  }
  pred_bn_ov[i] <- pred_bn_ov[i]/5
  pred_hc_ov[i] <- pred_hc_ov[i]/5
  pred_csbhc_ov[i] <- pred_csbhc_ov[i]/5
}

value <- c(loglik_hc,loglik_csbhc,loglik_bn,pred_hc,pred_csbhc,pred_bn,pred_hc_ov,pred_csbhc_ov,pred_bn_ov)
model <- rep(rep(c("ST_HC","ST_CSBHC","BN"),each =5),3)
score <- rep(c("Pred_LogLik","Acc_Life","Overall_Acc"), each = 15)
results <- data.frame(value,model,score)

ggplot(results, aes(model,value, color=model)) + geom_boxplot()+  facet_wrap( ~ score,scales = "free_y")
