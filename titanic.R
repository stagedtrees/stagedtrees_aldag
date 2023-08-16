options(download.file.method = "wininet")
devtools::install_github("stagedtrees/stagedtrees")
library(stagedtrees)
library(bnlearn)
library(qgraph)
library(igraph)
palette("Tableau")

## Required Functions

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



## LEARN BN (Figure 1)
titanic.df <- as.data.frame(Titanic)
titanic.df <- titanic.df[rep(row.names(titanic.df), titanic.df$Freq), 1:4]
plot(hc(titanic.df))
-2*BIC(hc(titanic.df),titanic.df)

## LEARN GENERIC STAGED TREE (Figure 2)
tree_hc <- stages_hc(full(Titanic,order = c("Class","Sex","Survived","Age"),join_unobserved = FALSE))
plot(tree_hc)
BIC(tree_hc)

## STAGED TREE REPRESENTATION OF TITANIC BN (Figure 3)
plot(as_sevt(bn.fit(hc(titanic.df),titanic.df)))

## ALDAG OF GENERIC STAGED TREE (Figure 4)
as_parentslist(tree_hc)
albn(tree_hc)

## CSBHC STAGED TREE (Figure 5)
tree_csbhc <- stages_csbhc(full(Titanic,order = c("Class","Sex","Survived","Age"),join_unobserved = FALSE))
plot(tree_csbhc)
albn(tree_csbhc)
