### preprocess raw data files from UCI

discretize <- function(x, l = 2){
  if (is.factor(x)) return(x)
  breaks <- quantile(x, probs = seq(0, 1, length.out = l + 1), na.rm = TRUE)
  cut(x, breaks = breaks, include.lowest = TRUE)
}

dir.create("dataset_ready", showWarnings = FALSE)

## abalone 
abalone <- read.table("datasets/abalone.data", sep = ",", stringsAsFactors = TRUE)
abalone <- na.exclude(abalone)
abalone <- data.frame(lapply(abalone, discretize))
saveRDS(abalone, "dataset_ready/abalone.rds")

## breast-cancer
bc <- read.table("datasets/breast-cancer.data", sep = ",", stringsAsFactors = TRUE, 
                 na.strings = "?")
bc <- na.exclude(bc)
## make V7 a factor
bc$V7 <- as.factor(bc$V7)

## regroup levels in V4
levels(bc$V4) <- list("0-19" = c("0-4", "5-9", "10-14", "15-19"),
                      "20-39" = c("20-24", "25-29", "30-34", "35-39"),
                      "40-54" = c("40-44", "45-49", "50-54"))
saveRDS(bc, "dataset_ready/breastcancer.rds")

## credit approval

cc <- read.table("datasets/crx.data", sep = ",", stringsAsFactors = TRUE, 
                 na.strings = "?")
cc <- na.exclude(cc)
cc$V11 <- cut(cc$V11, breaks = c(0,1,4,Inf), include.lowest = TRUE)
cc <- data.frame(lapply(cc, discretize))
cc <- cc[,-c(6,7)]
saveRDS(cc, "dataset_ready/creditapproval.rds")

## house votes

hv <- read.table("datasets/house-votes-84.data", sep = ",", na.strings = "?",
                 stringsAsFactors = TRUE)
hv <- na.exclude(hv)
saveRDS(hv, "dataset_ready/housevotes.rds")

## indian liver
 
inliv <- read.csv("datasets/Indian Liver Patient Dataset (ILPD).csv", 
                  header = FALSE, stringsAsFactors = TRUE)
inliv <- na.exclude(inliv)
inliv$V11 <- as.factor(inliv$V11)
inliv <- data.frame(lapply(inliv, discretize))

inliv <- inliv[,c(11, 1:10)]
saveRDS(inliv, "dataset_ready/indianliver.rds")

## nursery 

nurs <- read.table("datasets/nursery.data", sep = ",", stringsAsFactors = TRUE)
saveRDS(rev(nurs), "dataset_ready/nursery.rds")

## tic-tac-toe

ttt <- read.table("datasets/tic-tac-toe.data", sep = ",", stringsAsFactors = TRUE)
saveRDS(rev(ttt), "dataset_ready/tic-tac-toe.rds")

## PhdArticles 

data("PhDArticles")
saveRDS(PhDArticles, "dataset_ready/phdarticles.rds")

## Titanic
data("Titanic")
titanic <- as.data.frame(Titanic)
titanic <- titanic[rep(row.names(titanic), titanic$Freq), 1:4]

saveRDS(titanic, "dataset_ready/titanic.rds")

## 
