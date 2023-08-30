# install last stagedtrees from github
#remotes::install_github("stagedtrees/stagedtrees")
library(stagedtrees)

## load required packages, all available from CRAN
library(bnlearn)
library(xtable)
library(caret)
library(pbapply)
iscliavailable <- require(cli)  ## for progress bar and messages


extract_labels <- function(aldag){
  cc_pp <- sum(sapply(aldag, function(x) sum(x$context %in% x$partial)))

    
    cc <- sum(sapply(aldag, function(x) sum(length(x$context)))) - cc_pp
  pp <- sum(sapply(aldag, function(x) sum(length(x$partial)))) - cc_pp
  
  ll <- sum(sapply(aldag, function(x) sum(length(x$local))))
  all <- sum(sapply(aldag, function(x) length(x$parents)))
  tt <- all - (ll + cc + pp + cc_pp)
  return(c(cc_pp = cc_pp, cc = cc, pp = pp, ll = ll, tt = tt, all = all))
}



datasets <- list.files("dataset_ready/", include.dirs = FALSE)
results <- matrix(NA, nrow = length(datasets), ncol = 15 + 18, dimnames = list(datasets, 
                                                                          c("variables", "atomic events", 
                                                                            "time_dag.elapsed", "time_bhc.elapsed", "time_csbhc.elapsed",
                                                                            "time_aldag", 
                                                                            "bic.dag", "bic.bhc", "bic.csbhc",
                                                                            "ll_test.dag", "ll_test.bhc", "ll_test.csbhc",
                                                                            "pvals.dag", "pvals.bhc", "pvals.csbhc", ## 15
                                                                            "edges_bhc.cc_pp", "edges_bhc.cc", "edges_bhc.pp", "edges_bhc.ll", "edges_bhc.tt", "edges_bhc.all",
                                                                            "edges_csbhc.cc_pp", "edges_csbhc.cc", "edges_csbhc.pp", "edges_csbhc.ll", "edges_csbhc.tt", "edges_csbhc.all",
                                                                            "edges_dag.cc_pp", "edges_dag.cc", "edges_dag.pp", "edges_dag.ll", "edges_dag.tt", "edges_dag.all"
                                                                            )))


resultsMed <- results
set.seed(1)
if (iscliavailable){
  cli::cli_progress_bar(total = length(datasets))
}
for(d in datasets) {
  d.name <- substr(d, 1, nchar(d)-4)
  if (iscliavailable){
    cli::cli_progress_output("working on {d}...")
  }
  
  data_all <- readRDS(file.path("dataset_ready",d))
  
  folds <- createFolds(data_all[,1], k = 20, returnTrain = TRUE)
  
  atomic_events <- 1
  for(i in 1:NCOL(data_all)) {
    atomic_events <- atomic_events * length(levels(data_all[, i]))
  }
  
  results[d, "variables"] <- NCOL(data_all)
  resultsMed[d, "variables"] <- NCOL(data_all)
  results[d, "atomic events"] <- atomic_events
  resultsMed[d, "atomic events"] <- atomic_events
  
  
  
  res <- pbapply::pblapply(folds, function(trainix){
    data <- data_all[trainix,]
    test <-  data_all[-trainix,]
    time_dag <- system.time(dag <- tabu(data, tabu = 10))[3]
    bn <- bn.fit(dag, data)
    order <- node.ordering(bn)
    dag_edges <- NROW(dag$arcs)
    dag_st <- bn |> as_sevt(order = order) |> sevt_fit(data = data, lambda = 0)
    
    time_bhc <- system.time(tree1 <- dag_st |> stages_bhc())[3]
    time_csbhc <-  system.time(tree2 <- dag_st |> stages_csbhc())[3]
    time_aldag1 <- system.time(aldag1 <- as_parentslist(tree1, silent = TRUE))[3]
    time_aldag2 <- system.time(aldag2 <- as_parentslist(tree2, silent = TRUE))[3]
    time_aldag <- mean(c(time_aldag1, time_aldag2))
                               
    alls <- list(dag = dag_st, bhc = tree1, csbhc = tree2)
    
    bics <- sapply(alls, BIC)
    ll_test <- sapply(alls, function(mm) sum(prob(mm, test, log = TRUE)))
    
    #pvals <- sapply(alls, function(mm) lr_test(mm, dag_st)$`Pr(>Chisq)`[2])
    
    edges_bhc  <- extract_labels(aldag1)
    edges_csbhc <- extract_labels(aldag2)
    edges_dag <- extract_labels(as_parentslist(dag_st))
    
    return(c(time_dag = time_dag,
             time_bhc = time_bhc,
             time_csbhc = time_csbhc,
             time_aldag = time_aldag,
             bic = bics,
             ll_test = ll_test,
             #pvals = pvals,
             edges_bhc = edges_bhc,
             edges_csbhc = edges_csbhc,
             edges_dag = edges_dag))
  }, cl = 5)
  
  resm <- colMeans(t(as.data.frame(res))) 
  results[d, names(resm)] <- resm
  resultsMed[d, names(resm)] <- apply(t(as.data.frame(res)), 2, median)
  
  if (iscliavailable){
    cli::cli_progress_update(1)
  }
  warnings()
}


### save both average and median results
saveRDS(results, "results.rds")
saveRDS(resultsMed, "resultsMed.rds")

### use median 
results <- readRDS("resultsMed.rds")

nice <- data.frame(matrix(NA, nrow = length(datasets), ncol = 1, 
                          dimnames = list(sapply(rownames(results), \(d) substr(d, 1, nchar(d)-4)), 
                                          c("variables"))), check.names = FALSE)

nice[, "variables"] <- results[, "variables"]


nice[, "Edg. ALDAG-BHC"] <- paste0("(",
                                  results[, "edges_bhc.tt"], ",",
                                  results[, "edges_bhc.cc"], ",",
                                  results[, "edges_bhc.pp"], ",",
                                  results[, "edges_bhc.cc_pp"], ",",
                                  results[, "edges_bhc.ll"],
                                  ")")

nice[, "Edg. ALDAG-CSBHC"] <- paste0("(",
                                   results[, "edges_csbhc.tt"], ",",
                                   results[, "edges_csbhc.cc"], ",",
                                   results[, "edges_csbhc.pp"], ",",
                                   results[, "edges_csbhc.cc_pp"], ",",
                                   results[, "edges_csbhc.ll"],
                                   ")")

tabnice <- xtable(nice, digits = c(0, 0, 0, 0), auto = TRUE)
print(tabnice, booktabs = TRUE, file = "results_nice.tex", scalebox = 0.8, floating = FALSE)


nice <- data.frame(matrix(NA, nrow = length(datasets), ncol = 1, 
                          dimnames = list(sapply(rownames(results), \(d) substr(d, 1, nchar(d)-4)), 
                                          c("BIC DAG"))), check.names = FALSE)
nice[, "BIC DAG"] <- results[, "bic.dag"] 
nice[, "BIC BHC"] <- results[, "bic.bhc"] 
nice[, "BIC CSBHC"] <- results[, "bic.csbhc"] 

nice[, "log-lik test DAG"] <- results[, "ll_test.dag"]
nice[, "log-lik test BHC"] <-  results[, "ll_test.bhc"] 
nice[, "log-lik test CSBHC"] <-  results[, "ll_test.csbhc"] 


tabnice <- xtable(nice, digits = c(2, 2, 2, 2, 2, 2, 2), auto = TRUE)
print(tabnice, booktabs = TRUE, file = "results_nice_2.tex", scalebox = 0.8, floating = FALSE)


nice_times <- results[, c("time_dag.elapsed", "time_bhc.elapsed", "time_csbhc.elapsed",
                          "time_aldag")]
colnames(nice_times) <- c("DAG", "BHC", "CSBHC", "ALDAG")
rownames(nice_times) <- sapply(rownames(nice_times), \(d) substr(d, 1, nchar(d)-4))
tabnicetimes <-  xtable(nice_times, auto = TRUE)
print(tabnicetimes, booktabs = TRUE, file = "results_nice_times.tex", scalebox = 0.8, floating = FALSE)


tabnicesmall <- xtable(cbind(nice[, c(1,2,4,5)], nice_times[,c(1,2,4)]),
                       digits = c(0, 0, 2, 2, 2, 2, 2, 2))
print(tabnicesmall, booktabs = TRUE, file = "results_nice_small.tex", scalebox = 0.8, floating = FALSE)


