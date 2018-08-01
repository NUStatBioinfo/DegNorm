####################################################################################
# loading required library
library(foreach)
library(doParallel)
####################################################################################
# loading NMF functions
source("~/NMF_functions.R")

####################################################################################
# loading data

load("~/GBM_coverage_9.Rdata")

counts = read.table("~/GBM_filtered_9.txt")
counts = counts[,c(1:3,7:9)]

####################################################################################
# Estimating library sizes / scaling factors

norm.factor = colSums(counts)/median(colSums(counts))
cl = makeCluster(12)
registerDoParallel(cl)

start.time = proc.time()
rho = foreach(a = row.names(counts), .combine = "rbind") %dopar% {
  ratioSVD(f = coverage[[a]][,c(1:3,7:9)], norm.factor = norm.factor)
}
stopCluster(cl)

norm.factor = colSums(counts[apply(rho,1,max) < 0.10,])
norm.factor = norm.factor/median(norm.factor)


####################################################################################
# Calculating ratios
for(iter in 1:5){
  
  rho = NULL
  cl = makeCluster(12)
  registerDoParallel(cl)
  
  rho = foreach(a = row.names(counts), .combine = "rbind", .errorhandling = 'remove') %dopar% {
    optiNMF(f = coverage[[a]][,c(1:3,7:9)], norm.factor = norm.factor)$rho
  }
  stopCluster(cl)
  
  row.names(rho) = row.names(counts)
  colnames(rho) = colnames(counts)
  
  rho[rho < 0.9] = 0.9
  adjusted = counts / (1 - rho)
  ratio = 1 - rho
  ratio[apply(ratio,1,min)==1,] = 
    matrix(rep(colSums(counts)/colSums(adjusted), sum(apply(ratio,1,min)==1)), byrow = TRUE, ncol = 6)
  adjusted = counts / ratio
  
  norm.factor = colSums(adjusted)
  norm.factor = norm.factor / median(norm.factor)
  
  print(norm.factor)
  
  adjusted = sweep(adjusted, 2, norm.factor, FUN = "/")

  write.table(adjusted, file = "~/Output/GBM_10V4_NDP.txt")
  write.table(colSums(adjusted), file = "~/Output/GBM_10V4_NDP_librarysize.txt", append = TRUE)
  write.table(rho, file = "~/Output/GBM_10V4_NDP_rho.txt")
  
}

