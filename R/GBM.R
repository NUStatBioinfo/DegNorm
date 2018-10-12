####################################################################################
# loading required library
library(foreach)
library(doParallel)
####################################################################################
# loading NMF functions
source("/projects/b1042/WangLab/binx/DegNorm/NMF_functions.R")

####################################################################################
# loading coverage and count data

load("/projects/b1042/WangLab/binx/DegNorm/GBM_coverage.Rdata")

counts = read.table("/projects/b1042/WangLab/binx/DegNorm/GBM_filtered_9.txt")

####################################################################################
# Estimating sequencing depth by SVD.

rho = NULL

cl = makeCluster(12)
registerDoParallel(cl)
rho = foreach(a = row.names(counts), .combine = "rbind") %dopar% {
  ratioSVD(f = coverage[[a]])
}
stopCluster(cl)

# calculate the factor based on genes with low DI scores from SVD output
print("Number of genes used in sequencing depth calculation:")
print(sum(apply(rho,1,max) < 0.1))

norm.factor = colSums(counts[apply(rho,1,max) < 0.1,])
norm.factor = norm.factor/median(norm.factor)
scale = norm.factor

# log the scale for sequencing depth from SVD
print("SVD sequencing depth scale:")
print(scale)

# scale the raw counts by sequencing depth
counts.weighted = sweep(counts, 2, norm.factor, FUN = "/")
scale = norm.factor

####################################################################################
# DegNorm pipeline:
for(iter in 1:5){
  
  rho = NULL
  cl = makeCluster(12)
  registerDoParallel(cl)
  
  rho = foreach(a = row.names(counts), .combine = "rbind") %dopar% {
    optiNMF(f = coverage[[a]], norm.factor = scale)$rho
  }
  stopCluster(cl)
  
  row.names(rho) = row.names(counts)
  colnames(rho) = colnames(counts)
  
  # cap DI scores at 0.9 to reduce maximum inflation
  rho[rho > 0.9] = 0.9
  
  # calculate the DegNorm adjusted counts after removing sequencing depth and degradation effect
  adjusted = counts.weighted / (1 - rho)
  
  # For genes which are filtered out and didn't go through the optiNMF function, assign sample-average DI score to them.
  gene.out = (apply(rho,1,max) == 0)
  DI.sample = 1 - colSums(counts.weighted)/colSums(adjusted)
  rho[gene.out,] = matrix(rep(DI.sample, sum(gene.out)), byrow = TRUE, ncol = 9)
  print("Number of genes filtered out:")
  print(sum(gene.out))

  # Re-calculate the DegNorm adjusted counts after considering those filtered-out genes
  adjusted = counts.weighted / (1 - rho)
  
  # adjust degradation for sequencing depth estimation
  norm.factor = colSums(adjusted)
  norm.factor = norm.factor / median(norm.factor)
  
  print("Degradation adjusted:")
  print(norm.factor)
  
  # update the read counts matrix by adjusting degradation effect to sequencing depth
  counts.weighted = sweep(counts.weighted, 2, norm.factor, FUN = "/")
  
  # adjust sequencing depth factor (scale) for the next iteration (as input to optiNMF)
  scale = scale * norm.factor
  print("Sequencing depth updates:")
  print(scale)
  
  # output DI scores and adjusted counts
  write.table(adjusted, file = "/projects/b1042/WangLab/binx/DegNorm/GBM.txt")
  write.table(rho, file = "/projects/b1042/WangLab/binx/DegNorm/GBM_DI.txt")
  
}
