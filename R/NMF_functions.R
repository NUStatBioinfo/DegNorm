################################################################################################################
# NMF-OA functions
# Non-negative matrix factorization with over-approximation.
# input: coverage matrix (f) with columns representing samples
# output: truncated SVD (svd), fitted coverage matrix (fitted), residuals (res), Lagrangian multipler (Lambda).
################################################################################################################

NMF <- function(f, loop = 100){
  
  Lambda = matrix(0, ncol = dim(f)[2], nrow = dim(f)[1])
  model.svd = svd(f, nv = 1, nu = 1)
  fitted = model.svd$d[1] * model.svd$u %*% t(model.svd$v)
  res = fitted - f
  
  for(i in 1:loop){
    
    Lambda = Lambda - 1/sqrt(loop) * res
    Lambda[which(Lambda < 0)] = 0
    
    model.svd = svd(f+Lambda, nv = 1, nu = 1)
    fitted = model.svd$d[1] * model.svd$u %*% t(model.svd$v)
    res = fitted - f
    
  }
  
  return(list(svd = model.svd, fitted = fitted, res = res, Lambda = Lambda))
}

################################################################################################################
# baseline selection process to improve NMF
# input: coverage matrix (f), normalization factor for sequencing depth (norm.factor)
# output: DI scores (rho), baseline status (convergence) and abundance level (K).
################################################################################################################

optiNMF <- function(f, norm.factor, loop=100){
  
  # initialize the output with default values
  output = list(rho = NULL, convergence = FALSE, K = NULL)
  num.sample = dim(f)[2]
  output$rho = rep(0, num.sample)
  output$ran_baseline_selection <- FALSE
  
  f = t(as.matrix(f))
  # normalize sequencing depth based on the norm.factor
  f = f/norm.factor
  
  # filter out bases where the coverage are all low across all samples
  filter = ((apply(f,2,max) > 0.1 * max(apply(f,2,max))) == 1)
  
  # if there are still more than 50 bp left after filtering, proceed; otherwise return default values
  if(sum(filter) < 50){
    return(output)
  }else{
    
    # keep the filtered coverage and calculate the gene transcript length after filtering
    f = f[, filter]
    gene.length = dim(f)[2]
    
    # keep the filtered coverage as f.check for further calculation of abundance values K
    f.check = f

    # if there is any sample without coverage scores after filtering, return the default values
    if(sum(rowSums(f) > 0) < num.sample){ 
      return(output)
    }

    # run NMF to the filtered coverage and obtain the fitted values and the Lagrangian multiplier
    NMF.output = NMF(f, loop=loop)
    fitted = NMF.output$fitted
    lambda = NMF.output$Lambda
    
    # calculate the DI scores and the ratio
    rho = 1 - rowSums(f) / (rowSums(fitted) + 1) # fitted sum coverages are added by 1 to avoid 0 in the denominator
    ratio = 1 - rho
    # print(paste('iter = 0 -- rho = ', rho))

    # decide whether to search for baseline regions of a gene.
    # if min(rho, na.rm = TRUE) > 0.2: skip baseline selection because there aren't samples with consistent degradation patterns;
    # if gene.length < 200: for genes with short transcript, baseline selection may not be necessary;
    # if median(ratio, na.rm = TRUE) > 1: this is imposed to exclude extreme cases where NMF-OA result doesn't converge.
    if( (gene.length < 200) || (median(ratio, na.rm = TRUE) > 1) || (min(rho, na.rm = TRUE) > 0.2) ){
      # use the results from NMF-OA instead
      K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
      fitted = NMF.output$fitted
      fitted[f.check>fitted] = f.check[f.check>fitted]
      # exclude extreme cases where the NMF result doesn't converge, return default values
      if(median(ratio, na.rm = TRUE) > 1){
        return(output)
      }
    }else{
      # Enter baseline selection loop.
      
      # divide the remaining transcript into 20 bins (with almost equal size)
      bin.size = ceiling(gene.length/20)
      bin.num = gene.length %/% bin.size + 1

      # Continue while degradation is high.
      itr <- 1
      while(max(rho, na.rm = TRUE) > 0.1){
        output$ran_baseline_selection <- TRUE
        
        # calculate relative residual from NMF output
        res.std = NMF.output$res / (f + 1) # add 1 to avoid 0 in the denominator
        dat.res = matrix(c(apply(res.std, 2, function(x) max(x^2)),
                           rep(NA, bin.size * bin.num - dim(f)[2])), nrow = bin.size)
        
        # bin average
        bin.mean = apply(dat.res, 2, mean, na.rm=TRUE)
        # print(paste('iter = ', itr, '-- ss_r =', bin.mean))
        
        # drop the bin with maximum average weighted squared relative residual
        drop.bin = which.max(bin.mean)
        
        # update the number of bins
        bin.num = bin.num - 1
        
        # update the coverage matrix and the corresponding Lagrangian multiplier
        drop.bp = c(max(bin.size*(drop.bin-1)+1,0):min(bin.size*drop.bin,dim(f)[2]))
        f = f[,-drop.bp]
        lambda = lambda[,-drop.bp]
        
        # run NMF again with updated coverage matrix
        NMF.output = NMF(f, loop=loop)
        fitted = NMF.output$fitted
        fitted[f>fitted] = f[f>fitted]
        
        # stop if the fitted values are all zero for any sample (extreme cases)
        if(min(rowSums(NMF.output$fitted)) == 0){
          break
        }
        
        # update the DI scores
        rho = 1 - rowSums(f) / (rowSums(fitted) + 1) 
        ratio = 1 - rho
        # print(paste('iter = ', itr, '-- rho = ', rho))
        
        # Stop while there are only 20% of the bins remaining (4 out of 20)
        # Stop while the remaining transcript length is less than 200 bp (imposed by Bin)
        if((bin.num == 4) || (dim(f)[2] < 200)){
          # print('Exiting for number of dropped bins or length of remaining gene.')
          break
        }
        itr <- itr + 1
      }
      
      # After baseline selection loop, determine whether the baseline is identified based on the NMF output from the last iteration
      # if max(rho, na.rm = TRUE) < 0.2: the maximum DI score is less than 0.2. 
      #                                  Degradation effect is almost removed in the remaining transcript and we conclude the baseline is found.
      if(max(rho, na.rm = TRUE) < 0.2){
      
        # label the convergence status as TRUE since the baseline is found
        output$convergence = TRUE
        # use the abundance level calculated from the baseline regions to estimate the fitted curves
        K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
        G = apply(f.check, 2, function(x,K) max(x/K), K=K)
        fitted = K %*% G
        # update the DI scores and the ratio based on the whole filtered transcript
        rho = 1 - rowSums(f.check) / (rowSums(fitted) + 1) 
        ratio = 1 - rho
        # print(paste('Post baseline, max(rho) < 0.2. -- rho = ', rho))

        # eliminate extreme cases where there's high degradation, but it's unlikely that
        # the degradation comes from "degradation," degradation is just from noise. Prevents overinflating due to noise.
        # Happens typically for larger genes with low sequencing depth.
        if(max(rho, na.rm = TRUE) > 0.9){
          NMF.output = NMF(f.check, loop=loop)      
          K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
          fitted = NMF.output$fitted
          fitted[f.check>fitted] = f.check[f.check>fitted]
          rho = 1 - rowSums(f.check) / (rowSums(fitted) + 1) 
          # print(paste('Post baseline, max(rho) < 0.2, then max(rho) > 0.9 -- rho = ', rho))
        }
        
      }else{ # when the baseline is not identified, use the NMF-OA results instead
        NMF.output = NMF(f.check, loop=loop)      
        K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
        fitted = NMF.output$fitted
        fitted[f.check>fitted] = f.check[f.check>fitted]
        rho = 1 - rowSums(f.check) / (rowSums(fitted) + 1)
        # print(paste('Post baseline, max(rho) > 0.2. -- rho = ', rho))
      }
    }
    
    output$K = K
    output$rho = rho

    return(output)
    
  }
}

################################################################################################################
# function to calulate SVD for the raw coverage
# used to identify genes with similar coverage pattern across all samples (minimal degradation)
# input: coverage matrix (f) with columns representing samples
# output: DI scores (rho)
################################################################################################################

ratioSVD <- function(f){
  f = t(as.matrix(f))
  model.svd = svd(f, nv = 1, nu = 1)
  fitted = model.svd$d[1] * model.svd$u %*% t(model.svd$v)
  fitted[fitted < f] = f[fitted < f]
  return(1-rowSums(f)/(rowSums(fitted)+1))
}

