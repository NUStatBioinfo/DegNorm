####################################################################################
# NMF functions
# Non-negative matrix factorization with over-approximation.
# input: coverage matrix (f) with columns representing samples
# output: truncated SVD (svd), fitted coverage matrix (fitted), residuals (res), Lagrangian multipler (Lambda).
NMF <- function(f, loop = 100){
  
  Lambda = matrix(0, ncol = dim(f)[2], nrow = dim(f)[1])
  model_svd = svd(f, nv = 1, nu = 1)
  fitted = model_svd$d[1] * model_svd$u %*% t(model_svd$v)
  res = fitted - f
  
  for(i in 1:loop){
    
    Lambda = Lambda - 1/sqrt(loop) * res
    Lambda[which(Lambda < 0)] = 0
    
    model_svd = svd(f+Lambda, nv = 1, nu = 1)
    fitted = model_svd$d[1] * model_svd$u %*% t(model_svd$v)
    res = fitted - f

  }
  
  return(list(svd = model_svd, fitted = fitted, res = res, Lambda = Lambda))
}


# baseline selection process to improve NMF
# input: coverage matrix (f), normalization factor for sequencing depth (norm.factor)
# output: DI scores (rho), baseline status (convergence) and abundance level (K).
optiNMF <- function(f, norm.factor){
  
  f = t(as.matrix(f))
  f = f/norm.factor
  
  filter = (apply(f,2,max) > 0.1 * max(apply(f,2,max))) == 1
  
  output = list(rho = NULL, convergence = FALSE, K = NULL)
  output$rho = rep(0, dim(f)[1])
  
  num.sample = dim(f)[1]
  
  if(sum(filter) > 50){
    
    f = f[, filter]
    
    if(sum(rowSums(f) > 0) < num.sample){
      return(output)
    }

    f.check = f
    gene.length = dim(f)[2]
    
    NMF.output = NMF(f)
    fitted = NMF.output$fitted
    lambda = NMF.output$Lambda
    
    ratio = rowSums(f) / (rowSums(fitted) + 1)
    rho = 1 - ratio
    bin.size = ceiling(dim(f.check)[2]/20)
    bin.num = dim(f.check)[2] %/% bin.size + 1
    bin.label = rep(c(1:bin.num), each = bin.size)

    # decide whether to search for baseline regions of a gene.
    # if max(ratio) < 0.8, skip baseline selection because there aren't samples with consistent degradation patterns.
    # if gene.length < 200, and avg fragment length is 200, then reads aren't going to be reliable.
    if((min(ratio, na.rm = TRUE) > 0.9) || (max(ratio, na.rm = TRUE) < 0.8) || gene.length < 200){
      K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
      fitted = NMF.output$fitted
      fitted[f.check>fitted] = f.check[f.check>fitted]
    }else{
      
      # Enter baseline selection loop.
      # Continue while degradation is high.
      while(min(ratio) < 0.9){
        
        res.std = NMF.output$res / (f + 1)
        dat.res = matrix(c(apply(res.std, 2, function(x) max(x^2)),
                           rep(0, bin.size * bin.num - dim(f)[2])), nrow = bin.size)
        
        bin.mean = apply(dat.res, 2, mean)
        
        drop = which.max(bin.mean)
        bin.num = bin.num - 1
        
        # print(drop)
        drop = c(max(bin.size*(drop-1),0):min(bin.size*drop,dim(f)[2]))
        f = f[,-drop]
        lambda = lambda[,-drop]
        
        NMF.output = NMF(f)

        fitted = NMF.output$fitted
        fitted[f>fitted] = f[f>fitted]
        
        if(min(rowSums(NMF.output$fitted)) == 0){
          break
        }
        ratio = rowSums(f) / (rowSums(NMF.output$fitted) + 1)

        # "This step stops if the maximum DI score obtained from theremaining bins is <=0.1or 70% bins" - 
        # The first criteria is cited from Online Methods; second criteria was imposed by Bin.
        if((bin.num == 6) || (dim(f)[2] < 200)){
          break
        }
      }
      
      # if we stopped during baseline selection and all samples are "kind of" degraded, just consider
      # it that you HAVE found baseline. 
      if(min(ratio) > 0.8){
        
        output$convergence = TRUE
        K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
        G = apply(f.check, 2, function(x,K) max(x/K), K=K)
        fitted = K %*% G
        
        # eliminate extreme cases where there's high degradation, but it's unlikely that
        # the degradation comes from "degradation," degradation is just from noise. Prevents overinflating due to noise.
        # Happens typically for larger genes with low sequencing depth.
        if(min(rowSums(f.check) / (rowSums(fitted)+1)) < 0.1){
          NMF.output = NMF(f.check)      
          K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
          fitted = NMF.output$fitted
          fitted[f.check>fitted] = f.check[f.check>fitted]
        }
        
      # If the gene didn't go through baseline selection, just use NMF-OA estimate.
      }else{
        NMF.output = NMF(f.check)      
        K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
        fitted = NMF.output$fitted
        fitted[f.check>fitted] = f.check[f.check>fitted]
      }
    }
    
    output$K = K
    output$rho = 1 - rowSums(f.check) / (rowSums(fitted)+1)
    return(output)
    
  }else{
    return(output)
  }
  
}

# for intializing DegNorm iterations, use for obtaining sequencing depth factors.
ratioSVD <- function(f, norm.factor){
  f = t(as.matrix(f))
  model_svd = svd(f, nv = 1, nu = 1)
  fitted = model_svd$d[1] * model_svd$u %*% t(model_svd$v)
  fitted[fitted < f] = f[fitted < f]
  return(1-rowSums(f)/(rowSums(fitted)+1))
}

