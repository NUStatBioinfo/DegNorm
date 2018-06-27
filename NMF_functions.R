####################################################################################
# NMF functions
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

optiNMF <- function(f, norm.factor, print.plot = FALSE){
  
  f = t(as.matrix(f))
  f = f/norm.factor
  
  filter = (apply(f,2,max) > 0.1 * max(apply(f,2,max))) == 1
  
  output = list(ratio = NULL, convergence = FALSE, K = NULL)
  output$ratio = rep(1, dim(f)[1])
  
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
    bin.size = ceiling(dim(f.check)[2]/20)
    bin.num = dim(f.check)[2] %/% bin.size + 1
    bin.label = rep(c(1:bin.num), each = bin.size)

    if((min(ratio, na.rm = TRUE) > 0.9) || (max(ratio, na.rm = TRUE) < 0.8) || gene.length < 200){
      K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
      fitted = NMF.output$fitted
      fitted[f.check>fitted] = f.check[f.check>fitted]
    }else{
      
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

        if((bin.num == 4) || (dim(f)[2] < 200)){
          break
        }
      }
      
      if(min(ratio) > 0.8){
        
        output$convergence = TRUE
        K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
        G = apply(f.check, 2, function(x,K) max(x/K), K=K)
        fitted = K %*% G
        
        if(min(rowSums(f.check) / (rowSums(fitted)+1)) < 0.1){
          NMF.output = NMF(f.check)      
          K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
          fitted = NMF.output$fitted
          fitted[f.check>fitted] = f.check[f.check>fitted]
        }
        
      }else{
        NMF.output = NMF(f.check)      
        K = abs(NMF.output$svd$u) * NMF.output$svd$d[1]
        fitted = NMF.output$fitted
        fitted[f.check>fitted] = f.check[f.check>fitted]
      }
    }
    
    output$K = K
    output$ratio = rowSums(f.check) / (rowSums(fitted)+1)
    return(output)
    
  }else{
    return(output)
  }
  
}

ratioSVD <- function(f, norm.factor){
  f = t(as.matrix(f))
  nobel = svd(f, nv = 1, nu = 1)
  fitted = nobel$d[1] * nobel$u %*% t(nobel$v)
  fitted[fitted < f] = f[fitted < f]
  return(rowSums(f)/(rowSums(fitted)+1))
}

