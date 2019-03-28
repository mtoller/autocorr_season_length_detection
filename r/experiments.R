repAcf <- function(data,reps = 10)
{
  source('fastAcf.R')
  a <- data
  
  for (i in 1:reps)
  {
    a <- acf.fft(a)
  }
  return(a)
}

euclid <- function(a,b=NULL)
{
  if (!is.null(b)) return(sqrt(sum((a-b)^2)))
  return(apply(a,2,function(i){
    return(apply(a,2,function(j){
      return(euclid(i,j))
    }))
  }))
}
cos.sim <- function(a,b) 
{
  return(t(a)%*%b/sqrt(sum(a^2)*sum(b^2)))
}

repAcfSeq <- function(data, func = mean,limit = 100, epsilon = 1e-01,...)
{
  a <- data
  old_a <- a
  result <- c()
  for (i in 1:limit)
  {
    a <- acf.fft(a)
    result <- c(result,func(a,...))
    #result <- c(result,func(a,preprocess=F))
    if (CID(a,old_a) < epsilon) break;
    old_a <- a
    
  }
  return(result)
}
deepSazed <- function(data)
{
  require(pracma)
  source('fastAcf.R')
  if (anyNA(data))
  {
    print('NA is not supported')
    return(1)
  }
  if (any(is.infinite(data)) || any(is.complex(data)))
  {
    print('complex numbers are not supported')
    return(1)
  }
  if (var(data) == 0)
  {
    print('y has no variance')
    return(1)
  }
  #data <- scale(as.ts(detrend(as.vector(data))))
  data <- as.ts(detrend(as.vector(data)))
  #ms <- sazed(data)
  ms <- csazed(data)
  #ms <- unique(c(ms,repAcfSeq(data,func = sazed,method='down',preprocess=T)))
  ms <- unique(c(ms,repAcfSeq(data,func = csazed)))
  result <- compressedSimilarity(data,ms)
  print(list(ms,result))
  cat(paste0('Confidence: ', max(result),'\n'))
  return(ms[which.max(result)])
}
compressedFindFrequency <- function(data)
{
  source('fastAcf.R')
  if (anyNA(data))
  {
    print('NA is not supported')
    return(1)
  }
  if (any(is.infinite(data)) || any(is.complex(data)))
  {
    print('complex numbers are not supported')
    return(1)
  }
  if (var(data) == 0)
  {
    print('y has no variance')
    return(1)
  }
  data <- scale(as.ts(detrend(as.vector(data))))
  p <- spec.pgram(c(data),plot = F)
  ms <- unique(round(1/p$freq[which(p$spec > median(p$spec)+2*sd(p$spec))]))
  result <- compressedSimilarity(data,ms)
  print(list(ms,result))
  result <- ms[which.max(result)]
  return(result)
}
compressedSimilarity <- function(data,ms,method=cor,...)
{
  n <- length(data)
  result <- sapply(ms,function(m){
    k <- floor(n/m)
    if (m <=2 | k <= 3)
    {
      return(0)
    }
    subs <- sapply(1:k,function(i){data[(((i-1)*m)+1):(i*m)]})
    #print(subs)
    similarity <- cor(subs)#applyDistanceToPairs(matrix(unlist(subs),nrow = m,ncol = k,byrow = F),method)
    if (similarity[1,1] == 0) #convert distance to similarity
    {
      similarity <- similarity/max(similarity)
      similarity <- 1 - similarity
    }
    return(mean(similarity,na.rm = T))
  })
  return(result)
}

mdl <- function(data)
{
  u <- scale(data)
  n <- length(u)
  M <- 1
  mdls <- sapply(2:(trunc(sqrt(n))),function(P){
    N <- floor(n/P)
    if (P > N) break
    x <- sapply(1:N,function(i){u[(((i-1)*P)+1):(i*P)]})
    z <- unlist(as.list(apply(x,1,fftwtools::fftw)))
    s <- z %*% Conj(t(z))
    sks <- lapply(1:(N),function(i)s[((i-1)*P+1):(i*P),((i-1)*P+1):(i*P)])
    return(N*P*M*(log(pi)+1) + M*sum(log(unlist(lapply(sks,negDet)))) + 0.5*N*P*P*log(M))
  })
  print(mdls)
  return(which.max(mdls)+1)
}
variability <- function(x,limit=50,j=3)
{
  N <- length(x)
  if (N < limit) return(1)
  tau <- 0
  v_hat <- rep(0,limit-1)
  for (n in 2:limit)
  {
    p <- floor(N/n)
    blocked <- sapply(0:(p-1),function(i)x[(1+(i*n)):(n+(i*n))])
    m_hat <- sapply(blocked,mean)
    v_hat[n-1] <- var(m_hat)
    #r_hat <- apply(blocked,1,var) %>% as.numeric
    #print(r_hat)
    #print(var(as.numeric(r_hat)))
    #v_hat[n-1] <- var(r_hat)
  }
  #plot.ts(v_hat)
  v_hat %>% spec.pgram() %>% {1/.$freq[min(order(.$spec,decreasing = T)[1:j])]} %>% return
}
altSazed <- function(data)
{
  require(signal)
  require(forecast)
  require(pracma)
  require(bspec)
  source('sazed.R')
  if (anyNA(y))
  {
    print('NA is not supported')
    return(1)
  }
  if (any(is.infinite(y)) || any(is.complex(y)))
  {
    print('complex numbers are not supported')
    return(1)
  }
  if (var(y) == 0)
  {
    print('y has no variance')
    return(1)
  }
  
  results <- c()
  results <- c(results,S(data))
  results <- c(results,Sa(data))
  results <- c(results,ze(data))
  results <- c(results,aze(data))
  results <- c(results,zed(data))
  results <- c(results,azed(data))
  
  results <- results[which(!is.infinite(results))];
  results <- results[which(!is.na(results))];
  #print(results)
  if (var(results) == 0)
  {
    return(results[1])
  }
  ms <- unique(results)
  result <- evaluateSLSuggestions(data,ms)
  print(list(ms,result))
  result <- ms[which.max(result)]
  return(result)
}

testAllSLSuggestions <- function(data)
{
  require(pracma)
  data <- scale(as.ts(detrend(as.vector(data))))
  n <- length(data)
  #if (n > 8000) return(1)
  ms <- 2:floor(n/3)
  #return(evaluateSLSuggestions(data,ms))
  result <- compressedSimilarity(data,ms)
  print(list(ms,result))
  cat(paste0('Confidence: ',max(result),'\n'))
  return(ms[which.max(result)])
}
evaluateSLSuggestions <- function(data,ms)
{
  #similarityMeasures <- list(cor)
  distanceMeasures <- list(euclid,CID)
  similarityMeasures <- append(similarityMeasures,lapply(distanceMeasures,function(i){
    return(function(a,b){return(1/(1+i(a,b)))})
  }))
  
  n <- length(data)
  result <- sapply(ms,function(m){
    k <- floor(n/m)
    if (m <=2 | k <= 3)
    {
      return(0)
    }
    subs <- sapply(1:k,function(i){data[(((i-1)*m)+1):(i*m)]})
    #print(subs)
    mat <- matrix(unlist(subs),nrow = m,ncol = k,byrow = F)
    results <- lapply(similarityMeasures,function(s)applyDistanceToPairs(mat,s))
    print(results)
    r <- results[[1]]
    for (i in 2:length(results)) r <- r %*%results[[i]]
    #print(r)
    return(mean(r))
    #return(mean(unlist(lapply(results,function(i)mean(i,na.rm = T)))))
  })
  print(list(ms,result))
  return(ms[which.max(result)])
}
twed <- function(A,B,lambda=0,nu=1,timesA=NULL,timesB=NULL)
{
  #Padding
  A <- c(0,A)
  B <- c(0,B)
  
  #initialize variables
  n <- length(A)
  m <- length(B)
  
  if (is.null(timesA)) timesA <- c(1:n)
  if (is.null(timesB)) timesB <- c(1:m)
  
  #initialize array for dynamic programming
  DP <- matrix(rep(0,n*m),nrow = n, ncol = m)
  DP[1,] <- Inf
  DP[,1] <- Inf
  DP[1,1] <- 0
  
  for (i in 2:n)
  {
    for (j in 2:m)
    {
      cost <- euclid(A[i],B[j])
      
      #Calculate and save cost of various operations
      C <- c()
      #deletion in A
      C <- append(C, DP[i-1,j] + euclid(A[i-1],A[i]) + nu * (timesA[i] - timesA[i-1]) + lambda)
      #deletion in B
      C <- append(C, DP[i,j-1] + euclid(B[j-1],B[j]) + nu * (timesB[j] - timesB[j-1]) + lambda)
      #Keep data points in both time series
      C <- append(C, DP[i-1,j-1] +  euclid(A[i],B[j]) + euclid(A[i-1],B[j-1]) +
                    nu * (abs(timesA[i] - timesB[j]) + abs(timesA[i-1] - timesB[j-1])))
      
      DP[i,j] <- min(C)
      
    }
  }
  #print(DP)
  return(DP[n,m])
}
dtwWrapper <- function(a,b)
{
  require(dtw)
  return(dtw(a,b)$distance)
}
pearsonWrapper <- function(a,b) return(cor(a,b))
spearmanWrapper <- function(a,b) return(cor(a,b,method = 'spearman'))
kendallWrapper <- function(a,b) return(cor(a,b,method = 'kendall'))
applyDistanceToPairs <- function(X,dis)
{
  return(apply(X,2,function(i){
    return(apply(X,2,function(j){
      return(dis(i,j))
    }))
  }))
}
CID <- function(A,B)
{
  CE_A <- sqrt(sum(diff(A)^2))
  CE_B <- sqrt(sum(diff(B)^2))
  return(sqrt(sum((A-B)^2)) * max(c(CE_A,CE_B)) / min(c(CE_A,CE_B)))
}
deepCorrelation <- function(data)
{
  spec_dens <- spec.ar(data,plot = F)
  candidate <- round(1/spec_dens$freq[which.max(spec_dens$spec)])
  scaling_factor <- min(round(candidate/3),round(sqrt(length(data))))
  #if (length(data) > decision_boundary)
  #{
  #  scaling_factor <- round(sqrt(length(data)))
  data <- downsample(data,window_size = scaling_factor)
  #}
  ms <- 1:(scaling_factor*3*2)
  result <- compressedSimilarity(data,ms)
  ms <- ms * scaling_factor
  print(list(ms,result))
  #par(mfrow=c(2,1))
  #plot.ts(data)
  #plot(result)
  return(ms[which.max(result)])
}
rollingWindow <- function(data,m,func)
{
  n <- length(data)
  subs <- sapply(1:(n-m+1),function(x){data[x:(x+m-1)]})
  return(apply(subs,2,func))
}

kerasAttempt <- function(train_x, train_y, test_x, test_y)
{
  require(keras) #expects constant input size. useless.
  
  train_x <- lapply(train_x, scale)
  
  model <- keras_model_sequential()
  model %>% layer_dense(units = 64, activation = "relu",
                        input_shape = dim(train_data)[2]) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1)
  
  compile(model, loss = 'mse')
  
}

exp1 <- function(y, nu=100)
{
  require(dplyr)
  #downsampling
  #+ makes data readable
  #- short frequencies get lost
  n <- length(y)
  eta <- 1
  z <- y
  if (n > 1000)
  {
    eta <- round(n/nu)
    print(eta)
    z <- y %>% downsample(eta) %>% scale
    result <- testAllSLSuggestions(z)
  }
  return(testAllSLSuggestions(z)*eta)
}
downsample = function(data, window_size=2)
{
  require(zoo);
  return(ts(as.ts(rollapply(zoo(data),width=window_size,by=window_size,FUN=mean)),frequency=1))
}
csazed <- function(x){
  require(sazedR)
  require(signal)
  require(forecast)
  require(pracma)
  require(bspec)
  require(dplyr)
  source('sazed.R')
  if (anyNA(x))
  {
    print('NA is not supported')
    return(1)
  }
  if (any(is.infinite(x)) || any(is.complex(x)))
  {
    print('complex numbers are not supported')
    return(1)
  }
  if (var(x) == 0)
  {
    print('y has no variance')
    return(1)
  }
  
  p <- 0
  r <- 1
  f <- 0
  repeat{
    
    results <- c()
    results <- c(results,S(x))
    results <- c(results,Sa(x))
    results <- c(results,ze(x))
    results <- c(results,aze(x))
    results <- c(results,zed(x))
    results <- c(results,azed(x))
    
    results <- results[which(!is.infinite(results))]
    results <- results[which(!is.na(results))]
    results <- results[which(results >= 2)]
    
    if (length(results) < 1) return(1)
    if (var(results) == 0) return(results[1])
    
    #certainties <- sapply(results,function(d){
    #  
    #  k <- floor(length(x)/d)
    #  if (d <=2 | k <= 3)
    #  {
    #    return(0)
    #  }
    #  subs <- sapply(1:k,function(i){x[(((i-1)*d)+1):(i*d)]})
    #  subs %>% cor() %>% min %>% return()
    #})
    #if (any(certainties == 1)) return(max(certainties[which(certainties == 1)]))
    
    #print(results)
    #print(certainties)
    tab <- table(results)
    majorities <- which(tab==max(tab)) %>% {names(tab)[.]} %>% sapply(as.numeric) %>% as.numeric
    #print(table(results))
    #print(majorities)
    l <- length(majorities)
    #print(l)
    if (l == 1){
      return(majorities[1])
    }
    #print(majorities)
    certainties <- sapply(majorities,function(d){
      
      k <- floor(length(x)/d)
      if (d <=2 | k <= 3)
      {
        return(0)
      }
      subs <- sapply(1:k,function(i){x[(((i-1)*d)+1):(i*d)]})
      subs %>% cor() %>% min %>% return()
    })
    certainties[which(is.na(certainties))] <- 0
    p <- max(certainties)
    r <- certainties %>% {min(which(. == max(.)))} %>% {majorities[.]}
    #result %>% print
    #result %>% max %>% print
    #result %>% {max(which(. == max(.)))} %>% {majorities[.]} %>% return
    cat(paste0("p is ",p,"\n"))
    if (p > 0) break
    #x <- acf.fft(x)
    x <- downsample(x,window_size = 2)
    f <- f+1
    #if (f > 50) break
  }
  #return(r)
  return(r*2^f)
}