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

euclid <- function(a,b)
{
  return(sqrt(sum((a-b)^2)))
}
cos.sim <- function(a,b) 
{
  return(t(a)%*%b/sqrt(sum(a^2)*sum(b^2)))
}   

repAcfSeq <- function(data, func = mean,limit = 10, epsilon = 1e-07)
{
  a <- data
  old_a <- a
  result <- c()
  for (i in 1:limit)
  {
    a <- acf.fft(a)
    result <- c(result,func(a))
    if (euclid(a,old_a) < epsilon) break;
    old_a <- a
    
  }
  return(result)
}
deepSazed <- function(data)
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
  ms <- sazedR::sazed(data)
  ms <- unique(c(ms,repAcfSeq(scale(data),func = sazedR::sazed)))
  result <- testSLSuggestion(data,ms)
  print(list(ms,result))
  result <- ms[which.max(result)]
  return(result)
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
  result <- testSLSuggestion(data,ms)
  print(list(ms,result))
  result <- ms[which.max(result)]
  return(result)
}
testSLSuggestion <- function(data,ms)
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
    return(mean(cor(matrix(unlist(subs),nrow = m,ncol = k,byrow = F)),na.rm = T))
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
variability <- function(data)
{
  N <- length(data)
  for (n in 2:(floor(N)/50))
  {
    i <- floor(N/n)
    subs <- t(sapply(0:(n-1),function(j)data[(1:i)*n+j]))
    print(subs)
  }
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
  result <- testSLSuggestion(data,ms)
  print(list(ms,result))
  result <- ms[which.max(result)]
  return(result)
}

testAllSLSuggestions <- function(data)
{
  require(pracma)
  data <- scale(as.ts(detrend(as.vector(data))))
  n <- length(data)
  if (n > 8000) return(1)
  ms <- 2:floor(n/3)
  result <- testSLSuggestion(data,ms)
  print(list(ms,result))
  return(ms[which.max(result)])
}