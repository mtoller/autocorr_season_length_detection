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
  ms <- sazedR::sazed(data)
  ms <- c(ms,repAcfSeq(scale(data),func = sazedR::sazed))
  if (var(ms) == 0)
  {
    return(ms[1])
  }
  result <- testSLSuggestion(data,ms)
  #print(list(ms,result))
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
    return(mean(cor(matrix(unlist(subs),nrow = m,ncol = k,byrow = F))))
  })
  return(result)
}