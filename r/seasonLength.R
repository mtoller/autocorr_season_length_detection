linReg = function(x, y, degree=1,skip=0){
  X = matrix(1,nrow(y),1)
  
  for(i in 1:degree)
  {
    X = cbind(X,x^i);
  }
  
  if (FALSE)#(is.null(skip) && log10(rcond(t(X) %*% X)) < -10)
  {
    print('error');
    X = as.integer(0);
    theta = as.integer(0);
    return(NULL);
  }
  theta = solve(t(X) %*% X) %*% t(X) %*% y;
  return(list('X'=X,'theta'=theta));
}

detrend1 = function(y, degree,x=NULL,skip=NULL)
{
  if (is.null(x))
  {
    x = matrix(1:nrow(y),nrow(y),1);
  }
  
  
  trend = linReg(x,y,degree,skip);
  if (is.null(trend))
  {
    print('trend is null')
    return(y);
  }
  print(mean(as.vector(trend$X %*% trend$theta)),max = 10);
  return(y - trend$X %*% trend$theta);
}

expandData = function(x,y)
{
  if (nrow(x) != nrow(y))
  {
    print('Dimension mismatch');
    return;
  }
  ex = matrix(seq(x[1,1],x[nrow(x),1],by=0.01));
  ey = interp1(as.vector(x),as.vector(y),ex);
  return(list('ex' = ex, 'ey' = ey));
}

findThreshold = function(delta, limit, force_return=0)
{
  quotients = delta[2:nrow(delta)] / delta[1:(nrow(delta)-1)];
  
  gamma = which(abs(quotients[2:length(quotients)] - quotients[1:(length(quotients)-1)]) > limit);
  
  if (length(which(gamma == 1)) == 0 && length(quotients) > 0)
  {
    gamma = append(1,gamma);
  }
  gamma = append(gamma, length(quotients));
  print('gamma');
  print(gamma);
  
  if (length(gamma) <= 2 && length(delta) == 2)
  {
    return(c(delta[1],delta[2]));
  }
  gamma_diff = gamma[2:length(gamma)] - gamma[1:(length(gamma)-1)];
  index = which(gamma_diff == max(gamma_diff));
  return(c(delta[gamma[index]+1], delta[gamma[index+1]+1]));
  
}

seasonLength = function(y, degree=1, butter1=2, butter2=0.001, force_return=0){
  library(pracma);
  library(signal);
  n = nrow(y);
  
  x = matrix(1:n,n,1);
  
  ret1 = linReg(x,y,1);
  X1 = ret1$X;
  theta1 = ret1$theta;

  ret2 = linReg(x,y,1);
  X2 = ret2$X;
  theta2 = ret2$theta;
  
  mse1 = mean((y-X1%*%theta1)^2);
  mse2 = mean((y-X2%*%theta2)^2);
  
  if (log(abs(mse1-mse2)) > exp(2))
  {
    degree=2;
  }
  if (n > 2 && n < 10000)
  {
    ret = expandData(x,y);
    x = as.matrix(ret$ex);
    y = as.matrix(ret$ey);
  }
  ne = nrow(y);
  
  coeffs = butter(butter1, butter2);
  b = coeffs$b;
  a = coeffs$a;
  
  ys = as.matrix(filter(b,a,y));
  
  no_filter = 0;
  stop = 0;
  while(stop != 1)
  {
    print('iter')
    stop = 1;
    print('ys');
    print(mean(detrend1(ys,degree,skip = 1)));
    #print(detrend1(ys,degree,x));
    
    yt = acf.fft(detrend1(ys,degree,skip = 1));
    
    yt = yt[2:ne];
    
    yt = as.matrix(append(yt, yt[ne-1]));
    
    params = linReg(x,yt,1);
    X = params$X;
    theta = params$theta;
    
    alpha = which(abs(yt - X %*% theta) < 0.001);
    alpha = as.matrix(alpha);
    print('alpha');
    print(nrow(alpha));

    if (nrow(alpha) > 1)
    {
      distances = alpha[2:nrow(alpha)]-alpha[1:(nrow(alpha)-1)];
      delta = distances[order(-distances)];
      delta = delta[which(delta > 1)];
      delta = sort(delta);
    }
    else
    {
      delta = 0;
    }
    if (nrow(as.matrix(delta)) < 6 && no_filter == 0)
    {
      stop = 0;
      no_filter = 1;
      ys = y;
    }
  }
  if (nrow(as.matrix(delta)) < 2)
  {
    print('delta too small');
    return(0);
  }
  print('delta');
  print(delta);
  threshold = findThreshold(as.matrix(delta), 0.05, force_return);
  delta = delta[which(delta >= threshold[1])];
  delta = delta[which(delta <= threshold[2])];
  print('return value');
  return_value = 2 * mean(quantile(delta,c(0.1,0.25,0.5,0.75,0.9)));
  
  return_value = return_value * abs(x[2]-x[1]);
  return_value = round(return_value);
  print(return_value);
  
  if (return_value >= 0.5 * n || return_value <= 1)
  {
      return(0);
  }
  return(return_value);
}