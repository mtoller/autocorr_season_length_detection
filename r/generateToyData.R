generateToyData <- function(N=10){
  set.seed(0815)
  n <- 3000
  toy_data <- list()
  toy_expected <- list()
  for(i in 1:N){
    r <- arima.sim(list(ar=c(sample(seq(0,1,0.001))[1])),n=3000,rand.gen = rcauchy)
    period <- sample(3:1000)[1]
    amplitude <- sample(3:1000)[1]
    noise_amplitude <- sample(3:amplitude)[1]
    offset <- sample(0:100)[1]
    if (sample(c(T,F)[1])){
      method <- sin
    }
    else{
      method <- cos
    }
    
      
    trendslope <- rnorm(1)
    series <- amplitude*method(((1:n)+offset)*2*pi/period) + (1:n)*trendslope + r * noise_amplitude
    #k <- sample(500:2500)[1]
    #series[k:length(series)] <- series[k:length(series)] + rnorm(1,sd=sd(series)*2)
    toy_data <- append(toy_data,list(c(series)))
    toy_expected <- append(toy_expected,period)
  }
  #for (test in data){
  #  plot.ts(test)
  #  Sys.sleep(1.5)
  #}
  
}