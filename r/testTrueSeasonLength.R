testTrueSeasonLength = function(datafile){
  #print(datafile);
  source('trueSeasonLength.R');
  y = read.table(datafile);
  if (length(y) == 2)
  {
    y = y[2];
  }
  else
  {
    y = y[1];
  }  
  y = as.ts(y);
  return(trueSeasonLength(y));
}