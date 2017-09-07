fastMa = function(y,order){
  y = c(y,rep.int(0,order));
  n = length(y);
  csum = cumsum(y);
  ma = (csum[(order+1):n] - csum[1:(n-order)])/order;
  return(ma);
}