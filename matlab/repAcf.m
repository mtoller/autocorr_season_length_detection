function [x] = repAcf(x, reps=10)
  pkg load tsa;
  x = zscore(detrend(x));
  if reps < 1
    return
  endif
  n = length(x);
  for i=1:reps
    x = x - mean(x);
    x = [1;acorf(x',n-1)'];
  endfor
  return
endfunction