function [a] = callSeasonLength(filename='temp')
  ts = load(filename);
  m = seasonLength(ts);
  save('-ascii','temp','m');
  a = m;
endfunction
