function [a] = callSeasonLength(filename='../r/temp')
  warning('off', 'all');
  ts = load(filename);
  m = seasonLength(ts);
  save('-ascii','temp','m');
  a = m;
endfunction
