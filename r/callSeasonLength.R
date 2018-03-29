callSeasonLength = function(y)
{
  write(y,file = 'temp',ncolumns = 1);
  system('cd ../matlab;octave --silent --no-gui --eval "callSeasonLength()"');
  system('cd ../r');
  return(read.ts('temp')[1]);
}