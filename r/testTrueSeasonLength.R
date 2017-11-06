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

testAll = function()
{
  sum = 0;
  expected = c(88,57,3,6,5,25,214,900,10,1000,20,15,190,107,52,1000,20,35,50000,3);
  sum = sum + testSuite(expected,"../../test1/t",1);
  
  expected = c(10000,10,96,120,19,180,1332,45,19,20,140,18,30,1500,740,20,1590,10,29,5);
  sum = sum + testSuite(expected,"../../test2/t",2);
  
  expected = matrix(
    c(20,60,-1,
      3,10,-1,
      50,90,-1,
      60,300,-1,
      70,700,-1,
      40,3498,-1,
      8,32,-1,
      40,75,-1,
      400,1000,-1,
      9,20,-1,
      5,980,-1,
      10,64,-1,
      78,100,-1,
      6,12,18,
      3,7,-1,
      9,24,-1,
      7,70,700,
      400,800,-1,
      6,20,-1,
      500,1400,2000)
    ,nrow = 3,ncol = 20);
  expected = t(expected);
  sum = sum + testMultiple(expected,"../../test3/t",3);
  
  expected = c(128,128,128,128,128,156,156,156,156,156,37,37,37,37,37,8,8,8,8,8);
  sum = sum + testSuite(expected,"../../test4/t",4);

  expected = c(155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155,155);
  sum = sum + testSuite(expected,"../../test5/t",5);
  
  expected = c(4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32786,65536);
  sum = sum + testSuite(expected,"../../test6/t",6);

  expected = c(0,0,0,0,0,0,0,0,0,0)
  sum = sum + testSuite(expected,"../../test7/t",7);
  
  expected = c(12,12,12,4,12,6,6,6,12,6,4,12,12,12,12,12,3,3,12,12)
  sum = sum + testSuite(expected,"../../testR/t",8);
  
  expected = c(12,12,12,12,12,12,12,130,12,12,12,12,12,12,130,12,12,12,12,12)
  sum = sum + testSuite(expected,"../../testC/t",9);
  
  print(paste('Total: ', sum, ' out of 165'));
  print('original seasonLength passes 122');
  print('findFrequency passes 83');
}

testMultiple = function(expected,name,number)
{
  n = nrow(expected);
  passed = 0;
  results = c(1:n);
  bin = c(1:n);
  which_one = c(1:n);
  print(paste("Test suite ",number,":",sep=''));
  for (i in 1:n)
  {
    result = testTrueSeasonLength(paste(name,i,sep='',collapse = ' '));
    results[i] = result;
    for (j in 1:ncol(expected))
    {
      if ((result >= (expected[i,j]*0.8))  && (result <= (expected[i,j]*1.2)))
      {
        passed = passed + 1;
        bin[i] = 1;
        which_one[i] = j;
        break;
      }
      else if (j == ncol(expected))
      {
        bin[i] = 0;
      }
    }
  }
  matrix(append(c(1:n),which_one),nrow=n,ncol=2)
  print(paste(passed," out of ",n,sep='',collapse = ' '));
  print(t(matrix(c(bin,results,expected[which_one]),nrow=n,ncol=3)));
  #print(paste(c("Passed:   ",bin),sep='',collapse = ' '));
  #print(paste(c("Results:  ",results),sep='',collapse = ' '));
  #print(paste(c("Expected: ", expected),sep = '',collapse = ' '));
  print("");
  return(passed);
}

testSuite = function(expected,name,number)
{
  n = length(expected);
  passed = 0;
  results = c(1:n);
  bin = c(1:n);
  print(paste("Test suite ",number,":",sep=''));
  for (i in 1:n)
  {
    result = testTrueSeasonLength(paste(name,i,sep='',collapse = ' '));
    results[i] = result;
    if ((result >= (expected[i]*0.8))  && (result <= (expected[i]*1.2)))
    {
      passed = passed + 1;
      bin[i] = 1;
    }
    else
    {
      bin[i] = 0;
    }
  }
  
  print(paste(passed," out of ",n,sep='',collapse = ' '));
  print(t(matrix(c(bin,results,expected),nrow=n,ncol=3)));
  #print(paste(c("Passed:   ",bin),sep='',collapse = ' '));
  #print(paste(c("Results:  ",results),sep='',collapse = ' '));
  #print(paste(c("Expected: ", expected),sep = '',collapse = ' '));
  print("");
  return(passed);
}