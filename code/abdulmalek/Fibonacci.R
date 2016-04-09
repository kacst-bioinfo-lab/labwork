Fibonacci<-function(x){
  for(i in seq(along=(x))){
    
   if (x[i] == 0){
     x[i]== 0 
   }
    if (x[i]== 1){
      x[i]== 1
    }
    if (x[i] >= 2){
     x[i]=((x[i-1])+(x[i-2]))
      }
  }
  return(x)
  }
