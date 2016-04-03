#Fibonacci_Code
This function will compute the Fibonacci Code for any given namber by applying Fibonacci code conditions. First, we need to define the function name and arguments required.

`Fibonacci<-function(x=x) {}`

######The following code create the Fibonacci rules or conditions that we need to apply on the numerical variable x:

`for(i in seq(along=(x))) {`


`if (x[i] <= 2) {`

` x[i]=1`

`}`

` else {`

`x[i]=((x[i]-1)+(x[i]-2))`

`}`

`}`

######Then, we can ask the function to return the computed values of the variable x using the follwoing line:

`return(x)`

`}`



    
  
