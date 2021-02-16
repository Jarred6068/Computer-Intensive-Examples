
data=NULL
set.seed(23)
y=runif(500,0,2)
u=runif(500,0,1)
for ( i in 1:500){
  if(y[i]<=2*u[i]){
    
  data[i] = u[i]
    
  }
  
}