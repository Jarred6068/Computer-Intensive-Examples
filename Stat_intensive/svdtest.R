
n=50
a1=runif(n,0,10)
a2=runif(n,0,10)
a3=runif(n,0,10)
a4=runif(n,0,10)

A=cbind.data.frame(a1,a2,a3,a4)
A=as.matrix(A)

b1=round(runif(1, -5, 0.5))
b2=round(runif(1, -3, 8))
b3=round(runif(1, 1, 7))
b4=round(runif(1, 1, 15))

xi=c(b1,b2,b3,b4)
errors=rnorm(n, 0, 1)

y=A%*%xi+errors

data=cbind.data.frame(y, A)

rm(a1,a2,a3,a4, b1,b2,b3,b4)

obj=reg.coef(data)
svd.obj=svd(data[,-1])

print(xi)
print(obj)


#testing svd concepts
idx=dim(A)
n=idx[2]
m=idx[1]

U.svd=svd(A)$u
Sig.svd=diag(svd(A)$d)
V.svd=svd(A)$v

ATA=t(A)%*%A
AAT=A%*%t(A)


U=eigen(AAT)$vectors[, 1:min(n,m)]
