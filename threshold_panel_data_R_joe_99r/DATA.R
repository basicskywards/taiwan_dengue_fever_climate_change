##  DATA.R  ##    

invest <- read.table("invest.dat")
t <- 15
nt <- nrow(invest)
n <- nt/t

i <- invest[,1]       # investment/assets                 
q <- invest[,2]       # Tobin's Q                      
c <- invest[,3]       # cash-flow/assets                  
d <- invest[,4]       # debt/assets                         

max_lag <- 1
rhat1 <- 0.0157
rhat2 <- 0.5362

thresh <- d
tt <- t-max_lag

lag_value <- function(x,n,t,lagn,max_lag){
          y <- matrix(c(0),nrow=n,ncol=t)
          for (i in 1:n) {
              y[i,]<-x[(1+(i-1)*t):(t*i)]
          }  
          y <- y[,(1+max_lag-lagn):(t-lagn)]
          y <- matrix(t(y),nrow=nrow(y)*ncol(y),ncol=1)
          y
}

i0 <- lag_value(i,n,t,0,max_lag)
q1 <- lag_value(q,n,t,1,max_lag)
c1 <- lag_value(c,n,t,1,max_lag)
d1 <- lag_value(d,n,t,1,max_lag)

xx <- cbind(i0,q1,c1,d1)
nnt <- nrow(xx)

for (j in 1:4){
      xx[,j] <- sort(xx[,j])
}

qn1<-round(nnt/4)
qn2<-round(nnt/2)
qn3<-round(nnt*0.75)
x0 <- as.matrix(xx[1,])
x1 <- as.matrix(xx[qn1,])
x2 <- as.matrix(xx[qn2,])
x3 <- as.matrix(xx[qn3,])
x4 <- as.matrix(xx[nnt,])



e1 <- (d1 <= rhat1)
e2 <- (d1 <= rhat2) - e1
e3 <- 1 - e1 - e2
f1 <- matrix(c(0),nrow=n,ncol=tt)
f2 <- matrix(c(0),nrow=n,ncol=tt)
f3 <- matrix(c(0),nrow=n,ncol=tt)
for (i in 1:n){
    f1[i,]<-e1[(1+(i-1)*tt):(tt*i)]
    f2[i,]<-e2[(1+(i-1)*tt):(tt*i)]
    f3[i,]<-e3[(1+(i-1)*tt):(tt*i)]
}
g1 <- colMeans(f1)
g2 <- colMeans(f2)
g3 <- colMeans(f3)

g <- cbind(g1,g2,g3)
g <- round(g*100)

for (i in 1:1){
cat ("Full Sample Summary Statistics", "\n")
cat ("\n")
x0 <- format(x0, digits = 4)
x1 <- format(x1, digits = 4)
x2 <- format(x2, digits = 4)
x3 <- format(x3, digits = 4)
x4 <- format(x4, digits = 4)
for (j in 1:4){
cat (x0[j]," ",x1[j]," ",x2[j]," ",x3[j]," ",x4[j], "\n")
}
cat ("\n")
st <- seq(1974,1987,by=1)
cat ("Thresholds", rhat1, rhat2, "\n")
cat ("\n")
cat ("Percentage of Firms in Three Regimes, By Year", "\n")
for (j in 1:tt){
    cat("  ",st[j],"     ",g[j,1],"     ",g[j,2],"     ",g[j,3],"\n") 
}
cat ("\n")
cat ("\n")
}
save.image(file = "data.RData")





