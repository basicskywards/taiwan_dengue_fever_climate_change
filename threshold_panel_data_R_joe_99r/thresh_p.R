##  thresh_p.R
##  
##  This is a R program file.
##  It replicates the estimation, testing and graphs reported in
##  "Threshold Effects in Non-Dynamic Panels:
##  Estimation, Testing and Inference"
##  
##  For questions, please contact
##  
##  Bruce E. Hansen
##  Department of Economics
##  Social Science Building
##  University of Wisconsin
##  Madison, WI 53706-1393
##  bhansen@wisc.edu
##  http://www.ssc.wisc.edu/~bhansen/
##  
##  
##  This program file loads the dataset "invest.dat".
##  It creates the output file "thresh.RData"
##  
########################################################################

invest <- read.table("invest.dat")
t <- 15
nt <- nrow(invest)
n <- nt/t

i <- invest[,1]    # investment/assets                 
q <- invest[,2]    # Tobin's Q                      
c <- invest[,3]    # cash-flow/assets                  
d <- invest[,4]    # debt/assets     


qn <- 400          # number of quantiles to examine                
conf_lev <- 0.95   # confidence level for threshold                
vgraph <- 1        # set to 1 to graph likelihood ratios           
boot_1 <- 300      # number of replications, 0 for no bootstrap, single (300) 
boot_2 <- 300      # number of replications, 0 for no bootstrap, double (300) 
boot_3 <- 300      # number of replications, 0 for no bootstrap, triple (300) 
trim_1 <- 0.01     # percentage to trim before search, single      
trim_2 <- 0.01     # percentage to trim before search, double      
trim_3 <- 0.05     # percentage to trim before search, triple      

max_lag <- 1
tt <- t-max_lag
ty <- n*(t-max_lag-1)

lag_v <- function(x,lagn){
      yl <- matrix(c(0),nrow=n,ncol=t)
      for (i in 1:n) {
          yl[i,]<-x[(1+(i-1)*t):(t*i)]
      }  
      yl <- yl[,(1+max_lag-lagn):(t-lagn)]
      out <- matrix(t(yl),nrow=nrow(yl)*ncol(yl),ncol=1)
      out
}

tr <- function(y){
   yf <- matrix(c(0),nrow=n,ncol=tt)
   for (i in 1:n) {
       yf[i,]<-y[(1+(i-1)*tt):(tt*i)]
   }
   yfm <- yf- colMeans(t(yf))
   yfm <- yfm[,1:(tt-1)]
   out <- matrix(t(yfm),nrow=nrow(yfm)*ncol(yfm),ncol=1)
   out
}

y  <- lag_v(i,0)  
cf <- lag_v(c,1) 
q1 <- lag_v(q,1)
d1 <- lag_v(d,1)   # set to threshold variable 
yt <- tr(y)
ct <- tr(cf)

x <- cbind(q1,(q1^2),(q1^3),d1,(q1*d1))
k <- ncol(x)
xt <- matrix(c(0),nrow=nrow(yt),ncol=k) 
for (j in 1:k) xt[,j]=tr(x[,j])
thresh <- d1
dd <- unique(thresh)
dd <- as.matrix(sort(dd))
qnt1 <- qn*trim_1
sq <- as.matrix(seq(trim_1,trim_1+(1/qn)*(qn-2*qnt1),by=1/qn))   
qq1 <- as.matrix(dd[floor(sq*nrow(dd))])
qn1 <- nrow(qq1)
cc <- -2*log(1-sqrt(conf_lev))

sse_calc <- function(y,x){
         e <- y-x%*%qr.solve(x,y)
         out <- t(e)%*%e 
         out
} 

thr_sse <- function(y,q,r){
        nq <- nrow(q)
        sse <- matrix(c(0),nq,1)
        for (qi in 1:nq){
            if (r[1]==0) {rr <- q[qi] 
            }else{ rr <- rbind(r,q[qi])}
            rr <- as.matrix(sort(rr))
            xx <- cbind(xt,ct)
            for (j in 1:nrow(rr)){
                d <- (thresh < rr[j])
                xx <- cbind(xx,tr(cf*d))
            }
            sse[qi] <- sse_calc(y,xx)
        }
        sse
}

r_est <- function(y,r,trim){
      if (max(r)==0){
          qq <- qq1;
          rr <- 0;
      }else{rr <- as.matrix(sort(r))
            i <- as.matrix(seq(1,qn1,by=1))
            nn <- colSums(qq1%*%matrix(1,1,nrow(rr))< matrix(1,nrow(qq1),1)%*%t(rr))
            nn <- as.matrix(nn)
            qnt <- qn*trim          
            ii1 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn+qnt))
            ii2 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn-qnt))
            ii <- (ii1-ii2)%*%matrix(1,nrow(rr),1)
            qq <- as.matrix(qq1[ii!=1])
      } 
      sse <- thr_sse(y,qq,rr)
      rihat <- which.min(sse)
      list(sse_b=sse[rihat],rhat_b=qq[rihat])  
}

model <- function(r,trim,rep,it){       
      if (max(r)==0){
          qq <- qq1;
          rr <- 0;
      }else{rr <- as.matrix(sort(r))
            i <- as.matrix(seq(1,qn1,by=1))
            nn <- colSums(qq1%*%matrix(1,1,nrow(rr))< matrix(1,nrow(qq1),1)%*%t(rr))
            nn <- as.matrix(nn)
            qnt <- qn*trim          
            ii1 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn+qnt))
            ii2 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn-qnt))
            ii <- (ii1-ii2)%*%matrix(1,nrow(rr),1)
            qq <- as.matrix(qq1[ii!=1])
      }            
      sse <- thr_sse(yt,qq,rr)
      rihat <- which.min(sse)
      rhat <- qq[rihat]
      sse1 <- sse[rihat]
      lr <- (sse/sse1 - 1)*ty
      rhats <- as.matrix(qq[(lr < cc)])
      if (vgraph==1){
          if (it==0){ 
              titname=rbind("Figure 1","Confidence Interval Construction in Single Threshold Model")
              xname="Threshold Parameter"}
          if (it==1){ 
              titname=rbind("Figure 3","Confidence Interval Construction in Double Threshold Model")
              xname="First Threshold Parameter"}
          if (it==2){ 
              titname=rbind("Figure 2","Confidence Interval Construction in Double Threshold Model")
              xname="Second Threshold Parameter"}
          if (it==3){ 
              titname="Confidence Interval Construction in Triple Threshold Model"
              xname="Thrid Threshold Parameter"}
          yname="Likelihood Ratio"
          x11()
          xxlim <- range(qq)
          yylim <- range(rbind(lr,cc))
          plot(qq,lr,lty=1,col=1,xlim=xxlim,ylim=yylim,type="l",ann=0)
          lines(qq,matrix((1),nrow=nrow(qq),ncol=1)*cc,lty=2,col=2) 
          title(main=titname,xlab=xname,ylab=yname)      
      }     
      if (max(r) != 0){
          cat ("Fixed Thresholds       ", t(rr), "\n")
          rrr <- sort(rbind(rr,rhat))
      }else{ rrr <- rhat }
      rrr <- as.matrix(rrr)
      cat ("Threshold Estimate     ", rhat, "\n")
      cat ("Confidence Region      ", cbind(min(rhats),max(rhats)), "\n")
      cat ("Sum of Squared Errors  ", sse1, "\n")
      cat ("Trimming Percentage    ", trim, "\n")
      cat ("\n")
      cat ("\n")  
      nr <- nrow(rrr)
      xx <- xt
      dd <- matrix((0),nrow=nrow(thresh),ncol=nr)
      for (j in 1:nr){
          dd[,j] <- (thresh < rrr[j])
          d <- dd[,j]
          if (j>1) d <- d - dd[,(j-1)]
          xx <- cbind(xx,tr(cf*d))
      }
      d <- 1-dd[,nr]
      xx <- cbind(xx,tr(cf*d))
      xxi <- solve(t(xx)%*%xx)
      beta <- xxi%*%(t(xx)%*%yt)
      e <- yt - xx%*%beta
      xxe <- xx*(e%*%matrix((1),nrow=1,ncol=ncol(xx)))
      xxe <- t(xxe)%*%xxe
      sehet <- as.matrix(sqrt(diag(xxi%*%xxe%*%xxi)))
      sehomo <- as.matrix(sqrt(diag(xxi*as.vector((t(e)%*%e))/(ty-n-ncol(xx)))))
      beta <- cbind(beta,sehomo,sehet)
      cat ("Thresholds", "\n")
      cat (t(rrr), "\n")
      cat ("\n")
      cat ("Regime-independent Coefficients, standard errors, het standard errors", "\n")
      beta <- format(beta, digits = 4, scientific = FALSE)
      for (j in 1:k) cat (beta[j,], "\n")
      cat ("\n")
      cat ("Regime-dependent Coefficients, standard errors, het standard errors", "\n")
      for (j in (k+1):(k+nr+1)) cat (beta[j,], "\n")
      cat ("\n")
      cat ("\n")
      if (rep > 0){
          xx <- cbind(xt,ct)
          if (max(rr) != 0){
              for (j in 1:nrow(rr)) xx <- cbind(xx,tr(cf*(thresh < rr[j])))
          }
          yp <- xx%*%qr.solve(xx,yt)
          e <- yt-yp
          sse0 <- t(e)%*%e
          lrt <- (sse0/sse1-1)*ty
          cat ("LR Test for threshold effect  ", lrt, "\n")
          cat ("\n")
          cat ("\n")
          stats <- matrix(c(0),nrow=rep,ncol=1) 
          for (j in 1:rep){
              eb <- matrix(c(0),nrow=n,ncol=(tt-1))
              for (i in 1:n) {
                  eb[i,]<-e[(1+(i-1)*(tt-1)):((tt-1)*i)]
              }
              yeb <- t(eb[ceiling(runif(n)*n),])
              yb <- yp + matrix(yeb,nrow=nrow(yeb)*ncol(yeb),ncol=1)       
              sse0 <- sse_calc(yb,cbind(xt,ct))
              out <- r_est(yb,0,trim)
              sse1 <- out$sse_b
              rhat_b <- out$rhat_b
              rrr <- rhat_b
              if (max(r) != 0){
                  for (jj in 1:length(r)){
                      sse0 <- sse1
                      out <- r_est(yb,rrr,trim)
                      sse1 <- out$sse_b
                      rhat_b <- out$rhat_b
                      rrr <- rbind(rrr,rhat_b)
                  }
              }
              lrt_b <- (sse0/sse1-1)*ty
              stats[j] <- lrt_b
              cat ("Bootstrap Replication ", j, lrt_b, "\n")
          }
          cat ("\n")
          cat ("\n")
          stats <- as.matrix(sort(stats))
          crits <- as.matrix(stats[ceiling(rbind(.90,.95,.99)*rep)])
          cat ("Number of Bootstrap replications   ", rep, "\n")
          cat ("Bootstrap p-value                  ", mean(stats > as.vector(lrt)), "\n")
          cat ("Critical Values   ", crits[1], crits[2], crits[3], "\n")
          cat ("\n")
          cat ("\n")
      }
      rhat
}

for (i in 1:1){
cat ("Number of Firms        ", n, "\n")
cat ("Number of Years used   ", tt, "\n")
cat ("Total Observations     ", ty, "\n")
cat ("Number of Quantiles    ", qn, "\n")
cat ("Confidence Level       ", conf_lev, "\n")
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")
cat ("Zero Threshold Model", "\n")
sse0 <- sse_calc(yt,cbind(xt,ct))
cat ("Sum of Squared Errors                   ", sse0, "\n")
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")

cat ("Single Threshold Model",  "\n")
cat ("\n")
rhat1 <- model(0,trim_1,boot_1,0)
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")

cat ("Double Threshold Model", "\n")
cat ("Trimming Percentage    ", trim_2, "\n")
cat ("\n")
cat ("First Iteration", "\n")
rhat2 <- model(rhat1,trim_2,boot_2,2)
cat ("Second Iteration", "\n")
rhat1 <- model(rhat2,trim_2,0,1)
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")

cat ("Triple Threshold Model", "\n")
cat ("Trimming Percentage    ", trim_3, "\n")
cat ("\n")
rhat3 <- model(rbind(rhat1,rhat2),trim_3,boot_3,3)
cat ("\n")
cat ("\n")
cat ("*******************************************************", "\n")
cat ("\n")
cat ("\n")
}

save.image(file = "threshout.RData")


