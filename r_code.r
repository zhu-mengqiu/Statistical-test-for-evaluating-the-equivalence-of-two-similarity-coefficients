#simulate 2 lists of z-scores
sim.z <- function(n, pi, mu, sigma){
  g=3
  data <- sample(1:g, n, replace=TRUE, prob=pi)
  data.x <- data
  data.y <- data
  for(i in 1:g){
    l <- length(data[data==i]) 
    if(l>0){ 
      data.x[data==i] <- rnorm(l, mean=mu[i], sd=sigma[i])
      data.y[data==i] <- rnorm(l, mean=mu[i], sd=sigma[i])
    }
  }
  return(data.frame(cbind(data.x, data.y)))
}

#calculate the empirical power of two-sample test comparing independent similarity coefficients for a sequence of k
powerest<-function(K,B,n,pi,mu_1,mu_2,sigma_1,sigma_2,c,sim.z){
  z_hat<-FBM(K, B)
  tmp<-foreach (s = 1:B, .combine = 'c') %:%
    foreach (k = 1:K, .combine = 'c') %dopar% {
      zdata <- sim.z(n,pi, mu_1,sigma_1)
      z1<-zdata$data.x
      z2<-zdata$data.y
      z_and_l<-zdata[which(zdata$data.x>c[1] & zdata$data.y >c[1]),]
      z_and_s<-zdata[which(zdata$data.x<c[2] & zdata$data.y<c[2]),]
      z_or_l<-zdata[which(zdata$data.x>c[1] | zdata$data.y >c[1]),]
      z_or_s<-zdata[which(zdata$data.x<c[2] | zdata$data.y<c[2]),]
      z_cross1<-zdata[which(zdata$data.x<c[2] & zdata$data.y>c[1]),]
      z_cross2<-zdata[which(zdata$data.x>c[1] & zdata$data.y<c[2]),]
      p_x_hat1<-(dim(z_and_l)[1]+dim(z_and_s)[1])/n
      p_z_hat1<-(dim(z_or_l)[1]+dim(z_or_s)[1]-dim(z_cross1)[1]
                 -dim(z_cross2)[1]-dim(z_and_l)[1]-dim(z_and_s)[1])/n
      
      
      zdata <- sim.z(n,pi, mu_2,sigma_2)
      z1<-zdata$data.x
      z2<-zdata$data.y
      z_and_l<-zdata[which(zdata$data.x>c[1] & zdata$data.y >c[1]),]
      z_and_s<-zdata[which(zdata$data.x<c[2] & zdata$data.y<c[2]),]
      z_or_l<-zdata[which(zdata$data.x>c[1] | zdata$data.y >c[1]),]
      z_or_s<-zdata[which(zdata$data.x<c[2] | zdata$data.y<c[2]),]
      z_cross1<-zdata[which(zdata$data.x<c[2] & zdata$data.y>c[1]),]
      z_cross2<-zdata[which(zdata$data.x>c[1] & zdata$data.y<c[2]),]
      p_x_hat2<-(dim(z_and_l)[1]+dim(z_and_s)[1])/n
      p_z_hat2<-(dim(z_or_l)[1]+dim(z_or_s)[1]-dim(z_cross1)[1]
                 -dim(z_cross2)[1]-dim(z_and_l)[1]-dim(z_and_s)[1])/n
      
      theta_hat<-k*p_x_hat1/(k*p_x_hat1+p_z_hat1)
      -k*p_x_hat2/(k*p_x_hat2+p_z_hat2)
      variance_hat<-(k^2*p_z_hat1*p_x_hat1*(p_z_hat1+p_x_hat1))
      /(n*(k*p_x_hat1+p_z_hat1)^4)
      +(k^2*p_z_hat2*p_x_hat2*(p_z_hat2+p_x_hat2))
      /(n*(k*p_x_hat2+p_z_hat2)^4)
      sd_hat<-sqrt(variance_hat)
      z_hat[k,s]<-theta_hat/sd_hat
      
      NULL
    }
  power<-rep(0,K)
  for (i in 1:K){
    power[i]<- length(which(abs(z_hat[i,])>qnorm(1-0.05/2)))/B
  }
  return(power)
}

#calculate the theoritical power of two-sample test comparing independent similarity coefficients for a sequence of k
powertheory <- function(K,mu_1,mu_all,sigma_1,sigma_all,z_alpha,c,pi){
  mu=mu_1
  sigma=sigma_1
  p_x_theory1=0
  p_z_theory1=0
  for (i in 1:3){
    p_x<-(pnorm(c[2],mu[i],sigma[i]))^2
    +(1-pnorm(c[1],mu[i],sigma[i]))^2
    p_xi<-pnorm(c[2],mu[i],sigma[i])*(1-pnorm(c[1],mu[i],sigma[i]))*2
    p_z<-2*pnorm(c[2],mu[i],sigma[i])
    +2*(1-pnorm(c[1],mu[i],sigma[i]))-2*p_xi
    p_x_theory1<-p_x_theory1+pi[i]*p_x
    p_z_theory1<-p_z_theory1+pi[i]*p_z
  }
  
  power_theory<-matrix(0,K,length(mu_all))
  for (j in 1:length(mu_all)){
    mu=c(0,-mu_all[j],mu_all[j]) 
    sigma=rep(sigma_all[j],3)
    p_x_theory2=0
    p_z_theory2=0
    for (i in 1:3){
      p_x<-(pnorm(c[2],mu[i],sigma[i]))^2
      +(1-pnorm(c[1],mu[i],sigma[i]))^2
      p_xi<-pnorm(c[2],mu[i],sigma[i])
      *(1-pnorm(c[1],mu[i],sigma[i]))*2
      p_z<-2*pnorm(c[2],mu[i],sigma[i])
      +2*(1-pnorm(c[1],mu[i],sigma[i]))-2*p_xi
      p_x_theory2<-p_x_theory2+pi[i]*p_x
      p_z_theory2<-p_z_theory2+pi[i]*p_z
    }
    
    for (k in 1:K){
      theta <- k*p_x_theory1/(k*p_x_theory1+p_z_theory1)-k*p_x_theory2/(k*p_x_theory2+p_z_theory2)
      variance_hat <- (k^2*p_z_theory1*p_x_theory1*(p_z_theory1+p_x_theory1))/(n*(k*p_x_theory1+p_z_theory1)^4)+(k^2*p_z_theory2*p_x_theory2*(p_z_theory2+p_x_theory2))/(n*(k*p_x_theory2+p_z_theory2)^4)
      sd_hat<-sqrt(variance_hat)
      phi_1<-pnorm(z_alpha-theta/sd_hat)
      phi_2<-pnorm(-z_alpha-theta/sd_hat)
      power_theory[k,j]<-1-phi_1+phi_2
    }
  }
  return(power_theory)
}

#simulate 3 lists of z-scores
sim.z <- function(n, pi, mu1, sigma1, mu2, sigma2){
  g=3
  data <- sample(1:g, n, replace=TRUE, prob=pi)
  data.x <- data
  data.y <- data
  data.z <- data
  for(i in 1:g){
    l <- length(data[data==i]) 
    if(l>0){ 
      data.x[data==i] <- rnorm(l, mean=mu1[i], sd=sigma1[i]) 
      data.y[data==i] <- rnorm(l, mean=mu1[i], sd=sigma1[i])
      data.z[data==i] <- rnorm(l, mean=mu2[i], sd=sigma2[i])
    }
  }
  return(data.frame(cbind(data.x, data.y, data.z)))
}

#calculate the empirical power of two-sample test comparing dependent similarity coefficients for a sequence of k
powerest<-function(K,B,n,pi,mu1,mu2,sigma1,sigma2,c,sim.z){
  z_hat<-FBM(K, B)
  tmp<-foreach (s = 1:B, .combine = 'c') %:%
    foreach (k = 1:K, .combine = 'c') %dopar% {
      zdata <- sim.z(n,pi, mu1, sigma1, mu2, sigma2)
      z1<-zdata$data.x
      z2<-zdata$data.y
      z3 <- zdata$data.z
      
      z_1 <- zdata[which( ((z1>c[1] & z2 >c[1])|(z1<c[2] & z2<c[2]))
                          &((z1>c[1] & z3 >c[1])|(z1<c[2] & z3<c[2]))  ),]
      z_2 <- zdata[which( ((z1>c[1] & z2 <c[1])| (z1<c[2] & z2>c[2]) | ( z1<c[1] & z1>c[2] & z2>c[1])| ( z1<c[1] & z1>c[2] & z2<c[2]))
                          &((z1>c[1] & z3 <c[1])| (z1<c[2] & z3>c[2]) | ( z1<c[1] & z1>c[2] & z3>c[1])| ( z1<c[1] & z1>c[2] & z3<c[2]))  ),]
      z_3 <- zdata[which( (z1<c[1] & z1>c[2] & z2<c[1] & z2>c[2])
                          &(z1<c[1] & z1>c[2] & z3<c[1] & z3>c[2])  ),]
      z_4 <- zdata[which( ((z1>c[1] & z2 >c[1])|(z1<c[2] & z2<c[2]))
                          &((z1>c[1] & z3 <c[1])| (z1<c[2] & z3>c[2]) | ( z1<c[1] & z1>c[2] & z3>c[1])| ( z1<c[1] & z1>c[2] & z3<c[2]))  ),]
      z_5 <- zdata[which( ((z1>c[1] & z3 >c[1])|(z1<c[2] & z3<c[2]))
                          &((z1>c[1] & z2 <c[1])| (z1<c[2] & z2>c[2]) | ( z1<c[1] & z1>c[2] & z2>c[1])| ( z1<c[1] & z1>c[2] & z2<c[2]))  ),]
      z_6 <- zdata[which( (z1<c[1] & z1>c[2] & z2<c[1] & z2>c[2])
                          &((z1>c[1] & z3 <c[1])| (z1<c[2] & z3>c[2]) | ( z1<c[1] & z1>c[2] & z3>c[1])| ( z1<c[1] & z1>c[2] & z3<c[2]))  ),]
      z_7 <- zdata[which( (z1<c[1] & z1>c[2] & z3<c[1] & z3>c[2])
                          &((z1>c[1] & z2 <c[1])| (z1<c[2] & z2>c[2]) | ( z1<c[1] & z1>c[2] & z2>c[1])| ( z1<c[1] & z1>c[2] & z2<c[2]))  ),]
      p.hat1<-(dim(z_1)[1])/length(z1)
      p.hat2<-(dim(z_2)[1])/length(z1)
      p.hat3<-(dim(z_3)[1])/length(z1)
      p.hat4<-(dim(z_4)[1])/length(z1)
      p.hat5<-(dim(z_5)[1])/length(z1)
      p.hat6<-(dim(z_6)[1])/length(z1)
      p.hat7<-(dim(z_7)[1])/length(z1)
      
      
      theta_hat <- ((k*p.hat1+k*p.hat4)/(k*p.hat1+k*p.hat4+p.hat5+p.hat7+p.hat2)-(k*p.hat1+k*p.hat5)/(k*p.hat1+k*p.hat5+p.hat4+p.hat6+p.hat2))
      variance_hat<-(k^2*(p.hat1+p.hat4)*(p.hat5+p.hat7+p.hat2)*(p.hat1+p.hat4+p.hat5+p.hat7+p.hat2))/(n*(k*p.hat1+k*p.hat4+p.hat5+p.hat7+p.hat2)^4)+
        (k^2*(p.hat1+p.hat5)*(p.hat4+p.hat6+p.hat2)*(p.hat1+p.hat5+p.hat4+p.hat6+p.hat2))/(n*(k*p.hat1+k*p.hat5+p.hat4+p.hat6+p.hat2)^4)-
        (2*k^2*((p.hat1*p.hat2-p.hat4*p.hat5)*(p.hat1+p.hat2+p.hat4+p.hat5+p.hat6+p.hat7)+p.hat1*p.hat6*p.hat7))/(n*(k*p.hat1+k*p.hat4+p.hat5+p.hat7+p.hat2)^2*(k*p.hat1+k*p.hat5+p.hat4+p.hat6+p.hat2)^2)
      sd_hat<-sqrt(variance_hat)
      z_hat[k,s]<-theta_hat/sd_hat
      
    }
  
  power<-rep(0,K)
  for (i in 1:K){
    power[i]<- length(which(abs(z_hat[i,])>qnorm(1-0.05/2)))/B
  }
  return(power)
}

#calculate the theoritical power of two-sample test comparing dependent similarity coefficients for a sequence of k
powertheory <- function(K,mu1,mu_all,sigma1,sigma_all,z_alpha,c,pi){
  
  power_theory<-matrix(0,length(mu_all),K)
  for (j in 1:length(mu_all)){
    
    mu2=c(0,-mu_all[j],mu_all[j]) 
    sigma2=rep(sigma_all[j],3)
    
    p1<-0
    p2<-0
    p3<-0
    p4<-0
    p5<-0
    p6<-0
    p7<-0
    
    for (i in 1:3){
      p1.temp <- (pnorm(c[2],mu1[i],sigma1[i]))^2*(pnorm(c[2],mu2[i],sigma2[i]))^2 + 
        (1-pnorm(c[1],mu1[i],sigma1[i]))^2*(1-pnorm(c[1],mu2[i],sigma2[i]))
      p1 <- p1 + p1.temp*pi[i]
      p3.temp <- (pnorm(c[1],mu1[i],sigma1[i]) - pnorm(c[2],mu1[i],sigma1[i]))^2*(pnorm(c[1],mu2[i],sigma2[i]) - pnorm(c[2],mu2[i],sigma2[i]))
      p3 <- p3 + p3.temp*pi[i]
      p4.temp <- (1-pnorm(c[1],mu1[i],sigma1[i]))^2*(pnorm(c[1],mu2[i],sigma2[i]) - pnorm(c[2],mu2[i],sigma2[i])) + 
        (1-pnorm(c[1],mu1[i],sigma1[i]))^2*(pnorm(c[2],mu2[i],sigma2[i])) + 
        (pnorm(c[2],mu1[i],sigma1[i]))^2*(pnorm(c[1],mu2[i],sigma2[i]) - pnorm(c[2],mu2[i],sigma2[i])) + 
        (pnorm(c[2],mu1[i],sigma1[i]))^2*(1-pnorm(c[1],mu2[i],sigma2[i]))
      p4 <- p4 + p4.temp*pi[i]
      p5.temp <-  (1-pnorm(c[1],mu1[i],sigma1[i]))*(pnorm(c[2],mu1[i],sigma1[i]))*(1-pnorm(c[1],mu2[i],sigma2[i])) + 
        (1-pnorm(c[1],mu1[i],sigma1[i]))*(pnorm(c[2],mu1[i],sigma1[i]))*(pnorm(c[2],mu2[i],sigma2[i])) + 
        (pnorm(c[1],mu1[i],sigma1[i]) - pnorm(c[2],mu1[i],sigma1[i]))*(pnorm(c[2],mu1[i],sigma1[i]))*(pnorm(c[2],mu2[i],sigma2[i])) + 
        (pnorm(c[1],mu1[i],sigma1[i]) - pnorm(c[2],mu1[i],sigma1[i]))*(1-pnorm(c[1],mu1[i],sigma1[i]))*(1-pnorm(c[1],mu2[i],sigma2[i]))
      p5 <- p5 + p5.temp*pi[i]
      p6.temp <- (pnorm(c[1],mu1[i],sigma1[i]) - pnorm(c[2],mu1[i],sigma1[i]))^2*(pnorm(c[2],mu2[i],sigma2[i])) +
        (pnorm(c[1],mu1[i],sigma1[i]) - pnorm(c[2],mu1[i],sigma1[i]))^2*(1-pnorm(c[1],mu2[i],sigma2[i]))
      p6 <- p6 + p6.temp*pi[i]
      p7.temp <- (pnorm(c[1],mu1[i],sigma1[i]) - pnorm(c[2],mu1[i],sigma1[i]))*(pnorm(c[2],mu1[i],sigma1[i]))*(pnorm(c[1],mu2[i],sigma2[i]) - pnorm(c[2],mu2[i],sigma2[i])) + 
        (pnorm(c[1],mu1[i],sigma1[i]) - pnorm(c[2],mu1[i],sigma1[i]))*(1-pnorm(c[1],mu1[i],sigma1[i]))*(pnorm(c[1],mu2[i],sigma2[i]) - pnorm(c[2],mu2[i],sigma2[i]))
      p7 <- p7 + p7.temp*pi[i]
      p2.temp <- 1-p1.temp-p3.temp-p4.temp-p5.temp-p6.temp-p7.temp
      p2 <- p2 + p2.temp*pi[i]
      
    }
    
    for (k in 1:K){
      theta_hat <- ((k*p1+k*p4)/(k*p1+k*p4+p5+p7+p2)-(k*p1+k*p5)/(k*p1+k*p5+p4+p6+p2))
      variance_hat<-(k^2*(p1+p4)*(p5+p7+p2)*(p1+p4+p5+p7+p2))/(n*(k*p1+k*p4+p5+p7+p2)^4)+
        (k^2*(p1+p5)*(p4+p6+p2)*(p1+p5+p4+p6+p2))/(n*(k*p1+k*p5+p4+p6+p2)^4)-
        (2*k^2*((p1*p2-p4*p5)*(p1+p2+p4+p5+p6+p7)+p1*p6*p7))/(n*(k*p1+k*p4+p5+p7+p2)^2*(k*p1+k*p5+p4+p6+p2)^2)
      sd_hat<-sqrt(variance_hat)
      phi_1<-pnorm(z_alpha-theta_hat/sd_hat)
      phi_2<-pnorm(-z_alpha-theta_hat/sd_hat)
      power_theory[j,k]<-1-phi_1+phi_2
    }
  }
  return(power_theory)
}

#calculate test statistic for two sample test comparing independent similarity coefficients
similarity_coefficient.two_sample_test <- function(z1_1,z2_1,z1_2,z2_2,c,k,n1,n2){
  zdata_1<-data.frame(cbind(z1_1, z2_1))
  z.and.l_1<-zdata_1[which(z1_1>c[1] & z2_1 >c[1]),]
  z.and.s_1<-zdata_1[which(z1_1<c[2] & z2_1<c[2]),]
  z.or.l_1<-zdata_1[which(z1_1>c[1] | z2_1 >c[1]),]
  z.or.s_1<-zdata_1[which(z1_1<c[2] | z2_1<c[2]),]
  z.cross1_1<-zdata_1[which(z1_1<c[2] & z2_1>c[1]),]
  z.cross2_1<-zdata_1[which(z1_1>c[1] & z2_1<c[2]),]
  r.u.hat_1<-(dim(z.and.l_1)[1]+dim(z.and.s_1)[1])/length(z1_1)
  r.v.hat_1<-(dim(z.or.l_1)[1]+dim(z.or.s_1)[1]-dim(z.cross1_1)[1]-dim(z.cross2_1)[1]-dim(z.and.l_1)[1]-dim(z.and.s_1)[1])/length(z1_1)
  
  zdata_2<-data.frame(cbind(z1_2, z2_2))
  z.and.l_2<-zdata_2[which(z1_2>c[1] & z2_2 >c[1]),]
  z.and.s_2<-zdata_2[which(z1_2<c[2] & z2_2<c[2]),]
  z.or.l_2<-zdata_2[which(z1_2>c[1] | z2_2 >c[1]),]
  z.or.s_2<-zdata_2[which(z1_2<c[2] | z2_2<c[2]),]
  z.cross1_2<-zdata_2[which(z1_2<c[2] & z2_2>c[1]),]
  z.cross2_2<-zdata_2[which(z1_2>c[1] & z2_2<c[2]),]
  r.u.hat_2<-(dim(z.and.l_2)[1]+dim(z.and.s_2)[1])/length(z1_2)
  r.v.hat_2<-(dim(z.or.l_2)[1]+dim(z.or.s_2)[1]-dim(z.cross1_2)[1]-dim(z.cross2_2)[1]-dim(z.and.l_2)[1]-dim(z.and.s_2)[1])/length(z1_2)
  
  theta.hat <- k*r.u.hat_1/(k*r.u.hat_1+r.v.hat_1)-k*r.u.hat_2/(k*r.u.hat_2+r.v.hat_2)
  variance.hat_1 <- (k^2*r.v.hat_1*r.u.hat_1*(r.v.hat_1+r.u.hat_1))/(n1*(k*r.u.hat_1+r.v.hat_1)^4)
  variance.hat_2 <- (k^2*r.v.hat_2*r.u.hat_2*(r.v.hat_2+r.u.hat_2))/(n2*(k*r.u.hat_2+r.v.hat_2)^4)
  test_statistic <- theta.hat/sqrt(variance.hat_1 + variance.hat_2)
  p.value <- (1 - pnorm(abs(test_statistic)))*2
  
  return(cbind(test_statistic, p.value))
}

#calculate test statistic for two sample test comparing dependent similarity coefficients
similarity_coefficient.two_sample_test_dependent <- function(z1,z2,z3,c,k,n){
  zdata<-data.frame(cbind(z1, z2,z3))
  
  z_1 <- zdata[which( ((z1>c[1] & z2 >c[1])|(z1<c[2] & z2<c[2]))
                      &((z1>c[1] & z3 >c[1])|(z1<c[2] & z3<c[2]))  ),]
  z_2 <- zdata[which( ((z1>c[1] & z2 <c[1])| (z1<c[2] & z2>c[2]) | ( z1<c[1] & z1>c[2] & z2>c[1])| ( z1<c[1] & z1>c[2] & z2<c[2]))
                      &((z1>c[1] & z3 <c[1])| (z1<c[2] & z3>c[2]) | ( z1<c[1] & z1>c[2] & z3>c[1])| ( z1<c[1] & z1>c[2] & z3<c[2]))  ),]
  z_3 <- zdata[which( (z1<c[1] & z1>c[2] & z2<c[1] & z2>c[2])
                      &(z1<c[1] & z1>c[2] & z3<c[1] & z3>c[2])  ),]
  z_4 <- zdata[which( ((z1>c[1] & z2 >c[1])|(z1<c[2] & z2<c[2]))
                      &((z1>c[1] & z3 <c[1])| (z1<c[2] & z3>c[2]) | ( z1<c[1] & z1>c[2] & z3>c[1])| ( z1<c[1] & z1>c[2] & z3<c[2]))  ),]
  z_5 <- zdata[which( ((z1>c[1] & z3 >c[1])|(z1<c[2] & z3<c[2]))
                      &((z1>c[1] & z2 <c[1])| (z1<c[2] & z2>c[2]) | ( z1<c[1] & z1>c[2] & z2>c[1])| ( z1<c[1] & z1>c[2] & z2<c[2]))  ),]
  z_6 <- zdata[which( (z1<c[1] & z1>c[2] & z2<c[1] & z2>c[2])
                      &((z1>c[1] & z3 <c[1])| (z1<c[2] & z3>c[2]) | ( z1<c[1] & z1>c[2] & z3>c[1])| ( z1<c[1] & z1>c[2] & z3<c[2]))  ),]
  z_7 <- zdata[which( (z1<c[1] & z1>c[2] & z3<c[1] & z3>c[2])
                      &((z1>c[1] & z2 <c[1])| (z1<c[2] & z2>c[2]) | ( z1<c[1] & z1>c[2] & z2>c[1])| ( z1<c[1] & z1>c[2] & z2<c[2]))  ),]
  
  p.hat1<-(dim(z_1)[1])/length(z1)
  p.hat2<-(dim(z_2)[1])/length(z1)
  p.hat3<-(dim(z_3)[1])/length(z1)
  p.hat4<-(dim(z_4)[1])/length(z1)
  p.hat5<-(dim(z_5)[1])/length(z1)
  p.hat6<-(dim(z_6)[1])/length(z1)
  p.hat7<-(dim(z_7)[1])/length(z1)
  
  theta_hat <- ((k*p.hat1+k*p.hat4)/(k*p.hat1+k*p.hat4+p.hat5+p.hat7+p.hat2)-(k*p.hat1+k*p.hat5)/(k*p.hat1+k*p.hat5+p.hat4+p.hat6+p.hat2))
  variance_hat<-(k^2*(p.hat1+p.hat4)*(p.hat5+p.hat7+p.hat2)*(p.hat1+p.hat4+p.hat5+p.hat7+p.hat2))/(n*(k*p.hat1+k*p.hat4+p.hat5+p.hat7+p.hat2)^4)+
    (k^2*(p.hat1+p.hat5)*(p.hat4+p.hat6+p.hat2)*(p.hat1+p.hat5+p.hat4+p.hat6+p.hat2))/(n*(k*p.hat1+k*p.hat5+p.hat4+p.hat6+p.hat2)^4)-
    (2*k^2*((p.hat1*p.hat2-p.hat4*p.hat5)*(p.hat1+p.hat2+p.hat4+p.hat5+p.hat6+p.hat7)+p.hat1*p.hat6*p.hat7))/(n*(k*p.hat1+k*p.hat4+p.hat5+p.hat7+p.hat2)^2*(k*p.hat1+k*p.hat5+p.hat4+p.hat6+p.hat2)^2)
  sd_hat<-sqrt(variance_hat)
  test_statistic <- theta_hat/sd_hat
  p.value <- (1 - pnorm(abs(test_statistic)))*2
  
  return(cbind(test_statistic, p.value))
}

### Running samples ###
K=10
pi=c(0.6,0.2,0.2)
mu=c(0,-3,3) 
mu_1<-mu
mu_all<-c(3,3,3)
sigma_1=c(1,1,1)
sigma_all<-seq(0.94,0.98,by=0.02)
c<-c(1.96,-1.96)
powertheory(K,mu_1,mu_all,sigma_1,sigma_all,z_alpha,c,pi)


K=10
n=20000 
B=100000
pi=c(0.6,0.2,0.2)
c<-c(1.96,-1.96)
mu_1=c(0,-3,3) 
mu_2=c(0,-3,3)
sigma_1 <- c(1,1,1)
sigma_2 <- c(0.94,0.94,0.94)
library(bigstatsr)
library(foreach)
library(doParallel)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
powerest(K,B,n,pi,mu_1,mu_2,sigma_1,sigma_2,c,sim.z)
parallel::stopCluster(cl)


#TCGA data are available at https://gdac.broadinstitute.org
microarray<-read.delim("BRCA.medianexp.txt",header = FALSE)
RNAseq<-read.delim("BRCA.rnaseqv2__illuminahiseq_rnaseqv2
                     __unc_edu__Level_3__RSEM_genes_normalized
                     __data.data.txt",header = FALSE)
#case type name
RNAseqNT<-read.delim("RNAseq NT.txt",header = FALSE)
RNAseqNT<-RNAseqNT[RNAseqNT$V5%in%c("RSEM_genes"),]
RNAseqNT<-as.matrix(RNAseqNT)
RNAseqTM<-read.delim("RNAseq TM.txt",header = FALSE)
RNAseqTM<-as.matrix(RNAseqTM)
microarrayNT<-read.delim("microarray NT.txt",header = FALSE)
microarrayNT<-as.matrix(microarrayNT)
microarrayTM<-read.delim("microarray TM.txt",header = FALSE)
microarrayTM<-as.matrix(microarrayTM)


#seperate case
`%notin%` <- Negate(`%in%`)

microarray<-as.matrix(microarray)
colnames(microarray)<-microarray[1,]
microarray_NT<-microarray[names(microarray) %in% c(microarrayNT[,1])]
microarray_NT<-microarray_NT[c(-1,-2),]
microarray_TM<-microarray[names(microarray) %in% c(microarrayTM[,1])]
microarray_TM<-microarray_TM[c(-1,-2),]
microarray_TP<-microarray[names(microarray) %notin% c(microarrayTM[,1],microarrayNT[,1])]
microarray_TP<-microarray_TP[,-1]
microarray_TP<-microarray_TP[c(-1,-2),]


RNAseq1<-as.matrix(RNAseq)
colnames(RNAseq)<-RNAseq1[1,]
RNAseq_NT<-RNAseq[names(RNAseq) %in% c(RNAseqNT[,1])]
RNAseq_NT<-RNAseq_NT[c(-1,-2),]
RNAseq_TM<-RNAseq[names(RNAseq) %in% c(RNAseqTM[,1])]
RNAseq_TM<-RNAseq_TM[c(-1,-2),]
RNAseq_TP<-RNAseq[names(RNAseq)%notin%
                      c(RNAseqTM[,1],RNAseqNT[,1])]
RNAseq_TP<-RNAseq_TP[,-1]
RNAseq_TP<-RNAseq_TP[c(-1,-2),]

all_zvalues_sim<-function(RNAseq_NT_1,RNAseq_TP_1,
                          RNAseq_NT_2,RNAseq_TP_2,
                          microarray_NT_1,microarray_TP_1,
                          microarray_NT_2,microarray_TP_2){
  pvalue.statistic.RNAseq_1<-foreach (i = 1:dim(RNAseq_NT_1)[1],.combine = 'c') %dopar%{
    A <- wilcox.test(as.numeric(as.matrix(RNAseq_NT_1[i,])),
                     as.numeric(as.matrix(RNAseq_TP_1[i,])),alternative = "greater")
    rbind(A$p.value,A$statistic)
  }
  pvalue.statistic.RNAseq_2<-foreach (i = 1:dim(RNAseq_NT_2)[1],.combine = 'c') %dopar%{
    B <- wilcox.test(as.numeric(as.matrix(RNAseq_NT_2[i,])),
                     as.numeric(as.matrix(RNAseq_TP_2[i,])),alternative = "greater")
    rbind(B$p.value,B$statistic)
  }
  pvalue.statistic.microarray_1<-foreach (i = 1:dim(microarray_NT_1)[1],.combine = 'c') %dopar%{
    C <- wilcox.test(as.numeric(as.matrix(microarray_NT_1[i,])),
                     as.numeric(as.matrix(microarray_TP_1[i,])),alternative = "greater")
    rbind(C$p.value,C$statistic)
  }
  pvalue.statistic.microarray_2<-foreach (i = 1:dim(microarray_NT_2)[1],.combine = 'c') %dopar%{
    D <- wilcox.test(as.numeric(as.matrix(microarray_NT_2[i,])),
                     as.numeric(as.matrix(microarray_TP_2[i,])),alternative = "greater")
    rbind(D$p.value,D$statistic)
  }
  
  return(cbind(pvalue.statistic.RNAseq_1,pvalue.statistic.RNAseq_2,
               pvalue.statistic.microarray_1,pvalue.statistic.microarray_2 ))
}


