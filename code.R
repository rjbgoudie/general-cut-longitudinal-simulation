library(ggplot2)

simulation<-function(time,num_rep,init_N,rate,beta,sigma,X){
  observation<-rep(0,time*num_rep)
  dim(observation)<-c(time,num_rep)
  N_update <- rep(0,time-1)
  for(j in 1:(time-1)){
    N_update[j] <- beta*rate[j]-2  #here 2 is bias 
  }
  for(k in 1:num_rep){
    for(i in 1:time){
      if(i==1){
        observation[i,k] <- rnorm(1, mean = sum(X[i,k,]*c(init_N,1,rate[i])),sd=2)
      }else{
        observation[i,k] <- rnorm(1, mean = sum(X[i,k,]*c(N_update[i-1],1,rate[i])),sd=2)
      }
    } 
  }
  re<-list(observation=observation,N_real=c(init_N,N_update))
  return(re)
}


rate <- 10*sin(seq(1,100))
time<-length(rate)
num_rep <- 100

init_N <- 0
beta <- 1
sigma <- 3
sigma_model <- 0.001

X <- array(0,dim=c(time,num_rep,3))
X1 <- matrix(0,nrow = time, ncol = num_rep)
X2 <- matrix(0,nrow = (time+1), ncol = num_rep)
for(i in 1:time){
  for(k in 1:num_rep){
    X[i,k,] <- c(rnorm(1,0,1),1,rnorm(1,0,1))
    X1[i,k] <- X[i,k,3]
    X2[i,k] <- X[i,k,1]
  }
}

X2[(time+1),] <- rep(0,num_rep)



################# Full Bayesian Bias ############################

Gamma <- matrix(0,nrow = 2*time, ncol = 2*time)

for(i in 1:(time-1)){
  Gamma[(2*i-1):(2*i),(2*i+1):(2*(i+1))] <- matrix(c(0,0,0,sum(X1[i+1,]*X2[i+1,])))
}

Gamma <- Gamma + t(Gamma)

for(i in 1:time){
  Gamma[(2*i-1):(2*i),(2*i-1):(2*i)] <- matrix(c(num_rep+1,sum(X1[i,]),sum(X1[i,]),1+sum(X1[i,]*X1[i,])+sum(X2[i+1,]*X2[i+1,])))
}

Gamma_inverse <- solve(Gamma)

XZ <- rep(0, 2*time)
for(i in 1:time){
  if(i == 1){
    XZ[1] <- 0
    XZ[2] <- sum(X2[2,]*X2[2,])
  }
  if(i == time){
    XZ[2*time-1] <- sum(X2[time,])
    XZ[2*time] <- sum(X1[time,]*X2[time,])
  }
  if(i>1 & i<time){
    XZ[2*i-1] <- sum(X2[i,])
    XZ[2*i] <- sum(X1[i,]*X2[i,])+sum(X2[i+1,]*X2[i+1,])
  }
}

Full_Bayesian_bias <- Gamma_inverse %*% XZ

Full_Bayesian_variance <- diag(Gamma_inverse)

# Coefficient bias
hist(Full_Bayesian_bias[seq(2,2*time,2)])
# Intercept bias
hist(Full_Bayesian_bias[seq(1,2*time,2)])

################# Multiple Cut Bias ############################

Multiple_Cut_bias <- rep(0,2*time)

Multiple_Cut_variance <- rep(0,2*time)

# Coefficient bias
for(i in 2:time){
  Multiple_Cut_bias[2*i] <- ((num_rep+1)*sum(X1[i,]*X2[i,])-sum(X1[i,])*sum(X2[i,]))/((num_rep+1)*(sum(X1[i,]*X1[i,])+1)-(sum(X1[i,]))^2)
  Multiple_Cut_variance[2*i] <- (num_rep+1)/((num_rep+1)*(sum(X1[i,]*X1[i,])+1)-(sum(X1[i,]))^2)
}

hist(Multiple_Cut_bias[seq(2,2*time,2)])

# intercept bias
for(i in 2:time){
  Multiple_Cut_bias[2*i-1] <- (sum(X1[i,]*X1[i,])*sum(X2[i,])+sum(X2[i,])-sum(X1[i,])*sum(X1[i,]*X2[i,]))/((num_rep+1)*(sum(X1[i,]*X1[i,])+1)-(sum(X1[i,]))^2)
  Multiple_Cut_variance[2*i-1] <- (sum(X1[i,]*X1[i,])+1)/((num_rep+1)*(sum(X1[i,]*X1[i,])+1)-(sum(X1[i,]))^2)
}

hist(Multiple_Cut_bias[seq(1,2*time,2)])



#################### Model ##################################

simulated_data<-simulation(time,num_rep,init_N,rate,beta,sigma,X)
observation<- simulated_data$observation
N_real<- simulated_data$N_real
#observation1<-simulation(time,init_N,rate,beta,sigma_model)
plot(observation[,1])



#Bayesian(Full)
Lik<-function(observation,init_N,inter,rate){
  fun_mean <- function(x){
    return(sum(x*c(init_N,inter[1],rate[1])))
  }
  result <- sum(dnorm(observation[1,], mean = apply(X[1,,],MARGIN = 1,fun_mean),sd=2,log=T))
  for(k in 1:(time-1)){
    N_est <- beta*rate[k]
    fun_mean <- function(x){
      return(sum(x*c(N_est,inter[k+1],rate[k+1])))
    }
    result <- result + sum(dnorm(observation[k+1,], mean = apply(X[k+1,,],MARGIN = 1,fun_mean),sd=2,log=T))
  }
  
  result <- result + sum(dnorm(rate,0,1,log = T))
  
  return(result)
}

rate_old<-rep(0.5,time)

inter_old <- rep(0.5,time)

N<-10000
pro_sd <- 0.1
pro_sd2 <- 0.1

Lik_old <- Lik(observation,init_N,inter_old,rate_old)

store_rate <- rep(0,time*N)
dim(store_rate) <- c(time,N)

store_inter <- rep(0,time*N)
dim(store_inter) <- c(time,N)

for(i in 1:N){
  rate_new <- rnorm(length(rate_old), mean = rate_old, sd = pro_sd)
  inter_new <- rnorm(length(inter_old), mean = inter_old, sd = pro_sd2)
  Lik_new <- Lik(observation,init_N,inter_new,rate_new)
  acct_rate <- min(1, exp(Lik_new + sum(log(dnorm(rate_old, mean=rate_new,sd = pro_sd)))  +  sum(log(dnorm(inter_old, mean=inter_new,sd = pro_sd2)))- Lik_old -sum(log(dnorm(rate_new, mean=rate_old,sd = pro_sd))) -  sum(log(dnorm(inter_new, mean=inter_old,sd = pro_sd2)))))
  random_value <- runif(1)
  if(random_value < acct_rate){
    rate_old <- rate_new
    inter_old <- inter_new
    Lik_old <- Lik_new
  }else{
    rate_old <- rate_old
    inter_old <- inter_old
    Lik_old <- Lik_old
  }
  store_rate[,i] <- rate_old
  store_inter[,i] <- inter_old
  print(i)
  if(i%%100==0){
    plot(store_rate[1,1:i])
  }
}

rate_est2<-apply(store_rate[1:100,5000:10000],1,mean)
sd2<-apply(store_rate,1,sd)

inter_est2<-apply(store_inter,1,mean)
sdd2<-apply(store_inter,1,sd)
#pro_sd <- sd2/4
#pro_sd2 <- sdd2/4
boxplot(rate-rate_est2)

#Bayesian(Cut)
rate_old<-rep(0.5,time)
inter_old<-rep(0.5,time)
N<-10000
pro_sd <- 0.1



Lik_old<-rep(0,time)


store_rate <- rep(0,time*N)
dim(store_rate) <- c(time,N)

store_inter <- rep(0,time*N)
dim(store_inter) <- c(time,N)

for(i in 1:N){
  rate_new <-  rnorm(length(rate_old), mean = rate_old, sd = pro_sd)
  inter_new <-  rnorm(length(inter_old), mean = inter_old, sd = pro_sd)
  for(k in 1:time){
    if(k==1){
      fun_mean <- function(x){
        return(sum(x*c(init_N,inter_old[1],rate_old[1])))
      }
      Lik_old[k] <- sum(dnorm(observation[1,], mean = apply(X[1,,],MARGIN = 1,fun_mean),sd=2,log=T)) + dnorm(rate_old[k],0,1,log = T)
      fun_mean_new <- function(x){
        return(sum(x*c(init_N,inter_new[k],rate_new[k])))
      }
      Lik_new <- sum(dnorm(observation[1,], mean = apply(X[1,,],MARGIN = 1,fun_mean_new),sd=2,log=T)) + dnorm(rate_new[k],0,1,log = T)
      acct_rate <- min(1, exp(Lik_new + dnorm(rate_old[k], mean=rate_new[k],sd = pro_sd,log = T) + dnorm(inter_old[k], mean=inter_new[k],sd = pro_sd,log = T) - Lik_old[k] - dnorm(rate_new[k], mean=rate_old[k],sd = pro_sd,log = T)- dnorm(inter_new[k], mean=inter_old[k],sd = pro_sd,log = T)))
      random_value <- runif(1)
      rate_old[k] <- rate_new[k]*sign(random_value < acct_rate) + rate_old[k]*sign(random_value >= acct_rate)
      inter_old[k] <- inter_new[k]*sign(random_value < acct_rate) + inter_old[k]*sign(random_value >= acct_rate)
      Lik_old[k] <- Lik_new*sign(random_value < acct_rate) + Lik_old[k]*sign(random_value >= acct_rate)
    }else{
      N_est <- beta*rate_old[k-1]           #conditional on previous sample 'rate_old[k-1]'
      fun_mean <- function(x){
        return(sum(x*c(N_est,inter_old[k],rate_old[k])))
      }
      Lik_old[k] <- sum(dnorm(observation[k,], mean = apply(X[k,,],MARGIN = 1,fun_mean),sd=2,log=T)) + dnorm(rate_old[k],0,1,log = T)
      fun_mean_new <- function(x){
        return(sum(x*c(N_est,inter_new[k],rate_new[k])))
      }
      Lik_new <- sum(dnorm(observation[k,], mean = apply(X[k,,],MARGIN = 1,fun_mean_new),sd=2,log=T)) + dnorm(rate_new[k],0,1,log = T)
      acct_rate <- min(1, exp(Lik_new + dnorm(rate_old[k], mean=rate_new[k],sd = pro_sd,log = T) + dnorm(inter_old[k], mean=inter_new[k],sd = pro_sd,log = T)  - Lik_old[k] - dnorm(rate_new[k], mean=rate_old[k],sd = pro_sd,log = T) - dnorm(inter_new[k], mean=inter_old[k],sd = pro_sd,log = T)))
      random_value <- runif(1)
      rate_old[k] <- rate_new[k]*sign(random_value < acct_rate) + rate_old[k]*sign(random_value >= acct_rate)
      inter_old[k] <- inter_new[k]*sign(random_value < acct_rate) + inter_old[k]*sign(random_value >= acct_rate)
      Lik_old[k] <- Lik_new*sign(random_value < acct_rate) + Lik_old[k]*sign(random_value >= acct_rate)
    }
  }
  store_rate[,i] <- rate_old
  store_inter[,i] <- inter_old
  print(i)
  if(i%%100==0){
    plot(store_rate[1,1:i])
  }
}

rate_est3<-apply(store_rate[1:100,5000:10000],1,mean)
sd3<-apply(store_rate,1,sd)
inter_est3<-apply(store_inter,1,mean)
boxplot(rate-rate_est2,rate-rate_est3)

result<-data.frame(true_value=rep(rate,2),estimates=c(rate_est2,rate_est3),estimates_bias=(c(rate_est2,rate_est3)-rep(rate,2))/c(sd2,sd3),method=factor(c(rep('Full Bayesian',time),rep('Multiple Cut',time))))

p <- ggplot(result, aes(true_value, estimates)) + geom_point(aes(colour = method),size = 2) +
  geom_abline(intercept = 0, slope = 1)

p <- ggplot(result, aes(x=method, y=estimates_bias,fill=method)) + 
  geom_boxplot() + geom_hline(yintercept=0, linetype="dashed", color = "red")




