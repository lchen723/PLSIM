starting_time <- Sys.time()

method <- "naive"  # c("naive", "simexboost", "simex")
mmm <- "M2"  # c("M1", "M2")
distri <- "unif"  # c("normal", "unif")
sigmae_value <- 0.5 # c(0.15, 0.35, 0.50)
n <- 400 # c(400, 800)

library(mvtnorm)
library(numDeriv)
library(MASS)
options(scipen = 999)

G_true = NULL
esti_A_all = NULL
esti_B_all = NULL
esti_G_all = NULL
MSE_all = NULL

for(i in 1:100){
  #### Part 1: Data Generation ####
  #### (1-1) generate X, Z ####
  #n<-400
  X_vars_px<-5; Z_vars_pz<-6;      # Change px, pz  here!!! #
  tot_vars_p<-X_vars_px+Z_vars_pz  # generate X,Z together 
  # total:[(px+pz)*n]
  
  mu_X = rep(0,tot_vars_p)
  Sigma_X = diag(1,tot_vars_p)
  #total_vars = mvrnorm(n, mu_X, Sigma_X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  total_vars = matrix(runif(tot_vars_p*n, -1, 1),nrow=n)
  total_vars = t(total_vars)
  
  # separate X,Z by sampling from total_vars
  X_id<-sample(1:tot_vars_p,X_vars_px)
  X<-total_vars[1:X_vars_px,]  # X:(px*n)
  Z<-total_vars[(1+X_vars_px):tot_vars_p,] # Z:(pz*n)
  X <- scale(X); Z <- scale(Z)
  
  #### (1-2) specify beta, gamma, function g, epsilon~N(0,1) ####
  
  beta <- c(rep(0,2),rep(0,X_vars_px-2))
  gamma <- c(1,rep(0,Z_vars_pz-1))
  gamma <- gamma / sqrt(sum(gamma^2))
  epsilon<-matrix(rnorm(n,0,1),ncol=1)
  
  #### (1-3) generate Y from (2) ####
  Y<-t(X)%*%beta+sin(t(Z)%*%gamma)+epsilon
  Y <- scale(Y)
  
  #### (1-4) create sigmae ####
  # sigmae <- diag(sigmae_value,tot_vars_p)
  sigmae <- matrix(rep(NA,tot_vars_p^2),ncol=tot_vars_p)
  for(i in 1:tot_vars_p){
    for(j in 1:tot_vars_p){
      sigmae[i,j] <- sigmae_value^(abs(i - j) + 1)
    }
  }
  
  #### (1-5) generate e (p*n matrix) from sigmae  ####
  e<-t(rmvnorm(n=n, mean=rep(0,tot_vars_p), sigma = sigmae))
  
  #### (1-6) create (X_star,Z_star) by (X,Z)+e ####
  XZ_star<-rbind(X,Z)+e
  
  
  #### Part 2:Estimation ####
  
  #### Stage 1: Simulation ####
  
  #### Stage 2: Estimation ####
  
  #### (2-2-1) define # of knots ####
  # define J(# of) knots (J can be chosen between 5 to 10)
  J=5
  ## compute knots by quantiles change by diff. eta_star
  
  
  #### (2-2-2) create alpha=(1,0,1,....,0,0) ####
  m<-3
  lambda<-1 # a penalty parameter
  
  #### (2-2-3) q_star (all k, zi 3?FG-n-poa?FDX!LO) ####
  P <- diag(c(rep(0.1,m),rep(1,J))) # penalty on last J 
  # positive semidefinite matrix
  
  ##################
  #### boosting ####
  ##################
  
  #    alpha_boost <- matrix(rep(1,m+J),ncol=1)
  beta_boost <- matrix(rep(0,X_vars_px),ncol=1)
  gamma_boost <- matrix(rep(0,Z_vars_pz),ncol=1)
  
  alpha_zeta = NULL
  beta_zeta = NULL
  gamma_zeta = NULL
  
  
  beta_boost <- matrix(rep(0,X_vars_px),ncol=1)
  gamma_boost <- matrix(rep(0,Z_vars_pz),ncol=1)
  working_data = XZ_star
  
  for(t in 1:70) {
    eta_star_n<-t(working_data[-c(1:X_vars_px),])%*%matrix(gamma_boost,ncol=1)
    # rho<-rep(NA,J) # compute knots by quantiles (by diff. eta_star)
    # for (j in 1:J){
    #   rho[j]<-quantile(eta_star_n,(1/(J+1))*j)
    # }
    rho <- seq(0, 1, length.out = J+1)
    
    B_eta_star <- function(eta_star){
      box <- matrix(rep(NA,(m+J)),ncol=(m+J))
      for(i in 1:m){
        box[,i]<-eta_star^(i)
      }
      for(j in 1:J){
        box[,m+j]<-(eta_star-rho[j])^m
      }
      return(t(box))
    }
    
    #### alpha boost  ####
    sum_alpha_first <- matrix(rep(0,(m+J)^2),nrow=(m+J))
    for(i in 1:n){
      x_i_star_k_zi <- matrix(working_data[1:X_vars_px,i],ncol=1)
      z_i_star_k_zi <- matrix(working_data[-c(1:X_vars_px),i],ncol=1)
      eta_star<-t(z_i_star_k_zi)%*%matrix(gamma_boost,ncol=1)
      alpha_first <- B_eta_star(eta_star)%*%t(B_eta_star(eta_star))
      sum_alpha_first <- sum_alpha_first + alpha_first
    }
    xxx <- sum_alpha_first  - lambda * P
    solve(xxx) 
    
    
    #### alpha boost ####
    sum_alpha_last <- matrix(rep(0,(m+J)),nrow=(m+J))
    for(i in 1:n){
      x_i_star_k_zi <- matrix(working_data[1:X_vars_px,i],ncol=1)
      z_i_star_k_zi <- matrix(working_data[-c(1:X_vars_px),i],ncol=1)
      eta_star<-t(z_i_star_k_zi)%*%matrix(gamma_boost,ncol=1)
      alpha_last <- B_eta_star(eta_star) * as.numeric((Y[i]-t(x_i_star_k_zi)%*%beta_boost))
      sum_alpha_last <- sum_alpha_last + alpha_last
    }
    
    
    #### alpha boost ####
    alpha_boost <- matrix((solve(xxx) %*% sum_alpha_last),ncol=1)
    
    
    #### beta boost  ####
    sum_beta_last <- matrix(rep(0,X_vars_px),nrow=X_vars_px)
    for(i in 1:n){
      x_i_star_k_zi <- matrix(working_data[1:X_vars_px,i],ncol=1)
      z_i_star_k_zi <- matrix(working_data[-c(1:X_vars_px),i],ncol=1)
      eta_star<-t(z_i_star_k_zi)%*%matrix(gamma_boost,ncol=1)
      beta_last <- 2*x_i_star_k_zi %*% as.numeric((Y[i]-t(x_i_star_k_zi) %*% beta_boost - t(B_eta_star(eta_star)) %*% alpha_boost))
      sum_beta_last <- sum_beta_last + beta_last
    }
    
    
    #### beta boost ####
    beta_candidate <- sum_beta_last
    #      xi_beta <- 0.0001   ### !??U?FDy?A3]cw
    xi_beta <- 0.00058   ### !??U?FDy?A3]cw
    #      if(length(which((beta_candidate<0))) == 5){
    #      id_beta = which(abs(beta_candidate) >= .9*max(abs(beta_candidate)))} 
    #      if(length(which((beta_candidate<0))) != 5){
    id_beta = which((beta_candidate) >= 1*max((beta_candidate))) # }
    beta_boost[id_beta] = beta_boost[id_beta] + xi_beta * abs(beta_candidate[id_beta])
    #beta_boost[which(beta_boost < 0.2)] = 0
    
    #### gamma boost ?e?FDb ####    # ?FX?!PL?A?Oscaler
    sum_gamma_first <- matrix(rep(0,Z_vars_pz),nrow=Z_vars_pz)
    for(i in 1:n){
      z_i_star_k_zi <- matrix(working_data[-c(1:X_vars_px),i],ncol=1)
      eta_star<-t(z_i_star_k_zi)%*%matrix(gamma_boost,ncol=1)
      gamma_first <- 2*z_i_star_k_zi %*% t(jacobian(B_eta_star,eta_star))%*%alpha_boost * as.numeric((Y[i]-t(x_i_star_k_zi) %*% beta_boost - t(B_eta_star(eta_star)) %*% alpha_boost))
      sum_gamma_first <- sum_gamma_first + gamma_first
    }
    
    
    #### gamma boost ####
    gamma_candidate <- sum_gamma_first
    #      xi_gamma <- 0.0001    ### !??U?FDy?A3]cw
    xi_gamma <- 0.0004    ### !??U?FDy?A3]cw
    if(length(which((gamma_candidate<0))) == 6){
      id_gamma = which((gamma_candidate) >= 1*max((gamma_candidate)))}      
    if(length(which((gamma_candidate<0))) != 6){
      id_gamma = which((gamma_candidate) >= 1*max((gamma_candidate))) }
    gamma_boost[id_gamma] = gamma_boost[id_gamma] + xi_gamma * (abs(gamma_candidate[id_gamma]))
    gamma_boost <- gamma_boost / sqrt(sum(gamma_boost^2))
    
    
    
  }
  
  gamma_boost <- gamma_boost / sqrt(sum(gamma_boost^2))
  #alpha_boost
  #beta_boost
  #gamma_boost
  
  #eta_n<-t(Z)%*%matrix(gamma,ncol=1)
  eta_n<-t(working_data[-c(1:X_vars_px),])%*%matrix(gamma_boost,ncol=1)
  f_variables <- smooth.spline(x=eta_n,y=Y,cv=TRUE,all.knots=c(0,0.2,0.4,0.6,0.8,1))
  
  alpha_boost = fitted(f_variables)
  
  
  
  
  
  esti_A_all = rbind(esti_A_all, alpha_boost) 
  esti_B_all = rbind(esti_B_all, as.vector(beta_boost)) 
  esti_G_all = rbind(esti_G_all, as.vector(gamma_boost)) 
  G_true = rbind(G_true, sort(sin(t(Z)%*%gamma)))
  
  #### MSE ####
  MSE <- mean ( ( Y - ( t(X) %*% beta_boost + alpha_boost ) )^2 )
  MSE_all <- cbind( MSE_all, MSE )
}



colMeans(esti_A_all) -> estA
colMeans(esti_B_all) -> estB
colMeans(esti_G_all) -> estG

plot(colMeans(G_true))
points(sort(colMeans(esti_A_all)),col=2)
points(sort(estA),col=3)


######################################################################################

#### M2 version #### 

L1beta <- sum(abs(beta - estB)) # L1 norm of beta
L2beta <- sum( (beta - estB)^2 ) # L2 norm of beta
spe_beta <- sum(which(ifelse(estB < 0.1, 0, estB)==0)>0)/(X_vars_px-0) # specificity of beta M1
sen_beta <- "NULL" # sensitivity of beta M1

L1gamma <- sum( abs(gamma - estG) ) # L1 norm of gamma
L2gamma<- sum( (gamma - estG)^2 ) # L2 norm of gamma
spe_gamma <- sum(which(ifelse(estG < 0.1, 0, estG)==0)>1)/(Z_vars_pz-1) # specificity of gamma M1
sen_gamma <- sum(which(ifelse(estG < 0.1, 0, estG)!=0)<=1)/1# sensitivity of gamma M1

Inf_func <- max( abs(estA - sin(t(Z)%*%gamma)), na.rm=TRUE) # infinity norm of function g
L1_func <- sum( abs( estA - sin(t(Z)%*%gamma) ), na.rm=TRUE) # L1 norm of function g

meanMSE <- mean(MSE_all, na.rm=TRUE)

output <- data.frame(matrix(c(L1beta, L2beta, spe_beta, sen_beta, L1gamma, L2gamma, spe_gamma, sen_gamma, Inf_func, L1_func, meanMSE),nrow=1))
colnames(output) <- c("L1beta", "L2beta", "spe_beta", "sen_beta", "L1gamma", "L2gamma", "spe_gamma", "sen_gamma", "Inf_func", "L1_func", "meanMSE")
output

######################################################################################

#### save .RData & boxplot ####
#setwd("/Users/jcwu/Documents/paper/SIMEX invited revise/save RData")
setwd("C:/Users/User/Documents/RJ") # if remote PC

write.csv(output, paste(method, "_", mmm, "_", distri, "_", sigmae_value, "_", n, ".csv", sep=""))
save.image(file=paste(method, "_", mmm, "_", distri, "_", sigmae_value, "_", n, ".RData", sep=""))

#### store the boxplot ####
png(paste(method, "_", mmm, "_", distri, "_", sigmae_value, "_", n, ".png", sep=""), width = 800, height = 600) # open png device 
boxplot(as.numeric(MSE_all)) # plot
dev.off() # close the device 

ending_time <- Sys.time()
(ending_time - starting_time)/10