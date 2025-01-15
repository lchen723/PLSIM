starting_time <- Sys.time()

method <- "simexboost"  # c("naive", "simexboost", "simex")
mmm <- "real"  # c("M1", "M2")
distri <- "data"  # c("normal", "unif")
sigmae_value <- 0.50 # c(0.15, 0.35, 0.50)
#n <- 800 # c(400, 800)


library(mvtnorm)
library(numDeriv)
library(MASS)
options(scipen = 999)

G_true = NULL
esti_A_all = NULL
esti_B_all = NULL
esti_G_all = NULL
MSE_all = NULL


#setwd("/Users/jcwu/Documents/paper/SIMEX real data analysis a¡E¡¦c??/a¡E¡¦c??a??c??e3?a¡V?")
setwd("C:/Users/User/Documents/RJ") # if remote PC
data <- read.csv("all_diff_end.csv")
data <- data[,-1]
data <- scale(data)
# View(data)

y_col <- which(colnames(data)=="Points")
x_col <- c(which(colnames(data)=="Field_goal_2pt_attempts"),
           which(colnames(data)=="Field_goal_2pt_pct"),
           which(colnames(data)=="Field_goal_3pt_attempts"), 
           which(colnames(data)=="Field_goal_3pt_pct"), 
           which(colnames(data)=="Free_throw_pct")
)
z_col <- c(which(colnames(data)=="Minutes_played_OT"),
           which(colnames(data)=="Offensive_rebounds"),
           which(colnames(data)=="Defensive_rebounds"), 
           which(colnames(data)=="Assists"), 
           which(colnames(data)=="Blocks"), 
           which(colnames(data)=="Turnovers"), 
           which(colnames(data)=="Personal_fouls")
) 

Y <- data[,y_col]
X <- data[,x_col]	
Z <- data[,z_col]

X_vars_px<-length(x_col) 
Z_vars_pz<-length(z_col)
tot_vars_p<-X_vars_px+Z_vars_pz
n <- length(Y)

#### (1-4) create sigmae ####
# sigmae <- diag(sigmae_value,tot_vars_p)
sigmae <- matrix(rep(NA,tot_vars_p^2),ncol=tot_vars_p)
for(i in 1:tot_vars_p){
  for(j in 1:tot_vars_p){
    sigmae[i,j] <- sigmae_value^(abs(i - j) + 1)
  }
}
#### Part 2:Estimation ####
#### Stage 1: Simulation ####
#### (2-1-1) create K, M, xi_M ####
K<-500        # given positive integer (can change)
M<-7          # positive integer (can change)
xi_M<-1       # prespecified positive number (can change)
zi<-seq(0, xi_M, length = M+1)
zi = zi[-(M+1)]

#### Stage 2: Estimation ####
#### (2-2-1) define # of knots ####
# define J(# of) knots (J can be chosen between 5 to 10)
J=5
## compute knots by quantiles change by diff. eta_star
#### (2-2-2) create alpha=(1,0,1,....,0,0) ####
m<-3
lambda<-1 # a penalty parameter
#### (2-2-3) q_star (all k, zi ???n?p???X??) ####
P <- diag(c(rep(0.1,m),rep(1,J))) # penalty on last J 
# positive semidefinite matrix


#### raw data content ####
XZ_star<-rbind(t(X),t(Z))

#### (20220919 add) working data has i, k, zeta ####
#### (2-1-2) create working data ####
V<-array(NA,dim=c(K,n,tot_vars_p)) # k,i,p
for(k in 1:K){
  V[k,,]<-t(rmvnorm(n=n, mean=rep(0,tot_vars_p), sigma = sigmae))
}


###########################################################################################


#### bootstrap ####
tt <- 100
for(i in 1:tt){
  Y <- data[,y_col]
  X <- data[,x_col]	
  Z <- data[,z_col]
  ### bootstrap ###
  bootstrap_id <- sample(1:n, n, replace=T)
  Y <- Y[bootstrap_id]
  X <- X[bootstrap_id, ]
  Z <- Z[bootstrap_id, ]
  ### bootstrap ###
  
  XZ_star<-rbind(t(X),t(Z))
  
  #### (20220919 add) working data has i, k, zeta ####
  #### (2-1-2) create working data ####
  V<-array(NA,dim=c(K,n,tot_vars_p)) # k,i,p
  for(k in 1:K){
    V[k,,]<-t(rmvnorm(n=n, mean=rep(0,tot_vars_p), sigma = sigmae))
  }
  
  ##################
  #### boosting ####
  ##################
  
  #    alpha_boost <- matrix(rep(1,m+J),ncol=1)
  beta_boost <- matrix(rep(0,X_vars_px),ncol=1)
  gamma_boost <- matrix(rep(0,Z_vars_pz),ncol=1)
  
  alpha_zeta = NULL
  beta_zeta = NULL
  gamma_zeta = NULL
  G_zeta = NULL
  
  for(zi_m in 0:(M-1)) {
    alpha_boost_all <- NULL 
    beta_boost_all <- NULL
    gamma_boost_all <- NULL
    G_i_k_zi <- NULL
    for(k in 1:K) {
      beta_boost <- matrix(rep(0,X_vars_px),ncol=1)
      gamma_boost <- matrix(rep(0,Z_vars_pz),ncol=1)
      working_data = XZ_star + sqrt(zi[zi_m+1]) * t(rmvnorm(n=n, mean=rep(0,tot_vars_p), sigma = sigmae))
      
      for(t in 1:40) {
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
        
        #### alpha boost ?e?b ####
        sum_alpha_first <- matrix(rep(0,(m+J)^2),nrow=(m+J))
        for(i in 1:n){
          x_i_star_k_zi <- matrix(working_data[1:X_vars_px,i],ncol=1)
          z_i_star_k_zi <- matrix(working_data[-c(1:X_vars_px),i],ncol=1)
          eta_star<-t(z_i_star_k_zi)%*%matrix(gamma_boost,ncol=1)
          alpha_first <- B_eta_star(eta_star)%*%t(B_eta_star(eta_star))
          sum_alpha_first <- sum_alpha_first + alpha_first
        }
        xxx <- sum_alpha_first  - lambda * P
        solve(xxx) ### Error in solve.default(xxx) : ?t?I-p???W?O?W?S??: ???I¡Ó?????=7.76664e-18 ###  ### ?Y¡MN¡LM ##
        
        
        #### alpha boost ???b ####
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
        #      xi_beta <- 0.0001   ### ?U¢D??A3]?w
        xi_beta <- 0.0003   ### ?U¢D??A3]?w
        #      if(length(which((beta_candidate<0))) == 5){
        #      id_beta = which(abs(beta_candidate) >= .9*max(abs(beta_candidate)))} 
        #      if(length(which((beta_candidate<0))) != 5){
        id_beta = which((beta_candidate) >= 0.1*max((beta_candidate))) # }
        beta_boost[id_beta] = beta_boost[id_beta] + xi_beta * abs(beta_candidate[id_beta])
        #beta_boost[which(beta_boost < 0.2)] = 0
        
        #### gamma boost ####    
        sum_gamma_first <- matrix(rep(0,Z_vars_pz),nrow=Z_vars_pz)
        for(i in 1:n){
          z_i_star_k_zi <- matrix(working_data[-c(1:X_vars_px),i],ncol=1)
          eta_star<-t(z_i_star_k_zi)%*%matrix(gamma_boost,ncol=1)
          gamma_first <- 2*z_i_star_k_zi %*% t(jacobian(B_eta_star,eta_star))%*%alpha_boost * as.numeric((Y[i]-t(x_i_star_k_zi) %*% beta_boost - t(B_eta_star(eta_star)) %*% alpha_boost))
          sum_gamma_first <- sum_gamma_first + gamma_first
        }
        
        
        #### gamma boost ####
        gamma_candidate <- sum_gamma_first
        #      xi_gamma <- 0.0001    ### ?U¢D??A3]?w
        xi_gamma <- 0.0003    ### ?U¢D??A3]?w
        if(length(which((gamma_candidate<0))) == 6){
          id_gamma = which(abs(gamma_candidate) >= 0.1*max(abs(gamma_candidate)))}      
        if(length(which((gamma_candidate<0))) != 6){
          id_gamma = which(abs(gamma_candidate) >= 0.1*max(abs(gamma_candidate))) }
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
      
      alpha_boost_all <- cbind(alpha_boost_all,alpha_boost) 
      beta_boost_all <- cbind(beta_boost_all,beta_boost)
      gamma_boost_all <- cbind(gamma_boost_all,gamma_boost)
      
    }
    
    alpha_zeta = rbind(alpha_zeta,rowMeans(alpha_boost_all))
    beta_zeta = rbind(beta_zeta,rowMeans(beta_boost_all))
    gamma_zeta = rbind(gamma_zeta,((rowMeans(gamma_boost_all) / sqrt(sum(rowMeans(gamma_boost_all)^2))) ))
    #gamma_zeta[which(abs(gamma_zeta)<0.1)] = 0
  }
  
  dB = dim(beta_zeta)[1]
  dG = dim(gamma_zeta)[1]
  
  #ZZ = cbind(1,zi)
  ZZ = cbind(1,zi,zi^2)
  alpha_zeta[which(is.na(alpha_zeta))] = 0
  A = (solve(t(ZZ)%*%ZZ) %*% t(ZZ))[,1:dB] %*% (alpha_zeta[1:dB,])
  B = (solve(t(ZZ)%*%ZZ) %*% t(ZZ))[,1:dB] %*% (beta_zeta[1:dB,])
  G = (solve(t(ZZ)%*%ZZ) %*% t(ZZ))[,1:dG] %*% (gamma_zeta[1:dB,])
  esti_A = A[1,] - A[2,] + A[3,]
  esti_B = colSums(abs(B))
  esti_B[which(abs(esti_B)<0.1)]=0
  esti_G = colSums((G)) / sqrt( sum(colSums((G))^2)  )
  esti_G[which(abs(esti_G)<0.1)]=0 
  
  esti_A_all = rbind(esti_A_all, sort(esti_A)) 
  esti_B_all = rbind(esti_B_all, esti_B) 
  esti_G_all = rbind(esti_G_all, esti_G) 
  
  #### MSE ####
  MSE <- mean ( ( Y - ( X %*% beta_boost + alpha_boost ) )^2 )
  MSE_all <- cbind( MSE_all, MSE )
}  


colMeans(esti_A_all) -> estA
colMeans(esti_B_all) -> estB
colMeans(esti_G_all) -> estG

meanMSE <- mean(MSE_all, na.rm=TRUE)
mean(MSE_all)

#### save .RData & boxplot ####
#setwd("/Users/jcwu/Documents/paper/SIMEX invited revise/save RData")
setwd("C:/Users/User/Documents/RJ") # if remote PC

save.image(file=paste(method, "_", mmm, "_", distri, "_", sigmae_value, "_", n, ".RData", sep=""))

#### store the boxplot ####
png(paste(method, "_", mmm, "_", distri, "_", sigmae_value, ".png", sep=""), width = 800, height = 600) # open png device 
boxplot(as.numeric(MSE_all)) # plot
dev.off() # close the device 


ending_time <- Sys.time()
(ending_time - starting_time)/10
