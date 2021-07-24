library(evd)
library(Matrix)
library(tidyverse)
library(knitr)
library(INLA)
library(mapdata)
library(data.table)
library(rgdal)

# A few helper functions
source("helperFun.R")

#### Get the data ready
# Read in the data
# Get catchment descriptors and flow data
load("catchmentDesc.RData")  # catchmentDesc
load("flowDataSample.RData") # flowData

# Number of stations
n_st <- nrow(catchmentDesc)
stations <- unique(catchmentDesc$Station)

###################### ML step ############################################
# Constants
xi00 <- 0
alp3 <- 0.8 # c_phi
alpha <- 4
beta <- 4
sig3 = (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
b3 <- -sig3*(log(-log(1-0.5^alp3))) # a_phi
# Likelihood function to minimize
fn <- function(theta,y,t) {
  n <- length(t)
  t_e <- 1975
  delta_0 <- 1/125
  gamma <- theta[4]
  log_delta_i <- log(2*delta_0) + 2*gamma/delta_0 - log(1 + exp(2*gamma/delta_0))
  delta_i <- exp(log_delta_i)-delta_0
  mu_it <- (1 + delta_i*(t-t_e))*exp(theta[1])
  
  sigma_i <- exp(theta[2] + theta[1])
  if(sigma_i<=0){
    return(100000)
  }
  sig3 <- (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1))
  b3 <- -sig3*(log(-log(1-0.5^alp3)))
  xitheta <- (1 - exp(-exp((theta[3]-b3)/sig3)))^(1/alp3) - 0.5

  sum_loglik <- 0
  for(i in 1:n){
    sum_loglik <- sum_loglik + dgev(y[i], loc = mu_it[i], scale = sigma_i, shape = xitheta, log = TRUE)
  }
  res = -sum_loglik - 
    ((alpha - alp3)*log(xitheta + 0.5) + (beta-1)*log(0.5 - xitheta) + (theta[3]-b3)/sig3 - exp((theta[3]-b3)/sig3)  )  + 
    0.5*gamma^2/((0.5*delta_0)^2)
  return(res)
}

####################### ML step ##############################################
dt_mles <- data.frame()

for(i in 1:n_st){
  station <- stations[i]
  currentS <- flowData %>% filter(Station == station)
  # This is done to get some initial guess at the location, scale and shape parameters
  GEV_fit <- fgev(currentS$Flow)
  mu0 <- GEV_fit$estimate[1]
  sigma0 <- GEV_fit$estimate[2]
  xi0 <- GEV_fit$estimate[3]
  if (xi0 > 0) {
    xi0 <- min(xi0,0.45)
  } else {
    xi0 <- max(xi0,-0.45)
  }
  theta0 <- c(log(mu0),log(sigma0)-log(mu0),b3 + sig3*log(-log(1 - (xi0+0.5)^alp3)), 0)
  GEV_fit <- nlm(fn, theta <- theta0, y <- currentS$Flow, t <- currentS$year, hessian = T)
  
  S_d <- try(solve(GEV_fit$hessian),silent=T)
  
  if(i %% 10 == 0) print(i)
  
  res <- data.frame(
    Station = station,
    psi = GEV_fit$estimate[1], 
    tau = GEV_fit$estimate[2], 
    kappa = GEV_fit$estimate[3], 
    gamma = GEV_fit$estimate[4], 
    v_p = S_d[1,1],
    v_p_t = S_d[1,2],
    v_p_k = S_d[1,3],
    v_p_g = S_d[1,4],
    v_t = S_d[2,2],
    v_t_k = S_d[2,3],
    v_t_g = S_d[2,4],
    v_k = S_d[3,3],
    v_k_g = S_d[3,4],
    v_g = S_d[4,4]
  )
  dt_mles <- rbind(res, dt_mles)
}

################################## Make covariates ready #####################################################################
names <- c(
  "Int.",
  "log(Area)",
  "log(SAAR)",
  "log(FARL)",
  "BFIHOST^2",
  "log(FPEXT)",
  "log(URBEXT+1)",
  "log(DPLBAR)",
  "log(DPSBAR)",
  "log(LDP)",
  "log(SPRHOST)",
  "log(ASPBAR)",
  "log(ALTBAR)",
  "log(ASPVAR)",
  "log(PROPWET)")
ncols <- length(names)
covariates <- matrix(c(
  rep(1,n_st), # Intercept
  transform(catchmentDesc$Area, "AREA"),
  transform(catchmentDesc$SAAR, "SAAR"),
  transform(catchmentDesc$FARL, "FARL"),
  transform(catchmentDesc$BFIHOST, "BFIHOST"),
  transform(catchmentDesc$FPEXT, "FPEXT"),
  transform(catchmentDesc$URBEXT2000, "URBEXT"),
  transform(catchmentDesc$DPLBAR, "DPLBAR"),
  transform(catchmentDesc$DPSBAR, "DPSBAR"),
  transform(catchmentDesc$LDP, "LDP"),
  transform(catchmentDesc$SPRHOST, "SPRHOST"),
  transform(catchmentDesc$ASPBAR, "ASPBAR"),
  transform(catchmentDesc$ALTBAR, "ALTBAR"),
  transform(catchmentDesc$ASPVAR, "ASPVAR"),
  transform(catchmentDesc$PROPWET, "PROPWET")
), ncol = 15, nrow = n_st
)
colnames(covariates) <- names

################################ Make X matrices for psi, tau and xi #########################################################
cov_names_psi <- c("Int.","log(Area)", "log(SAAR)", "log(FARL)", "BFIHOST^2")
cov_names_tau <- c("Int.","log(Area)", "log(SAAR)", "log(FARL)", "log(FPEXT)", "log(URBEXT+1)")
cov_names_kappa <- c("Int.", "log(FPEXT)")
cov_names_gamma <- c("Int.", "log(PROPWET)")
X_psi <- covariates[, c("Int.","log(Area)", "log(SAAR)", "log(FARL)", "BFIHOST^2")]
X_tau <- covariates[, c("Int.","log(Area)", "log(SAAR)", "log(FARL)", "log(FPEXT)", "log(URBEXT+1)")]
X_kappa <- as.matrix(covariates[, c("Int.", "log(FPEXT)")])
X_gamma <- as.matrix(covariates[, c("Int.", "log(PROPWET)")])

############################### Make mesh and A matrices for spatial component ###############################################
coords <- cbind(catchmentDesc$long, catchmentDesc$lat)
mesh <- inla.mesh.2d(
  loc=coords,
  offset = 0.08,
  max.edge = 0.07,
  cutoff = 0.005)

# To plot mesh
# plot(mesh)

# Make A matrix - The same for psi and tau
A_mat <- inla.spde.make.A(mesh, loc=coords)
A_psi <- A_mat
A_tau <- A_mat

################# Fit INLA models, to have good starting proposal distribution ###############
# Add location to MLE data
dt_mles <- dt_mles %>% 
  left_join(catchmentDesc %>% select(Northing, Easting, Station))
mdl_sep_psi <- fit_model_psi(data = dt_mles, desc = X_psi)
mdl_sep_tau <- fit_model_tau(data = dt_mles, desc = X_tau)
mdl_sep_xi <- fit_model_kappa(data = dt_mles, desc = X_kappa)
mdl_sep_gamma <- fit_model_gamma(data = dt_mles, desc = X_gamma)

# The spde object is only used for getting the precision matrix for a given range and sd
# so priors do not matter here
d_spde <- inla.spde2.pcmatern(mesh, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

# Make the Z matrix
Z <- bdiag(cbind(X_psi,A_psi), cbind(X_tau, A_tau), X_kappa, X_gamma)

N_psi <- dim(X_psi)[2]
N_tau <- dim(X_tau)[2]
N_xi <- dim(X_kappa)[2]
N_gamma <- dim(X_gamma)[2]
N_colsA <- dim(A_psi)[2]

################################################ Helper functions ########################################################

get_mode_mean_sd_for_hyperparameter <- function(mdl, fun, hyperparameter){
  x0 <- mdl$marginals.hyperpar[[hyperparameter]]
  E <- inla.emarginal(function(x) c(fun(x), fun(x)^2), x0)
  sd <- sqrt(E[2] - E[1]^2)
  mean <- E[1]
  
  mt_tmp <- inla.tmarginal(fun, x0)
  dat_tmp <- inla.smarginal(mt_tmp)
  mode <- dat_tmp$x[which(dat_tmp$y == max(dat_tmp$y))]
  return(list(sd = sd, mean = mean, mode = mode))
}
log_f_x_given_theta <- function(kappa, Q_x){
  res <- determinant(Q_x, logarithm = T)
  res <- .5*res$modulus[1]*res$sign
  return(res)
}
log_f_x_given_eta_hat_and_theta <- function(kappa, Q_x_given_etaHat, N, Z, eta_hat, Q_etay, K){
  #Q_x_given_etaHat <- Q_x + (t(B) %*% Q_etay) %*% B
  #Q_x_given_etaHat <- Matrix(Q_x_given_etaHat, sparse = T)
  mu_x_given_etaHat <- solve(a = Cholesky(Q_x_given_etaHat),b = b)
  determ <- determinant(Q_x_given_etaHat, logarithm = T)
  res <- -.5*t(mu_x_given_etaHat)%*%(Q_x_given_etaHat%*%mu_x_given_etaHat)
  res <- res + .5*determ$modulus[1]*determ$sign
  return(res)
}
calc_Q_x_given_etaHat <- function(Q_x, N, Q_etay){
  Q_x_given_etaHat <- Q_x
  Q_x_given_etaHat[1:(4*N), 1:(4*N)] <- Q_x[1:(4*N), 1:(4*N)] + Q_etay
  return(Q_x_given_etaHat)
}
calc_Q_x <- function(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi,Q_beta_gamma, N, Z, d_spde){
  
  Q_u_psi <- makeQ_u(s = exp(kappa$u_psi), rho = exp(kappa$v_psi), d_spde)
  Q_u_tau <- makeQ_u(s = exp(kappa$u_tau), rho = exp(kappa$v_tau), d_spde)
  
  Q_nu <- bdiag(Q_beta_psi,Q_u_psi,Q_beta_tau,Q_u_tau,Q_beta_xi,Q_beta_gamma)
  
  Q_epsilon <- bdiag(Diagonal(N,exp(kappa$kappa_psi)), Diagonal(N,exp(kappa$kappa_tau)), 
                     Diagonal(N,exp(kappa$kappa_xi)),Diagonal(N,exp(kappa$kappa_gamma)))
  K <- dim(Q_nu)[1]
  Q_x[1:(4*N),  1:(4*N)] <- Q_epsilon
  Q_x[1:K + 4*N,1:(4*N)] <- -t(Z)%*%Q_epsilon
  Q_x[1:(4*N),  1:K + 4*N] <- - (Q_epsilon %*% Z)
  Q_x[1:K + 4*N, 1:K + 4*N] <- Q_nu + t(Z)%*%Q_epsilon%*%Z
  return(Q_x)
}

################################## Prior matrices for covariate coefficients ###################################

betaVarPrior <- 10000
K <- 2*N_colsA + N_psi +N_tau +N_xi + N_gamma
Q_beta_psi <- Diagonal(N_psi, 1/betaVarPrior)
Q_beta_tau <- Diagonal(N_tau, 1/betaVarPrior)
Q_beta_xi <- Diagonal(N_xi, 1/betaVarPrior)
Q_beta_gamma <- Diagonal(N_gamma, 1/betaVarPrior)

################# get variance and mode for hyper parameters for proposal distribution ###############################

tau_psi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_psi$mdl, function(x) log(x), "Precision for idx")
tau_tau_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_tau$mdl, function(x) log(x), "Precision for idx")
tau_xi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_xi$mdl, function(x) log(x), "Precision for idx")
tau_gamma_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_gamma$mdl, function(x) log(x), "Precision for idx")
rho_psi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_psi$mdl, function(x) log(x), "Range for s")
rho_tau_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_tau$mdl, function(x) log(x), "Range for s")
sigma_psi_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_psi$mdl, function(x) log(x), "Stdev for s")
sigma_tau_0 <- get_mode_mean_sd_for_hyperparameter(mdl_sep_tau$mdl, function(x) log(x), "Stdev for s")

# Number of stations
N <- nrow(dt_mles)
eta_hat <- matrix(c(dt_mles$psi, dt_mles$tau, dt_mles$kappa, dt_mles$gamma), nrow = 4 * N)

###### Make sigma_eta_y
sigma_eta_y <- makeSigma_etay_w_gamma(dt = dt_mles)
Q_etay <- solve(sigma_eta_y)

# Mean and covariance matrix of the proposal distribution in the first loop
kappa_0 <- c(tau_psi_0$mode, rho_psi_0$mode, sigma_psi_0$mode,
             tau_tau_0$mode, rho_tau_0$mode, sigma_tau_0$mode, 
             tau_xi_0$mode, tau_gamma_0$mode) %>% as.matrix()
Sigma_kappa_0 <- diag(c(tau_psi_0$sd^2*2, rho_psi_0$sd^2, sigma_psi_0$sd^2,
                        tau_tau_0$sd^2*2, rho_tau_0$sd^2, sigma_tau_0$sd^2, 
                        tau_xi_0$sd^2/2, tau_gamma_0$sd^2/2))*0.2
kappa_0 <- apply(chain,MARGIN =  1, mean) %>% as.matrix()
Sigma_kappa_0 <- (chain %>% t() %>% data.frame() %>% distinct() %>% cov())*.5
tmp <- read.table("../../../Downloads/grein_copy/chain_oos42.csv") %>% as.matrix()
Sigma_kappa_0 <- cov(t(tmp))*0.4
library(MASS)
library(SparseM)
#################################### Model fitting - Hyperparameters (One chain) ##################
N_samples <- 1000
Q_x <- Matrix(0,nrow = 4*N+K, ncol = 4*N+K)
B <- Matrix(0,nrow = 4*N, ncol = 4*N+K)
B[1:(4*N),1:(4*N)] <- diag(1,4*N)
b <- t(B)%*%(Q_etay%*%eta_hat)
kappa_k <- mvrnorm(mu = kappa_0, Sigma = Sigma_kappa_0)
kappa_mat <- matrix(NA, nrow = length(kappa_k), ncol = N_samples)
kappa_k_is_kappa_star <- F
accept <- 0
for(i in 1:N_samples){
  kappa_star <- mvrnorm(mu = kappa_k, Sigma = Sigma_kappa_0)
  k_st <- list(
    kappa_psi = kappa_star[1],
    v_psi = kappa_star[2],
    u_psi = kappa_star[3],
    kappa_tau = kappa_star[4],
    v_tau = kappa_star[5],
    u_tau = kappa_star[6],
    kappa_xi = kappa_star[7],
    kappa_gamma = kappa_star[8]
  )
  k_k <- list(
    kappa_psi = kappa_k[1],
    v_psi = kappa_k[2],
    u_psi = kappa_k[3],
    kappa_tau = kappa_k[4],
    v_tau = kappa_k[5],
    u_tau = kappa_k[6],
    kappa_xi = kappa_k[7],
    kappa_gamma = kappa_k[8]
  )
  if(i == 1){
    Q_x <- calc_Q_x(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi, Q_beta_gamma, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    l_x_etaHat_k_st <- log_f_x_given_eta_hat_and_theta(kappa = k_st, Q_x_given_etaHat = Q_x_given_etaHat, 
                                                       N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_st <- log_f_x_given_theta(kappa = k_st, Q_x)
    l_prior_k_st <- log_prior_logTheta(k_st)
    
    Q_x <- calc_Q_x(kappa = k_k, Q_beta_psi, Q_beta_tau, Q_beta_xi, Q_beta_gamma, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    
    l_x_etaHat_k_k <- log_f_x_given_eta_hat_and_theta(kappa = k_k, Q_x_given_etaHat = Q_x_given_etaHat,
                                                      N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_k <- log_f_x_given_theta(kappa = k_k, Q_x)
    l_prior_k_k <- log_prior_logTheta(kappa = k_k)
  }else if(kappa_k_is_kappa_star){
    Q_x <- calc_Q_x(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi, Q_beta_gamma, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    
    l_x_etaHat_k_k <- l_x_etaHat_k_st
    l_x_k_k <- l_x_k_st
    l_prior_k_k <- l_prior_k_st
    
    l_x_etaHat_k_st <- log_f_x_given_eta_hat_and_theta(kappa = k_st, Q_x_given_etaHat = Q_x_given_etaHat, 
                                                       N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_st <- log_f_x_given_theta(kappa = k_st, Q_x)
    l_prior_k_st <- log_prior_logTheta(k_st)
  }else{
    Q_x <- calc_Q_x(kappa = k_st, Q_beta_psi, Q_beta_tau, Q_beta_xi, Q_beta_gamma, N, Z, d_spde)
    Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x, N = N, Q_etay = Q_etay)
    
    l_x_etaHat_k_st <- log_f_x_given_eta_hat_and_theta(kappa = k_st, Q_x_given_etaHat = Q_x_given_etaHat,
                                                       N = N, Z = Z, eta_hat = eta_hat, Q_etay = Q_etay,K = K)
    l_x_k_st <- log_f_x_given_theta(kappa = k_st, Q_x)
    l_prior_k_st <- log_prior_logTheta(k_st)
  }
  
  r <- l_prior_k_st-l_prior_k_k+l_x_k_st-l_x_k_k+l_x_etaHat_k_k-l_x_etaHat_k_st
  
  if(as.numeric(r) > log(runif(1,0,1))){
    kappa_k <- kappa_star
    kappa_k_is_kappa_star = T
    accept <- accept + 1
  }else{
    kappa_k <- kappa_k
    kappa_k_is_kappa_star = F
  }
  kappa_mat[,i] <- kappa_k
}
accept/N_samples
# Take first 20% as burn-in.
chain <- kappa_mat[,seq(.5*N_samples,N_samples,1)]
chain %>% t() %>% data.frame() %>% unique()
sigma_psi <- sqrt(1/exp(chain[1,]))
sigma_tau <- sqrt(1/exp(chain[4,]))
sigma_xi <- sqrt(1/exp(chain[7,]))
sigma_gamma <- sqrt(1/exp(chain[8,]))
r_psi <- exp(chain[2,])
s_psi <- exp(chain[3,])
r_tau <- exp(chain[5,])
s_tau <- exp(chain[6,])

# cov(kappa_mat %>% t())
data.frame(val = sigma_psi, indx = 1:length(sigma_psi)) %>% ggplot(aes(indx,val)) + geom_line()
data.frame(val = s_psi, indx = 1:length(s_psi)) %>% ggplot(aes(indx,val)) + geom_line()
data.frame(val = sigma_gamma, indx = 1:length(sigma_gamma)) %>% ggplot(aes(indx,val)) + geom_line()
data.frame(val = r_tau, indx = 1:length(r_tau)) %>% ggplot(aes(indx,val)) + geom_line()
data.frame(val = s_tau, indx = 1:length(s_tau)) %>% ggplot(aes(indx,val)) + geom_line()
############################ ############################ Spatial hyper for psi ############################ ############################ 
psi_hyper_sp <- data.frame(
  rho_psi = r_psi,
  s_psi = s_psi
)

N_dim_dens <- 100

kd <- with(psi_hyper_sp, MASS::kde2d(s_psi, rho_psi, n = 100))

tmpDraslDt <- data.table()
for(i in 1:N_dim_dens){
  tmpDraslDt <- rbind(tmpDraslDt, data.table(z = kd$z[i,], x = rep(kd$x[i],N_dim_dens), y=kd$y))
}

p1 <- ggplot(tmpDraslDt, aes(x, y, z = z))+ geom_raster(aes(fill = z)) + geom_contour(col = "white") + xlim(c(0.22,0.4)) + ylim(c(NA,0.08)) + 
  labs(
    x = expression(s[psi]),
    y = expression(rho[psi]),
    fill = "Density"
  ) + scale_fill_gradientn(colours = terrain.colors(10)) + theme_bw()

long_hyper_psi_sp <- gather(psi_hyper_sp, key = Type, value = Value)
long_hyper_psi_sp$Type <- factor(long_hyper_psi_sp$Type)
levels(long_hyper_psi_sp$Type) <- c("rho[psi]", "s[psi]")
p2 <- ggplot(long_hyper_psi_sp, aes(Value)) + geom_density(fill = "black", alpha = .4) + facet_wrap(~Type, scales = "free", labeller = label_parsed) + theme_bw() + labs(x = "", y = "Density")

multiplot(p1, p2)

############################ ############################ Spatial hyper for tau ############################ ############################ 
tau_hyper_sp <- data.table(
  rho_tau = r_tau,
  s_tau = s_tau
)

N_dim_dens <- 100

kd <- with(tau_hyper_sp, MASS::kde2d(s_tau, rho_tau, n = 100))

tmpDraslDt2 <- data.table()
for(i in 1:N_dim_dens){
  tmpDraslDt2 <- rbind(tmpDraslDt2, data.table(z = kd$z[i,], x = rep(kd$x[i],N_dim_dens), y=kd$y))
}

p1 <- ggplot(tmpDraslDt2, aes(x, y, z = z))+ geom_raster(aes(fill = z)) + geom_contour(col = "white")  + ylim(c(NA,0.16)) + xlim(c(NA,0.25)) + 
  labs(
    x = expression(s[tau]),
    y = expression(rho[tau]),
    fill = "Density"
  ) + scale_fill_gradientn(colours = terrain.colors(10)) + theme_bw()

long_hyper_tau_sp <- gather(tau_hyper_sp, key = Type, value = Value)
long_hyper_tau_sp$Type <- factor(long_hyper_tau_sp$Type)
levels(long_hyper_tau_sp$Type) <- c("rho[tau]", "s[tau]")
p2 <- ggplot(long_hyper_tau_sp, aes(Value)) + geom_density(fill = "black", alpha = .4) + facet_wrap(~Type, scales = "free", labeller = label_parsed) + theme_bw() + labs(x = "", y = "Density")

multiplot(p1, p2)

##################################################  iid random effects ############################ ############################ 

model_error_hyper <- data.table(
  sigma_psi = sigma_psi,
  sigma_tau = sigma_tau,
  sigma_xi = sigma_xi,
  sigma_gamma = sigma_gamma
)

long_model_error_hyper <- gather(model_error_hyper, key = Type, value = Value)
long_model_error_hyper$Type <- factor(long_model_error_hyper$Type, levels = c("sigma_psi","sigma_tau","sigma_xi","sigma_gamma"))
levels(long_model_error_hyper$Type) <- c("sigma[psi*epsilon]","sigma[tau*epsilon]","sigma[phi*epsilon]","sigma[gamma*epsilon]")

ggplot(long_model_error_hyper %>% filter(!(Type == "sigma[gamma*epsilon]" & Value > 0.0015)), aes(Value)) + geom_density(fill = "black", alpha = .3) + 
  facet_wrap(~Type, scales = "free", ncol = 2, labeller = label_parsed) + 
  labs(x="",y="Density") + theme_bw()


#################################### Latent parameter posterior samples ####################################
N_x_loops <- dim(chain)[2]
x_mat <- matrix(NA, nrow = K+4*N, ncol = N_x_loops)
kappa_old <- rep(0,8)
start <- proc.time()
Q_x <- Matrix(0,nrow = 4*N+K, ncol = 4*N+K)
B <- Matrix(0,nrow = 4*N, ncol = 4*N+K)
B[1:(4*N),1:(4*N)] <- Diagonal(4*N, 1)
b <- t(B)%*%(Q_etay%*%eta_hat)
for(i in 1:N_x_loops){
  k <- list(
    kappa_psi = chain[1,i],
    v_psi = chain[2,i],
    u_psi = chain[3,i],
    kappa_tau = chain[4,i],
    v_tau = chain[5,i],
    u_tau = chain[6,i],
    kappa_xi = chain[7,i],
    kappa_gamma = chain[8,i]
  )
  Q_x <- calc_Q_x(kappa = k, Q_beta_psi = Q_beta_psi, Q_beta_tau = Q_beta_tau, 
                  Q_beta_xi = Q_beta_xi, Q_beta_gamma = Q_beta_gamma, N = N, Z = Z, d_spde = d_spde)
  Q_x_given_etaHat <- calc_Q_x_given_etaHat(Q_x = Q_x,N = N, Q_etay = Q_etay)
  mu_x_given_etaHat <- solve(a = Cholesky(Q_x_given_etaHat),b = b)
  x <- inla.qsample(n = 1, mu = mu_x_given_etaHat, Q = Q_x_given_etaHat)
  x_mat[,i] <- x
  if(i %% 100 == 0){
    print(i)
  }
}


##############################################  covariate coefficients ##############################################

nu <- x_mat[(4*N+1):nrow(x_mat),]
beta_psi <- nu[1:N_psi,]
u_psi <- nu[1:N_colsA + N_psi,]
beta_tau <- nu[1:N_tau + N_psi + N_colsA,]
u_tau <- nu[1:N_colsA + N_tau + N_psi + N_colsA,]
beta_xi <- nu[1:N_xi + N_colsA + N_tau + N_psi + N_colsA,]
beta_gamma <- nu[1:N_gamma + N_xi + N_colsA + N_tau + N_psi + N_colsA,]

apply(beta_psi, 1, mean)
apply(beta_psi, 1, sd)
apply(beta_psi, 1, quantile, probs = 0.025)
apply(beta_psi, 1, quantile, probs = 0.5)
apply(beta_psi, 1, quantile, probs = 0.975)

apply(beta_tau, 1, mean)
apply(beta_tau, 1, sd)
apply(beta_tau, 1, quantile, probs = 0.025)
apply(beta_tau, 1, quantile, probs = 0.5)
apply(beta_tau, 1, quantile, probs = 0.975)

apply(beta_xi, 1, mean)
apply(beta_xi, 1, sd)
apply(beta_xi, 1, quantile, probs = 0.025)
apply(beta_xi, 1, quantile, probs = 0.5)
apply(beta_xi, 1, quantile, probs = 0.975)

apply(beta_gamma, 1, mean)
apply(beta_gamma, 1, sd)
apply(beta_gamma, 1, quantile, probs = 0.025)
apply(beta_gamma, 1, quantile, probs = 0.5)
apply(beta_gamma, 1, quantile, probs = 0.975)

