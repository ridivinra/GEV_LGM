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
  # Bæta við normal prior sem er með sd = hálft delda 0, mu = 0
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
hess <- data.frame()
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
Sigma_kappa_0 <- diag(c(tau_psi_0$sd^2, rho_psi_0$sd^2, sigma_psi_0$sd^2,
                        tau_tau_0$sd^2, rho_tau_0$sd^2, sigma_tau_0$sd^2, 
                        tau_xi_0$sd^2, tau_gamma_0$sd^2/3)
)*0.5

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

# Take first 20% as burn-in.
chain <- kappa_mat[,seq(.2*loops,loops,1)]

sigma_psi <- sqrt(1/exp(chain[1,]))
sigma_tau <- sqrt(1/exp(chain[4,]))
sigma_xi <- sqrt(1/exp(chain[7,]))
sigma_gamma <- sqrt(1/exp(chain[8,]))
r_psi <- exp(chain[2,])
s_psi <- exp(chain[3,])
r_tau <- exp(chain[5,])
s_tau <- exp(chain[6,])

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

### Hyperparameters ### ### ### ### ### ### ### ### ### ### ### 
sigma_psi <- sqrt(1/exp(chain[1,]))
sigma_tau <- sqrt(1/exp(chain[4,]))
sigma_xi <- sqrt(1/exp(chain[7,]))
sigma_gamma <- sqrt(1/exp(chain[8,]))
mean(sigma_psi)
sd(sigma_psi)
quantile(sigma_psi, c(.025,.5,.975))
mean(sigma_tau)
sd(sigma_tau)
quantile(sigma_tau, c(.025,.5,.975))
mean(sigma_xi)
sd(sigma_xi)
quantile(sigma_xi, c(.025,.5,.975))
mean(sigma_gamma)
sd(sigma_gamma)
quantile(sigma_gamma, c(.025,.5,.975))

r_psi <- exp(chain[2,])
mean(r_psi)
sd(r_psi)
quantile(r_psi, c(.025,.5,.975))
ggplot(data.frame(x = r_psi), aes(x)) + geom_density()

s_psi <- exp(chain[3,])
mean(s_psi)
sd(s_psi)
quantile(s_psi, c(.025,.5,.975))
ggplot(data.frame(x = s_psi), aes(x)) + geom_density()

r_tau <- exp(chain[5,])
mean(r_tau)
sd(r_tau)
quantile(r_tau, c(.025,.5,.975))
ggplot(data.frame(x = r_psi), aes(x)) + geom_density()

s_tau <- exp(chain[6,])
mean(s_tau)
sd(s_tau)
quantile(s_tau, c(.025,.5,.975))
ggplot(data.frame(x = s_psi), aes(x)) + geom_density()

###### Hyperparamter plots ###### Hyperparamter plots ###### Hyperparamter plots ###### Hyperparamter plots ########
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

##################################################  Model error parameters ############################ ############################ 

model_error_hyper <- data.table(
  sigma_psi = sigma_psi,
  sigma_tau = sigma_tau,
  sigma_xi = sigma_xi,
  sigma_gamma = sigma_gamma
)

long_model_error_hyper <- gather(model_error_hyper, key = Type, value = Value)
long_model_error_hyper$Type <- factor(long_model_error_hyper$Type, levels = c("sigma_gamma","sigma_psi","sigma_tau","sigma_xi"))
levels(long_model_error_hyper$Type) <- c("sigma[gamma*epsilon]","sigma[psi*epsilon]","sigma[tau*epsilon]","sigma[phi*epsilon]")

ggplot(long_model_error_hyper, aes(Value)) + geom_density(fill = "black", alpha = .3) + 
  facet_wrap(~Type, scales = "free", ncol = 1, labeller = label_parsed) + 
  labs(x="",y="Density") + theme_bw()

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

##############################################  covariate coefficients ##############################################
##### Covariate table comparing median with quartiles ################################################
X_psi <- covariates[, c("Int.","log(Area)", "log(SAAR)", "log(FARL)", "BFIHOST^2")]
X_tau <- covariates[, c("Int.","log(Area)", "log(SAAR)", "log(FARL)", "log(FPEXT)", "log(URBEXT+1)")]
X_kappa <- covariates[, c("Int.", "log(FPEXT)")]

beta_psi_mean <- apply(beta_psi[-1,], 1, mean)
beta_tau_mean <- apply(beta_tau[-1,], 1, mean)
beta_sigma_mean <- c(
  beta_psi_mean[1]+beta_tau_mean[1],
  beta_psi_mean[2]+beta_tau_mean[2],
  beta_psi_mean[3]+beta_tau_mean[3],
  beta_tau_mean[4],
  beta_tau_mean[5]
)
dt_table_fixed_med_quartile <- data.table(
  desc = c(cov_names_psi[-1], cov_names_tau[-1]), 
  Mean = c(beta_psi_mean, beta_sigma_mean),
  Median = c(
    quantile(data$Area,.5),
    quantile(data$SAAR,.5),
    quantile(data$FARL,.5),
    quantile(exp(data$BFIHOST^2),.5),
    quantile(data$Area,.5),
    quantile(data$SAAR,.5),
    quantile(data$FARL,.5),
    quantile(data$FPEXT,.5),
    quantile(data$URBEXT2000 + 1,.5)
  ),
  firstQ = c(
    quantile(data$Area,.25),
    quantile(data$SAAR,.25),
    quantile(data$FARL,.25),
    quantile(exp(data$BFIHOST^2),.25),
    quantile(data$Area,.25),
    quantile(data$SAAR,.25),
    quantile(data$FARL,.25),
    quantile(data$FPEXT,.25),
    quantile(data$URBEXT2000 + 1,.25)
  ),
  thirdQ = c(
    quantile(data$Area,.75),
    quantile(data$SAAR,.75),
    quantile(data$FARL,.75),
    quantile(exp(data$BFIHOST^2),.75),
    quantile(data$Area,.75),
    quantile(data$SAAR,.75),
    quantile(data$FARL,.75),
    quantile(data$FPEXT,.75),
    quantile(data$URBEXT2000 + 1,.75)
  )
)

table_final_quart_med <- dt_table_fixed_med_quartile %>% 
  mutate(third_med = (thirdQ/Median)^Mean, first_med = (firstQ/Median)^Mean) %>% 
  dplyr::select(desc, third_med, first_med)

beta_xi_mean <- apply(beta_xi, 1, mean)

xi_trans <- function(phi){
  xi00 <- 0
  alp3 <- 0.8 # c_phi
  alpha <- 4
  beta <- 4
  sig3 = (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
  b3 <- -sig3*(log(-log(1-0.5^alp3)))
  xihat <- ((1-(exp(-exp((phi - b3)/sig3))))^(1/alp3) - 0.5)
  return(xihat)
}

phi_vals <- data.frame(med = xi_trans(beta_xi_mean[1] + beta_xi_mean[2]*median(log(data$FPEXT))), 
                       Q1 = xi_trans(beta_xi_mean[1] + beta_xi_mean[2]*quantile(log(data$FPEXT),.25)), 
                       Q3 = xi_trans(beta_xi_mean[1] + beta_xi_mean[2]*quantile(log(data$FPEXT),.75)))

phi_vals[3]-phi_vals[1]
phi_vals[2]-phi_vals[1]

beta_gamma_mean <- apply(beta_gamma, 1, mean)

gamma_vals <- data.frame(med = d(beta_gamma_mean[1] + beta_gamma_mean[2]*median(log(data$PROPWET))), 
                         Q1 = d(beta_gamma_mean[1] + beta_gamma_mean[2]*quantile(log(data$PROPWET),.25)), 
                         Q3 = d(beta_gamma_mean[1] + beta_gamma_mean[2]*quantile(log(data$PROPWET),.75)))

gamma_vals[3]-gamma_vals[1]
gamma_vals[2]-gamma_vals[1]


a_upsi <- apply(A_mat %*% u_psi, 1, mean)
a_utau <- apply(A_mat %*% u_tau, 1, mean)
a_usigma <- apply(A_mat %*% u_tau + A_mat %*% u_psi, 1, mean)

exp(quantile(a_upsi, probs = 0.75))/exp(median(a_upsi))
exp(quantile(a_upsi, probs = 0.25))/exp(median(a_upsi))

exp(quantile(a_utau, probs = 0.75))/exp(median(a_utau))
exp(quantile(a_utau, probs = 0.25))/exp(median(a_utau))

exp(quantile(a_usigma, probs = 0.75))/exp(median(a_usigma))
exp(quantile(a_usigma, probs = 0.25))/exp(median(a_usigma))

## Do the same as ^ but for quantiles.
eta <- x_mat[1:(4*554),]
psi <- eta[1:554,]
tau <- eta[1:554 + 554,]
phi <- eta[1:554 + 2*554,]
gamma <- eta[1:554 + 3*554,]
dim(gamma)
xi_trans <- function(phi){
  xi00 <- 0
  alp3 <- 0.8 # c_phi
  alpha <- 4
  beta <- 4
  sig3 = (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
  b3 <- -sig3*(log(-log(1-0.5^alp3)))
  xihat <- ((1-(exp(-exp((phi - b3)/sig3))))^(1/alp3) - 0.5)
  return(xihat)
}


mu <- exp(psi)
sigma <- exp(psi + tau)
xif <- xi_trans(phi)

mu_p <- apply(mu, 1, mean)
sigma_p <- apply(sigma, 1, mean)
xif_p <- apply(xif, 1, mean)

med_mu <- median(mu_p)
med_sigma <- median(sigma_p)
med_xi <- median(xif_p)

quantileGEV <- function(mu,sigma,xi,p){
  q <- mu + sigma*((-log(p))^(-xi)-1)/xi
  return(q)
}

med_q <- quantileGEV(med_mu,med_sigma,med_xi,.99)

# Area
area_q <- data.frame(
  upper = quantileGEV(2.093*med_mu, 2*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(0.488*med_mu, 0.510*med_sigma,med_xi, p = .99)/med_q
)
# SAAR
saar_q <- data.frame(
  upper = quantileGEV(1.548*med_mu, 1.294*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(.673*med_mu, .792*med_sigma,med_xi, p = .99)/med_q
)
#farl
farl_q <- data.frame(
  upper = quantileGEV(1.056*med_mu, 1.04*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(.891*med_mu, 0.921*med_sigma,med_xi, p = .99)/med_q
)
#BFIHOST
bfihost_q <- data.frame(
  upper = quantileGEV(0.749*med_mu, 1*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(1.187*med_mu, 1*med_sigma,med_xi, p = .99)/med_q
)
#URBEXT
urbext_q <- data.frame(
  upper = quantileGEV(1*med_mu, .984*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(1*med_mu, 1.004*med_sigma,med_xi, p = .99)/med_q
)
#fpext
fpext_q <- data.frame(
  upper = quantileGEV(1*med_mu, .952*med_sigma,med_xi-0.013, p = .99)/med_q,
  lower = quantileGEV(1*med_mu, 1.054*med_sigma,med_xi+0.014, p = .99)/med_q
)
#propwet
propwet_q <- data.frame(
  upper = quantileGEV(1*med_mu, 1*med_sigma,med_xi-0.013, p = .99)/med_q,
  lower = quantileGEV(1*med_mu, 1*med_sigma,med_xi+0.014, p = .99)/med_q
)

#spatialpsi
spPSI_q <- data.frame(
  upper = quantileGEV(1.176*med_mu, 1.176*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(0.835*med_mu, 0.835*med_sigma,med_xi, p = .99)/med_q
)

#spatialpsitau
spPSITAU_q <- data.frame(
  upper = quantileGEV(1*med_mu, 1.094*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(1*med_mu, 0.896*med_sigma,med_xi, p = .99)/med_q
)
#spatialtau
spTAU_q <- data.frame(
  upper = quantileGEV(1*med_mu, 1.111*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(1*med_mu, 0.908*med_sigma,med_xi, p = .99)/med_q
)

apply(beta_gamma,1,mean)

# Unstructured randomeffect
sigma_psi_m <- mean(sigma_psi)
qnorm(p = c(0.25,0.75), sd = sigma_psi_m)
exp(0.1666197)/exp(0)
exp(-0.1666197)/exp(0)
uepsi_q <- data.frame(
  upper = quantileGEV(1.176796*med_mu, 1*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(0.8497647*med_mu, 1*med_sigma,med_xi, p = .99)/med_q
)
sigma_tau_m <- mean(sigma_tau)
qnorm(p = c(0.25,0.75), sd = sigma_tau_m)
exp(0.0897694)/exp(0)
exp(-0.0897694)/exp(0)
uepsi_q <- data.frame(
  upper = quantileGEV(1*med_mu, 1.093922*med_sigma,med_xi, p = .99)/med_q,
  lower = quantileGEV(1*med_mu, 0.914142*med_sigma,med_xi, p = .99)/med_q
)
sigma_xi_m <- mean(sigma_xi)
qnorm(p = c(0.25,0.75), sd = sigma_xi_m)
phi_vals2 <- data.frame(med = xi_trans(beta_xi_mean[1] + beta_xi_mean[2]*median(log(data$FPEXT))), 
                        Q1 = xi_trans(beta_xi_mean[1] + beta_xi_mean[2]*median(log(data$FPEXT)) + 0.04574787), 
                        Q3 = xi_trans(beta_xi_mean[1] + beta_xi_mean[2]*median(log(data$FPEXT)) -0.04574787 ))
phi_vals2[2]- phi_vals2[1]
phi_vals2[3]- phi_vals2[1]
uepsi_q <- data.frame(
  upper = quantileGEV(1*med_mu, 1*med_sigma,med_xi + 0.04664247, p = .99)/med_q,
  lower = quantileGEV(1*med_mu, 1*med_sigma,med_xi-0.04511683, p = .99)/med_q
)
sigma_gamma_m <- mean(sigma_gamma)
qnorm(p = c(0.25,0.75), sd = sigma_gamma_m)
phi_vals2 <- data.frame(med = d(beta_gamma_mean[1] + beta_gamma_mean[2]*median(log(data$PROPWET))), 
                        Q1 = d(beta_gamma_mean[1] + beta_gamma_mean[2]*median(log(data$PROPWET)) - 0.0003503697), 
                        Q3 = d(beta_gamma_mean[1] + beta_gamma_mean[2]*median(log(data$PROPWET)) +0.0003503697 ))
beta_gamma_mean = apply(beta_gamma,1,mean)
phi_vals2[2]- phi_vals2[1]
phi_vals2[3]- phi_vals2[1]
####### CDF and quantile figures ######################################################################
tmp_countYears <- dt_mles %>% left_join(allFlow %>% unique(), by = "Station") %>% group_by(Station) %>% summarise(n=n())
tmp_countYears[which(tmp_countYears$n > 60),]
xif_p[which(tmp_countYears$n > 60)]
dt_mles[which(round(xif_p,10) == 0.0910233676),]$Station
dt_mles[which(round(xif_p,8) == -0.09461199),]$Station
dt_mles[which(round(xif_p,10) == 0.0126497601),]$Station
dt_mles[which(round(xif_p,10) == 0.0837216170),]$Station
rand_stations <- c("44008", "37020", "54020", "25011")
rndIndx <- which(data$Station %in% rand_stations)
quantiles <- array(dim = c(99,dim(x_mat)[2],4))
ps <- seq(0.01,.99,0.01)
d <- function(gamma){
  delta_0 <- 1/125
  res <- 2*delta_0*exp( 2*gamma/delta_0 )/(1 + exp(2*gamma/delta_0)) - delta_0
  return(res)
}
xi_trans <- function(phi){
  xi00 <- 0
  alp3 <- 0.8 # c_phi
  alpha <- 4
  beta <- 4
  sig3 = (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
  b3 <- -sig3*(log(-log(1-0.5^alp3)))
  xihat <- ((1-(exp(-exp((phi - b3)/sig3))))^(1/alp3) - 0.5)
  return(xihat)
}
for(i in 1:length(rndIndx)){
  j <- rndIndx[i]
  psi_tmp_calc <- x_mat[j,]
  tau_tmp_calc <- x_mat[j + N,]
  xi_tmp_calc <- x_mat[j + 2*N,]
  for(k in 1:length(ps)){
    p <- ps[k]
    quantiles[k,,i] <- quantileGEV(psi_tmp_calc, tau_tmp_calc, xi_trans(xi_tmp_calc), p = p)
  }
}

quantileGEV <- function(psi,tau,xi,p){
  mu <- exp(psi)
  sigma <- exp(tau + psi)
  q <- mu + sigma*((-log(p))^(-xi)-1)/xi
  return(q)
}
### Vantar hliðrun hér skv time-trend!
delta_mean <- apply(d(gamma), 1, mean)
mu_p <- apply(mu, 1, mean)
med_mu_r <- mu_p[rndIndx]
time_trend <- data.frame(Station = data[rndIndx,]$Station, delta = delta_mean[rndIndx], mu = med_mu_r)
sorted_data <- unique(filter(allAmax,Station %in% data[rndIndx,]$Station)) %>% arrange(Station, Flow) 
sorted_data <- sorted_data %>% left_join(time_trend, by = "Station") %>% 
  mutate(Flow = Flow - mu*delta*(WaterYear-1975))

sorted_dataNs <- (sorted_data %>% group_by(Station) %>% summarise(n = n()))$n
sordet_data_sampleStations <- mutate(sorted_data, 
                                     indx = sorted_dataNs %>% sapply(function(x){return(1:x)}) %>% unlist(),
                                     n = sorted_dataNs %>% sapply(function(x){return(rep(c(x),x))}) %>% unlist() ) %>%
  mutate(prob = indx/(n+1))

stations <- data.table(
  mean = c(apply(quantiles[,,1], 1, mean),
           apply(quantiles[,,2], 1, mean),
           apply(quantiles[,,3], 1, mean),
           apply(quantiles[,,4], 1, mean)),
  q_025 = c(apply(quantiles[,,1], 1, quantile, probs = .025),
            apply(quantiles[,,2], 1, quantile, probs = .025),
            apply(quantiles[,,3], 1, quantile, probs = .025),
            apply(quantiles[,,4], 1, quantile, probs = .025)),
  q_975 = c(apply(quantiles[,,1], 1, quantile, probs = .975),
            apply(quantiles[,,2], 1, quantile, probs = .975),
            apply(quantiles[,,3], 1, quantile, probs = .975),
            apply(quantiles[,,4], 1, quantile, probs = .975)),
  Prob = rep(ps,4),
  T_rate = rep(1/(1-ps),4),
  Stations = rep(paste("Station", data[rndIndx,]$Station), each = length(ps))
)
sordet_data_sampleStations <- sordet_data_sampleStations %>% mutate(Stations = paste("Station",Station))

p_quantile1 <- ggplot(stations %>% dplyr::filter(Stations %in% c("Station 54020", "Station 37020", "Station 25011")),y="log") + geom_line(aes(x = T_rate, y = mean)) + 
  geom_line(mapping = aes(x = T_rate, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = T_rate, y = q_975), col = "blue",linetype = "dashed") + 
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 54020", "Station 37020", "Station 25011")),
            aes(x = 1/(1-Mean),y=Q_high), col="red",linetype="dashed") +
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 54020", "Station 37020", "Station 25011")),
            aes(x = 1/(1-Mean),y=Q_low), col="red",linetype="dashed") +
  geom_point(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 54020", "Station 37020", "Station 25011")),aes(x = 1/(1-Mean),y=Flow), size = .6) +
  facet_wrap(~Stations, scales = "free") + labs(y = expression(paste("Flow [",m^3,"/s]")), x="Return period [years]") + theme_bw() + 
  scale_x_log10()

p_quantile1_1 <- ggplot(stations %>% filter(Stations %in% c("Station 54020")),y="log") + geom_line(aes(x = T_rate, y = mean)) + 
  geom_line(mapping = aes(x = T_rate, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = T_rate, y = q_975), col = "blue",linetype = "dashed") + 
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 54020")),
            aes(x = 1/(1-Mean),y=Q_high), col="red",linetype="dashed") +
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 54020")),
            aes(x = 1/(1-Mean),y=Q_low), col="red",linetype="dashed") +
  geom_point(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 54020")),aes(x = 1/(1-Mean),y=Flow), size = .6) +
  facet_wrap(~Stations, scales = "free") + labs(y = expression(paste("Flow [",m^3,"/s]")), x="") + theme_bw() + 
  scale_x_log10()

p_quantile1_2 <- ggplot(stations %>% filter(Stations %in% c("Station 37020")),y="log") + geom_line(aes(x = T_rate, y = mean)) + 
  geom_line(mapping = aes(x = T_rate, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = T_rate, y = q_975), col = "blue",linetype = "dashed") + 
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 37020")),
            aes(x = 1/(1-Mean),y=Q_high), col="red",linetype="dashed") +
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 37020")),
            aes(x = 1/(1-Mean),y=Q_low), col="red",linetype="dashed") +
  geom_point(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 37020")),aes(x = 1/(1-Mean),y=Flow), size = .6) +
  facet_wrap(~Stations, scales = "free") + labs(y = "",x="Return period [years]") + theme_bw() + 
  scale_x_log10()

p_quantile1_3 <- ggplot(stations %>% filter(Stations %in% c("Station 25011")),y="log") + geom_line(aes(x = T_rate, y = mean)) + 
  geom_line(mapping = aes(x = T_rate, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = T_rate, y = q_975), col = "blue",linetype = "dashed") + 
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 25011")),
            aes(x = 1/(1-Mean),y=Q_high), col="red",linetype="dashed") +
  geom_line(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 25011")),
            aes(x = 1/(1-Mean),y=Q_low), col="red",linetype="dashed") +
  geom_point(data = sordet_data_sampleStations %>% filter(Stations %in% c("Station 25011")),aes(x = 1/(1-Mean),y=Flow), size = .6) +
  facet_wrap(~Stations, scales = "free") + labs(y = "",x = "") + theme_bw() + 
  scale_x_log10()

multiplot(plotlist = list(p_quantile1_1,p_quantile1_2,p_quantile1_3),cols = 3)

################################################### CDF !!! ####################################################################

cdfs_samples <- array(dim = c(100,dim(x_mat)[2],4))

cdfGEV <- function(psi,tau,xi,x){
  mu <- exp(psi)
  sigma <- exp(tau + psi)
  t <- (1 + xi * ( (x-mu)/sigma) )^(-1/xi)
  res <- exp(-t)
}

for(i in 1:length(rndIndx)){
  j <- rndIndx[i]
  flows <- unique(filter(allAmax, Station == data[j,]$Station))
  psi_tmp_calc <- x_mat[j,]
  tau_tmp_calc <- x_mat[j + N,]
  xi_tmp_calc <- x_mat[j + 2*N,]
  ps <- seq(min(flows$Flow), max(flows$Flow), length.out = 100)
  for(k in 1:length(ps)){
    p <- ps[k]
    cdfs_samples[k,,i] <- cdfGEV(psi = psi_tmp_calc, tau = tau_tmp_calc, xi = xi_trans(xi_tmp_calc), x = p)
  }
}

smp_st_flow <- unique(filter(allAmax, Station %in% data[rndIndx,]$Station)) %>% group_by(Station) %>% summarise(max = max(Flow), min = min(Flow))

cdf_for_plot <- data.table(
  mean = c(apply(cdfs_samples[,,4], 1, median, na.rm = T),
           apply(cdfs_samples[,,3], 1, median, na.rm = T),
           apply(cdfs_samples[,,2], 1, median, na.rm = T),
           apply(cdfs_samples[,,1], 1, median, na.rm = T)),
  q_025 = c(apply(cdfs_samples[,,4], 1, quantile, probs = .025, na.rm = T),
            apply(cdfs_samples[,,3], 1, quantile, probs = .025, na.rm = T),
            apply(cdfs_samples[,,2], 1, quantile, probs = .025, na.rm = T),
            apply(cdfs_samples[,,1], 1, quantile, probs = .025, na.rm = T)),
  q_975 = c(apply(cdfs_samples[,,4], 1, quantile, probs = .975, na.rm = T),
            apply(cdfs_samples[,,3], 1, quantile, probs = .975, na.rm = T),
            apply(cdfs_samples[,,2], 1, quantile, probs = .975, na.rm = T),
            apply(cdfs_samples[,,1], 1, quantile, probs = .975, na.rm = T)),
  Prob = c(seq(smp_st_flow$min[1], smp_st_flow$max[1], length.out = 100),
           seq(smp_st_flow$min[2], smp_st_flow$max[2], length.out = 100),
           seq(smp_st_flow$min[3], smp_st_flow$max[3], length.out = 100),
           seq(smp_st_flow$min[4], smp_st_flow$max[4], length.out = 100)),
  Stations = rep(paste("Station", sort(data[rndIndx,]$Station)), each = length(ps))
)
allAmax$Stations = paste("Station", allAmax$Station)

p_cdfPlot <- ggplot(cdf_for_plot %>% filter(Stations %in% c("Station 12001", "Station 27034", "Station 33005"))) + 
  geom_line(mapping = aes(x = Prob, y = mean)) + 
  geom_line(mapping = aes(x = Prob, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = Prob, y = q_975), col = "blue",linetype = "dashed") + 
  stat_ecdf(data = filter(allAmax, Stations %in% c("Station 12001", "Station 27034", "Station 33005")) ,mapping = aes(Flow)) + 
  facet_wrap(~Stations, scales = "free") + labs(y = "Cdf", x=expression(paste("Flow [",m^3,"/s]"))) + theme_bw()

p_cdfPlot_1 <- ggplot(cdf_for_plot %>% filter(Stations %in% c("Station 54020"))) + 
  geom_line(mapping = aes(x = Prob, y = mean)) + 
  geom_line(mapping = aes(x = Prob, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = Prob, y = q_975), col = "blue",linetype = "dashed") + 
  stat_ecdf(data = filter(allAmax, Stations %in% c("Station 54020")) ,mapping = aes(Flow)) + 
  facet_wrap(~Stations, scales = "free") + labs(y = "Cdf", x="") + theme_bw()

p_cdfPlot_2 <- ggplot(cdf_for_plot %>% filter(Stations %in% c("Station 37020"))) + 
  geom_line(mapping = aes(x = Prob, y = mean)) + 
  geom_line(mapping = aes(x = Prob, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = Prob, y = q_975), col = "blue",linetype = "dashed") + 
  stat_ecdf(data = filter(allAmax, Stations %in% c("Station 37020")) ,mapping = aes(Flow)) + 
  facet_wrap(~Stations, scales = "free") + labs(y = "", x=expression(paste("Flow [",m^3,"/s]"))) + theme_bw()

p_cdfPlot_3 <- ggplot(cdf_for_plot %>% filter(Stations %in% c("Station 25011"))) + 
  geom_line(mapping = aes(x = Prob, y = mean)) + 
  geom_line(mapping = aes(x = Prob, y = q_025),col = "blue",linetype = "dashed") + 
  geom_line(mapping = aes(x = Prob, y = q_975), col = "blue",linetype = "dashed") + 
  stat_ecdf(data = filter(allAmax, Stations %in% c("Station 25011")) ,mapping = aes(Flow)) + 
  facet_wrap(~Stations, scales = "free") + labs(y = "", x="") + theme_bw()


######################################## Quantile prediction ########################################################################
rand_stations <- c("12001", "27034", "33005", "45816")
sorted_data <- unique(filter(allAmax,Station %in% rand_stations)) %>% arrange(Station, Flow)  %>% unique()
sorted_dataNs <- (sorted_data %>% group_by(Station) %>% summarise(n = n()))$n
sordet_data_sampleStations <- mutate(sorted_data, 
                                     indx = sorted_dataNs %>% sapply(function(x){return(1:x)}) %>% unlist(),
                                     n = sorted_dataNs %>% sapply(function(x){return(rep(c(x),x))}) %>% unlist() ) %>%
  mutate(prob = indx/(n+1))
sordet_data_sampleStations$Low <- NA
sordet_data_sampleStations$Mean <- NA
sordet_data_sampleStations$High <- NA
for(i in 1:dim(sordet_data_sampleStations)[1]){
  a_tmp <- sordet_data_sampleStations$indx[i]
  b_tmp <- sordet_data_sampleStations$n[i] + 1 - sordet_data_sampleStations$indx[i]
  sordet_data_sampleStations$Low[i] <- qbeta(0.025, a_tmp, b_tmp)
  sordet_data_sampleStations$Mean[i] <- a_tmp/(a_tmp+b_tmp)
  sordet_data_sampleStations$High[i] <- qbeta(0.975, a_tmp, b_tmp)
}

sordet_data_sampleStations <- sordet_data_sampleStations %>%
  mutate(
    Q_low = NA,
    Q_mean = NA,
    Q_high = NA
  )
quantileGEV2 <- function(mu,sigma,xi,p){
  q <- mu + sigma*( (-log(p))^(-xi) - 1 )/xi
  return(q)
}


st_indx <- 0
for(i in 1:dim(sordet_data_sampleStations)[1]){
  if(st_indx == 0){
    st_indx <- which(data$Station == sordet_data_sampleStations$Station[i])
    mu_t <- exp(x_mat[st_indx,])
    sigma_t <- exp(x_mat[st_indx+N,]+x_mat[st_indx,])
    xi_t <- xi_trans(x_mat[st_indx+2*N,])
  }else if(st_indx != which(data$Station == sordet_data_sampleStations$Station[i])){
    st_indx <- which(data$Station == sordet_data_sampleStations$Station[i])
    mu_t <- exp(x_mat[st_indx,])
    sigma_t <- exp(x_mat[st_indx+N,]+x_mat[st_indx,])
    xi_t <- xi_trans(x_mat[st_indx+2*N,])
  }
  p_low <- sordet_data_sampleStations$Low[i]
  p_mean <- sordet_data_sampleStations$Mean[i]
  p_high <- sordet_data_sampleStations$High[i]
  sordet_data_sampleStations$Q_low[i] = mean(quantileGEV2(mu = mu_t, sigma = sigma_t, xi = xi_t,p = p_low))
  sordet_data_sampleStations$Q_mean[i] = mean(quantileGEV2(mu = mu_t, sigma = sigma_t, xi = xi_t,p = p_mean))
  sordet_data_sampleStations$Q_high[i] = mean(quantileGEV2(mu = mu_t, sigma = sigma_t, xi = xi_t,p = p_high))
}
sordet_data_sampleStations$Stations
p_quantilF <- ggplot(sordet_data_sampleStations %>% filter(Stations %in% c("Station 12001", "Station 27034", "Station 33005")), 
                     aes(x = 1/(1-Mean), y=Q_mean)) + 
  geom_line() + geom_line(aes(y=Q_low), col="blue",linetype="dashed") + 
  geom_line(aes(y=Q_high), col="blue",linetype="dashed") +
  facet_wrap(~Stations, scales = "free") + theme_bw() + geom_point(aes(y=Flow), size = .6) + 
  scale_x_log10(limits = c(1,100),breaks = c(5,10,20,40,100)) +
  labs(x = "Return period [years]", y=expression(paste("Flow [",m^3,"/s]")))

multiplot(plotlist = list( p_quantile1_1,p_cdfPlot_1,p_quantile1_2,p_cdfPlot_2,p_quantile1_3,p_cdfPlot_3), cols = 3)

library(egg)
ggarrange(p_quantile1_1,p_quantile1_2,p_quantile1_3,p_cdfPlot_1,p_cdfPlot_2,p_cdfPlot_3, ncol=3)
######################################## Fit models without spatial component ##################################################

get_mean_sd_for_hyperparameter <- function(mdl, fun, hyperparameter){
  x0 <- mdl$marginals.hyperpar[[hyperparameter]]
  E <- inla.emarginal(function(x) c(fun(x), fun(x)^2), x0)
  sd <- sqrt(E[2] - E[1]^2)
  mean <- E[1]
  
  return(list(sd = sd, mean = mean))
}
dim(data)

mdl_whitoutSp_psi <- fit_model_psi3(data = data, desc = X_psi)
mdl_whitoutSp_tau <- fit_model_tau3(data = data, desc = X_tau)

mdl_whitoutSp_psi$mdl$summary.hyperpar

get_mean_sd_for_hyperparameter(mdl = mdl_whitoutSp_psi$mdl, fun = function(x){1/sqrt(x)}, hyperparameter = 1)

get_mean_sd_for_hyperparameter(mdl = mdl_whitoutSp_tau$mdl, fun = function(x){1/sqrt(x)}, hyperparameter = 1)


############################################ Fit model for out of sample stations ################################################


#############################################  ATH ÞETTA!!!!1!!!!!!!!!!! #########################################################
d <- function(gamma){
  delta_0 <- 1/125
  res <- 2*delta_0*exp( 2*gamma/delda_0 )/(1 + exp(2*gamma/delda_0)) - delta_0
  return(res)
}
d_and <- function(delta){
  
}
xi_trans <- function(phi){
  xi00 <- 0
  alp3 <- 0.8 # c_phi
  alpha <- 4
  beta <- 4
  sig3 = (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
  b3 <- -sig3*(log(-log(1-0.5^alp3)))
  xihat <- ((1-(exp(-exp((phi - b3)/sig3))))^(1/alp3) - 0.5)
  return(xihat)
}
# Compare eta with MLE
eta <- x_mat[1:(4*N),]
psi <- eta[1:554,]
tau <- eta[1:554 + 554,]
phi <- eta[1:554 + 2*554,]
gamma <- eta[1:554 + 3*554,]
mu <- exp(psi)
sigma <- exp(psi + tau)
xi <- xi_trans(phi)
delta <- d(gamma)

plot_bayes_dt <- data.frame(
  values = c(
    apply(psi, 1, mean),
    apply(tau, 1, mean),
    apply(phi, 1, mean),
    apply(gamma, 1, mean),
    apply(mu, 1, mean),
    apply(sigma, 1, mean),
    apply(xi, 1, mean),
    apply(delta, 1, mean)
  ),
  type = factor(rep(c("psi","tau","phi", "gamma","mu","sigma","xi", "Delta"), each = N), 
                levels = c("psi","tau","phi", "gamma","mu","sigma","xi", "Delta"))
)

plot_bayes_dt %>% ggplot(aes(values)) + geom_histogram(bins = 25) + 
  facet_wrap(~type, scales = "free",labeller = label_parsed,ncol = 4) + 
  labs(x = "Model estimates", y = "Frequency") + theme_bw()

plot_MLE_dt <- data.frame(
  values = c(
    dt_mles$psi,
    dt_mles$tau,
    dt_mles$kappa,
    dt_mles$gamma,
    exp(dt_mles$psi),
    exp(dt_mles$psi + dt_mles$tau),
    xi_trans(dt_mles$kappa),
    d(dt_mles$gamma)
  ),
  type = factor(rep(c("hat(psi)","hat(tau)","hat(phi)","hat(gamma)","hat(mu)","hat(sigma)","hat(xi)","hat(Delta)"), each = N), 
                levels = c("hat(psi)","hat(tau)","hat(phi)","hat(gamma)","hat(mu)","hat(sigma)","hat(xi)","hat(Delta)"))
)

dummy <- data.frame(values = c(-.6,.6,-.6,.6),
                    type = factor(rep(c("hat(phi)","hat(xi)"), each = 2), levels = c("hat(psi)","hat(tau)","hat(phi)","hat(gamma)","hat(mu)","hat(sigma)","hat(xi)","hat(Delta)")), stringsAsFactors=FALSE)

plot_MLE_dt %>% ggplot(aes(values)) + 
  geom_blank(data = dummy) + geom_histogram(bins = 25)  + 
  facet_wrap(~type, scales = "free",labeller = label_parsed, ncol = 4) + 
  labs(x = "ML estimates", y = "Frequency") + theme_bw() 


plot_mle_bayes_dt <- plot_MLE_dt %>% rename(ML_values = values) %>% dplyr::select(ML_values) %>% cbind(plot_bayes_dt)
psi_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "psi") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "psi") %>% select(ML_values))))
tau_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "tau") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "tau") %>% select(ML_values))))
phi_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "phi") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "phi") %>% select(ML_values))))
gamma_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "gamma") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "gamma") %>% select(ML_values))))
mu_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "mu") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "mu") %>% select(ML_values))))
sigma_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "sigma") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "sigma") %>% select(ML_values))))
xi_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "xi") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "xi") %>% select(ML_values))))
delta_range <- 
  range(c(range(plot_mle_bayes_dt %>% filter(type == "Delta") %>% select(values)),
          range(plot_mle_bayes_dt %>% filter(type == "Delta") %>% select(ML_values))))

dummy_lims <- data.frame(
  values = c(psi_range, tau_range, phi_range,gamma_range, mu_range, sigma_range, xi_range, delta_range),
  ML_values = c(psi_range, tau_range, phi_range,gamma_range, mu_range, sigma_range, xi_range,delta_range),
  type = factor(rep(c("psi","tau","phi","gamma","mu","sigma","xi","Delta"), each = 2), levels = c("psi","tau","phi","gamma","mu","sigma","xi","Delta"))
)

plot_mle_bayes_dt %>% ggplot(aes(values, ML_values)) + geom_blank(data = dummy_lims) + geom_abline(slope = 1, intercept = 0, col ="red") +
  geom_point(size = 0.5) + facet_wrap(~type, scales = "free",labeller = label_parsed, ncol = 4) + theme_bw() + labs(x = "Model estimates", y = "ML estimates")

test_plot_ptp <- data.frame(values = c(
  apply(pred_psi, 1, mean),data$psi,
  apply(pred_tau,1,mean),data$tau,
  apply(pred_xi,1,mean),data$kappa), parameter = rep(c("psi","tau","phi"), each = 2*N), type = rep(rep(c("Posterior estimates","ML estimates"),3), each = N))

psi_p <- ggplot(test_plot_ptp %>% filter(parameter == "psi"), aes(values)) + geom_histogram() + facet_wrap(~type,nrow = 2) + theme_bw() + labs(title = "psi")
tau_p <- ggplot(test_plot_ptp %>% filter(parameter == "tau"), aes(values)) + geom_histogram() + facet_wrap(~type,nrow = 2) + theme_bw() + labs(title = "tau")
phi_p <- ggplot(test_plot_ptp %>% filter(parameter == "phi"), aes(values)) + geom_histogram() + facet_wrap(~type,nrow = 2) + theme_bw() + labs(title = "phi")

multiplot(plotlist = list(psi_p,tau_p,phi_p), cols = 3)

test_plot_msx <- data.frame(values = c(
  apply(exp(pred_psi), 1, mean),data$mu,
  apply(exp(pred_tau + pred_psi),1,mean),data$sigma,
  apply(xi_trans(pred_xi),1,mean),data$xi), parameter = rep(c("mu","sigma","xi"), each = 2*N), type = rep(rep(c("Posterior estimates","ML estimates"),3), each = N))

mu_p <- ggplot(test_plot_msx %>% filter(parameter == "mu"), aes(values)) + geom_histogram() + facet_wrap(~type,nrow = 2) + theme_bw() + labs(title = "mu")
sigma_p <- ggplot(test_plot_msx %>% filter(parameter == "sigma"), aes(values)) + geom_histogram() + facet_wrap(~type,nrow = 2) + theme_bw() + labs(title = "sigma")
xi_p <- ggplot(test_plot_msx %>% filter(parameter == "xi"), aes(values)) + geom_histogram() + facet_wrap(~type,nrow = 2) + theme_bw() + labs(title = "xi")

multiplot(plotlist = list(mu_p,sigma_p,xi_p), cols = 3)

####### Forward selection mynd #####

dt_f_newest <- read.table("forwardSel.csv")
dt_f_newest$parameter <- factor(dt_f_newest$parameter, levels = c("psi","tau","phi","gamma"))
dt_plot_f_p <- dt_f_newest %>% group_by(parameter) %>% mutate(min_error = min(model_error)) %>%
  mutate(perc = model_error / min_error)
ggplot(dt_f_newest, aes(x = ind - 1, y = model_error, col = type)) + geom_point() + geom_line() + 
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) + scale_color_grey(name = "") +  
  scale_x_continuous(breaks = seq(0,15,1), name = "Number of covariates") + scale_y_continuous(name = "Test error") + theme_bw()

ggplot(dt_plot_f_p, aes(x = ind - 1, y = perc, col = type)) + geom_point() + geom_line() + 
  geom_abline(slope = 0,intercept = 1.01, col = "red", alpha = 0.4) + 
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) + scale_color_grey(name = "") +  
  scale_x_continuous(breaks = seq(0,15,1), name = "Number of covariates") + 
  scale_y_continuous(name = "RMSE over min RMSE") + theme_bw()

# Spatial figures
library(fields)
local.plot.field = function(field, mesh, col, xlim=c(-0.1,1), ylim=c(-0.1,1), ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  proj = inla.mesh.projector(mesh, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 plotting grid 
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), col = col,
             xlim = xlim, ylim = ylim, ...)  
}

local.plot.field = function(field, mesh, xlim=c(-0.1,1), ylim=c(-0.1,1), ...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  proj = inla.mesh.projector(mesh, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 plotting grid 
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj),
             xlim = xlim, ylim = ylim, ...)  
}

u_psi_m <- apply(u_psi, 1, mean)
u_tau_m <- apply(u_tau, 1, mean)

library('RColorBrewer')

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))   # make colors
r <- rf(64)

colorTable<- designer.colors(20, c( "#5ab4ac","white", "#d8b365"), 
                             x = c(min(u_psi_m), 0, max(u_psi_m)) / (max(u_psi_m) -min(u_psi_m)))
colorTable<- designer.colors(100, c( "#e66101","#fdb863","#b2abd2", "#5e3c99"), 
                             x = c(min(u_psi_m),quantile(u_psi_m, 0.25), 0, max(u_psi_m)) / (max(u_psi_m) -min(u_psi_m)))
colorTable2<- designer.colors(100, c( "#e66101","#fdb863","#b2abd2", "#5e3c99"), 
                              x = c(min(u_tau_m),quantile(u_tau_m, 0.25), 0, max(u_tau_m)) / (max(u_tau_m) -min(u_tau_m)))

max(u_psi_m)
min(u_psi_m)
max(u_tau_m)
min(u_tau_m)
length(u_psi_m)
u_psi_m
local.plot.field(field = u_psi_m, mesh)
points(map_scaled$long, map_scaled$lat, pch = 20, cex=.2)
points(coords[,1],coords[,2], pch = 20, col = "black", cex = .3)

local.plot.field(field = u_tau_m, mesh)
points(map_scaled$long, map_scaled$lat, pch = 20, cex=.2)
points(coords[,1],coords[,2], pch = 20, col = "black", cex = .3)

local.plot.field(field = u_psi_m, mesh, col = colorTable)
points(map_scaled$long, map_scaled$lat, pch = 20, cex=.2)
points(coords[,1],coords[,2], pch = 20, col = "black", cex = .3)

local.plot.field(field = u_tau_m, mesh, col = colorTable2)
points(map_scaled$long, map_scaled$lat, pch = 20, cex=.2)
points(coords[,1],coords[,2], pch = 20, col = "black", cex = .3)


################################ Anderson-darling test ################################################################################
library(goftest)

cdfGEV2 <- function(data, m, s, xi, d, t){
  t <- ( 1 + xi*( (data-(m+d*0.008*(t-1975)))/s ))^(-1/xi)
  res <- exp(-t)
  return(res)
}

fnB <- function(theta,y,t) {
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
  # Bæta við normal prior sem er með sd = hálft delda 0, mu = 0
  sum_loglik <- 0
  for(i in 1:n){
    sum_loglik <- sum_loglik + dgev(y[i], loc = mu_it[i], scale = sigma_i, shape = xitheta, log = TRUE)
  }
  res = -sum_loglik - 
    ((alpha - alp3)*log(xitheta + 0.5) + (beta-1)*log(0.5 - xitheta) + (theta[3]-b3)/sig3 - exp((theta[3]-b3)/sig3)  )  + 
    0.5*gamma^2/((0.5*delta_0)^2)
  return(res)
}

xi_trans <- function(phi){
  xi00 <- 0
  alp3 <- 0.8 # c_phi
  alpha <- 4
  beta <- 4
  sig3 = (log(1 - (xi00+0.5)^alp3))*(1 - (xi00+0.5)^alp3)*(-(alp3)^(-1)*(xi00+0.5)^(-alp3+1)) # b_phi
  b3 <- -sig3*(log(-log(1-0.5^alp3)))
  xihat <- ((1-(exp(-exp((phi - b3)/sig3))))^(1/alp3) - 0.5)
  return(xihat)
}

eta <- x_mat[1:(3*N),]
data$pB <- apply(exp(eta[1:N,]),1,mean)
data$tB <- apply(exp(eta[1:N,]+eta[1:N + N,]),1,mean)
data$xB <- apply(xi_trans(eta[1:N + 2*N,]), 1, mean)

nu <- x_mat[-c(1:(3*N)),]
beta_psi <- nu[1:N_psi,]
u_psi <- nu[1:N_colsA + N_psi,]
beta_tau <- nu[1:N_tau + N_colsA + N_psi,]
u_tau <- nu[1:N_colsA + N_tau + N_colsA + N_psi,]
beta_xi <- nu[1:N_xi + N_colsA + N_tau + N_colsA + N_psi,]

pred_psi <- X_psi %*% beta_psi + A_mat %*% u_psi 
pred_tau <- X_tau %*% beta_tau + A_mat %*% u_tau
pred_xi <- X_kappa %*% beta_xi 
pred_gamma <- X_gamma %*% beta_gamma

dt_mles$pBB <- apply(exp(pred_psi),1,mean)
dt_mles$tBB <- apply(exp(pred_psi+pred_tau),1,mean)
dt_mles$xBB <- apply(xi_trans(pred_xi),1,mean)
dt_mles$gBB <- apply(d(pred_gamma),1,mean)

p_vals_ad <- data.table()
adTests <- list()
for(i in 1:N){
  item <- dt_mles[i,]
  currentF <- allAmax %>% filter(Station == item$Station) %>% dplyr::select(Flow, WaterYear) %>% unique()
  adTestB <- tryCatch({
    ad.test(x = as.matrix(currentF$Flow), null = "cdfGEV2", m = item$pBB, s = item$tBB, xi = item$xBB, d = item$gBB, t=as.matrix(currentF$WaterYear))
  }, error = function(e){return(0)})
  adTestMle <- ad.test(x = as.matrix(currentF$Flow), null = "cdfGEV2", m = exp(item$psi), s = exp(item$tau)*exp(item$psi), 
                       xi = xi_trans(item$kappa), d = d(item$gamma), t = as.matrix(currentF$WaterYear))
  adTests[[i]] <- adTestB
  #ksB <- ks.test(currentF$Flow,y = "cdfGEV2", m = item$pB, s = item$tB, xi = item$xB)
  #ksMLE <- ks.test(currentF$Flow,y = "cdfGEV2", m = item$mu, s = item$sigma, xi = item$xi)
  if(typeof(adTestB) == "double"){
    p_vals_ad <- rbind(p_vals_ad, data.table(p_B = adTestB,
                                             #p_B_ks = ksB$p.value,
                                             p_MLE = adTestMle$p.value
                                             #p_MLE_ks = ksMLE$p.value
    ))
  }else{
    p_vals_ad <- rbind(p_vals_ad, data.table(p_B = adTestB$p.value ,
                                             #p_B_ks = ksB$p.value,
                                             p_MLE = adTestMle$p.value
                                             #p_MLE_ks = ksMLE$p.value
    ))
  }
}
ggplot(p_vals_ad, aes(p_B_ks)) + geom_histogram(breaks = seq(0,1,0.05), col = "black") + 
  labs(y="Number of stations",x="P-value") + theme_bw() + scale_fill_discrete(guide=FALSE) 

ggplot(p_vals_ad, aes(p_MLE)) + geom_histogram(breaks = seq(0,1,0.05), col = "black") + 
  labs(y="Number of stations",x="P-value") + theme_bw() + scale_fill_discrete(guide=FALSE) 

ggplot(p_vals_ad, aes(p_B)) + geom_histogram()

p_vals_ad$n <- dplyr::select(inner_join(unique(allAmax) %>% group_by(Station) %>% summarise(n=n()), fd2),n)

ggplot(p_vals_ad, aes(n,p_B, col = signif)) + geom_point()

sum(p_vals_ad$p_MLE<0.05)


