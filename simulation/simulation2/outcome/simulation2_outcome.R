rm(list=ls())

# Import functions
source(sub("simulation/simulation2/outcome$", "functions/FI_functions.R", getwd()))
source(sub("simulation/simulation2/outcome$", "functions/PPMA_functions.R", getwd()))

# Import libraries
library(doParallel)
library(doRNG)
library(bestglm)
library(ggplot2)

# common setups
n=1000; sig_level = 0.1; B=1000; M=100
mu_x1 <- 0.5; mu_x2 <- 0.5; sigma_x1 <- sqrt(0.5); sigma_x2 <- sqrt(0.5); sigma_y <- sigma_x1; mu_y <- mu_x1
identifiability <- 1

# outcome_model <- "normal"
# outcome_model <- "gamma"
# outcome_model <- "mixture"
# outcome_model <- "log_normal"
# outcome_model <- "heterogeneous"
# 
# rho <- 0.5
# rho <- 0.8
# 
# phi <- c(0, 2, 0, 0)
# phi <- c(0, 1, 0, 1)
# phi <- c(0, 0, 0, 2)

# Simulation function
simulation2_outcome <- function(outcome_model, rho, phi) {
  alp <- phi[-length(phi)]; beta <- phi[length(phi)]
  sigma_e <- sqrt(1-rho^2)*sigma_y
  
  kappa1 <- sqrt(rho^2*sigma_y^2/(sigma_x1^2+identifiability^2*sigma_x2^2))
  kappa2 <- identifiability*kappa1
  kappa2/kappa1
  
  sim <- function(it) {
    if (outcome_model=="normal") {
      kappa0 <- mu_y - kappa1*mu_x1 - kappa2*mu_x2
      kap <- c(kappa0, kappa1, kappa2)
      x1 <- rnorm(n, mean = mu_x1, sd = sigma_x1)
      x2 <- rnorm(n, mean = mu_x2, sd = sigma_x2)
      design_out <- cbind(1, x1, x2)
      mu <- design_out %*% kap
      e <- rnorm(n=n, mean = 0, sd=sigma_e)
    } else if (outcome_model=="gamma") {
      kappa0 <- mu_y - kappa1*mu_x1 - kappa2*mu_x2
      kap <- c(kappa0, kappa1, kappa2)
      theta_x1 <- sigma_x1^2/mu_x1; k_x1 <- mu_x1^2/sigma_x1^2
      theta_x2 <- sigma_x2^2/mu_x2; k_x2 <- mu_x2^2/sigma_x2^2
      x1 <- rgamma(n=n, shape=k_x1, scale=theta_x1)
      x2 <- rgamma(n=n, shape=k_x2, scale=theta_x2)
      design_out <- cbind(1, x1, x2)
      mu <- design_out %*% kap
      e <- rnorm(n=n, mean = 0, sd=sigma_e)
    } else if (outcome_model=="mixture") {
      kappa0 <- mu_y - kappa1*mu_x1 - kappa2*mu_x2
      kap <- c(kappa0, kappa1, kappa2)
      lambda <- 0.5
      r <- 3
      s_e <- sigma_e/sqrt(r)
      mu_e <- sqrt(sigma_e^2 - s_e^2)
      x1 <- rnorm(n, mean = mu_x1, sd = sigma_x1)
      x2 <- rnorm(n, mean = mu_x2, sd = sigma_x2)
      design_out <- cbind(1, x1, x2)
      mu <- design_out %*% kap
      mix <- rbinom(n=n, size=1, prob=lambda)
      e1 <- rnorm(n=n, mean=mu_e, sd=s_e); e2 <- rnorm(n=n, mean=-mu_e, sd=s_e)
      e <- mix*e1 + (1-mix)*e2
    } else if (outcome_model=="log_normal") {
      kappa0 <- mu_y - kappa1*mu_x1 - kappa2*mu_x2 - 1
      kap <- c(kappa0, kappa1, kappa2)
      theta_e <- log((1-rho^2)*sigma_y^2 + 1)
      x1 <- rnorm(n, mean = mu_x1, sd = sigma_x1)
      x2 <- rnorm(n, mean = mu_x2, sd = sigma_x2)
      design_out <- cbind(1, x1, x2)
      mu <- design_out %*% kap
      e <- rlnorm(n=n, meanlog=-theta_e/2, sdlog=sqrt(theta_e))
    } else if (outcome_model=="heterogeneous") {
      kappa0 <- mu_y - kappa1*mu_x1 - kappa2*mu_x2
      kap <- c(kappa0, kappa1, kappa2)
      x1 <- rnorm(n, mean = mu_x1, sd = sigma_x1)
      x2 <- rnorm(n, mean = mu_x2, sd = sigma_x2)
      design_out <- cbind(1, x1, x2)
      mu <- design_out %*% kap
      if (rho==0.5) {
        a <- 0.2
      } else if (rho==0.8) {
        a <- 0.1
      }
      e <- rnorm(n=n, mean=0, sd=sqrt(a*(1+x1^2)))
    }
    y <- mu + e
    design_res <- cbind(1, x1, x2, y)
    lin_res <- design_res %*% phi
    prob <- exp(lin_res)/(1+exp(lin_res))
    R <- rbinom(n=n, size=1, prob=prob); n1 <- sum(R); n0 <- n-n1
    sum(R)/n
    mean(y)
    mean(x1)
    mean(x2)
    var(y)
    var(x1)
    var(x2)
    cor(y, mu)
    hist(e, breaks=100)
    hist(x1, breaks=100)
    hist(x2, breaks=100)
    
    data <- data.frame(x1=x1, x2=x2, R=R, y=y)
    full <- mean(y)
    
    # CC analysis
    CC <- mean(data[data$R==1, ]$y)
    var_CC <- var(data[data$R==1, ]$y)/nrow(data[data$R==1, ])
    lower.bound_CC <- CC - sqrt(var_CC) * qnorm(0.975)
    upper.bound_CC <- CC + sqrt(var_CC) * qnorm(0.975)
    CI_CC <- c(lower.bound_CC, upper.bound_CC)
    
    prop <- setdiff(colnames(data), c("y", "R"))
    data_response <- data; data_response$y <- NULL; data_response <- data_response %>% rename(y=R)
    MARmodel1 <- glm(y~1, data=data_response, family=binomial(link="logit"))
    MARmodel2 <- glm(y~x1, data=data_response, family=binomial(link="logit"))
    MARmodel3 <- glm(y~x2, data=data_response, family=binomial(link="logit"))
    MARmodel4 <- glm(y~x1+x2, data=data_response, family=binomial(link="logit"))
    BestModels <- list(MARmodel1, MARmodel2, MARmodel3, MARmodel4)
    BestModel <- BestModels[[which.min(unlist(lapply(BestModels, BIC)))]]
    bestglm_out <- list(BestModel=BestModel)
    res_MAR <- colnames(bestglm_out$BestModel$model)[-1]
    VE_MAR_out <- VE_MAR(data=data, XX_names = res_MAR, Y_name = "y", R_name = "R")
    
    res_xy <- list(NULL, c("x1"), c("x2"), c("x1", "x2"))
    FI_out_MNAR_list <- vector(mode="list", length = length(res_xy))
    VE_FI_MNAR_list <- vector(mode="list", length = length(res_xy))
    BIC_MNAR <- vector(length = length(res_xy))
    p_values <- rep(1, length(res_xy))
    
    for (j in 1:length(res_xy)) {
      success <- FALSE
      attempt <- 1
      max_attempts <- 5  # Maximum number of attempts to retry
      while (!success && attempt <= max_attempts) {
        tryCatch({
          # Execute FI function
          FI_out_MNAR_list[[j]] <- FI(pop=data, prop=prop, response=c(res_xy[[j]], "y"), phi_curr=rep(0, length(c(res_xy[[j]], "y"))+1), M=M, mixture=TRUE)
          BIC_MNAR[j] <- FI_out_MNAR_list[[j]]$BIC_full
          
          # Execute VE_NMAR function
          VE_FI_MNAR_list[[j]] <- VE_NMAR(FI_out_MNAR_list[[j]])
          p_values[j] <- VE_FI_MNAR_list[[j]]$summary_table["y", "Pr(>|z|)"]
          
          # If no error occurs, set success to TRUE
          success <- TRUE
        }, error = function(e) {
          attempt <<- attempt + 1
        })
      }
    }
    
    best_MNAR <- which.min(BIC_MNAR)
    if (p_values[best_MNAR] < sig_level) {
      New_out <- FI_out_MNAR_list[[best_MNAR]]$mu.hat
      selected_model <- paste0("NMAR", best_MNAR, collapse = "")
    } else {
      New_out <- VE_MAR_out$mean_PS_MAR
      selected_model <- "MAR"
    }
    
    ppm_fit <- lm(y~x1+x2, data=data[R==1, ])
    x_ppm <- predict(ppm_fit, newdata=data)
    y_ppm <- y; y_ppm[R==0] <- NA
    PPM_0 <- mle(x_ppm, y_ppm, 0)
    PPM_1 <- mle(x_ppm, y_ppm, 1)
    PPM_inf <- mle(x_ppm, y_ppm, Inf)
  
    if (identical(phi, c(0, 2, 0, 0)) | identical(phi, c(0, 1, 0, 1))) {
      FI_c <- FI_out_MNAR_list[[2]]
    } else if (identical(phi, c(0, 0, 0, 2))) {
      FI_c <- FI_out_MNAR_list[[1]]
    }
    mu_estimate <- c(CC, PPM_0$muY, PPM_1$muY, PPM_inf$muY, VE_MAR_out$mean_PS_MAR, FI_c$mu.hat, New_out)
    names(mu_estimate) <- c("CC", "PPM_0", "PPM_1", "PPM_inf", "MAR", "FI_c", "FI_m")
    # var_mu_estimate <- c(var_CC, PPM_0$muYvar, PPM_1$muYvar, PPM_inf$muYvar, VE_MAR_out$variance_mu, VE_FI_MNAR_list[[2]]$variance_mu)
    # names(var_mu_estimate) <- c("CC", "PPM_0", "PPM_1", "PPM_inf", "MAR", "FI_c")
    
    # beta_estimate <- c(FI_out_MNAR_list[[2]]$phi.hat[length(FI_out_MNAR_list[[2]]$phi.hat)])
    # var_beta_estimate <- c(diag(VE_FI_MNAR_list[[2]]$variance_phi)[length(diag(VE_FI_MNAR_list[[2]]$variance_phi))])
     
    # out <- list(mu_estimate=mu_estimate,
    #             var_mu_estimate=var_mu_estimate,
    #             beta_estimate=beta_estimate,
    #             var_beta_estimate=var_beta_estimate,
    #             p_values=p_values,
    #             selected_model=selected_model,
    #             full=full)
    out <- list(mu_estimate=mu_estimate,
                p_values=p_values,
                selected_model=selected_model,
                full=full)
    return(out)
  }
  
  # parallel computation
  t<-proc.time()
  stopimplicitCluster2 <- function() {
    options <- doParallel:::.options
    if (exists(".revoDoParCluster", where = options) &&
        !is.null(get(".revoDoParCluster", envir = options))) {
      stopCluster(get(".revoDoParCluster", envir = options))
      remove(".revoDoParCluster", envir = options)
    }
  }
  
  registerDoParallel(max(detectCores()-1,1))
  registerDoRNG(1)
  Sys.time()
  
  ptime <- system.time({
    output1 <- foreach(it = 1:B, .export=ls(envir=parent.frame()), .packages=loadedNamespaces(), .errorhandling = "pass") %dopar% {
      return(sim(it))
      invisible(gc())
      invisible(gc())
    }
  })[3]
  stopimplicitCluster2()
  tictoc <- proc.time() - t 
  
  return(output1)
}

tic <- proc.time()

# Normal Outcome Model
normal_rho0.5_MAR <- simulation2_outcome(outcome_model="normal", phi=c(0, 2, 0, 0), rho=0.5)
normal_rho0.5_NMAR1 <- simulation2_outcome(outcome_model="normal", phi=c(0, 1, 0, 1), rho=0.5)
normal_rho0.5_NMARinf <- simulation2_outcome(outcome_model="normal", phi=c(0, 0, 0, 2), rho=0.5)
save(normal_rho0.5_MAR, file = "normal_rho0.5_MAR.Rdata")
load(file = "normal_rho0.5_MAR.Rdata")
save(normal_rho0.5_NMAR1, file = "normal_rho0.5_NMAR1.Rdata")
load(file = "normal_rho0.5_NMAR1.Rdata")
save(normal_rho0.5_NMARinf, file = "normal_rho0.5_NMARinf.Rdata")
load(file = "normal_rho0.5_NMARinf.Rdata")
see_result(normal_rho0.5_MAR)
see_result(normal_rho0.5_NMAR1)
see_result(normal_rho0.5_NMARinf)

normal_rho0.8_MAR <- simulation2_outcome(outcome_model="normal", phi=c(0, 2, 0, 0), rho=0.8)
normal_rho0.8_NMAR1 <- simulation2_outcome(outcome_model="normal", phi=c(0, 1, 0, 1), rho=0.8)
normal_rho0.8_NMARinf <- simulation2_outcome(outcome_model="normal", phi=c(0, 0, 0, 2), rho=0.8)
save(normal_rho0.8_MAR, file = "normal_rho0.8_MAR.Rdata")
load(file = "normal_rho0.8_MAR.Rdata")
save(normal_rho0.8_NMAR1, file = "normal_rho0.8_NMAR1.Rdata")
load(file = "normal_rho0.8_NMAR1.Rdata")
save(normal_rho0.8_NMARinf, file = "normal_rho0.8_NMARinf.Rdata")
load(file = "normal_rho0.8_NMARinf.Rdata")
see_result(normal_rho0.8_MAR)
see_result(normal_rho0.8_NMAR1)
see_result(normal_rho0.8_NMARinf)


# Gamma Outcome Model
gamma_rho0.5_MAR <- simulation2_outcome(outcome_model="gamma", phi=c(0, 2, 0, 0), rho=0.5)
gamma_rho0.5_NMAR1 <- simulation2_outcome(outcome_model="gamma", phi=c(0, 1, 0, 1), rho=0.5)
gamma_rho0.5_NMARinf <- simulation2_outcome(outcome_model="gamma", phi=c(0, 0, 0, 2), rho=0.5)
save(gamma_rho0.5_MAR, file = "gamma_rho0.5_MAR.Rdata")
load(file = "gamma_rho0.5_MAR.Rdata")
save(gamma_rho0.5_NMAR1, file = "gamma_rho0.5_NMAR1.Rdata")
load(file = "gamma_rho0.5_NMAR1.Rdata")
save(gamma_rho0.5_NMARinf, file = "gamma_rho0.5_NMARinf.Rdata")
load(file = "gamma_rho0.5_NMARinf.Rdata")
see_result(gamma_rho0.5_MAR)
see_result(gamma_rho0.5_NMAR1)
see_result(gamma_rho0.5_NMARinf)

gamma_rho0.8_MAR <- simulation2_outcome(outcome_model="gamma", phi=c(0, 2, 0, 0), rho=0.8)
gamma_rho0.8_NMAR1 <- simulation2_outcome(outcome_model="gamma", phi=c(0, 1, 0, 1), rho=0.8)
gamma_rho0.8_NMARinf <- simulation2_outcome(outcome_model="gamma", phi=c(0, 0, 0, 2), rho=0.8)
save(gamma_rho0.8_MAR, file = "gamma_rho0.8_MAR.Rdata")
load(file = "gamma_rho0.8_MAR.Rdata")
save(gamma_rho0.8_NMAR1, file = "gamma_rho0.8_NMAR1.Rdata")
load(file = "gamma_rho0.8_NMAR1.Rdata")
save(gamma_rho0.8_NMARinf, file = "gamma_rho0.8_NMARinf.Rdata")
load(file = "gamma_rho0.8_NMARinf.Rdata")
see_result(gamma_rho0.8_MAR)
see_result(gamma_rho0.8_NMAR1)
see_result(gamma_rho0.8_NMARinf)


# Mixture Outcome Model
mixture_rho0.5_MAR <- simulation2_outcome(outcome_model="mixture", phi=c(0, 2, 0, 0), rho=0.5)
mixture_rho0.5_NMAR1 <- simulation2_outcome(outcome_model="mixture", phi=c(0, 1, 0, 1), rho=0.5)
mixture_rho0.5_NMARinf <- simulation2_outcome(outcome_model="mixture", phi=c(0, 0, 0, 2), rho=0.5)
save(mixture_rho0.5_MAR, file = "mixture_rho0.5_MAR.Rdata")
load(file = "mixture_rho0.5_MAR.Rdata")
save(mixture_rho0.5_NMAR1, file = "mixture_rho0.5_NMAR1.Rdata")
load(file = "mixture_rho0.5_NMAR1.Rdata")
save(mixture_rho0.5_NMARinf, file = "mixture_rho0.5_NMARinf.Rdata")
load(file = "mixture_rho0.5_NMARinf.Rdata")
see_result(mixture_rho0.5_MAR)
see_result(mixture_rho0.5_NMAR1)
see_result(mixture_rho0.5_NMARinf)

mixture_rho0.8_MAR <- simulation2_outcome(outcome_model="mixture", phi=c(0, 2, 0, 0), rho=0.8)
mixture_rho0.8_NMAR1 <- simulation2_outcome(outcome_model="mixture", phi=c(0, 1, 0, 1), rho=0.8)
mixture_rho0.8_NMARinf <- simulation2_outcome(outcome_model="mixture", phi=c(0, 0, 0, 2), rho=0.8)
save(mixture_rho0.8_MAR, file = "mixture_rho0.8_MAR.Rdata")
load(file = "mixture_rho0.8_MAR.Rdata")
save(mixture_rho0.8_NMAR1, file = "mixture_rho0.8_NMAR1.Rdata")
load(file = "mixture_rho0.8_NMAR1.Rdata")
save(mixture_rho0.8_NMARinf, file = "mixture_rho0.8_NMARinf.Rdata")
load(file = "mixture_rho0.8_NMARinf.Rdata")
see_result(mixture_rho0.8_MAR)
see_result(mixture_rho0.8_NMAR1)
see_result(mixture_rho0.8_NMARinf)


# Log-normal Outcome Model
log_normal_rho0.5_MAR <- simulation2_outcome(outcome_model="log_normal", phi=c(0, 2, 0, 0), rho=0.5)
log_normal_rho0.5_NMAR1 <- simulation2_outcome(outcome_model="log_normal", phi=c(0, 1, 0, 1), rho=0.5)
log_normal_rho0.5_NMARinf <- simulation2_outcome(outcome_model="log_normal", phi=c(0, 0, 0, 2), rho=0.5)
save(log_normal_rho0.5_MAR, file = "log_normal_rho0.5_MAR.Rdata")
load(file = "log_normal_rho0.5_MAR.Rdata")
save(log_normal_rho0.5_NMAR1, file = "log_normal_rho0.5_NMAR1.Rdata")
load(file = "log_normal_rho0.5_NMAR1.Rdata")
save(log_normal_rho0.5_NMARinf, file = "log_normal_rho0.5_NMARinf.Rdata")
load(file = "log_normal_rho0.5_NMARinf.Rdata")
see_result(log_normal_rho0.5_MAR)
see_result(log_normal_rho0.5_NMAR1)
see_result(log_normal_rho0.5_NMARinf)

log_normal_rho0.8_MAR <- simulation2_outcome(outcome_model="log_normal", phi=c(0, 2, 0, 0), rho=0.8)
log_normal_rho0.8_NMAR1 <- simulation2_outcome(outcome_model="log_normal", phi=c(0, 1, 0, 1), rho=0.8)
log_normal_rho0.8_NMARinf <- simulation2_outcome(outcome_model="log_normal", phi=c(0, 0, 0, 2), rho=0.8)
save(log_normal_rho0.8_MAR, file = "log_normal_rho0.8_MAR.Rdata")
load(file = "log_normal_rho0.8_MAR.Rdata")
save(log_normal_rho0.8_NMAR1, file = "log_normal_rho0.8_NMAR1.Rdata")
load(file = "log_normal_rho0.8_NMAR1.Rdata")
save(log_normal_rho0.8_NMARinf, file = "log_normal_rho0.8_NMARinf.Rdata")
load(file = "log_normal_rho0.8_NMARinf.Rdata")
see_result(log_normal_rho0.8_MAR)
see_result(log_normal_rho0.8_NMAR1)
see_result(log_normal_rho0.8_NMARinf)


# Heterogeneous Outcome Model
heterogeneous_rho0.5_MAR <- simulation2_outcome(outcome_model="heterogeneous", phi=c(0, 2, 0, 0), rho=0.5)
heterogeneous_rho0.5_NMAR1 <- simulation2_outcome(outcome_model="heterogeneous", phi=c(0, 1, 0, 1), rho=0.5)
heterogeneous_rho0.5_NMARinf <- simulation2_outcome(outcome_model="heterogeneous", phi=c(0, 0, 0, 2), rho=0.5)
save(heterogeneous_rho0.5_MAR, file = "heterogeneous_rho0.5_MAR.Rdata")
load(file = "heterogeneous_rho0.5_MAR.Rdata")
save(heterogeneous_rho0.5_NMAR1, file = "heterogeneous_rho0.5_NMAR1.Rdata")
load(file = "heterogeneous_rho0.5_NMAR1.Rdata")
save(heterogeneous_rho0.5_NMARinf, file = "heterogeneous_rho0.5_NMARinf.Rdata")
load(file = "heterogeneous_rho0.5_NMARinf.Rdata")
see_result(heterogeneous_rho0.5_MAR)
see_result(heterogeneous_rho0.5_NMAR1)
see_result(heterogeneous_rho0.5_NMARinf)

heterogeneous_rho0.8_MAR <- simulation2_outcome(outcome_model="heterogeneous", phi=c(0, 2, 0, 0), rho=0.8)
heterogeneous_rho0.8_NMAR1 <- simulation2_outcome(outcome_model="heterogeneous", phi=c(0, 1, 0, 1), rho=0.8)
heterogeneous_rho0.8_NMARinf <- simulation2_outcome(outcome_model="heterogeneous", phi=c(0, 0, 0, 2), rho=0.8)
save(heterogeneous_rho0.8_MAR, file = "heterogeneous_rho0.8_MAR.Rdata")
load(file = "heterogeneous_rho0.8_MAR.Rdata")
save(heterogeneous_rho0.8_NMAR1, file = "heterogeneous_rho0.8_NMAR1.Rdata")
load(file = "heterogeneous_rho0.8_NMAR1.Rdata")
save(heterogeneous_rho0.8_NMARinf, file = "heterogeneous_rho0.8_NMARinf.Rdata")
load(file = "heterogeneous_rho0.8_NMARinf.Rdata")
see_result(heterogeneous_rho0.8_MAR)
see_result(heterogeneous_rho0.8_NMAR1)
see_result(heterogeneous_rho0.8_NMARinf)

toc <- proc.time()
toc - tic

#### Output Organization ####
see_result <- function(output1) {
  mu <- mean(unlist(do.call(rbind, lapply(output1, FUN = function(x) x$full))))
  mu_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$mu_estimate))
  
  point_summary <- function(true, point_estimate) {
    RB <- apply((point_estimate - true)/true*100, MARGIN = 2, FUN = mean)
    # SE <- apply(point_estimate, MARGIN = 2, FUN = sd)
    RMSE <- sqrt(apply((point_estimate - true)^2, MARGIN = 2, FUN = mean))
    return(data.frame(RB=round(RB, 1), RMSE=round(RMSE, 3)))
  }
  
  point_sum_mu <- point_summary(mu, mu_estimates)
  model_selection_result <- summary(as.factor(unlist(lapply(output1, FUN = function(x) x$selected_model))))
  
  return(list(point_sum_mu=point_sum_mu, model_selection_result=model_selection_result))
}

point_summary <- function(output1) {
  mu <- mean(unlist(do.call(rbind, lapply(output1, FUN = function(x) x$full))))
  mu_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$mu_estimate))
  RB <- apply((mu_estimates - mu)/mu*100, MARGIN = 2, FUN = mean)
  RMSE <- sqrt(apply((mu_estimates - mu)^2, MARGIN = 2, FUN = mean))
  return(data.frame(RB=round(RB, 1), RMSE=round(RMSE, 3)))
}

organize_result <- function(output1, output2, output3) {
  return(cbind(point_summary(output1), point_summary(output2), point_summary(output3)))
}

outcome_result <- rbind(cbind(organize_result(normal_rho0.5_MAR, normal_rho0.5_NMAR1, normal_rho0.5_NMARinf), organize_result(normal_rho0.8_MAR, normal_rho0.8_NMAR1, normal_rho0.8_NMARinf)),
                        cbind(organize_result(gamma_rho0.5_MAR, gamma_rho0.5_NMAR1, gamma_rho0.5_NMARinf), organize_result(gamma_rho0.8_MAR, gamma_rho0.8_NMAR1, gamma_rho0.8_NMARinf)),
                        cbind(organize_result(log_normal_rho0.5_MAR, log_normal_rho0.5_NMAR1, log_normal_rho0.5_NMARinf), organize_result(log_normal_rho0.8_MAR, log_normal_rho0.8_NMAR1, log_normal_rho0.8_NMARinf)),
                        cbind(organize_result(heterogeneous_rho0.5_MAR, heterogeneous_rho0.5_NMAR1, heterogeneous_rho0.5_NMARinf), organize_result(heterogeneous_rho0.8_MAR, heterogeneous_rho0.8_NMAR1, heterogeneous_rho0.8_NMARinf)))
write.csv(outcome_result, file = "outcome_result.csv")


# point estimation performance with model selection
result_m <- cbind(rbind(normal_rho0.5$sim_summary_mu_real, gamma_rho0.5$sim_summary_mu_real, mixture_rho0.5$sim_summary_mu_real, log_normal_rho0.5$sim_summary_mu_real, heterogeneous_rho0.5$sim_summary_mu_real),
                  rbind(normal_rho0.8$sim_summary_mu_real, gamma_rho0.8$sim_summary_mu_real, mixture_rho0.8$sim_summary_mu_real, log_normal_rho0.8$sim_summary_mu_real, heterogeneous_rho0.8$sim_summary_mu_real))
write.csv(result_m, file = "result_m.csv")

# model selection result
model_selection_result <- bind_rows(
  normal_rho0.5 = normal_rho0.5$selected_model,
  gamma_rho0.5 = gamma_rho0.5$selected_model,
  mixture_rho0.5 = mixture_rho0.5$selected_model,
  lognormal_rho0.5 = log_normal_rho0.5$selected_model,
  hetero_rho0.5 = heterogeneous_rho0.5$selected_model,
  normal_rho0.8 = normal_rho0.8$selected_model,
  gamma_rho0.8 = gamma_rho0.8$selected_model,
  mixture_rho0.8 = mixture_rho0.8$selected_model,
  lognormal_rho0.8 = log_normal_rho0.8$selected_model,
  hetero_rho0.8 = heterogeneous_rho0.8$selected_model,
  .id = "source"
)
model_selection_result[is.na(model_selection_result)] <- 0
model_selection_result <- model_selection_result %>% mutate(total=apply(model_selection_result[, -1], MARGIN=1, FUN=sum))
write.csv(model_selection_result, file = "model_selection_result.csv")

# point estimation performance with correctly specified response model
rho0.5_mu <- rbind(normal_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   gamma_rho0.5$sim_summary_mu[c("CC", "PPM_inf", "FI_c"), ],
                   mixture_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   log_normal_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   heterogeneous_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ])
rho0.8_mu <- rbind(normal_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   gamma_rho0.8$sim_summary_mu[c("CC", "PPM_inf", "FI_c"), ],
                   mixture_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   log_normal_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   heterogeneous_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ])
result_c_mu <- cbind(rho0.5_mu, rho0.8_mu)
result_c <- rbind(result_c_mu)
write.csv(result_c, file = "result_c.csv")
