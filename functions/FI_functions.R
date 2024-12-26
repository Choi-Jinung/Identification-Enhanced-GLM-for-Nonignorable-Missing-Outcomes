library(dplyr)
library(nleqslv)
library(flexmix)

# Arguments
## 1. pop: Full sample. Response variable should be denoted by "y".
## 2. prop: predictor variables of the proposal distribution
## 3. response: predictor variables to be included in the response model. "y" should be included explicitly.
## 4. M: imputation size
## 5. phi_curr: initial value of the propensity score model parameters
## 6. binary: If true, the response variable is binary and logistic outcome model will be fitted.
FI <- function(pop, prop, response, M, phi_curr, mixture=F) {
  pop$ID <- 1:nrow(pop)
  res <- pop[pop$R==1, ]; nonres <- pop[pop$R==0, ]; nonres$y <- NULL
  n_res <- nrow(res); n_nonres <- nrow(nonres); n_pop <- nrow(pop)
  prop_formula <- as.formula(paste("y", "~", paste(prop, collapse = "+")))
  
  binary <- all(res[, "y"] %in% c(0, 1))
  if (binary==F) {
    if (mixture==F) {
      # Fit the respondents' data with a Gaussian regression model 
      normal_fit <- lm(prop_formula, data=res)
      normal_fit <- step(normal_fit)
      best_model <- normal_fit
      best_model_name <- "normal"
      prop <- colnames(best_model$model)[-1]
      kappa <- best_model$coefficients; sig2 <- sum(best_model$residuals^2)/best_model$df.residual  
      gamma.hat <- c(kappa, sqrt(sig2))
      
      # Generate imputed values
      nonres_aug <- nonres[rep(1:n_nonres, each=M), ]
      fitted_values <- predict(best_model, newdata = nonres_aug)
      nonres_aug$y <- rnorm(n=n_nonres*M, mean=fitted_values, sd=sqrt(sig2))
    } else {
      # Fit the respondents' data with a Mixture Gaussian regression model 
      kk <- 2 # number of mixture components
      mixture_fit <- flexmix(prop_formula, data = res, k = kk)
      best_model <- mixture_fit
      best_model_name <- "mixture_normal"
      if (ncol(parameters(mixture_fit))==2) {
        kappa1 <- parameters(best_model)[1:(nrow(parameters(best_model))-1), 1]
        kappa2 <- parameters(best_model)[1:(nrow(parameters(best_model))-1), 2]
        sigma1 <- parameters(best_model)[nrow(parameters(best_model)), 1]
        sigma2 <- parameters(best_model)[nrow(parameters(best_model)), 2]
        w <- best_model@prior
      } else if (ncol(parameters(mixture_fit))==1) {
        kappa1 <- parameters(best_model)[1:(nrow(parameters(best_model))-1), 1]
        sigma1 <- parameters(best_model)[nrow(parameters(best_model)), 1]
        kappa2 <- kappa1
        sigma2 <- sigma1
        w <- c(0.5, 0.5)
      }
      mix <- t(rmultinom(n=n_nonres*M, size=1, w))
      gamma.hat <- c(w[1], kappa1, sigma1, kappa2, sigma2)
      
      # Generate imputed values
      nonres_aug <- nonres[rep(1:n_nonres, each=M), ]
      nonres_aug$y <- 0
      fitted_value1 <- as.matrix(cbind(1, nonres_aug[, prop]))%*%kappa1
      fitted_value2 <- as.matrix(cbind(1, nonres_aug[, prop]))%*%kappa2
      fitted_values <- list(fitted_value1, fitted_value2)
      sigmas <- c(sigma1, sigma2)
      for (k in 1:kk) {
        nonres_aug$y <- nonres_aug$y + mix[, k]*rnorm(n=n_nonres*M, mean=fitted_values[[k]], sd=sigmas[k])
      }
    }
  } else {
    # Fit the respondents' data with a Logistic regression model 
    logistic_fit <- glm(prop_formula, data=res, family = binomial(link="logit"))
    logistic_fit <- step(logistic_fit)
    # logistic_fit <- step(logistic_fit, k=log(nrow(res)))
    best_model <- logistic_fit
    best_model_name <- "logistic"
    prop <- colnames(best_model$model)[-1]
    gamma.hat <- best_model$coefficients
    
    # Generate imputed values
    # For the binary case, imputed values are 1 and 0.
    M <- 2  
    nonres_aug <- nonres[rep(1:n_nonres, each=M), ]
    nonres_aug$y <- rep(c(1,0), n_nonres) 
  }
  
  # Solve the mean score equation
  est.function <- function(phi) {
    alp <- phi[-length(phi)]; beta <- phi[length(phi)]
    if (binary == F) {
      nonres_aug$frac.w <- exp(-nonres_aug[, "y"]*beta)
      nonres_aug$frac.w <- nonres_aug %>% group_by(ID) %>% reframe(frac.w = frac.w/sum(frac.w)) %>% .[[2]]
    } else {
      w1 <- exp(as.matrix(cbind(1, nonres[, prop])) %*% gamma.hat - beta) / (1 + exp(as.matrix(cbind(1, nonres[, prop])) %*% gamma.hat - beta))
      w0 <- 1 - w1
      nonres_aug$frac.w <- NA
      nonres_aug$frac.w[seq(from=1, to=n_nonres*M, by=2)] <- w1
      nonres_aug$frac.w[seq(from=2, to=n_nonres*M, by=2)] <- w0
    }
    z.mat_res <- as.matrix(cbind(1, res[, response]))
    pi_res <- exp(as.matrix(cbind(1, res[, response]))%*%phi)/(1+exp(as.matrix(cbind(1, res[, response]))%*%phi))
    RHS1 <- apply(sweep(z.mat_res, MARGIN = 1, STATS = pi_res, FUN = "*"), MARGIN = 2, FUN = sum)
    pi_nonres <- exp(as.matrix(cbind(1, nonres_aug[, response]))%*%phi)/(1+exp(as.matrix(cbind(1, nonres_aug[, response]))%*%phi))
    z.mat_nonres <- cbind(1, nonres_aug[, response])
    LHS <- apply(z.mat_res, MARGIN = 2, FUN = sum)
    RHS2 <- apply(sweep(z.mat_nonres, MARGIN = 1, STATS = nonres_aug$frac.w*pi_nonres, FUN = "*"), MARGIN = 2, FUN = sum)
    RHS <- RHS1 + RHS2
    return(LHS-RHS)
  }
  nleqslv_result <- nleqslv(phi_curr, est.function); termcd <- nleqslv_result$termcd
  phi <- nleqslv_result$x; alp <- phi[-length(phi)]; beta <- phi[length(phi)]
  
  if (binary == F) {
    nonres_aug$frac.w <- exp(-nonres_aug[, "y"]*beta)
    nonres_aug$frac.w <- nonres_aug %>% group_by(ID) %>% reframe(frac.w = frac.w/sum(frac.w)) %>% .[[2]]
  } else {
    w1 <- exp(as.matrix(cbind(1, nonres[, prop])) %*% gamma.hat - beta) / (1 + exp(as.matrix(cbind(1, nonres[, prop])) %*% gamma.hat - beta))
    w0 <- 1 - w1
    nonres_aug$frac.w <- NA
    nonres_aug$frac.w[seq(from=1, to=n_nonres*M, by=2)] <- w1
    nonres_aug$frac.w[seq(from=2, to=n_nonres*M, by=2)] <- w0
  }
  
  full_obs_lik_func <- function(phi, data, best_model_name, response) {
    alp <- phi[-length(phi)]; beta <- phi[length(phi)]
    response_x <- response[-length(response)]
    h_alp <- as.matrix(cbind(1, data[, response_x])) %*% alp
    
    if (best_model_name=="normal") {
      mean_structure <- as.matrix(cbind(1, data[, prop])) %*% kappa
      mgf <- exp(-mean_structure*beta + sig2*(beta^2)/2)
    } else if (best_model_name=="mixture_normal") {
      mean_structure1 <- as.matrix(cbind(1, data[, prop])) %*% kappa1
      mean_structure2 <- as.matrix(cbind(1, data[, prop])) %*% kappa2
      mgf <- w[1]*exp(-mean_structure1*beta + (sigma1^2)*(beta^2)/2) + w[2]*exp(-mean_structure2*beta + (sigma2^2)*(beta^2)/2)
    } else if (best_model_name=="logistic") {
      lin_pred <- as.matrix(cbind(1, data[, prop])) %*% gamma.hat
      mgf <- (1 + exp(lin_pred - beta)) / (1 + exp(lin_pred))
    }
    prob <- 1/(1 + exp(-h_alp)*mgf)
    obs_lik <- sum(log(prob[data$R==1])) + sum(log(1-prob[data$R==0]))
    return(obs_lik)
  }
  full_obs_lik <- full_obs_lik_func(phi, pop, best_model_name, response)
  AIC_full <- -2*full_obs_lik + 2*(length(response)+1)
  BIC_full <- -2*full_obs_lik + (length(response)+1)*log(n_pop)
  
  res$frac.w <- 1
  data_aug <- rbind(res, nonres_aug)
  res_pi <- res
  res_pi$pi <- exp(as.matrix(cbind(1, res[, response]))%*%phi)/(1+exp(as.matrix(cbind(1, res[, response]))%*%phi))
  
  # Estimating function for the target parameter
  mu.ps <- weighted.mean(res_pi$y, w=1/res_pi$pi)
  mu.mi <- (weighted.mean(nonres_aug$y, w=nonres_aug$frac.w)*n_nonres + sum(res$y))/n_pop
  
  if (mixture==F) {
    outcome_model <- colnames(best_model$model)
  } else {
    outcome_model <- c("y", prop)
  }
  output1 <- list(gamma.hat=gamma.hat, phi.hat=phi, mu.hat=mu.ps, mu.mi=mu.mi, 
                  AIC_full=AIC_full, BIC_full=BIC_full, 
                  data=pop, data_aug=data_aug, res_pi=res_pi, nonres_aug=nonres_aug, nonres=nonres,
                  outcome_model = outcome_model, best_model_name=best_model_name, prop=prop, response_model = c("R", response),
                  best_model = best_model,
                  M=M,
                  termcd=termcd)
  return(output1)
}

# Arguments
## 1. FI_result: Output list from the "FI" function.
VE_NMAR <- function(FI_result) {
  best_model_name <- FI_result$best_model_name
  if (best_model_name=="normal") {
    gamma <- FI_result$gamma.hat[-length(FI_result$gamma.hat)]; sigma <- FI_result$gamma.hat[length(FI_result$gamma.hat)] 
    phi <- FI_result$phi.hat; alpha <- phi[-length(phi)]; beta <- phi[length(phi)] 
    mu <- FI_result$mu.hat
    n <- nrow(FI_result$data); n1 <- sum(FI_result$data$R); n0 <- n-n1
    outcome_model <- FI_result$outcome_model; response_model <- FI_result$response_model
    data_aug <- as.matrix(FI_result$data_aug)
    data1 <- data_aug[data_aug[, "R"]==1,  ]
    data0_aug <- data_aug[data_aug[, "R"]==0,  ]
    M <- FI_result$M
    
    # score function for the respondents outcome model
    score_1 <- function(data) {
      y <- data[, "y"]; design <- cbind(1, data[, outcome_model[-1]]); yhat <- design%*%gamma
      resid <- y - yhat
      dkappa <- sweep(design, MARGIN = 1, resid/sigma^2, FUN = "*")
      dsigma_sq <- -1/(2*sigma^2) + resid^2/(2*sigma^4)
      return(cbind(dkappa, dsigma_sq))
    }
  } else if (best_model_name=="mixture_normal") {
    gamma <- FI_result$gamma.hat
    lambda <- gamma[1]; gamma_1_2 <- gamma[-1]  
    gamma1 <- gamma_1_2[1:(length(gamma_1_2)/2)]
    gamma2 <- gamma_1_2[(length(gamma_1_2)/2+1):length(gamma_1_2)]
    kappa1 <- gamma1[-length(gamma1)]; sigma1 <- gamma1[length(gamma1)]
    kappa2 <- gamma2[-length(gamma2)]; sigma2 <- gamma2[length(gamma2)]
    phi <- FI_result$phi.hat; alpha <- phi[-length(phi)]; beta <- phi[length(phi)] 
    mu <- FI_result$mu.hat
    n <- nrow(FI_result$data); n1 <- sum(FI_result$data$R); n0 <- n-n1
    outcome_model <- FI_result$outcome_model; response_model <- FI_result$response_model
    data_aug <- as.matrix(FI_result$data_aug)
    data1 <- data_aug[data_aug[, "R"]==1,  ]
    data0_aug <- data_aug[data_aug[, "R"]==0,  ]
    M <- FI_result$M
    
    # score function for the respondents outcome model
    score_1 <- function(data) {
      y <- data[, "y"]; design <- cbind(1, data[, outcome_model[-1]])
      mu1 <- design%*%kappa1; mu2 <- design%*%kappa2
      resid1 <- y - mu1 ; resid2 <- y - mu2
      f1 <- dnorm(y, mean=mu1, sd=sigma1)
      f2 <- dnorm(y, mean=mu2, sd=sigma2)
      dlambda <- (f1 - f2)/(lambda*f1 + (1-lambda)*f2) 
      dkappa1 <- sweep(design, MARGIN = 1, lambda*f1/(lambda*f1 + (1-lambda)*f2)*resid1/(sigma1^2), FUN = "*")
      dsigma1_sq <- lambda*f1/(lambda*f1 + (1-lambda)*f2)*(-1/(2*sigma1^2) + resid1^2/(2*sigma1^4))
      dkappa2 <- sweep(design, MARGIN = 1, (1-lambda)*f2/(lambda*f1 + (1-lambda)*f2)*resid2/(sigma2^2), FUN = "*")
      dsigma2_sq <- (1-lambda)*f2/(lambda*f1 + (1-lambda)*f2)*(-1/(2*sigma2^2) + resid2^2/(2*sigma2^4))
      return(cbind(dlambda, dkappa1, dsigma1_sq, dkappa2, dsigma2_sq))
    }
  } else if (best_model_name=="logistic") {
    gamma <- FI_result$gamma.hat
    phi <- FI_result$phi.hat; alpha <- phi[-length(phi)]; beta <- phi[length(phi)] 
    mu <- FI_result$mu.hat
    n <- nrow(FI_result$data); n1 <- sum(FI_result$data$R); n0 <- n-n1
    outcome_model <- FI_result$outcome_model; response_model <- FI_result$response_model
    data_aug <- as.matrix(FI_result$data_aug)
    data1 <- data_aug[data_aug[, "R"]==1,  ]
    data0_aug <- data_aug[data_aug[, "R"]==0,  ]
    M <- 2
    
    # score function for the respondents outcome model
    score_1 <- function(data) {
      y <- data[, "y"]; design <- cbind(1, data[, outcome_model[-1]])
      yhat <- 1/(1+exp(-design%*%gamma))
      resid <- y-yhat
      return(sweep(design, MARGIN = 1, resid, FUN = "*"))
    }
  }
  
  # response probability
  pi <- function(data) {
    lin_pred <- cbind(1, data[, response_model[-1]]) %*% phi
    return(1/(1+exp(-lin_pred)))
  }
  
  # score function for the response model
  score_p <- function(data) {
    R <- data[, "R"]; design <- cbind(1, data[, response_model[-1]])
    resid <- R - pi(data)
    return(sweep(design, MARGIN = 1, resid, FUN = "*"))
  }
  
  # estimating function for the mean of Y
  u_function <- function(data) {
    y <- data[, "y"]
    return(y-mu)
  }
  
  # score function for the respondent outcome model
  score_1_data1 <- score_1(data1)
  I11 <- t(score_1_data1)%*%score_1_data1/n
  
  # score function for the response model
  s0 <- score_p(data0_aug)
  s0 <- as.data.frame(cbind(s0, data0_aug[, c("frac.w", "ID")]))
  colnames(s0)[1:length(c("intercept", response_model[-1]))] <- c("intercept", response_model[-1])
  s0.bar <- s0 %>% group_by(ID) %>% summarize(across(all_of(c("intercept", response_model[-1])), ~sum(. * frac.w))); s0.bar$ID <- NULL; s0.bar <- as.matrix(s0.bar)
  s0$frac.w <- NULL; s0$ID <- NULL
  s0_diff <- s0 - s0.bar[rep(1:n0, each=M), ]
  score_1_data0aug <- score_1(data0_aug)
  I21 <- (-1/n)*t(sweep(s0_diff, MARGIN = 1, data0_aug[, "frac.w"], FUN = "*"))%*%score_1_data0aug
  
  z0 <- as.data.frame(cbind(1, data0_aug)) %>% rename(intercept = V1)
  z0.bar <- z0 %>% group_by(ID) %>% summarize(across(all_of(c("intercept", response_model[-1])), ~sum(. * frac.w))); z0.bar$ID <- NULL; z0.bar <- as.matrix(z0.bar)
  
  I22 <- (-1/n)*t(s0.bar)%*%z0.bar
  
  score_p_data1 <- score_p(data1)
  J <- score_p_data1 - score_1_data1%*%t(I21%*%solve(I11))
  V1 <- (t(J)%*%J + t(s0.bar)%*%s0.bar)/n
  
  pi_data1 <- pi(data1)
  tau <- (-1/n)*sum(1/pi_data1)
  kappa <- I21%*%solve(I11)
  
  u_data1 <- u_function(data1)
  
  B <- t(apply(sweep(score_p_data1, MARGIN = 1, u_data1/pi_data1, FUN = "*"), MARGIN = 2, FUN = sum))%*%solve(I22)/n
  L <- u_data1/pi_data1 - (score_p_data1 - score_1_data1%*%t(kappa))%*%t(B)
  V2 <- ( sum(L^2) + sum((B%*%t(s0.bar))^2) )/n
  
  variance_phi <- solve(I22)%*%V1%*%t(solve(I22))/n
  variance_mu <- solve(tau)%*%V2%*%t(solve(tau))/n
  lower.bound_phi <- phi - sqrt(diag(variance_phi)) * qnorm(0.975)
  upper.bound_phi <- phi + sqrt(diag(variance_phi)) * qnorm(0.975)
  CI_phi <- cbind(lower.bound_phi, upper.bound_phi)
  lower.bound_mu <- mu - sqrt(diag(variance_mu)) * qnorm(0.975)
  upper.bound_mu <- mu + sqrt(diag(variance_mu)) * qnorm(0.975)
  CI_mu <- cbind(lower.bound_mu, upper.bound_mu)
  
  Estimate <- FI_result$phi.hat
  Std.Error <- sqrt(diag(variance_phi))
  z_value <- Estimate/Std.Error
  calculate_p_value <- function(z) {
    if (z>0) {
      p_value <- 2*(1-pnorm(z))
    } else {
      p_value <- 2*pnorm(z)
    }
    return(p_value)
  }
  p_value <- sapply(z_value, calculate_p_value)
  summary_table <- data.frame(Estimate, Std.Error, z_value, p_value)
  colnames(summary_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  VE_result <- list(variance_phi = variance_phi, variance_mu = variance_mu,
                    CI_phi = CI_phi, CI_mu = CI_mu, summary_table=summary_table)
  return(VE_result)
}

VE_MAR <- function(data, XX_names, Y_name, R_name) {
  n <- nrow(data)
  XX <- data[, XX_names, drop=F]
  Y <- data[, Y_name]
  R <- data[, R_name]
  data <- data.frame(XX, Y=Y, R=R)
  data_obs <- (data[data$R==1, ])[, -ncol(data)]
  
  # Response model fitting
  if (length(XX_names)!=0) {
    res_formula <- as.formula(paste("R", "~", paste(XX_names, collapse = "+")))
  } else {
    res_formula <- R~1
  }
  est_res <- glm(res_formula, data=data, family = "binomial")
  e_alp_MAR <- est_res$coefficients
  IPW_MAR <- 1/est_res$fitted.values
  mean_PS_MAR <- sum(IPW_MAR[R==1]*Y[R==1])/sum(IPW_MAR[R==1])
  
  # variance estimation of the PS estimator under MAR assumption
  XX_res <- as.matrix(cbind(1, XX))
  p2 <- ncol(XX_res)
  A11 <- -as.matrix(mean(R*IPW_MAR))
  B11 <- mean(R*(IPW_MAR^2)*(Y-mean_PS_MAR)^2)
  
  coeff <- (R/exp(XX_res%*%e_alp_MAR))*(Y-mean_PS_MAR)
  B12 <- apply(matrix(rep(coeff, p2), ncol=p2)*XX_res, FUN=mean, MARGIN=2)
  
  coeff2 <- R/(1+exp(XX_res%*%e_alp_MAR))
  XX_sqrt <- matrix(rep(sqrt(coeff2), p2), ncol=p2)*XX_res
  B22 <- t(XX_sqrt)%*%XX_sqrt/n
  
  A21 <- matrix(0, nrow=nrow(B22), ncol=ncol(A11))
  A <- rbind(cbind(A11, t(B12)),
             cbind(A21, B22))
  B <- rbind(cbind(B11, t(B12)),
             cbind(t(t(B12)), B22))
  
  cov_mat <- solve(A) %*% B %*% t(solve(A)) / n
  vari_mean_MAR <- diag(cov_mat)[1]
  vari_alp_MAR <- diag(cov_mat)[-1]
  
  lower.bound_alp <- e_alp_MAR - sqrt(vari_alp_MAR) * qnorm(0.975)
  upper.bound_alp <- e_alp_MAR + sqrt(vari_alp_MAR) * qnorm(0.975)
  lower.bound_mu <- mean_PS_MAR - sqrt(vari_mean_MAR) * qnorm(0.975)
  upper.bound_mu <- mean_PS_MAR + sqrt(vari_mean_MAR) * qnorm(0.975)
  
  CI_alp_MAR <- cbind(lower.bound_alp, upper.bound_alp)
  CI_mean_MAR <- cbind(lower.bound_mu, upper.bound_mu)
  
  Estimate <- c(e_alp_MAR)
  Std.Error <- sqrt(diag(cov_mat)[-1])
  z_value <- Estimate/Std.Error
  calculate_p_value <- function(z) {
    if (z>0) {
      p_value <- 2*(1-pnorm(z))
    } else {
      p_value <- 2*pnorm(z)
    }
    return(p_value)
  }
  p_value <- sapply(z_value, calculate_p_value)
  summary_table <- data.frame(Estimate, Std.Error, z_value, p_value)
  colnames(summary_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  VE_result <- list(variance_phi = vari_alp_MAR, variance_mu = vari_mean_MAR,
                    CI_phi = CI_alp_MAR, CI_mu = CI_mean_MAR, summary_table=summary_table,
                    mean_PS_MAR = mean_PS_MAR)
  
  return(VE_result)
}
