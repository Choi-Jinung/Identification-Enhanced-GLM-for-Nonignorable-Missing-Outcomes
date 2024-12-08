getwd()
setwd(getwd())
rm(list=ls())

# Import functions
source(sub("simulation/simulation2/outcome$", "functions/FI_functions.R", getwd()))
source(sub("simulation/simulation2/outcome$", "functions/PPMA_functions.R", getwd()))

# Import libraries
library(doParallel)
library(doRNG)
library(bestglm)
library(ggplot2)
library(dplyr)

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
    # estimators <- c("CC", "PPM_0", "PPM_1", "PPM_inf", "MAR", 
    #                 "FI1", "FI2", "FI3", "FI4")
    # mu_estimate <- c(CC, PPM_0$muY, PPM_1$muY, PPM_inf$muY, VE_MAR_out$mean_PS_MAR, 
    #                  FI_out_MNAR_list[[1]]$mu.hat, FI_out_MNAR_list[[2]]$mu.hat, FI_out_MNAR_list[[3]]$mu.hat, FI_out_MNAR_list[[4]]$mu.hat)
    # var_mu_estimate <- c(var_CC, PPM_0$muYvar, PPM_1$muYvar, PPM_inf$muYvar, VE_MAR_out$variance_mu,
    #                      VE_FI_MNAR_list[[1]]$variance_mu, VE_FI_MNAR_list[[2]]$variance_mu, VE_FI_MNAR_list[[3]]$variance_mu, VE_FI_MNAR_list[[4]]$variance_mu)
    estimators <- c("PPM_0", "PPM_1", "PPM_inf", "MAR", 
                    "FI1", "FI2", "FI3", "FI4")
    mu_estimate <- c(PPM_0$muY, PPM_1$muY, PPM_inf$muY, VE_MAR_out$mean_PS_MAR, 
                     FI_out_MNAR_list[[1]]$mu.hat, FI_out_MNAR_list[[2]]$mu.hat, FI_out_MNAR_list[[3]]$mu.hat, FI_out_MNAR_list[[4]]$mu.hat)
    var_mu_estimate <- c(PPM_0$muYvar, PPM_1$muYvar, PPM_inf$muYvar, VE_MAR_out$variance_mu,
                         VE_FI_MNAR_list[[1]]$variance_mu, VE_FI_MNAR_list[[2]]$variance_mu, VE_FI_MNAR_list[[3]]$variance_mu, VE_FI_MNAR_list[[4]]$variance_mu)
    lower_ci <- mu_estimate - qnorm(0.975)*sqrt(var_mu_estimate)
    upper_ci <-  mu_estimate + qnorm(0.975)*sqrt(var_mu_estimate)
    est_result <- data.frame(
      Estimator = estimators,
      point_estimate = mu_estimate,
      var_estimate = var_mu_estimate,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      ci_length = round(upper_ci - lower_ci, digits = 2)
    )
    out <- list(est_result=est_result, BIC_MNAR=BIC_MNAR, p_values=p_values, selected_model=selected_model)
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

# Normal Outcome Model
normal_rho0.5_NMAR1 <- simulation2_outcome(outcome_model="normal", phi=c(0, 1, 0, 1), rho=0.5)
normal_rho0.8_NMAR1 <- simulation2_outcome(outcome_model="normal", phi=c(0, 1, 0, 1), rho=0.8)

normal_rho0.5_NMAR1_mat <- do.call(rbind, lapply(normal_rho0.5_NMAR1, FUN = function(x) x$est_result$ci_length))
colnames(normal_rho0.5_NMAR1_mat) <- c("PPM_0", "PPM_1", "PPM_inf", "MAR", "FI1", "FI2", "FI3", "FI4")
apply(normal_rho0.5_NMAR1_mat, MARGIN = 2, FUN = median)
apply(normal_rho0.5_NMAR1_mat, MARGIN = 2, FUN = sd)

normal_rho0.8_NMAR1_mat <- do.call(rbind, lapply(normal_rho0.8_NMAR1, FUN = function(x) x$est_result$ci_length))
colnames(normal_rho0.8_NMAR1_mat) <- c("PPM_0", "PPM_1", "PPM_inf", "MAR", "FI1", "FI2", "FI3", "FI4")
apply(normal_rho0.8_NMAR1_mat, MARGIN = 2, FUN = median)
apply(normal_rho0.8_NMAR1_mat, MARGIN = 2, FUN = sd)

median_length_mat <- rbind(apply(normal_rho0.5_NMAR1_mat, MARGIN = 2, FUN = median),
                           apply(normal_rho0.8_NMAR1_mat, MARGIN = 2, FUN = median))
write.csv(median_length_mat, "median_length_mat.csv")

# Visualization
b <- 9
data_plot <- rbind(normal_rho0.5_NMAR1[[b]]$est_result, normal_rho0.8_NMAR1[[b]]$est_result)
data_plot <- data_plot %>% mutate(case=rep(1:2, each=8))
data_plot <- data_plot %>% mutate(Method=rep(c(rep("PPM", 3), rep("FI", 5)), 2))
data_plot$Estimator <- factor(data_plot$Estimator, levels = c("CC", "PPM_0", "PPM_1", "PPM_inf", "MAR", "FI1", "FI2", "FI3", "FI4"),
                              labels = c("CC", "PPM[0]", "PPM[1]", "PPM[\"\\U221E\"]", "MAR", "FI[1]", "FI[2]", "FI[3]", "FI[4]"))
custom_colors <- c("PPM" = "gray45", "FI" = "black")

true_value <- 0.5
ggplot(data_plot, aes(x = Estimator, y = point_estimate, color=Method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, linewidth=0.8) +
  geom_text(aes(label = round(lower_ci, 2), y = lower_ci - 0.02), size = 3) +
  geom_text(aes(label = round(upper_ci, 2), y = upper_ci + 0.02), size = 3) +
  geom_text(aes(label = round(point_estimate, 2), y = point_estimate, x = Estimator), 
            nudge_x = 0.35, size = 3, color = "blue") +  
  geom_hline(aes(yintercept = true_value), linetype = "dashed", color = "red") +
  facet_wrap(~ case) +
  labs(x = "Estimator",
       y = "y") +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  scale_color_manual(values = custom_colors) +
  theme(strip.text = element_blank(), legend.position = "none")

normal_rho0.5_NMAR1[[b]]$selected_model
normal_rho0.8_NMAR1[[b]]$selected_model
