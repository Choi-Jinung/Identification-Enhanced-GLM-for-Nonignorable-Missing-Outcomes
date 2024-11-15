rm(list=ls())

# Import functions
source(sub("simulation/simulation1/binary$", "functions/FI_functions.R", getwd()))
source(sub("simulation/simulation1/binary$", "functions/binaryPPMA_functions.R", getwd()))

# Import libraries
library(doParallel)
library(doRNG)
library(bestglm)
library(ggplot2)
library(msm)

# common setups
n=1000; sig_level = 0.1; B=5000
mu_x1 <- 0.5; mu_x2 <- 0.5; sigma_x1 <- sqrt(0.25); sigma_x2 <- sqrt(0.25)
phi <- c(0, 1, 1); alp <- phi[-length(phi)]; beta <- phi[length(phi)]; rho <- 0.8

# outcome_link <- "probit"
# outcome_link <- "logit"
# identifiability <- 1.5
# identifiability <- 0.5
# identifiability <- 0.1
# identifiability <- 0

# Simulation function
simulation1_binary <- function(outcome_link="logit", identifiability) {
  if (outcome_link=="logit") {
    sigma_e <- sqrt(3.3)
  } else if (outcome_link=="probit") {
    sigma_e <- sqrt(1)
  }
  
  kappa1 <- sqrt(rho^2/(1-rho^2)*sigma_e^2/(sigma_x1^2+identifiability^2*sigma_x2^2))
  kappa2 <- identifiability*kappa1 
  kappa2/kappa1
  kappa0 <- -kappa1*mu_x1 - kappa2*mu_x2
  kap <- c(kappa0, kappa1, kappa2)
  
  sim <- function(it) {
    if (outcome_link=="logit") {
      e <- rlogis(n=n, location=0, scale=1)
    } else if (outcome_link=="probit") {
      e <- rnorm(n=n, mean=0, sd=1)
    }
    x1 <- rnorm(n=n, mean=mu_x1, sd=sigma_x1); x2 <- rnorm(n=n, mean=mu_x2, sd=sigma_x2)
    design_out <- cbind(1, x1, x2)
    proxy <- design_out %*% kap # proxy variable
    u <- proxy + e # latent variable
    cor(proxy, u)
    y <- (u > 0)*1 
    mean(y)
    mean(x1)
    mean(x2)
    var(y)
    var(x1)
    var(x2)
    glm(y~x1+x2, family = binomial(link=outcome_link))
    design_res <- cbind(1, x1, y)
    prob <- exp(design_res %*% phi) / (1+exp(design_res %*% phi))
    sum(prob)/n
    R <- rbinom(n=n, size=1, prob=prob)
    data <- data.frame(x1=x1, x2=x2, y=y, R=R)
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
      # Execute FI function
      FI_out_MNAR_list[[j]] <- FI(pop=data, prop=prop, response=c(res_xy[[j]], "y"), phi_curr=rep(0, length(c(res_xy[[j]], "y"))+1), M=2)
      BIC_MNAR[j] <- FI_out_MNAR_list[[j]]$BIC_full
      
      # Execute VE_NMAR function
      try({
        VE_FI_MNAR_list[[j]] <- VE_NMAR(FI_out_MNAR_list[[j]])
        p_values[j] <- VE_FI_MNAR_list[[j]]$summary_table["y", "Pr(>|z|)"]
      })
    }
    
    best_MNAR <- which.min(BIC_MNAR)
    if (p_values[best_MNAR] < sig_level) {
      New_out <- FI_out_MNAR_list[[best_MNAR]]$mu.hat
      selected_model <- paste0("NMAR", best_MNAR, collapse = "")
    } else {
      New_out <- VE_MAR_out$mean_PS_MAR
      selected_model <- "MAR"
    }
    
    # Unmodified PPM
    ppm_fit <- glm(y~x1+x2, data=data[R==1, ], family=binomial(link="probit"))
    x_ppm <- predict(ppm_fit, newdata=data)
    y_ppm <- y; y_ppm[R==0] <- NA
    data_ppm <- data.frame(x_ppm=x_ppm, y_ppm=y_ppm)
    PPM_0 <- mleFull(x_ppm, y_ppm, 0)
    PPM_1 <- mleFull(x_ppm, y_ppm, 0.5)
    PPM_inf <- mleFull(x_ppm, y_ppm, 1)
    # # Modified PPM
    # mPPM_0 <- mle2step(x_ppm, y_ppm, 0)
    # mPPM_1 <- mle2step(x_ppm, y_ppm, 0.5)
    # mPPM_inf <- mle2step(x_ppm, y_ppm, 1)
    # # Bootstrap variance estimation
    # n_boot <- 1000
    # boot_estimates <- matrix(nrow=n_boot, ncol=3)
    # colnames(boot_estimates) <- c("mPPM_0", "mPPM_1", "mPPM_inf")
    # for (i in 1:n_boot) {
    #   index_boot <- sample(1:nrow(data_ppm), size=nrow(data_ppm), replace=TRUE)
    #   boot_sample <- data_ppm[index_boot, ]
    #   x_boot <- boot_sample$x_ppm; y_boot <- boot_sample$y_ppm
    #   boot_estimates[i, "mPPM_0"] <- mle2step(x_boot, y_boot, 0)$muY
    #   boot_estimates[i, "mPPM_1"] <- mle2step(x_boot, y_boot, 0.5)$muY
    #   boot_estimates[i, "mPPM_inf"] <- mle2step(x_boot, y_boot, 1)$muY
    # }
    # mPPM_0$muYvar <- var(boot_estimates[, "mPPM_0"])
    # mPPM_1$muYvar <- var(boot_estimates[, "mPPM_1"])
    # mPPM_inf$muYvar <- var(boot_estimates[, "mPPM_inf"])
    
    mu_estimate <- c(CC, PPM_0$muY, PPM_1$muY, PPM_inf$muY, VE_MAR_out$mean_PS_MAR, FI_out_MNAR_list[[2]]$mu.hat, New_out)
    names(mu_estimate) <- c("CC", "PPM_0", "PPM_1", "PPM_inf", "MAR", "FI_c", "FI_m")
    var_mu_estimate <- c(var_CC, PPM_0$muYvar, PPM_1$muYvar, PPM_inf$muYvar, VE_MAR_out$variance_mu, VE_FI_MNAR_list[[2]]$variance_mu)
    names(var_mu_estimate) <- c("CC", "PPM_0", "PPM_1", "PPM_inf", "MAR", "FI_c")
    
    beta_estimate <- c(FI_out_MNAR_list[[2]]$phi.hat[length(FI_out_MNAR_list[[2]]$phi.hat)])
    var_beta_estimate <- c(diag(VE_FI_MNAR_list[[2]]$variance_phi)[length(diag(VE_FI_MNAR_list[[2]]$variance_phi))])
    
    out <- list(mu_estimate=mu_estimate,
                var_mu_estimate=var_mu_estimate,
                beta_estimate=beta_estimate,
                var_beta_estimate=var_beta_estimate,
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

see_result <- function(output1) {
  mu <- mean(unlist(do.call(rbind, lapply(output1, FUN = function(x) x$full))))
  mu_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$mu_estimate))
  var_mu_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$var_mu_estimate))
  beta_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$beta_estimate))
  var_beta_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$var_beta_estimate))
  
  point_summary <- function(true, point_estimate) {
    RB <- apply((point_estimate - true)/true*100, MARGIN = 2, FUN = mean)
    SE <- apply(point_estimate, MARGIN = 2, FUN = sd)
    RMSE <- sqrt(apply((point_estimate - true)^2, MARGIN = 2, FUN = mean))
    return(data.frame(RB=round(RB, 1), SE=round(SE, 3), RMSE=round(RMSE, 3)))
  }
  var_summary <- function(true, point_estimate, var_estimate) {
    Var <- apply(point_estimate, MARGIN = 2, FUN = var)
    Var_mat <- matrix(rep(Var, B), nrow=B, byrow=T)
    RB_var <- apply((var_estimate - Var_mat)/Var_mat*100, MARGIN = 2, FUN = mean)
    lower <- point_estimate - qnorm(0.975)*sqrt(var_estimate)
    upper <- point_estimate + qnorm(0.975)*sqrt(var_estimate)
    CR <- apply(lower < true & true < upper, MARGIN = 2, FUN = mean)*100
    return(data.frame(RB=round(RB_var, 1), CR=round(CR, 1)))
  }
  point_sum_mu <- point_summary(mu, mu_estimates)
  piont_sum_beta <- point_summary(beta, beta_estimates)
  var_sum_mu <- var_summary(mu, mu_estimates[, colnames(var_mu_estimates)], var_mu_estimates)
  var_sum_beta <- var_summary(beta, beta_estimates, var_beta_estimates)
  point_sum <- rbind(point_sum_mu, piont_sum_beta)
  var_sum <- rbind(var_sum_mu, var_sum_beta)
  
  model_selection_result <- summary(as.factor(unlist(lapply(output1, FUN = function(x) x$selected_model))))
  
  return(list(point_sum=point_sum, var_sum=var_sum, model_selection_result=model_selection_result))
}
# strong1_out <- simulation1_binary(identifiability = 2)
# # save(strong1_out, file = "strong1_out.Rdata")
# see_result(strong1_out)
# 
# strong2_out <- simulation1_binary(identifiability = 1.5)
# # save(strong2_out, file = "strong2_out.Rdata")
# see_result(strong2_out)

strong3_out <- simulation1_binary(identifiability = 1)
save(strong3_out, file = "strong3_out.Rdata")
see_result(strong3_out)

moderate_out <- simulation1_binary(identifiability = 0.6)
save(moderate_out, file = "moderate_out.Rdata")
see_result(moderate_out)

weak_out <- simulation1_binary(identifiability = 0.3)
save(weak_out, file = "weak_out.Rdata")
see_result(weak_out)

unidentifiable_out <- simulation1_binary(identifiability = 0)
save(unidentifiable_out, file = "unidentifiable_out.Rdata")
see_result(unidentifiable_out)


#### Output Organization ####
# model selection result
model_selection_result <- bind_rows(
  strong = see_result(strong3_out)$model_selection_result,
  moderate = see_result(moderate_out)$model_selection_result,
  weak = see_result(weak_out)$model_selection_result,
  unidentifiable = see_result(unidentifiable_out)$model_selection_result,
  .id = "source"
)
model_selection_result[is.na(model_selection_result)] <- 0
model_selection_result
write.csv(model_selection_result, file = "model_selection_result.csv")

# point estimation performance with model selection
point_result <- cbind(see_result(strong3_out)$point_sum,
                      see_result(moderate_out)$point_sum,
                      see_result(weak_out)$point_sum,
                      see_result(unidentifiable_out)$point_sum)
point_result
write.csv(point_result, file = "point_result.csv")

# variance estimation performance
var_result <- cbind(see_result(strong3_out)$var_sum,
                    see_result(moderate_out)$var_sum,
                    see_result(weak_out)$var_sum,
                    see_result(unidentifiable_out)$var_sum)
var_result
write.csv(var_result, file = "var_result.csv")

# box plot of p-value
p_val_strong <- do.call(rbind, lapply(strong3_out, FUN = function(x) x$p_value))[, 2]
p_val_moderate <- do.call(rbind, lapply(moderate_out, FUN = function(x) x$p_value))[, 2]
p_val_weak <- do.call(rbind, lapply(weak_out, FUN = function(x) x$p_value))[, 2]
p_val_unidentifiable <- do.call(rbind, lapply(unidentifiable_out, FUN = function(x) x$p_value))[, 2]
p_value_data <- data.frame(ratio = rep(c("1.0", "0.6", "0.3", "0.0"), each=B),
                           p_value = c(p_val_strong, p_val_moderate, p_val_weak, p_val_unidentifiable))
p_value_data$ratio <- factor(p_value_data$ratio, levels = c("1.0", "0.6", "0.3", "0.0"))

ggplot(p_value_data, aes(x = ratio, y = p_value)) +
  geom_boxplot() +
  labs(x = "Ratio",
       y = "p-value") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
  theme_minimal()
p_value_data_binary <- p_value_data
save(p_value_data_binary, file = "p_value_data_binary.Rdata")
