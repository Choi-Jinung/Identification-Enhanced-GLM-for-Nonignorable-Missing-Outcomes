rm(list=ls())

# Import functions
source(sub("simulation/simulation1/continuous$", "functions/FI_functions.R", getwd()))
source(sub("simulation/simulation1/continuous$", "functions/PPMA_functions.R", getwd()))

# Import libraries
library(doParallel)
library(doRNG)
library(bestglm)
library(ggplot2)
library(dplyr)

# common setups
n=1000; sig_level = 0.1; B=5000; M=100
mu_x1 <- 0.5; mu_x2 <- 0.5; sigma_x1 <- sqrt(0.5); sigma_x2 <- sqrt(0.5); sigma_y <- sigma_x1; mu_y <- mu_x1
phi=c(0, 1, 0, 1); alp <- phi[-length(phi)]; beta <- phi[length(phi)]; rho <- 0.8

# identifiability <- 1
# identifiability <- 0.6
# identifiability <- 0.3
# identifiability <- 0

# Simulation function
simulation1_conti <- function(identifiability) {
  sigma_e <- sqrt(1-rho^2)*sigma_y
  kappa1 <- sqrt(rho^2*sigma_y^2/(sigma_x1^2+identifiability^2*sigma_x2^2))
  kappa2 <- identifiability*kappa1
  kappa2/kappa1
  round(kappa1, digits = 2)
  round(kappa2, digits = 2)
  kappa0 <- mu_y - kappa1*mu_x1 - kappa2*mu_x2
  kap <- c(kappa0, kappa1, kappa2)
  
  sim <- function(it) {
    x1 <- rnorm(n, mean = mu_x1, sd = sigma_x1)
    x2 <- rnorm(n, mean = mu_x2, sd = sigma_x2)
    design_out <- cbind(1, x1, x2)
    mu <- design_out %*% kap
    e <- rnorm(n=n, mean = 0, sd=sigma_e)
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
          FI_out_MNAR_list[[j]] <- FI(pop=data, prop=prop, response=c(res_xy[[j]], "y"), phi_curr=rep(0, length(c(res_xy[[j]], "y"))+1), M=M)
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

strong3_out <- simulation1_conti(identifiability = 1)
save(strong3_out, file = "strong3_out.Rdata")
see_result(strong3_out)

moderate_out <- simulation1_conti(identifiability = 0.6)
save(moderate_out, file = "moderate_out.Rdata")
see_result(moderate_out)

weak_out <- simulation1_conti(identifiability = 0.3)
save(weak_out, file = "weak_out.Rdata")
see_result(weak_out)

unidentifiable_out <- simulation1_conti(identifiability = 0)
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
load(file="p_value_data_binary.Rdata")

p_value_data_combined <- rbind(p_value_data, p_value_data_binary)
p_value_data_combined <- p_value_data_combined %>% mutate(Outcome = factor(rep(c("Continuous", "Binary"), each=B*4), levels=c("Continuous", "Binary")))

ggplot(p_value_data_combined, aes(x = ratio, y = p_value, fill = Outcome)) +
  geom_boxplot() +
  labs(x = expression(theta[2]/theta[1]),
       y = "p-value") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
  theme_minimal() + 
  facet_wrap(~ Outcome) + 
  scale_fill_brewer(palette = "Set2")