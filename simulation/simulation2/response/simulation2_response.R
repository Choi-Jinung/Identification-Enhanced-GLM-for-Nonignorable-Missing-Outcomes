rm(list=ls())

# Import functions
source(sub("simulation/simulation2/response$", "functions/FI_functions.R", getwd()))
source(sub("simulation/simulation2/response$", "functions/PPMA_functions.R", getwd()))

# Import libraries
library(doParallel)
library(doRNG)
library(bestglm)
library(ggplot2)

# Link functions
logit <- function(x) 1 / (1 + exp(-x))
probit <- function(x) pnorm(x)
loglog <- function(x) exp(-exp(-x))
cloglog <- function(x) 1 - exp(-exp(x))
inv_quad_logit <- function(x) 1.5 - (1 / (1 + exp(-x^2)))

# common setups
n=1000; sig_level = 0.1; B=1000; M=100
mu_x1 <- 0.5; mu_x2 <- 0.5; sigma_x1 <- sqrt(0.5); sigma_x2 <- sqrt(0.5); sigma_y <- sigma_x1; mu_y <- mu_x1
identifiability <- 1

# response_model <- "logit"
# response_model <- "probit"
# response_model <- "loglog"
# response_model <- "cloglog"
# response_model <- "inv_quad_logit"

# rho <- 0.5
# rho <- 0.8

# res_mch <- "MAR"
# res_mch <- "NMAR1"
# res_mch <- "NMARinf"

# Simulation function
simulation2_response <- function(response_model, res_mch, rho) {
  sigma_e <- sqrt(1-rho^2)*sigma_y
  
  kappa1 <- sqrt(rho^2*sigma_y^2/(sigma_x1^2+identifiability^2*sigma_x2^2))
  kappa2 <- identifiability*kappa1
  kappa2/kappa1
  kappa0 <- mu_y - kappa1*mu_x1 - kappa2*mu_x2
  kap <- c(kappa0, kappa1, kappa2)
  
  if (response_model=="logit") {
    link <- logit
    alp0 <- 0
  } else if (response_model=="probit") {
    link <- probit
    alp0 <- -0.25
  } else if (response_model=="loglog") {
    link <- loglog
    alp0 <- 0.3
  } else if (response_model=="cloglog") {
    link <- cloglog
    alp0 <- -0.7
  } else if (response_model=="inv_quad_logit") {
    link <- inv_quad_logit
    alp0 <- 0.4
  }
  
  if (res_mch=="MAR") {
    alp1 <- 2; alp2 <- 0; beta <- 0
  } else if (res_mch=="NMAR1") {
    alp1 <- 1; alp2 <- 0; beta <- 1
  } else if (res_mch=="NMARinf") {
    alp1 <- 0; alp2 <- 0; beta <- 2
  }
  alp <- c(alp0, alp1, alp2); phi <- c(alp, beta)
  
  sim <- function(it) {
    x1 <- rnorm(n, mean = mu_x1, sd = sigma_x1)
    x2 <- rnorm(n, mean = mu_x2, sd = sigma_x2)
    design_out <- cbind(1, x1, x2)
    mu <- design_out %*% kap
    e <- rnorm(n=n, mean = 0, sd=sigma_e)
    y <- mu + e
    design_res <- cbind(1, x1, x2, y)
    lin_res <- design_res %*% phi
    prob <- link(lin_res)
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
  
    if (res_mch=="MAR" | res_mch=="NMAR1") {
      FI_c <- FI_out_MNAR_list[[2]]
    } else if (res_mch=="NMARinf") {
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

see_result <- function(output1) {
  mu <- mean(unlist(do.call(rbind, lapply(output1, FUN = function(x) x$full))))
  mu_estimates <- do.call(rbind, lapply(output1, FUN = function(x) x$mu_estimate))
  
  point_summary <- function(true, point_estimate) {
    RB <- apply((point_estimate - true)/true*100, MARGIN = 2, FUN = mean)
    SE <- apply(point_estimate, MARGIN = 2, FUN = sd)
    RMSE <- sqrt(apply((point_estimate - true)^2, MARGIN = 2, FUN = mean))
    return(data.frame(RB=round(RB, 1), SE=round(SE, 3), RMSE=round(RMSE, 3)))
  }
  
  point_sum_mu <- point_summary(mu, mu_estimates)
  model_selection_result <- summary(as.factor(unlist(lapply(output1, FUN = function(x) x$selected_model))))
  
  return(list(point_sum_mu=point_sum_mu, model_selection_result=model_selection_result))
}

tic <- proc.time()

# Logit Response Model
logit_rho0.5_MAR <- simulation2_response(response_model="logit", rho=0.5, res_mch="MAR")
logit_rho0.5_NMAR1 <- simulation2_response(response_model="logit", rho=0.5, res_mch="NMAR1")
logit_rho0.5_NMARinf <- simulation2_response(response_model="logit", rho=0.5, res_mch="NMARinf")
save(logit_rho0.5_MAR, file = "logit_rho0.5_MAR.Rdata")
save(logit_rho0.5_NMAR1, file = "logit_rho0.5_NMAR1.Rdata")
save(logit_rho0.5_NMARinf, file = "logit_rho0.5_NMARinf.Rdata")
see_result(logit_rho0.5_MAR)
see_result(logit_rho0.5_NMAR1)
see_result(logit_rho0.5_NMARinf)

logit_rho0.8_MAR <- simulation2_response(response_model="logit", rho=0.8, res_mch="MAR")
logit_rho0.8_NMAR1 <- simulation2_response(response_model="logit", rho=0.8, res_mch="NMAR1")
logit_rho0.8_NMARinf <- simulation2_response(response_model="logit", rho=0.8, res_mch="NMARinf")
save(logit_rho0.8_MAR, file = "logit_rho0.8_MAR.Rdata")
save(logit_rho0.8_NMAR1, file = "logit_rho0.8_NMAR1.Rdata")
save(logit_rho0.8_NMARinf, file = "logit_rho0.8_NMARinf.Rdata")
see_result(logit_rho0.8_MAR)
see_result(logit_rho0.8_NMAR1)
see_result(logit_rho0.8_NMARinf)


# Probit Response Model
probit_rho0.5_MAR <- simulation2_response(response_model="probit", rho=0.5, res_mch="MAR")
probit_rho0.5_NMAR1 <- simulation2_response(response_model="probit", rho=0.5, res_mch="NMAR1")
probit_rho0.5_NMARinf <- simulation2_response(response_model="probit", rho=0.5, res_mch="NMARinf")
save(probit_rho0.5_MAR, file = "probit_rho0.5_MAR.Rdata")
save(probit_rho0.5_NMAR1, file = "probit_rho0.5_NMAR1.Rdata")
save(probit_rho0.5_NMARinf, file = "probit_rho0.5_NMARinf.Rdata")
see_result(probit_rho0.5_MAR)
see_result(probit_rho0.5_NMAR1)
see_result(probit_rho0.5_NMARinf)


probit_rho0.8_MAR <- simulation2_response(response_model="probit", rho=0.8, res_mch="MAR")
probit_rho0.8_NMAR1 <- simulation2_response(response_model="probit", rho=0.8, res_mch="NMAR1")
probit_rho0.8_NMARinf <- simulation2_response(response_model="probit", rho=0.8, res_mch="NMARinf")
save(probit_rho0.8_MAR, file = "probit_rho0.8_MAR.Rdata")
save(probit_rho0.8_NMAR1, file = "probit_rho0.8_NMAR1.Rdata")
save(probit_rho0.8_NMARinf, file = "probit_rho0.8_NMARinf.Rdata")
see_result(probit_rho0.8_MAR)
see_result(probit_rho0.8_NMAR1)
see_result(probit_rho0.8_NMARinf)


# Log-log Response Model
loglog_rho0.5_MAR <- simulation2_response(response_model="loglog", rho=0.5, res_mch="MAR")
loglog_rho0.5_NMAR1 <- simulation2_response(response_model="loglog", rho=0.5, res_mch="NMAR1")
loglog_rho0.5_NMARinf <- simulation2_response(response_model="loglog", rho=0.5, res_mch="NMARinf")
save(loglog_rho0.5_MAR, file = "loglog_rho0.5_MAR.Rdata")
save(loglog_rho0.5_NMAR1, file = "loglog_rho0.5_NMAR1.Rdata")
save(loglog_rho0.5_NMARinf, file = "loglog_rho0.5_NMARinf.Rdata")
see_result(loglog_rho0.5_MAR)
see_result(loglog_rho0.5_NMAR1)
see_result(loglog_rho0.5_NMARinf)

loglog_rho0.8_MAR <- simulation2_response(response_model="loglog", rho=0.8, res_mch="MAR")
loglog_rho0.8_NMAR1 <- simulation2_response(response_model="loglog", rho=0.8, res_mch="NMAR1")
loglog_rho0.8_NMARinf <- simulation2_response(response_model="loglog", rho=0.8, res_mch="NMARinf")
save(loglog_rho0.8_MAR, file = "loglog_rho0.8_MAR.Rdata")
save(loglog_rho0.8_NMAR1, file = "loglog_rho0.8_NMAR1.Rdata")
save(loglog_rho0.8_NMARinf, file = "loglog_rho0.8_NMARinf.Rdata")
see_result(loglog_rho0.8_MAR)
see_result(loglog_rho0.8_NMAR1)
see_result(loglog_rho0.8_NMARinf)


# C.Log-log Response Model
cloglog_rho0.5_MAR <- simulation2_response(response_model="cloglog", rho=0.5, res_mch="MAR")
cloglog_rho0.5_NMAR1 <- simulation2_response(response_model="cloglog", rho=0.5, res_mch="NMAR1")
cloglog_rho0.5_NMARinf <- simulation2_response(response_model="cloglog", rho=0.5, res_mch="NMARinf")
save(cloglog_rho0.5_MAR, file = "cloglog_rho0.5_MAR.Rdata")
save(cloglog_rho0.5_NMAR1, file = "cloglog_rho0.5_NMAR1.Rdata")
save(cloglog_rho0.5_NMARinf, file = "cloglog_rho0.5_NMARinf.Rdata")
see_result(cloglog_rho0.5_MAR)
see_result(cloglog_rho0.5_NMAR1)
see_result(cloglog_rho0.5_NMARinf)

cloglog_rho0.8_MAR <- simulation2_response(response_model="cloglog", rho=0.8, res_mch="MAR")
cloglog_rho0.8_NMAR1 <- simulation2_response(response_model="cloglog", rho=0.8, res_mch="NMAR1")
cloglog_rho0.8_NMARinf <- simulation2_response(response_model="cloglog", rho=0.8, res_mch="NMARinf")
save(cloglog_rho0.8_MAR, file = "cloglog_rho0.8_MAR.Rdata")
save(cloglog_rho0.8_NMAR1, file = "cloglog_rho0.8_NMAR1.Rdata")
save(cloglog_rho0.8_NMARinf, file = "cloglog_rho0.8_NMARinf.Rdata")
see_result(cloglog_rho0.8_MAR)
see_result(cloglog_rho0.8_NMAR1)
see_result(cloglog_rho0.8_NMARinf)


# Inverse Quad Logit Response Model
inv_quad_logit_rho0.5_MAR <- simulation2_response(response_model="inv_quad_logit", rho=0.5, res_mch="MAR")
inv_quad_logit_rho0.5_NMAR1 <- simulation2_response(response_model="inv_quad_logit", rho=0.5, res_mch="NMAR1")
inv_quad_logit_rho0.5_NMARinf <- simulation2_response(response_model="inv_quad_logit", rho=0.5, res_mch="NMARinf")
save(inv_quad_logit_rho0.5_MAR, file = "inv_quad_logit_rho0.5_MAR.Rdata")
save(inv_quad_logit_rho0.5_NMAR1, file = "inv_quad_logit_rho0.5_NMAR1.Rdata")
save(inv_quad_logit_rho0.5_NMARinf, file = "inv_quad_logit_rho0.5_NMARinf.Rdata")
see_result(inv_quad_logit_rho0.5_MAR)
see_result(inv_quad_logit_rho0.5_NMAR1)
see_result(inv_quad_logit_rho0.5_NMARinf)

inv_quad_logit_rho0.8_MAR <- simulation2_response(response_model="inv_quad_logit", rho=0.8, res_mch="MAR")
inv_quad_logit_rho0.8_NMAR1 <- simulation2_response(response_model="inv_quad_logit", rho=0.8, res_mch="NMAR1")
inv_quad_logit_rho0.8_NMARinf <- simulation2_response(response_model="inv_quad_logit", rho=0.8, res_mch="NMARinf")
save(inv_quad_logit_rho0.8_MAR, file = "inv_quad_logit_rho0.8_MAR.Rdata")
save(inv_quad_logit_rho0.8_NMAR1, file = "inv_quad_logit_rho0.8_NMAR1.Rdata")
save(inv_quad_logit_rho0.8_NMARinf, file = "inv_quad_logit_rho0.8_NMARinf.Rdata")
see_result(inv_quad_logit_rho0.8_MAR)
see_result(inv_quad_logit_rho0.8_NMAR1)
see_result(inv_quad_logit_rho0.8_NMARinf)


#### Output Organization ####
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

response_result <- rbind(cbind(organize_result(logit_rho0.5_MAR, logit_rho0.5_NMAR1, logit_rho0.5_NMARinf), organize_result(logit_rho0.8_MAR, logit_rho0.8_NMAR1, logit_rho0.8_NMARinf)),
                         cbind(organize_result(probit_rho0.5_MAR, probit_rho0.5_NMAR1, probit_rho0.5_NMARinf), organize_result(probit_rho0.8_MAR, probit_rho0.8_NMAR1, probit_rho0.8_NMARinf)),
                         cbind(organize_result(loglog_rho0.5_MAR, loglog_rho0.5_NMAR1, loglog_rho0.5_NMARinf), organize_result(loglog_rho0.8_MAR, loglog_rho0.8_NMAR1, loglog_rho0.8_NMARinf)),
                         cbind(organize_result(cloglog_rho0.5_MAR, cloglog_rho0.5_NMAR1, cloglog_rho0.5_NMARinf), organize_result(cloglog_rho0.8_MAR, cloglog_rho0.8_NMAR1, cloglog_rho0.8_NMARinf)))
write.csv(response_result, file = "response_result.csv")



# point estimation performance with model selection
result_m <- cbind(rbind(logit_rho0.5$sim_summary_mu_real, probit_rho0.5$sim_summary_mu_real, loglog_rho0.5$sim_summary_mu_real, cloglog_rho0.5$sim_summary_mu_real, inv_quad_logit_rho0.5$sim_summary_mu_real),
                  rbind(logit_rho0.8$sim_summary_mu_real, probit_rho0.8$sim_summary_mu_real, loglog_rho0.8$sim_summary_mu_real, cloglog_rho0.8$sim_summary_mu_real, inv_quad_logit_rho0.8$sim_summary_mu_real))
write.csv(result_m, file = "result_m.csv")

# model selection result
model_selection_result <- bind_rows(
  logit_rho0.5 = logit_rho0.5$selected_model,
  probit_rho0.5 = probit_rho0.5$selected_model,
  loglog_rho0.5 = loglog_rho0.5$selected_model,
  cloglog_rho0.5 = cloglog_rho0.5$selected_model,
  invquad_rho0.5 = inv_quad_logit_rho0.5$selected_model,
  logit_rho0.8 = logit_rho0.8$selected_model,
  probit_rho0.8 = probit_rho0.8$selected_model,
  loglog_rho0.8 = loglog_rho0.8$selected_model,
  cloglog_rho0.8 = cloglog_rho0.8$selected_model,
  invquad_rho0.8 = inv_quad_logit_rho0.8$selected_model,
  .id = "source"
)
model_selection_result[is.na(model_selection_result)] <- 0
model_selection_result <- model_selection_result %>% mutate(total=apply(model_selection_result[, -1], MARGIN=1, FUN=sum))
write.csv(model_selection_result, file = "model_selection_result.csv")


# point estimation performance with correctly specified response model
rho0.5_mu <- rbind(logit_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   probit_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   loglog_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   cloglog_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   inv_quad_logit_rho0.5$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ])
rho0.8_mu <- rbind(logit_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   probit_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   loglog_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   cloglog_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ],
                   inv_quad_logit_rho0.8$sim_summary_mu[c("CC", "PPM_1", "FI_c"), ])
result_c_mu <- cbind(rho0.5_mu, rho0.8_mu)
result_c <- rbind(result_c_mu)
write.csv(result_c, file = "result_c.csv")
