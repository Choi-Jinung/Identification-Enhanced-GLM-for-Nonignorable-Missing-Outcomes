rm(list=ls())

# Import functions
source(sub("real_data/SRI$", "functions/FI_functions.R", getwd()))
source(sub("real_data/SRI$", "functions/PPMA_functions.R", getwd()))

# Import libraries
library(doParallel)
library(doRNG)
library(bestglm)
library(ggplot2)
library(dplyr)

# Exclude households with incomes below the Korean legal minimum living cost and zero expenditures (658) 
load("pop.Rdata")
pop <- pop %>% filter(income>658) 
pop <- pop %>% filter(income>0 & exp1>0 & exp2>0)
pop <- pop %>% mutate(geo = as.factor(geo), gender = as.factor(gender), edu = as.factor(edu), age_g = as.factor(age_g), mrg = as.factor(mrg), stt = as.factor(stt))
pop <- pop %>% mutate(income = log(income), exp1 = log(exp1), exp2 = log(exp2), age=age/10)
pop <- pop %>% mutate(edu1 = case_when(edu==1~1, TRUE~0), edu2 = case_when(edu==2~1, TRUE~0), edu3 = case_when(edu==3~1, TRUE~0))
data <- pop

# Assume that households with male household heads are observed, and female household heads are missing
data <- data %>% mutate(R = case_when(gender==1~1, gender==2~0))
data <- data %>% rename(x1 = exp1, x2 = exp2, x3 = age, y = income) %>%
  mutate(x1.2 = x1^2, x2.2 = x2^2, x3.2=x3^2, x13 = x1*x3, x23 = x2*x3, x12 = x1*x2, x123 = x1*x2*x3) %>%
  dplyr::select(x1, x2, x3, x1.2, x2.2, x3.2, x12, x13, x23, x123, y, R)

Full <- mean(data$y)
round(Full, digits = 3)
CC <- mean(data[data$R==1, ]$y)
var_CC <- var(data[data$R==1, ]$y)/nrow(data[data$R==1, ])
lower.bound_CC <- CC - sqrt(var_CC) * qnorm(0.975)
upper.bound_CC <- CC + sqrt(var_CC) * qnorm(0.975)
CI_CC <- c(lower.bound_CC, upper.bound_CC)

prop <- setdiff(colnames(data), c("y", "R"))
data_response <- data; data_response$y <- NULL; data_response <- data_response %>% rename(y=R)
bestglm_out <- bestglm(data_response, family = binomial(link="logit"), TopModels = 5)
res_MAR <- colnames(bestglm_out$BestModel$model)[-1]
VE_MAR_out <- VE_MAR(data=data, XX_names = res_MAR, Y_name = "y", R_name = "R")

best_models <- bestglm_out$BestModels
best_models <- best_models[ , -ncol(best_models)]  
res_xy <- apply(best_models, 1, function(row) {
  names(row)[which(row)]
})
FI_out_MNAR_list <- vector(mode="list", length = length(res_xy))
VE_FI_MNAR_list <- vector(mode="list", length = length(res_xy))
BIC_MNAR <- vector(length = length(res_xy))
p_values <- rep(1, length(res_xy))

set.seed(20240817)
for (j in 1:length(res_xy)) {
  success <- FALSE
  attempt <- 1
  max_attempts <- 5  # Maximum number of attempts to retry
  while (!success && attempt <= max_attempts) {
    tryCatch({
      # Execute FI function
      FI_out_MNAR_list[[j]] <- FI(pop=data, prop=prop, response=c(res_xy[[j]], "y"), phi_curr=rep(0, length(c(res_xy[[j]], "y"))+1), M=200)
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
sig_level <- 0.1
if (p_values[best_MNAR] < sig_level) {
  FI_m_out <- FI_out_MNAR_list[[best_MNAR]]
  FI_m_mu <- FI_m_out$mu.hat
  FI_m_var_out <- VE_NMAR(FI_m_out)
  FI_m_var <- FI_m_var_out$variance_mu
  selected_model <- paste0("NMAR", best_MNAR, collapse = "")
} else {
  FI_m_mu <- VE_MAR_out$mean_PS_MAR
  FI_m_var <- VE_MAR_out$variance_mu
  selected_model <- "MAR"
}

FI_c_out <- FI(pop=data, prop=prop, response=c("x1", "x2", "x3", "y"), phi_curr=rep(0, length(c("x1", "x2", "x3", "y"))+1), M=200)
FI_c_mu <- FI_c_out$mu.hat
FI_c_var_out <- VE_NMAR(FI_c_out)
FI_c_var <- FI_c_var_out$variance_mu

ppm_fit <- FI_out_MNAR_list[[1]]$best_model
summary(ppm_fit)
x_ppm <- predict(ppm_fit, newdata=data)
y_ppm <- data$y; y_ppm[data$R==0] <- NA
cor(x_ppm[data$R==1], y_ppm[data$R==1])
PPM_0 <- mle(x_ppm, y_ppm, 0)
PPM_1 <- mle(x_ppm, y_ppm, 1)
PPM_inf <- mle(x_ppm, y_ppm, Inf)

estimators <- c("CC", "MAR", "PPM_0", "PPM_1", "PPM_inf", "FI_c", "FI_m")
mu_estimate <- c(CC, VE_MAR_out$mean_PS_MAR, PPM_0$muY, PPM_1$muY, PPM_inf$muY, FI_c_mu, FI_m_mu)
var_mu_estimate <- c(var_CC, VE_MAR_out$variance_mu, PPM_0$muYvar, PPM_1$muYvar, PPM_inf$muYvar, FI_c_var, FI_m_var)
lower_ci <- mu_estimate - qnorm(0.975)*sqrt(var_mu_estimate)
upper_ci <-  mu_estimate + qnorm(0.975)*sqrt(var_mu_estimate)
est_result <- data.frame(
  estimator = estimators,
  point_estimate = mu_estimate,
  var_estimate = var_mu_estimate,
  lower_ci = lower_ci,
  upper_ci = upper_ci
)
write.csv(est_result, file = "est_result.csv")

# Plotting
true_value <- mean(data$y)
est_result_plot <- est_result[-1, ]
est_result_plot <- est_result_plot %>% mutate(estimator=factor(estimator, levels=c("MAR", "PPM_0", "PPM_1", "PPM_inf", "FI_c", "FI_m")))
ggplot(est_result_plot, aes(x = estimator, y = point_estimate)) +
  geom_point(size = 3) +  # Plot point estimates
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +  # Plot error bars
  geom_hline(yintercept = true_value, linetype = "dashed", color = "red") +  # Add true value line
  labs(x = "Method",
       y = "y") +
  theme_minimal()  

