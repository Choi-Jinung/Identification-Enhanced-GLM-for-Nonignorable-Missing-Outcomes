rm(list=ls())

# Import functions
source(sub("real_data/ExitPoll$", "functions/FI_functions.R", getwd()))
source(sub("real_data/ExitPoll$", "functions/PPMA_functions.R", getwd()))

# Import libraries
library(bestglm)
library(ggplot2)
library(readxl)
library(dplyr)
library(patchwork)

true_vote_rate <- read_excel("true_vote_rate.xlsx", sheet = "true_vote_rate")
cands <- read_excel("true_vote_rate.xlsx", sheet = "candidates")
true_elected <- ifelse(true_vote_rate$true_value>0.5, cands$cand1, cands$cand2)
exitpoll_data <- read.csv("election19_Seoul.csv", header = T, fileEncoding = "euc-kr")
district_groups <- split(exitpoll_data, exitpoll_data$district)
district_groups[[39]] <- NULL
output_list <- vector("list", length = length(district_groups))
summary(as.factor(true_elected))

# district <- district_groups[[46]]; cand1 <- 1; cand2 <- 2; true_value <- true_vote_rate$true_value[46]
get_result_one_district <- function(district, cand1, cand2, true_value) {
  # Data Cleaning
  district <- district %>% filter(sex!=0 & age!=0) # exclude units with missing gender and age information
  # gender
  district <- district %>% mutate(Z = case_when(sex==1~1, TRUE~0) ) 
  # age
  district <- district %>% mutate(I_20 = case_when(age==1~1, TRUE~0),
                                  I_30 = case_when(age==2~1, TRUE~0),
                                  I_40 = case_when(age==3~1, TRUE~0),
                                  I_50 = case_when(age==4~1, TRUE~0),
                                  I_60 = case_when(age==5~1, TRUE~0) ) 
  # polynomial coding
  II <- cbind(district$I_20, district$I_30, district$I_40, district$I_50, district$I_60) %*% contr.poly(5)
  district <- district %>% mutate(L = II[, 1], Q = II[, 2], C = II[, 3], Fourth = II[, 4])
  district <- district %>% mutate(ZL = Z*L, ZQ = Z*Q, ZC = Z*C, ZFourth = Z*Fourth)
  
  # response variable
  # exclude units who did not vote for candidate 1 or 2
  district <- district %>% filter(cand==cand1 | cand==cand2 | cand==0)
  # response variable: cand1 -> 1, cand2 -> 0, 0 -> 4 (NA)
  district <- district %>% mutate(y = case_when(cand==cand1~1, cand==cand2~0, cand==0~4))
  
  # missing indicator
  district <- district %>% mutate(R = case_when(y==4~0, TRUE~1))
  data <- district %>% select(L, Q, C, Fourth, Z, ZL, ZQ, ZC, ZFourth, y, R)
  
  # CC analysis
  CC <- mean(data[data$R==1, ]$y)
  var_CC <- var(data[data$R==1, ]$y)/nrow(data[data$R==1, ])
  
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
    FI_m_mu <- FI_out_MNAR_list[[best_MNAR]]$mu.hat
    FI_m_var <- VE_FI_MNAR_list[[best_MNAR]]$variance_mu
    selected_model <- paste0("NMAR", best_MNAR, collapse = "")
  } else {
    FI_m_mu <- VE_MAR_out$mean_PS_MAR
    FI_m_var <- VE_MAR_out$variance_mu
    selected_model <- "MAR"
  }
  is_MNAR <- (selected_model!="MAR")
  FI_c_out <- FI(pop=data, prop=prop, response=c("L", "Z", "ZL", "y"), phi_curr=rep(0, length(c("L", "Z", "ZL", "y"))+1), M=2)
  FI_c_mu <- FI_c_out$mu.hat
  FI_c_var_out <- VE_NMAR(FI_c_out)
  FI_c_var <- FI_c_var_out$variance_mu
  
  estimators <- c("CC", "MAR", "FI_c", "FI_m")
  mu_estimate <- c(CC, VE_MAR_out$mean_PS_MAR, FI_c_mu,  FI_m_mu)
  bias <- mu_estimate - true_value
  var_mu_estimate <- c(var_CC, VE_MAR_out$variance_mu, FI_c_var, FI_m_var)
  lower_ci_95 <- mu_estimate - qnorm(0.975)*sqrt(var_mu_estimate)
  upper_ci_95 <-  mu_estimate + qnorm(0.975)*sqrt(var_mu_estimate)
  length_95 <- upper_ci_95 - lower_ci_95
  lower_ci_99 <- mu_estimate - qnorm(0.995)*sqrt(var_mu_estimate)
  upper_ci_99 <-  mu_estimate + qnorm(0.995)*sqrt(var_mu_estimate)
  length_99 <- upper_ci_99 - lower_ci_99
  est_result <- data.frame(
    estimator = estimators,
    point_estimate = mu_estimate,
    var_estimate = var_mu_estimate,
    bias = bias,
    lower_ci_95 = lower_ci_95,
    upper_ci_95 = upper_ci_95,
    length_95 = length_95,
    lower_ci_99 = lower_ci_99,
    upper_ci_99 = upper_ci_99,
    length_99 = length_99
  )
  boolean_to_cand <- c("TRUE" = cand1, "FALSE" = cand2)
  winner_predicted <- boolean_to_cand[as.character(mu_estimate > 0.5)]
  ci_contain_true_95 <- (lower_ci_95 < true_value & true_value < upper_ci_95)
  ci_contain_true_99 <- (lower_ci_99 < true_value & true_value < upper_ci_99)
  out <- list(VE_FI_MAR=VE_MAR_out, est_result=est_result,
              FI_out_MNAR_list=FI_out_MNAR_list, VE_FI_MNAR_list=VE_FI_MNAR_list, BIC_MNAR=BIC_MNAR, p_values=p_values,
              is_MNAR=is_MNAR, winner_predicted=winner_predicted, ci_contain_true_95=ci_contain_true_95, ci_contain_true_99=ci_contain_true_99)
  return(out)
}

sig_level <- 0.10
for (i in 1:length(district_groups)) {
  output_list[[i]] <- get_result_one_district(district = district_groups[[i]], 
                                              cand1 = cands$cand1[i], 
                                              cand2 = cands$cand2[i], 
                                              true_value = true_vote_rate$true_value[i])
}

true_seats <- summary(as.factor(true_elected))
true_seats_MNAR <- summary(as.factor(true_elected[district_is_MNAR]))
true_seats_MNAR <- c(true_seats_MNAR, 0)
names(true_seats) <- c("A", "B", "Others")
names(true_seats_MNAR) <- c("A", "B", "Others")  

winner_predicted <- do.call(rbind, lapply(output_list, FUN = function(x) x$winner_predicted))
colnames(winner_predicted) <- c("CC", "MAR", "FI_c", "FI_m")
winner_predicted_combined <- cbind(winner_predicted)
predicted_seats <- t(apply(winner_predicted_combined, MARGIN = 2, FUN = function(x) summary(as.factor(x))))
colnames(predicted_seats) <- c("A", "B", "Others")
rownames(predicted_seats) <- c("CC", "MAR", "FI_c", "FI_m")

winner_predicted_MNAR <- winner_predicted[district_is_MNAR, ]
predicted_seats_MNAR <- rbind(summary(as.factor(winner_predicted_MNAR[, "CC"])),
                              summary(as.factor(winner_predicted_MNAR[, "MAR"])),
                              c(summary(as.factor(winner_predicted_MNAR[, "FI_c"])), 0),
                              c(summary(as.factor(winner_predicted_MNAR[, "FI_m"])), 0))
colnames(predicted_seats_MNAR) <- c("A", "B", "Others")
rownames(predicted_seats_MNAR) <- c("CC", "MAR", "FI_c", "FI_m")

bias <- apply(do.call(rbind, lapply(output_list, FUN=function(x) x$est_result$bias)), MARGIN = 2, FUN = mean)
length_95 <- apply(do.call(rbind, lapply(output_list, FUN=function(x) x$est_result$length_95)), MARGIN = 2, FUN = mean)
length_99 <- apply(do.call(rbind, lapply(output_list, FUN=function(x) x$est_result$length_99)), MARGIN = 2, FUN = mean)

district_names <- names(district_groups)
district_names <- substr(district_names, nchar(district_names)-1, nchar(district_names))
district_is_MNAR <- unlist(lapply(output_list, FUN = function(x) x$is_MNAR))
MNAR_districts <- district_names[district_is_MNAR]

coverage_95 <- apply(do.call(rbind, lapply(output_list, FUN = function(x) x$ci_contain_true_95)), MARGIN = 2, FUN = sum)
coverage_99 <- apply(do.call(rbind, lapply(output_list, FUN = function(x) x$ci_contain_true_99)), MARGIN = 2, FUN = sum)

coverage_95_MNAR <- apply(do.call(rbind, lapply(output_list, FUN = function(x) x$ci_contain_true_95))[district_is_MNAR, ], MARGIN = 2, FUN = sum)
coverage_99_MNAR <- apply(do.call(rbind, lapply(output_list, FUN = function(x) x$ci_contain_true_99))[district_is_MNAR, ], MARGIN = 2, FUN = sum)

exitpoll_summary <- cbind(predicted_seats, predicted_seats_MNAR, bias, length_95, coverage_95, coverage_95_MNAR, length_99, coverage_99, coverage_99_MNAR)
write.csv(exitpoll_summary, file = "exitpoll_summary.csv")
true_seats
true_seats_MNAR
exitpoll_summary

data_plot <- data.frame(
  District = rep(paste("District", MNAR_districts), each = 4),
  Method = rep(c("CC", "MAR", "FI_c", "FI_m"), length(MNAR_districts)),
  Actual = rep(true_vote_rate$true_value[district_is_MNAR], each=4),
  Estimate = unlist(lapply(output_list[district_is_MNAR], FUN = function(x) x$est_result$point_estimate)),
  Lower_95 = unlist(lapply(output_list[district_is_MNAR], FUN = function(x) x$est_result$lower_ci_95)),
  Upper_95 = unlist(lapply(output_list[district_is_MNAR], FUN = function(x) x$est_result$upper_ci_95)),
  Lower_99 = unlist(lapply(output_list[district_is_MNAR], FUN = function(x) x$est_result$lower_ci_99)),
  Upper_99 = unlist(lapply(output_list[district_is_MNAR], FUN = function(x) x$est_result$upper_ci_99))
)
data_plot <- data_plot %>% mutate(Method=factor(Method, level=c("CC", "MAR", "FI_c", "FI_m")))
write.csv(data_plot, "data_plot.csv")

ggplot(data_plot, aes(x = Method, y = Estimate, color = District)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower_95, ymax = Upper_95), width = 0.2, linewidth=0.8) +
  geom_hline(aes(yintercept = Actual), linetype = "dashed", color = "black") +
  facet_wrap(~ District) +
  labs(x = "Method",
       y = "Vote percentage",
       color = "District") +
  theme_minimal()
