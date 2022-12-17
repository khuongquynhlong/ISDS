library(tidyverse)
library(magrittr)
library(boot)
library(stdReg)
library(survey)

#---------- Simulation set up
#===============================================================================

# U ∼ N (1, 1)
# X1 ∼ Bin(n, px1); with: logit(px1) = −2 + 0.25∗U
# X2 = 1 + 0.5∗U + εx2; with: εx2 ∼ N (0, 0.1)
# X3 ∼ N (1, 1)
# 
# Small confounding
# A ∼ Bin(n, pa); with: logit(pa) = −2 + 0.25∗X1 + 0.5∗X3
# Y ∼ Bin(n, py); with: logit(py) = −2 + A + 0.25∗X1
# 
# Medium confounding
# A ∼ Bin(n, pa); with: logit(pa) = −2 + 0.25∗X1 + 0.5∗X3
# Y ∼ Bin(n, py); with: logit(py) = −2 + A + 0.25∗X1 + 0.5∗X3
# 
# Large confounding
# A ∼ Bin(n, pa); with: logit(pa) = −2 + 0.25∗X1 + 0.5∗X3 + 0.5∗X2
# Y ∼ Bin(n, py); with: logit(py) = −2 + A + 0.25∗X1 + 0.75∗X3 + 0.75∗X2


genData <- function(n, sen = "small") {
    u    <- rnorm(n, 1, 1)
    x1   <- rbinom(n, size=1, prob = plogis(-2 + 0.25*u))
    x2   <- 1 + 0.5*u + rnorm(n, 0, 0.1)
    x3   <- rnorm(n, 1, 1)
    # Set up treatment A and outcome Y under different scenarios
    if (sen == "small") {
        A    <- rbinom(n, size=1, prob = plogis(-2 + x1 + 0.5*x3))
        Y_1  <- rbinom(n, size=1, prob = plogis(-2 + 1 + x1)) # Y(1)
        Y_0  <- rbinom(n, size=1, prob = plogis(-2 + 0 + x1)) # Y(0)
    } else if (sen == "medium") {
        A    <- rbinom(n, size=1, prob = plogis(-2 + x1 + 0.5*x3))
        Y_1  <- rbinom(n, size=1, prob = plogis(-2 + 1 + x1 + 0.5*x3)) # Y(1)
        Y_0  <- rbinom(n, size=1, prob = plogis(-2 + 0 + x1 + 0.5*x3)) # Y(0)
    } else {
        A    <- rbinom(n, size=1, prob = plogis(-2 + x1 + 0.5*x3 + 0.5*x2))
        Y_1  <- rbinom(n, size=1, prob = plogis(-2 + 1 + x1 + 0.75*x3 + 0.75*x2)) # Y(1)
        Y_0  <- rbinom(n, size=1, prob = plogis(-2 + 0 + x1 + 0.75*x3 + 0.75*x2)) # Y(0)
    }
    # Observed outcome
    Y <- Y_1*A + Y_0*(1 - A)
    
    # Return data.frame
    return(data.frame(x1, x2, x3, A, Y, Y_1, Y_0)) # not include U to indicate unobserved variable
}


#---------- Test bias for 3 scenarios (i.e., bias with no adjustment)
# Small confounding
set.seed(12345)
df_small  <- genData(n = 5000, sen = "small")
trueATE_s <- mean(df_small$Y_1 - df_small$Y_0);trueATE_s 
biasATE_s <- mean(df_small$Y[df_small$A==1]) - mean(df_small$Y[df_small$A==0]); biasATE_s
# biasATE_s - trueATE_s

# Medium confounding
df_medium <- genData(n = 5000, sen = "medium")
trueATE_m <- mean(df_medium$Y_1 - df_medium$Y_0);trueATE_m 
biasATE_m <- mean(df_medium$Y[df_small$A==1]) - mean(df_small$Y[df_medium$A==0]); biasATE_m
# biasATE_m - trueATE_m

# large confounding
df_large  <- genData(n = 5000, sen = "large")
trueATE_l <- mean(df_large$Y_1 - df_large$Y_0);trueATE_l 
trueRR_l  <-  mean(df_large$Y_1)/mean(df_large$Y_0);trueRR_l 
biasATE_l <- mean(df_large$Y[df_small$A==1]) - mean(df_small$Y[df_large$A==0]); biasATE_l
# biasATE_l - trueATE_l


biasATE_s - trueATE_s; biasATE_m - trueATE_m; biasATE_l - trueATE_l


#---------- Compare different methods
#===============================================================================
doSimulation <- function(n, sen = "small") {
    if (sen == "small") {
        data     <- genData(n = n, sen = "small")
        formula1 <- Y ~  x1 
        formula2 <- Y ~  A + x1
        formula3 <- A ~  x1
    } else if (sen == "medium") {
        data     <- genData(n = n, sen = "medium")
        formula1 <- Y ~  x1 + x3
        formula2 <- Y ~  A + x1 + x3
        formula3 <- A ~  x1 + x3
    } else if (sen == "large") {
        data     <- genData(n = n, sen = "large")
        formula1 <- Y ~  x1 + x2 + x3
        formula2 <- Y ~  A + x1 + x2 + x3
        formula3 <- A ~  x1 + x2 + x3
    } else if (sen == "large2") {
        data      <- genData(n = n, sen = "large")
        formula1  <- Y ~  x1 + x2 + x3     # Model correctly specified
        formula2  <- Y ~  A + x1 + x2 + x3 # Model correctly specified
        formula3  <- A ~  x1 + x2 + x3     # Model correctly specified
        formula1f <- Y ~  1     # For model misspecification
        formula2f <- Y ~  A     # For model misspecification
        formula3f <- A ~  1     # For model misspecification
    }
    
    #--- True ATE and RR
    trueATE  <- mean(data$Y_1 - data$Y_0)
    trueRR   <- mean(data$Y_1)/mean(data$Y_0)
    
    
    #--- No adjustment bias
    max_bias_ATE <- mean(data$Y[data$A==1]) - mean(data$Y[data$A==0])
    max_bias_RR  <- mean(data$Y[data$A==1])/mean(data$Y[data$A==0])
    
    
    #--- Regression adjustment
    # Create function for Regression adjustment (later for bootstrap CI)
    reg_adj_f <- function(data, indices) {
        dat   <- data[indices,]
        if (sen == "large2") {
            reg_mod <- glm(formula2f, family = poisson(link="log"), data = data) # Model misspecification
        } else {
            reg_mod <- glm(formula2, family = poisson(link="log"), data = data) # Model correctly specified
        }
        Y1    <- predict(reg_mod, newdata=data.frame(A = 1, df_large[, c("x1", "x2", "x3")]), type="response")
        Y0    <- predict(reg_mod, newdata=data.frame(A = 0, df_large[, c("x1", "x2", "x3")]), type="response")
        ATE   <- mean((Y1) - mean(Y0))
        RR    <- exp(coef(reg_mod)[[2]])
        return(c(ATE, RR))
    }
    
    reg_adj <- reg_adj_f(data, indices=1:nrow(data))
    
    #--- Non-parametric G-formula
    # Create function for G-formula (later for bootstrap CI)
    if (sen == "small") {
        nonpar_gform_f <- function(data, indices) {
            dat  <- data[indices,]
            pr.l <- prop.table(table(dat$x1))
            # ATE  <- ((mean(dat$Y[dat$A==1 & dat$x1==1]) - mean(dat$Y[dat$A==0 & dat$x1==1]))*pr.l[[2]]) + 
            #     (mean(dat$Y[dat$A==1 & dat$x1==0]) - mean(dat$Y[dat$A==0 & dat$x1==0]))*pr.l[[1]]
            Y1    <- mean(dat$Y[dat$A==1 & dat$x1==1])*pr.l[[2]] + mean(dat$Y[dat$A==1 & dat$x1==0])*pr.l[[1]]
            Y0    <- mean(dat$Y[dat$A==0 & dat$x1==1])*pr.l[[2]] + mean(dat$Y[dat$A==0 & dat$x1==0])*pr.l[[1]]
            ATE   <- Y1 - Y0
            RR    <- Y1/Y0
            return(c(ATE, RR))
        }
        nonpar_gform <- nonpar_gform_f(data, indices=1:nrow(data))
        # # bootstrap CI
        # nonpar_gform_boot <- boot.ci(boot(data, nonpar_gform, 1000), type = "perc", conf = 0.95) 
        # nonpar_gform_ATE_CI <- c(nonpar_gform_boot$percent[4], nonpar_gform_boot$percent[5])
    } else {
        nonpar_gform <- c(NA, NA)
    }
    
    
    
    #--- Parametric G-formula
    par_gform_f <- function(data, indices) {
        dat   <- data[indices,]
        if (sen == "large2") {
            glm1  <- glm(formula1f, family="binomial", data=dat[dat$A==1,]) # Model misspecification
            glm2  <- glm(formula1f, family="binomial", data=dat[dat$A==0,]) # Model misspecification
        } else {
            glm1  <- glm(formula1, family="binomial", data=dat[dat$A==1,])  # Model correctly specified
            glm2  <- glm(formula1, family="binomial", data=dat[dat$A==0,])  # Model correctly specified
        }
        Y1    <- predict(glm1, newdata=data.frame(A = 1, dat[, c("x1", "x2", "x3")]), type="response")
        Y0    <- predict(glm2, newdata=data.frame(A = 0, dat[, c("x1", "x2", "x3")]), type="response")
        ATE   <- mean((Y1) - mean(Y0))
        RR    <- mean((Y1)/mean(Y0))
        return(c(ATE, RR))
    }
    
    par_gform <- par_gform_f(data, indices=1:nrow(data))
    # # bootstrap CI
    # par_gform_boot <- boot.ci(boot(data, par_gform, 1000), type = "perc", conf = 0.95)
    # par_gform_ATE_CI <- c(par_gform_boot$percent[4], par_gform_boot$percent[5])
    
    
    
    #--- IPW
    iptw_f <- function(data, indices) {
        dat    <-  data[indices,]
        if (sen == "large2") {
            ps_mod <- glm(formula3f, data=dat, family="binomial") # Model misspecification
        } else {
            ps_mod <- glm(formula3, data=dat, family="binomial") # Model correctly specified
        }
        pscore <- ifelse(dat$A == 0, 1 - predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
        dat$w  <- 1/pscore
        # ATE    <- mean(dat$w*as.numeric(dat$A==1)*dat$Y) - mean(dat$w*as.numeric(dat$A==0)*dat$Y)
        Y1     <- mean(dat$w*as.numeric(dat$A==1)*dat$Y)
        Y0     <- mean(dat$w*as.numeric(dat$A==0)*dat$Y)
        ATE    <- Y1 - Y0
        RR     <- Y1/Y0
        return(c(ATE, RR))
    }
    
    IPW <- iptw_f(data, indices=1:nrow(data))
    # # bootstrap CI
    # IPW_boot <- boot.ci(boot(data, iptw, 1000), type = "perc", conf = 0.95)
    # IPW_ATE_CI <- c(IPW_boot$percent[4], IPW_boot$percent[5])
    
    
    
    #--- Double-robust methods 
    # For misspecification of treatment assignment (i.e., unbias estimate)
    double_robust_A <- function(data, indices) {
        dat    <-  data[indices,]
        # Pscore\
        if (sen == "large2") {
            ps_mod <- glm(formula3f, data=dat, family="binomial") # Model misspecification
        } else {
            ps_mod <- glm(formula3, data=dat, family="binomial") # Model correctly specified
        }
        pscore <- ifelse(dat$A == 0, 1 - predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
        dat$w  <- 1/pscore
        # g-formula + IPW (using Horvitz-Thompson estimator with sampling weights)
        glm1   <- svyglm(formula1, family="quasibinomial", design = svydesign(~ 1, weights = ~ dat$w[dat$A==1], data = dat[dat$A==1,]))
        glm2   <- svyglm(formula1, family="quasibinomial", design = svydesign(~ 1, weights = ~ dat$w[dat$A==0], data = dat[dat$A==0,]))
        Y1     <- predict(glm1, newdata = data.frame(A = 1, dat[, c("x1", "x2", "x3")]), type="response")
        Y0     <- predict(glm2, newdata = data.frame(A = 0, dat[, c("x1", "x2", "x3")]), type="response")
        ATE    <- mean((Y1) - mean(Y0))
        RR     <- mean((Y1)/mean(Y0))
        return(c(ATE, RR))
    }
    
    # For misspecification of outcome model (i.e., unbias estimate)
    double_robust_Y <- function(data, indices) {
        dat    <-  data[indices,]
        # Pscore
        ps_mod <- glm(formula3, data=dat, family="binomial") # Model correctly specified
        pscore <- ifelse(dat$A == 0, 1 - predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
        dat$w  <- 1/pscore
        # g-formula + IPW (using Horvitz-Thompson estimator with sampling weights)
        if (sen == "large2") {
            glm1   <- svyglm(formula1f, family="quasibinomial", 
                             design = svydesign(~ 1, weights = ~ dat$w[dat$A==1], data = dat[dat$A==1,])) # Model misspecification
            glm2   <- svyglm(formula1f, family="quasibinomial", 
                             design = svydesign(~ 1, weights = ~ dat$w[dat$A==0], data = dat[dat$A==0,])) # Model misspecification
            Y1     <- predict(glm1, newdata = data.frame(A = 1, dat[, c("x1", "x2", "x3")]), type="response")
            Y0     <- predict(glm2, newdata = data.frame(A = 0, dat[, c("x1", "x2", "x3")]), type="response")
            ATE    <- mean((Y1) - mean(Y0))
            RR     <- mean((Y1)/mean(Y0))
            return(c(ATE, RR))
        } else {
            return(c(NA, NA))
        }
    }
    # For misspecification of both treatment assignment and outcome models (i.e., Bias estimate)
    double_robust_AY <- function(data, indices) {
        dat    <-  data[indices,]
        if (sen == "large2") {
            # Pscore
            ps_mod <- glm(formula3f, data=dat, family="binomial") # Model misspecification
            pscore <- ifelse(dat$A == 0, 1 - predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
            dat$w  <- 1/pscore
            # g-formula + IPW (using Horvitz-Thompson estimator with sampling weights)
            glm1   <- svyglm(formula1f, family="quasibinomial", 
                             design = svydesign(~ 1, weights = ~ dat$w[dat$A==1], data = dat[dat$A==1,])) # Model misspecification
            glm2   <- svyglm(formula1f, family="quasibinomial", 
                             design = svydesign(~ 1, weights = ~ dat$w[dat$A==0], data = dat[dat$A==0,])) # Model misspecification
            Y1     <- predict(glm1, newdata = data.frame(A = 1, dat[, c("x1", "x2", "x3")]), type="response")
            Y0     <- predict(glm2, newdata = data.frame(A = 0, dat[, c("x1", "x2", "x3")]), type="response")
            ATE    <- mean((Y1) - mean(Y0))
            RR     <- mean((Y1)/mean(Y0))
            return(c(ATE, RR))
        } else {
            return(c(NA, NA))
        }
    }
    DB_A <- double_robust_A(data, indices=1:nrow(data))
    DB_Y <- double_robust_Y(data, indices=1:nrow(data))
    DB_AY <- double_robust_AY(data, indices=1:nrow(data))
    
    # # bootstrap CI
    # DB_boot <- boot.ci(boot(data, double_robust, 1000), type = "perc", conf = 0.95)
    # DB_ATE_CI <- c(DB_boot$percent[4], DB_boot$percent[5])
    
    # df <- data.frame(method = c("Non-parametric g-form", "Parametric g-form", "IPW", "Double-robust"),
    #                  est = c(nonpar_gform_ATE, par_gform_ATE, IPW_ATE, DB_ATE),
    #                  lb = c(nonpar_gform_ATE_CI[1], par_gform_ATE_CI[1], IPW_ATE_CI[1], DB_ATE_CI[1]),
    #                  ub = c(nonpar_gform_ATE_CI[2], par_gform_ATE_CI[2], IPW_ATE_CI[2], DB_ATE_CI[2]),
    #                  trueATE = rep(trueATE, 4))
    
    df <- data.frame(method = c("1.Regression adjustment", "2.Non-parametric g-form", 
                                "3.Parametric g-form", "4.IPW", "5.DB-unbias1",
                                "6.DB-unbias2", "7.DB-bias"),
                     trueATE = rep(trueATE, 7),
                     max_biasATE = rep(max_bias_ATE, 7),
                     estATE = c(reg_adj[1], nonpar_gform[1], par_gform[1], IPW[1], DB_A[1], 
                                DB_Y[1], DB_AY[1]),
                     trueRR = rep(trueRR, 7),
                     max_biasRR = rep(max_bias_ATE, 7),
                     estRR = c(reg_adj[2], nonpar_gform[2], par_gform[2], IPW[2], DB_A[2], 
                                DB_Y[2], DB_AY[2]))
    return(df)
}


doSimulation(n = 1000, sen = "large")


#---------- Run the loop to calculate average bias
result_tab <- data.frame(method = NA, trueATE = NA, max_biasATE = NA, estATE = NA, 
                         trueRR = NA, max_biasRR = NA, estRR = NA, scenario = NA, iter = NA)

set.seed(12345)
for (k in c("small", "medium", "large", "large2")) {
    for (i in 1: 1000) {
        df_out          <- doSimulation(n = 5000, sen = k)
        df_out$scenario <- k
        df_out$iter     <- i
        result_tab      <- rbind(result_tab, df_out)
    }
}

result_tab1 <- result_tab %>% slice(-1) %>% mutate(biasATE = abs(trueATE - estATE),
                                     max_biasATE = abs(trueATE - max_biasATE),
                                     biasRR = abs(trueRR/estRR),
                                     max_biasRR = abs(trueRR/max_biasRR))


result_tab %>% group_by(scenario, method) %>%
    summarise(trueATE_m = mean(trueATE),
              estATE_m = mean(estATE),
              biasATE_m = mean(biasATE),
              max_biasATE_m = mean(max_biasATE),
              trueRR_m = mean(trueRR),
              estRR_m = mean(estRR),
              biasRR_m = mean(biasRR),
              max_biasRR_m = mean(max_biasRR)) %>% ungroup() %>%
    mutate(scenario = factor(scenario, levels = c("small", "medium", "large", "large2"))) %>%
    knitr::kable()


#---------- Visualization

# result_tab %>% gather(-c(method, scenario, iter, bias), value = "value", key = "type") %>%
#     mutate(scenario = factor(scenario, levels = c("small", "medium", "large"))) %>%
#     ggplot(aes(x = value, color = type)) +
#     geom_density(alpha = 0.5, size = 1) +
#     scale_color_manual(labels = c("Estimate ATE", "True ATE"),
#                        values = c("#d7301f", "#2171b5")) +
#     facet_grid(scenario ~ method) +
#     labs(x = NULL, y = NULL, color = NULL) +
#     theme_bw() +
#     theme(
#         legend.position = "top",
#         legend.text = element_text(size = 14),
#         strip.text = element_text(face = "bold", size = 14)
#     ) -> result_plot
# 
# result_tab %>% gather(-c(method, scenario, iter, bias), value = "value", key = "type") %>%
#     na.omit() %>% 
#     mutate(scenario = factor(scenario, levels = c("small", "medium", "large"))) %>%
#     ggplot(aes(x = value, fill = type)) +
#     geom_histogram(position = position_dodge(), color = "white", alpha = 0.9) +
#     scale_fill_manual(labels = c("Estimate ATE", "True ATE"),
#                       values = c("#d7301f", "#2171b5")) +
#     facet_grid(scenario ~ method) +
#     labs(x = NULL, y = NULL, fill = NULL) +
#     theme_bw() +
#     theme(
#         legend.position = "top",
#         legend.text = element_text(size = 14),
#         strip.text = element_text(face = "bold", size = 14)
#     ) -> result_plot2
# 
# 
# png("Fig1_result_plot.png", units="in", width = 16, height = 9, res = 300)
# result_plot
# dev.off()
# 
# png("Fig2_result_plot.png", units="in", width = 16, height = 9, res = 300)
# result_plot2
# dev.off()











