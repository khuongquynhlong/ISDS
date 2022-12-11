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


genData <- function(n, sen = "small", seed) {
    set.seed(seed)
    u    <- rnorm(n, 1, 1)
    x1   <- rbinom(n, size=1, prob = plogis(-2 + 0.25*u))
    x2   <- 1 + 0.5*u + rnorm(n, 0, 0.1)
    x3   <- rnorm(n, 1, 1)
    # Set up treatment A and outcome Y under different scenarios
    if (sen == "small") {
        A    <- rbinom(n, size=1, prob = plogis(-2 + 0.25*x1 + 0.5*x3))
        Y_1  <- rbinom(n, size=1, prob = plogis(-2 + 1 + 0.25*x1)) # Y(1)
        Y_0  <- rbinom(n, size=1, prob = plogis(-2 + 0 + 0.25*x1)) # Y(0)
    } else if (sen == "medium") {
        A    <- rbinom(n, size=1, prob = plogis(-2 + 0.25*x1 + 0.5*x3))
        Y_1  <- rbinom(n, size=1, prob = plogis(-2 + 1 + 0.25*x1 + 0.5*x3)) # Y(1)
        Y_0  <- rbinom(n, size=1, prob = plogis(-2 + 0 + 0.25*x1 + 0.5*x3)) # Y(0)
    } else {
        A    <- rbinom(n, size=1, prob = plogis(-2 + 0.25*x1 + 0.5*x3 + 0.5*x2))
        Y_1  <- rbinom(n, size=1, prob = plogis(-2 + 1 + 0.25*x1 + 0.75*x3 + 0.75*x2)) # Y(1)
        Y_0  <- rbinom(n, size=1, prob = plogis(-2 + 0 + 0.25*x1 + 0.75*x3 + 0.75*x2)) # Y(0)
    }
    # Observed outcome
    Y <- Y_1*A + Y_0*(1 - A)
    
    # Return data.frame
    return(data.frame(x1, x2, x3, A, Y, Y_1, Y_0)) # not include U to indicate unobserved variable
}


#---------- Test bias for 3 scenarios (i.e., bias with no adjustment)
# Small confounding
df_small  <- genData(n = 5000, sen = "small", seed = 12345)
trueATE_s <- mean(df_small$Y_1 - df_small$Y_0);trueATE_s 
biasATE_s <- mean(df_small$Y[df_small$A==1]) - mean(df_small$Y[df_small$A==0]); biasATE_s
biasATE_s - trueATE_s

# Medium confounding
df_medium <- genData(n = 5000, sen = "medium", seed = 12345)
trueATE_m <- mean(df_medium$Y_1 - df_medium$Y_0);trueATE_m 
biasATE_m <- mean(df_medium$Y[df_small$A==1]) - mean(df_small$Y[df_medium$A==0]); biasATE_m
biasATE_m - trueATE_m

# large confounding
df_large  <- genData(n = 5000, sen = "large", seed = 12345)
trueATE_l <- mean(df_large$Y_1 - df_large$Y_0);trueATE_l 
trueRR_l  <-  mean(df_large$Y_1)/mean(df_large$Y_0);trueRR_l 
biasATE_l <- mean(df_large$Y[df_small$A==1]) - mean(df_small$Y[df_large$A==0]); biasATE_l
biasATE_l - trueATE_l


#---------- Compare different methods
#===============================================================================
doSimulation <- function(n, sen = "small", seed) {
    if (sen == "small") {
        data     <- genData(n = n, sen = "small", seed = seed)
        trueATE  <- mean(data$Y_1 - data$Y_0)
        formula1 <- Y ~  x1 
        formula2 <- Y ~  A + x1
        formula3 <- A ~  x1
    } else if (sen == "medium") {
        data     <- genData(n = n, sen = "medium", seed = seed)
        trueATE  <- mean(data$Y_1 - data$Y_0)
        formula1 <- Y ~  x1 + x3
        formula2 <- Y ~  A + x1 + x3
        formula3 <- A ~  x1 + x3
    } else {
        data     <- genData(n = n, sen = "large", seed = seed)
        trueATE  <- mean(data$Y_1 - data$Y_0)
        formula1 <- Y ~  x1 + x2 + x3
        formula2 <- Y ~  A + x1 + x2 + x3
        formula3 <- A ~  x1 + x2 + x3
    }
    
    #--- Regression adjustment
    reg_RR <- exp(glm(formula2, family = poisson(link="log"), data = data)$coef[[2]])
    
    
    
    #--- Non-parametric G-formula
    # Create function for G-formula (later for bootstrap CI)
    if (sen == "small") {
        nonpar_gform <- function(data, indices) {
            dat  <- data[indices,]
            pr.l <- prop.table(table(dat$x1))
            ATE  <- ((mean(dat$Y[dat$A==1 & dat$x1==1]) - mean(dat$Y[dat$A==0 & dat$x1==1]))*pr.l[[2]]) + 
                (mean(dat$Y[dat$A==1 & dat$x1==0]) - mean(dat$Y[dat$A==0 & dat$x1==0]))*pr.l[[1]]
            return(ATE)
        }
        
        nonpar_gform_ATE <- nonpar_gform(data, indices=1:nrow(data))
        # # bootstrap CI
        # nonpar_gform_boot <- boot.ci(boot(data, nonpar_gform, 1000), type = "perc", conf = 0.95) 
        # nonpar_gform_ATE_CI <- c(nonpar_gform_boot$percent[4], nonpar_gform_boot$percent[5])
    } else {
        nonpar_gform_ATE <- NA
        nonpar_gform_ATE_CI <- c(NA, NA)
    }
    
    
    
    #--- Parametric G-formula
    par_gform <- function(data, indices) {
        dat   <- data[indices,]
        glm1  <- glm(formula1, family="binomial", data=dat[dat$A==1,])
        glm2  <- glm(formula1, family="binomial", data=dat[dat$A==0,])
        Y1    <- predict(glm1, newdata=data.frame(A = 1, dat[, c("x1", "x2", "x3")]), type="response")
        Y0    <- predict(glm2, newdata=data.frame(A = 0, dat[, c("x1", "x2", "x3")]), type="response")
        ATE   <- mean((Y1) - mean(Y0))
        return(ATE)
    }
    
    par_gform_ATE <- par_gform(data, indices=1:nrow(data))
    # # bootstrap CI
    # par_gform_boot <- boot.ci(boot(data, par_gform, 1000), type = "perc", conf = 0.95)
    # par_gform_ATE_CI <- c(par_gform_boot$percent[4], par_gform_boot$percent[5])
    
    
    
    #--- IPW
    iptw <- function(data, indices) {
        dat    <-  data[indices,]
        ps_mod <- glm(formula3, data=dat, family="binomial")
        pscore <- ifelse(dat$A == 0, 1 - predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
        dat$w  <- 1/pscore
        ATE    <- mean(dat$w*as.numeric(dat$A==1)*dat$Y) - mean(dat$w*as.numeric(dat$A==0)*dat$Y)
        return(ATE)
    }
    
    IPW_ATE <- iptw(data, indices=1:nrow(data))
    # # bootstrap CI
    # IPW_boot <- boot.ci(boot(data, iptw, 1000), type = "perc", conf = 0.95)
    # IPW_ATE_CI <- c(IPW_boot$percent[4], IPW_boot$percent[5])
    
    
    
    #--- Double-robust methods
    double_robust <- function(data, indices) {
        dat    <-  data[indices,]
        # Pscore
        ps_mod <- glm(formula3, data=dat, family="binomial")
        pscore <- ifelse(dat$A == 0, 1 - predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
        dat$w  <- 1/pscore
        # g-formula + IPW (using Horvitz-Thompson estimator with sampling weights)
        glm1   <- svyglm(formula1, family="quasibinomial", design = svydesign(~ 1, weights = ~ dat$w[dat$A==1], data = dat[dat$A==1,]))
        glm2   <- svyglm(formula1, family="quasibinomial", design = svydesign(~ 1, weights = ~ dat$w[dat$A==0], data = dat[dat$A==0,]))
        Y1     <- predict(glm1, newdata = data.frame(A = 1, dat[, c("x1", "x2", "x3")]), type="response")
        Y0     <- predict(glm2, newdata = data.frame(A = 0, dat[, c("x1", "x2", "x3")]), type="response")
        ATE    <- mean(Y1 - Y0)
        return(ATE)
    }
    
    DB_ATE <- double_robust(data, indices=1:nrow(data))
    # # bootstrap CI
    # DB_boot <- boot.ci(boot(data, double_robust, 1000), type = "perc", conf = 0.95)
    # DB_ATE_CI <- c(DB_boot$percent[4], DB_boot$percent[5])
    
    # df <- data.frame(method = c("Non-parametric g-form", "Parametric g-form", "IPW", "Double-robust"),
    #                  est = c(nonpar_gform_ATE, par_gform_ATE, IPW_ATE, DB_ATE),
    #                  lb = c(nonpar_gform_ATE_CI[1], par_gform_ATE_CI[1], IPW_ATE_CI[1], DB_ATE_CI[1]),
    #                  ub = c(nonpar_gform_ATE_CI[2], par_gform_ATE_CI[2], IPW_ATE_CI[2], DB_ATE_CI[2]),
    #                  trueATE = rep(trueATE, 4))
    
    df <- data.frame(method = c("Non-parametric g-form", "Parametric g-form", "IPW", "Double-robust"),
                     estATE = c(nonpar_gform_ATE, par_gform_ATE, IPW_ATE, DB_ATE),
                     trueATE = rep(trueATE, 4))
    return(df)
}


#---------- Run the loop to calculate average bias
result_tab <- data.frame(method = NA, estATE = NA, trueATE = NA, scenario = NA, iter = NA)

for (k in c("small", "medium", "large")) {
    for (i in 1: 1000) {
        df_out          <- doSimulation(n = 1000, sen = k, seed = i)
        df_out$scenario <- k
        df_out$iter     <- i
        result_tab      <- rbind(result_tab, df_out)
    }
}

result_tab %<>% slice(-1) %>% mutate(bias = abs(trueATE - estATE)*100)

result_tab %>% group_by(scenario, method) %>%
    summarise(trueATE_m = mean(trueATE),
              estATE_m = mean(estATE),
              bias_m = mean(bias))


#---------- Visualization

result_tab %>% gather(-c(method, scenario, iter, bias), value = "value", key = "type") %>%
    mutate(scenario = factor(scenario, levels = c("small", "medium", "large"))) %>%
    ggplot(aes(x = value, color = type)) +
    geom_density(alpha = 0.5, size = 1) +
    scale_color_manual(labels = c("Estimate ATE", "True ATE"),
                       values = c("#d7301f", "#2171b5")) +
    facet_grid(scenario ~ method) +
    labs(x = NULL, y = NULL, color = NULL) +
    theme_bw() +
    theme(
        legend.position = "top",
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 14)
    ) -> result_plot

result_tab %>% gather(-c(method, scenario, iter, bias), value = "value", key = "type") %>%
    na.omit() %>% 
    mutate(scenario = factor(scenario, levels = c("small", "medium", "large"))) %>%
    ggplot(aes(x = value, fill = type)) +
    geom_histogram(position = position_dodge(), color = "white", alpha = 0.9) +
    scale_fill_manual(labels = c("Estimate ATE", "True ATE"),
                      values = c("#d7301f", "#2171b5")) +
    facet_grid(scenario ~ method) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_bw() +
    theme(
        legend.position = "top",
        legend.text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 14)
    ) -> result_plot2


png("Fig1_result_plot.png", units="in", width = 16, height = 9, res = 300)
result_plot
dev.off()

png("Fig2_result_plot.png", units="in", width = 16, height = 9, res = 300)
result_plot2
dev.off()











