library(survival)
library(survminer)
library(KMsurv)
library(MASS)

#load data
data("Melanoma")

#data structure
str(Melanoma)
head(Melanoma)
summary(Melanoma) 

# Create event indicator (1 = died from melanoma, 0 = censored)
Melanoma$event = ifelse(Melanoma$status == 1, 1, 0)
Melanoma$event

# Factor variables; sex=Male, Female  ulcer=Yes, No
Melanoma$sex = factor(Melanoma$sex, levels = c(1, 0), labels = c("Male", "Female"))
Melanoma$ulcer = factor(Melanoma$ulcer, levels = c(0, 1), labels = c("No", "Yes"))

table(Melanoma$status, Melanoma$event)

# Convert time from days to years for easy interpretation
Melanoma$time_years = Melanoma$time / 365.25 #365.25 accounts for leap year

# Create survival object
surv_obj = Surv(time = Melanoma$time_years, event = Melanoma$event)
surv_obj

# Fit Kaplan-Meier estimator for overall survival
KM_fit = survfit(surv_obj ~ 1, data = Melanoma)

# Summary of overall survival
print("Kaplan-Meier Estimate Overall Survival")
summary(KM_fit)

# Plot Kaplan-Meier overall survival curve
ggsurvplot(KM_fit,
           data = Melanoma,
           title = "Kaplan–Meier Survival Estimates for Melanoma Patients",
           xlab = "Time (years)",
           ylab = "Survival Probability",
           conf.int = TRUE,
           risk.table = "nrisk_cumevents",
           risk.table.col = "strata",
           ggtheme = theme_minimal(),
           palette = "purple")

# Median survival time
median_survival = surv_median(KM_fit)

if (is.na(median_survival$median)) {
  print("Median survival time: undefined")
} else {
  print(paste("Median survival time:", round(median_survival$median, 2), "years"))
}


# Fitting Kaplan-Meier by sex
KM_sex = survfit(surv_obj ~ sex, data = Melanoma)
KM_sex
print("Survival by sex")
summary(KM_sex)

# Plot survival curves by sex
ggsurvplot(KM_sex,
           data = Melanoma,
           title = "Survival Comparison: Males vs Females",
           xlab = "Time (years)",
           ylab = "Survival Probability",
           conf.int = TRUE,
           pval = TRUE,
           risk.table = "nrisk_cumevents",
           risk.table.col = "strata",
           tables.height = 0.30,
           legend.labs = c("Male", "Female"),
           ggtheme = theme_minimal(),
           palette = c("blue", "pink"))

logrank_test = survdiff(surv_obj ~ sex, data = Melanoma)
logrank_test

# Calculate p-value
p_value = 1 - pchisq(logrank_test$chisq, df = 1)
p_value

if (p_value > 0.05) {
  print("p > 0.05: Fail to reject the null hypothesis. There's no significant difference in survival between groups.")
} else {
  print("p <= 0.05: Reject the null hypothesis. There is a significant difference in survival between groups.")
}

# Median survival by sex
median_sex = surv_median(KM_sex)
print("Median survival times by sex:")
print(median_sex)


# Median survival by sex
median_sex <- surv_median(KM_sex)

print("Median survival times by sex:")

for (i in 1:nrow(median_sex)) {
  group <- median_sex$strata[i]
  med   <- median_sex$median[i]
  
  if (is.na(med)) {
    message <- paste(group, ": median survival time is undefined.")
  } else {
    message <- paste(group, ": median survival time is", round(med, 2), "years.")
  }
  
  print(message)
}

#Finding Survival Probabilities for 5-15years
summary(KM_sex, times = c(5,10,15))



# Check for missing values
sapply(Melanoma, function(x) sum(is.na(x)))

# Cox model with all variables
cox_all = coxph(surv_obj ~ sex + age + thickness + ulcer + year, data = Melanoma)
summary(cox_all)

# Model selection using stepwise selection
cox_step = step(cox_all, direction = "backward")
summary(cox_step)

#Finding the parsimonious model (based on clinical importance and statistical significance)
cox_test = coxph(surv_obj ~ thickness + ulcer + sex + age + year, data = Melanoma)
summary(cox_test)

# Check proportional hazards assumption for all variables
PH_test = cox.zph(cox_test)
PH_test

# Final parsimonious model (based on clinical importance and statistical significance)
cox_model = coxph(surv_obj ~ thickness + ulcer + sex, data = Melanoma)
summary(cox_model)

# Check proportional hazards assumption
PH_model = cox.zph(cox_model)
PH_model


# Plot Schoenfeld residuals
par(mfrow = c(2, 2))
plot(PH_test)

# Plot survival curves for the final model
ggsurvplot(survfit(cox_model), 
           data = Melanoma,
           title = "Survival Curves from Cox Model",
           palette = c("purple"),
           ggtheme = theme_minimal())

# Calculate risk scores and create risk groups
Melanoma$risk_score = predict(cox_model, type = "risk")
Melanoma$risk_score
Melanoma$risk_group = cut(Melanoma$risk_score, 
                          breaks = quantile(Melanoma$risk_score, c(0, 0.33, 0.67, 1)),
                          labels = c("Low", "Medium", "High"),
                          include.lowest = TRUE)
Melanoma$risk_group
# Plot survival by risk groups
KM_risk = survfit(surv_obj ~ risk_group, data = Melanoma)
ggsurvplot(KM_risk,
           data = Melanoma,
           title = "Survival by Risk Groups (Cox Model)",
           xlab = "Time (years)",
           ylab = "Survival Probability",
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("Low Risk", "Medium Risk", "High Risk"),
           ggtheme = theme_minimal(),
           tables.height = 0.3,
           palette = c("green", "blue", "red"))



