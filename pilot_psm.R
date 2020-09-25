## import necessary modules
library(Matching)
library(rms)
library(Epi)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
require('survminer')

## check working directory
setwd("C:/Users/user name/Desktop/V‹KŒ¤‹†_•¶Œ£’²¸/New_Research/")
getwd()

## import dataset for PSM
pilot_psm <- read.csv('calculate_ps.csv',
                       header = TRUE)
attach(pilot_psm)
head(pilot_psm)

## preprocessing binomial variances
pilot_psm$Drug_id <- as.factor(Drug_id)
pilot_psm$GCT <- as.factor(GCT)
pilot_psm$Pr_JS_binom <- as.factor(Pr_JS_binom)
pilot_psm$Country_EU <- as.factor(Country_EU)
pilot_psm$Country_JP <- as.factor(Country_JP)
pilot_psm$Country_US <- as.factor(Country_US)
pilot_psm$Orphan <- as.factor(Orphan)
pilot_psm$FIC <- as.factor(FIC)
pilot_psm$ATC_J <- as.factor(ATC_J)
pilot_psm$ATC_L <- as.factor(ATC_L)
pilot_psm$ATC_others <- as.factor(ATC_others)
pilot_psm$Post_2014 <- as.factor(Post_2014)
pilot_psm$Event <- as.factor(Event)

## check data types
sapply(pilot_psm, class)

## check the relarionship between GCT and Event (PMSE)
table(GCT, Event)

## treatment effect with PSM
# caluculate PS by Logistic Regression
model <- glm(GCT ~ Pr_JS_binom +
                   Country_JP +
                   Country_US +
                   Orphan +
                   FIC +
                   ATC_J +
                   ATC_L,
             family = binomial)
# fit model
ps <- model$fitted.values
# caluculate recommended caliper (0.2 * sd(logit(ps)))
logit = function(x) {
  log10(x / (1 - x))
}
r_caliper <- 0.2 * sd(logit(ps))
r_caliper

# caluculate c-statistics of model
model_lrm <- lrm(GCT ~ Pr_JS_binom +
                       Country_JP +
                       Country_US +
                       Orphan +
                       FIC +
                       ATC_J +
                       ATC_L)
model_lrm

# plot ROC curve
Epi::ROC(form = GCT ~ Pr_JS_binom +
                      Country_JP +
                      Country_US +
                      FIC +
                      ATC_J +
                      ATC_L,
         plot = 'ROC')

# matching by ps
mout <- Match(Y = Event,
              Tr = GCT,
              X = ps,
              estimand = 'ATE',
              caliper = r_caliper)
summary(mout)

# check match balance
match_balance <- MatchBalance(GCT ~ Pr_JS_binom +
                                    Country_JP +
                                    Country_US +
                                    FIC +
                                    ATC_J +
                                    ATC_L,
                              match.out = mout,
                              nboots = 1000)

fit <- survfit(Surv(TimetoEvent, Event) ~ GCT, data = pilot_psm)
ggsurvplot(fit,
           fun = function(x)1-x,
           size = 1.4,
           surv.plot.height = 1, 
           surv.scale = "percent",
           censor.shape="|",
           censor.size = 3,
           font.x = c(12, "black"), font.y = c(12, "black"),
           risk.table = TRUE,
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE,
           risk.table.col = "black",
           risk.table.fontsize = 4,
           risk.table.height = 0.25,
           legend.labs = c("No. at risk"),
           palette = "#2E9FDF",
           xlim = c(0,8), ylim = c(0.0,1.0),
           break.time.by = 1,
           ggtheme = theme_bw(),
           data = pilot_psm)



gct.cox <- coxph(Surv(TimetoEvent, Event) ~ GCT)
summary(gct.cox)

gct_ps.cox <- coxph(Surv(TimetoEvent, Event) ~ GCT + ps)
summary(gct_ps.cox)
