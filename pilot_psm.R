## import necessary modules
install.packages('Matching')
library(Matching)
install.packages('rms')
library(rms)
install.packages('Epi')
library(Epi)

## check working directory
setwd("C:/Users/Watanabe Yuta/Desktop/V‹KŒ¤‹†_•¶Œ£’²¸/New_Research/")
getwd()

## import dataset for PSM
pilot_psm <- read.csv('pilot_psm_dummy.csv',
                       header = TRUE)
attach(pilot_psm)
head(pilot_psm)

## preprocessing binomial variances
pilot_psm$Drug_id <- as.factor(Drug_id)
pilot_psm$GCT <- as.factor(GCT)
pilot_psm$Pr_JS_binom <- as.factor(Pr_JS_binom)
pilot_psm$Orphan <- as.factor(Orphan)
pilot_psm$Biologic <- as.factor(Biologic)
pilot_psm$FIC <- as.factor(FIC)
pilot_psm$ATC_J <- as.factor(ATC_J)
pilot_psm$ATC_L <- as.factor(ATC_L)
pilot_psm$ATC_others <- as.factor(ATC_others)
pilot_psm$Approval_year_binom <- as.factor(Approval_year_binom)
pilot_psm$Approval_year_2011 <- as.factor(Approval_year_2011)
pilot_psm$Approval_year_2012 <- as.factor(Approval_year_2012)  
pilot_psm$Approval_year_2013 <- as.factor(Approval_year_2013)
pilot_psm$Approval_year_2014 <- as.factor(Approval_year_2014)
pilot_psm$Approval_year_2015 <- as.factor(Approval_year_2015)
pilot_psm$Event <- as.factor(Event)

## check data types
sapply(pilot_psm, class)

## check the relarionship between GCT and Event (PMSE)
table(GCT, Event)

## treatment effect with PSM
# caluculate PS by Logistic Regression
model <- glm(GCT ~ Lag_log +
                   Pr_JS_binom +
                   FIC +
                   ATC_J +
                   ATC_L,
             family = binomial)
# fit model
ps <- model$fitted.values
head(ps)

# caluculate c-statistics of model
model_lrm <- lrm(GCT ~ Lag_log +
                       Pr_JS_binom +
                       FIC +
                       ATC_J +
                       ATC_L)
model_lrm

# plot ROC curve
Epi::ROC(form = GCT ~ Lag_log +
                      Pr_JS_binom +
                      FIC +
                      ATC_J +
                      ATC_L,
         plot = 'ROC')

# matching by ps
mout <- Match(Y = Event,
             Tr = GCT,
             X = ps,
             estimand = 'ATE',
             caliper = 0.05)
summary(mout)

# check match balance
match_balance <- MatchBalance(GCT ~ Lag_log +
                                    Pr_JS_binom +
                                    FIC +
                                    ATC_J +
                                    ATC_L,
                              match.out = mout,
                              nboots = 1000)

## treatment effect without PSM
te_without_psm <- glm(Event ~ GCT,
                      family = binomial('logit'))
summary(te_without_psm)
## treatment effect with PSM
te_with_psm <- glm(Event ~ GCT +
                           ps,
                   family = binomial('logit'))
summary(te_with_psm)