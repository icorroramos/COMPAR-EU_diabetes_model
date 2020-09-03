# Changes to formula: gender.f coded as factor (this is needed when 0-1 is not used), age scaled to avoid convergence problems
# Note df in the equation below is your DID_single_final
informal_care_pred <- glmer(formula = inf_care ~ factor(gender.f) + scale(age) + heartattack_new2.f + 
                              stroke_new2.f + CKD_new2.f + (1 | mergeid), data = df, family = binomial)
# These are the coefficients
informal_care_pred

# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: inf_care ~ factor(gender.f) + scale(age) + heartattack_new2.f +      stroke_new2.f + CKD_new2.f + (1 | mergeid)
# Data: df
# AIC       BIC    logLik  deviance  df.resid 
# 11372.833 11423.598 -5679.417 11358.833     10419 
# Random effects:
#   Groups  Name        Std.Dev.
# mergeid (Intercept) 1.244   
# Number of obs: 10426, groups:  mergeid, 5228
# Fixed Effects:
#   (Intercept)    factor(gender.f)2           scale(age)  heartattack_new2.f1       stroke_new2.f1          CKD_new2.f1  
# -1.8322               0.4369               0.6316               0.4805               0.5426               0.2713  

# To test the equation I predicted the value for one patient (I took the first one from df)
patient <- df[1,]
df[1,]

# These are the characteristics
# inf_care gender.f age heartattack_new2.f stroke_new2.f CKD_new2.f      mergeid
# 1        0        2  63                  0             0          0 AT-000674-01

# The equation results in a log odds which can be obtained from here:
log.oods.inf_care <- predict(informal_care_pred, patient)
log.oods.inf_care
# -2.228942 is the log odds for that particular patient to receive informal care in one week (if I understood correctly)

# Now the log odds are converted into probability with the following formula
p.inf_care <- exp(log.oods.inf_care)/(1+exp(log.oods.inf_care))
p.inf_care
# 0.09718139 is the probability for that particular patient to receive informal care in one week (if I understood correctly)

# Since the simulation model works with annual probabilities, we can use
rbinom(1,52,p.inf_care)
# to draw at random how many weeks per year that patient will receive informal care

