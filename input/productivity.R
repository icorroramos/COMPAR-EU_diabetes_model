
df <- read.csv("C:/Users/Isaac/Dropbox/COMPAR/panel_prod.csv", sep=",")

head(df)

df$age_scale <- scale(df$age)

library("lme4")


# Changes to formula: gender.f coded as factor (this is needed when 0-1 is not used), age scaled to avoid convergence problems
productivity_pred <- glmer(formula = unempl ~ factor(female) 
                           + age_scale 
                           + newheart 
                           + newstroke 
                           + newCHKD + (1 | mergeid),
                           data = df, family = binomial,
                           control = glmerControl(optimizer = c("bobyqa"))
                           )

# These are the coefficients
productivity_pred

#restart

ss <- getME(productivity_pred,c("theta","fixef"))
m2 <- update(productivity_pred,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))


#https://www.google.com/search?client=firefox-b-d&sxsrf=ALeKk03j-V2saD2jojqjjYQGH9nCwQZgzQ%3A1596197782203&ei=lgskX930C4G1kwXzr4foDA&q=%22Model+failed+to+converge+with+max%7Cgrad%7C%22&oq=%22Model+failed+to+converge+with+max%7Cgrad%7C%22&gs_lcp=CgZwc3ktYWIQAzICCAAyBAgAEB4yBAgAEB4yBAgAEB4yBAgAEB4yBggAEAUQHjIGCAAQBRAeOgoIIxCuAhCwAxAnUOLYiwFY4tiLAWC524sBaAFwAHgAgAFKiAFKkgEBMZgBAKABAqABAaoBB2d3cy13aXrAAQE&sclient=psy-ab&ved=0ahUKEwjdl7KOvPfqAhWB2qQKHfPXAc0Q4dUDCAs&uact=5
#https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
#https://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/
#https://github.com/lme4/lme4/issues/489
#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2018q3/027015.html 
#https://stats.stackexchange.com/questions/88960/lme4-glmer-problems-with-offset
#http://rstudio-pubs-static.s3.amazonaws.com/4951_47ef851e46ab4de68d9bab5dcfabbf82.html
#https://stackoverflow.com/questions/27990948/error-message-error-in-fnx-downdated-vtv-is-not-positive-definite
#https://github.com/lme4/lme4/issues/173
#http://shape-of-code.coding-guidelines.com/2013/07/31/i-made-a-mistake-please-dont-shoot-me/
#https://github.com/lme4/lme4/issues/483
#https://www.google.com/search?client=firefox-b-d&q=%22unable+to+evaluate+scaled+gradient%22
#https://stackoverflow.com/questions/53834754/scaling-predictors-in-lme4-glmer-doesnt-resolve-eigenvalue-warnings-neither-do




# To test the equation I predicted the value for one patient (I took the first one from df)
patient <- df[1,]
df[1,]

# These are the characteristics
# productivity gender.f age heartattack_new2.f stroke_new2.f CKD_new2.f      mergeid
# 1        0        2  63                  0             0          0 AT-000674-01

# The equation results in a log odds which can be obtained from here:
log.oods.productivity <- predict(productivity_pred, patient)
log.oods.productivity
# -2.228942 is the log odds for that particular patient to receive informal care in one week (if I understood correctly)

# Now the log odds are converted into probability with the following formula
p.productivity <- exp(log.oods.productivity)/(1+exp(log.oods.productivity))
p.productivity


rbinom(100,1,p.productivity)












# library(GLMMadaptive)
# fm <- mixed_model(unempl ~ factor(female), random = ~ 1 | mergeid, data = df, family = binomial())
# summary(fm)
# You can set the initial values in the mixed_model() function of my
# GLMMadaptive package using the initial_values argument. For more info,
# check: https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.html
