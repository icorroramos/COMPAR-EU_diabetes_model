# ##gen al health states
# 
# df <- expand.grid(MO=1:5,SC=1:5,UA=1:5,PD=1:5,AD=1:5)
# 
# ##apply EQ tariff
# 
# library(eq5d)
# 
# NLutil <- eq5d(scores=df, country = "Netherlands", version = "5L", type ="VT")
# 
# UKvanHout <- eq5d(scores=df, country = "UK", version = "5L", type ="CW")
# 
# UKdevlin <- eq5d(scores=df, country = "England", version = "5L", type ="VT")
# 
# 
# 
# #make transformation algorithms
# 
# Devlin_to_NL <- lm(NLutil ~ UKdevlin)
# 
# VanHout_to_NL <- lm(NLutil ~ UKvanHout)


####5L####

##gen all EQ5D health states
df      <- expand.grid(MO=1:5,SC=1:5,UA=1:5,PD=1:5,AD=1:5)
##apply EQ tariff
library(eq5d)
NL      <- eq5d(scores=df, country = "Netherlands", version = "5L", type ="VT")
UK      <- eq5d(scores=df, country = "UK", version = "5L", type ="CW")
GER     <- eq5d(scores=df, country = "Germany", version = "5L", type ="VT")
FR      <- eq5d(scores=df, country = "France", version = "5L", type ="VT")
SPAIN   <- eq5d(scores=df, country = "Spain", version = "5L", type ="VT")

# QUESTION TO MATTHIJS: why type CW for UK and VT for the others?

#make transformation algorithms
UK_NL   <- lm(NL ~ UK)
UK_GER  <- lm(GER ~ UK)
UK_FR   <- lm(FR ~ UK)
UK_SPAIN<- lm(SPAIN ~ UK)

#store coefficients
UK_to_NL_coeff    <- coefficients(UK_NL)
UK_to_GER_coeff   <- coefficients(UK_GER)
UK_to_FR_coeff    <- coefficients(UK_FR)
UK_to_SPAIN_coeff <- coefficients(UK_SPAIN)

#Enter observed value
observed_UK       <- 0.73

#transform value and report result
Dutch_observed    <- observed_UK*UK_to_NL_coeff[2]+UK_to_NL_coeff[1]
GER_observed      <- observed_UK*UK_to_GER_coeff[2]+UK_to_GER_coeff[1]
FR_observed       <- observed_UK*UK_to_FR_coeff[2]+UK_to_FR_coeff[1]
SPAIN_observed    <- observed_UK*UK_to_SPAIN_coeff[2]+UK_to_SPAIN_coeff[1]

Dutch_observed
GER_observed
FR_observed
SPAIN_observed


##########
### 3L ###
##########

# Beaudet et al. 2014,{#REF} suggests a reference set of utility values which is in line with NICE requirements.
# Therefore, even though not explicitly mentioned in the paper, a 3L version of EQ5D is assumed.

# The package eq5d is used for these calculations.
library(eq5d)

# First a data frame representing all EQ5D health states (all possible combinations between states and levels) is generated.
df_3L      <- expand.grid(MO=1:3,SC=1:3,UA=1:3,PD=1:3,AD=1:3)

# Country-specific tariff is then calculated using the eq5d function available in the R package.
NL_3L      <- eq5d(scores=df_3L, country = "Netherlands", version = "3L", type ="TTO")
UK_3L      <- eq5d(scores=df_3L, country = "UK", version = "3L", type ="TTO")
GER_3L     <- eq5d(scores=df_3L, country = "Germany", version = "3L", type ="TTO")
FR_3L      <- eq5d(scores=df_3L, country = "France", version = "3L", type ="TTO")
SPAIN_3L   <- eq5d(scores=df_3L, country = "Spain", version = "3L", type ="TTO")

# Then, transformation algorithms (linear regression) with UK values as predictors are applied. 
UK_NL_3L    <- lm(NL_3L ~ UK_3L)
UK_GER_3L   <- lm(GER_3L ~ UK_3L)
UK_FR_3L    <- lm(FR_3L ~ UK_3L)
UK_SPAIN_3L <- lm(SPAIN_3L ~ UK_3L)

# The coefficients of these algorithms will be used to predict country-specific utilities derived from UK values.
UK_to_NL_coeff_3L    <- coefficients(UK_NL_3L)
UK_to_GER_coeff_3L   <- coefficients(UK_GER_3L)
UK_to_FR_coeff_3L    <- coefficients(UK_FR_3L)
UK_to_SPAIN_coeff_3L <- coefficients(UK_SPAIN_3L)

# Enter observed UK values
uncomplicated_diabetes_UK_3L <- 0.785
minor_hypo_UK_3L <- uncomplicated_diabetes_UK_3L -0.014
major_hypo_UK_3L <- uncomplicated_diabetes_UK_3L -0.047
excess_BMI_UK_3L <- uncomplicated_diabetes_UK_3L -0.006
catarat_UK_3L <- uncomplicated_diabetes_UK_3L -0.016
mod_retinopathy_UK_3L <- uncomplicated_diabetes_UK_3L -0.040
macular_oedema_UK_3L <- uncomplicated_diabetes_UK_3L -0.040 
sev_retinopathy_UK_3L <- uncomplicated_diabetes_UK_3L -0.070
vision_loss_UK_3L <- uncomplicated_diabetes_UK_3L -0.074
proteinuria_UK_3L <- uncomplicated_diabetes_UK_3L -0.048
renal_transplant_UK_3L <- uncomplicated_diabetes_UK_3L -0.082
haemodialysis_UK_3L <- uncomplicated_diabetes_UK_3L -0.164
peri_dialysis_UK_3L <- uncomplicated_diabetes_UK_3L -0.204
pvd_UK_3L <- uncomplicated_diabetes_UK_3L -0.061
neuropathy_UK_3L <- uncomplicated_diabetes_UK_3L -0.084
ulcer_UK_3L <- uncomplicated_diabetes_UK_3L -0.170
amputation_UK_3L <- uncomplicated_diabetes_UK_3L -0.280
stroke_UK_3L <- uncomplicated_diabetes_UK_3L -0.164
mi_UK_3L <- uncomplicated_diabetes_UK_3L -0.055
ihd_UK_3L <- uncomplicated_diabetes_UK_3L -0.090
hf_UK_3L <- uncomplicated_diabetes_UK_3L -0.108

### ADD CONFIDENCE INTERVALS!

# Create a data.frame object with all utility values
UK.utils <- data.frame(
  UK_3L = c(uncomplicated_diabetes_UK_3L, minor_hypo_UK_3L, major_hypo_UK_3L, excess_BMI_UK_3L, catarat_UK_3L,
            mod_retinopathy_UK_3L, macular_oedema_UK_3L, sev_retinopathy_UK_3L, vision_loss_UK_3L, proteinuria_UK_3L,
            renal_transplant_UK_3L, haemodialysis_UK_3L, peri_dialysis_UK_3L, pvd_UK_3L, neuropathy_UK_3L,
            ulcer_UK_3L, amputation_UK_3L, stroke_UK_3L, mi_UK_3L, ihd_UK_3L, hf_UK_3L))

UK.disutils <- UK.utils - UK.utils[1,1]
UK.disutils

# Predict country-specific values
NL.utils <- predict(UK_NL_3L, newdata = UK.utils)
NL.utils

# Calculate transformed utility decrements
NL.disutils <- NL.utils - NL.utils[1]
NL.disutils

# Add other countries above when needed

# Report results
UK.utils.model <- c(UK.utils[1,1],UK.disutils[-1,1])
NL.utils.model <- round(c(NL.utils[1],NL.disutils[-1]),3)
#GER_observed_3L      <- uncomplicated_diabetes_UK_3L*UK_to_GER_coeff_3L[2]+UK_to_GER_coeff_3L[1]
#FR_observed_3L       <- uncomplicated_diabetes_UK_3L*UK_to_FR_coeff_3L[2]+UK_to_FR_coeff_3L[1]
#SPAIN_observed_3L    <- uncomplicated_diabetes_UK_3L*UK_to_SPAIN_coeff_3L[2]+UK_to_SPAIN_coeff_3L[1]

util_matrix <- data.frame(UK.utils.model,NL.utils.model)
colnames(util_matrix) <- c("UK","NL")

rownames(util_matrix) <- c("uncomplicated_diabetes", "minor_hypo", "major_hypo", "excess_BMI", "catarat",
                           "mod_retinopathy", "macular_oedema", "sev_retinopathy", "vision_loss", "proteinuria",
                           "renal_transplant", "haemodialysis", "peri_dialysis", "pvd", "neuropathy", "ulcer",
                           "amputation", "stroke", "mi", "ihd", "hf")
util_matrix




# Here step by step and attached the file where you just need to change the input parameter (Enter observed value) to get what you need. It works with negatives as well as positives, but if you are faced with disutilities we have to think it through again.
#
# DF is a dataframe with all possible values that patients can report on EQ5D-5L, i.e. 5^5th combinations
# 
# > head(df)
# 
# MO SC UA PD AD
# 
# 1  1  1  1  1  1
# 
# 2  2  1  1  1  1
# 
# 3  3  1  1  1  1
# 
# 4  4  1  1  1  1
# 
# 5  5  1  1  1  1
# 
# 6  1  2  1  1  1
# 
# 
# If your values are sourced from 5L, you can use the code I wrote. If they are sourced from 3L, df needs to be replaced with the following code:
#   
#   df <- expand.grid(MO=1:3,SC=1:3,UA=1:3,PD=1:3,AD=1:3)
# 
# 
# 
# After I have made a dataset (df) with all possible values of EQ-5D, I attach country specific utility values for all possible 5^5th combinations using the eq5d library. (Note, change to 3L if your values are sourced from a 3L valueset)
# 
# 
# 
# library(eq5d)
# 
# NL          <- eq5d(scores=df, country = "Netherlands", version = "5L", type ="VT")
# 
# UK         <- eq5d(scores=df, country = "UK", version = "5L", type ="CW")
# 
# GER       <- eq5d(scores=df, country = "Germany", version = "5L", type ="VT")
# 
# FR          <- eq5d(scores=df, country = "France", version = "5L", type ="VT")
# 
# SPAIN   <- eq5d(scores=df, country = "Spain", version = "5L", type ="VT")
# 
# 
# 
# Then we can estimate the linear models between the country value sets (we do not have Italy do you need it?)
# 
# 
# 
# UK_NL                  <- lm(NL ~ UK)
# 
# UK_GER               <- lm(GER ~ UK)
# 
# UK_FR                  <- lm(FR ~ UK)
# 
# UK_SPAIN          <- lm(SPAIN ~ UK)
# 
# 
# 
# Then indeed we have the transformations rather easily. First we store the coefficients:
#   
#   
#   
#   #store coefficients
#   
#   UK_to_NL_coeff              <- coefficients(UK_NL)
# 
# UK_to_GER_coeff           <- coefficients(UK_GER)
# 
# UK_to_FR_coeff              <- coefficients(UK_FR)
# 
# UK_to_SPAIN_coeff       <- coefficients(UK_SPAIN)
# 
# 
# 
# Then we apply it to our observed EQ-5D-5L value from the meta-analysis you described
# 
# 
# 
# #Enter observed value
# 
# observed_UK                    <- 0.73
# 
# 
# 
# #transform value and report results
# 
# Dutch_observed <- observed_UK*UK_to_NL_coeff[2]+UK_to_NL_coeff[1]
# 
# GER_observed <- observed_UK*UK_to_GER_coeff[2]+UK_to_GER_coeff[1]
# 
# FR_observed <- observed_UK*UK_to_FR_coeff[2]+UK_to_FR_coeff[1]
# 
# SPAIN_observed <- observed_UK*UK_to_SPAIN_coeff[2]+UK_to_SPAIN_coeff[1]
# 
# 
# 
# 
# 
# And report results (as you can see the value sets for Germany and France report much higher values which is reflected below).
# 
# 
# 
# 8)  > Dutch_observed
# 
# 9)         UK
# 
# 10)   0.7433108
# 
# 11)   > GER_observed
# 
# 12)          UK
# 
# 13)   0.9002068
# 
# 14)   > FR_observed
# 
# 15)          UK
# 
# 16)   0.9244409
# 
# 17)   > SPAIN_observed
# 
# 18)          UK
# 
# 19)   0.7918904
# 
# You can enter any value you want at step 6 and it will work. Just remember to check if you need transformations from 3_L or 5_L. I have the full code for 3L and 5L attached for you.
