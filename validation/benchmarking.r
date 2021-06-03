library(microbenchmark)
library(dplyr)
library(data.table)

# Selection

MI  <- 1
CHF <- 0
IHD <- 0

current_patient <- data.frame(MI,CHF,IHD)
risk_factors_macrovascular <- c("MI","CHF","IHD")

select1 <- function() {current_patient_macrovascular <- current_patient %>% select(risk_factors_macrovascular)}
select2 <- function() {current_patient_macrovascular <- current_patient[,risk_factors_macrovascular]}

microbenchmark(select1(),select2())

# Unit: microseconds
# expr    min      lq     mean  median      uq    max neval cld
# select1() 1364.3 1424.05 1557.868 1479.05 1642.10 2718.7   100   b
# select2()   10.1   11.65   26.306   13.25   18.55 1065.7   100  a 

# conclusion: avoid using %>%


# Binding rows

bind1 <- function(){
  current_patient <- bind_rows(current_patient, current_patient[risk_factors_macrovascular])
  current_patient
  
}

bind2 <- function(){
  current_patient <- rbindlist(current_patient, current_patient[risk_factors_macrovascular])
  current_patient
  
}

microbenchmark(bind1(),bind2())



