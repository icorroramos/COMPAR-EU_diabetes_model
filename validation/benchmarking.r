MI  <- 1
CHF <- 0
IHD <- 0

current_patient <- data.frame(MI,CHF,IHD)
risk_factors_macrovascular <- c("MI","CHF")

select1 <- function() {current_patient_macrovascular <- current_patient %>% select(risk_factors_macrovascular)}
select2 <- function() {current_patient_macrovascular <- current_patient[,risk_factors_macrovascular]}

microbenchmark(select1(),select2())