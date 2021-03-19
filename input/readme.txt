Input files to be saved here. 

UKPDS equations used for all countries. 

# Macrovascular complications: UKPDS paper ESM Table 4. Macrovascular complications include CHF, IHD, MI and stroke.
# The risk factors used to predict macrovascular complications are the following (notation used in the model, definition from 
# UKPDS paper ESM Table 2):

# 1. AFRO: 1 == afro-caribbean ethnicity, 0 == Otherwise. HR referent == Caucasian
# 2. AGE.DIAG: age in years at diagnosis of diabetes. HR per year increase in age at diagnosis
# 3. FEMALE: 1 == female, 0 == male. HR referent == male
# 4. INDIAN: 1 == asian-indian ethnicity, 0 == Otherwise. HR referent == Caucasian
# 5. ATFIB: 1 == attrial fibrillation, 0 == otherwise. Defined from Minnesota codes 831 and 833. HR referent == no ATFIB
# 6. BMI: Body mass index (m/kg^2) measured continuously. HR per unit increase in BMI
# 7. eGFR: Estimated glomerular filtration rate (ml/min/1.73m^2) from modification of diet in renal disease (MDRD) formula.
#          Continuously measured and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase 
# 8. eGFR60less: Same as eGFR. Continuous spline (knot at 60) and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase if < 60
# 9. HbA1c: percentage HbA1c measured continuously. HR per 1% increase in HbA1c
# 10. HDL: high density lipoprotein cholesterol (mmol/l) measured continuously and MULTIPLIED by 10. HR per 0.1 mg/dL increase in HDL
# 11. LDL: low density lipoprotein cholesterol (mmol/l) measured continuously and MULTIPLIED by 10. HR per 0.1 mg/dL increase in LDL
# 12. LDL35more: same as LDL. Continuous spline (knot at 3.5) and MULTIPLIED by 10. HR per 0.1 mg/dL increase in LDL if > 3.5 mg/dL
# 13. MMALB: presence of micro- or macro-albuminuria. 1 == urine albumin >= 50 mg/l, 0 == otherwise. HR referent == no albuminuria
# 14. PVD: 1 == peripheral vascular disease, 0 == otherwise. Defined from presence of intermittent claudication or ankle brachial
#          pressure index < 0.9. HR referent == no presence/evidence of PVD
# 15. SBP: systolic blood pressure (mm Hg) measured continuously and further DIVIDED by 10. HR per 10 mm Hg increase in SBP --> 10 or 100?
# 16. SMOKER: 1 == current smoker, 0 == otherwise. HR referent == not current smoker
# 17. WBC: white blood cell count measured continuously. HR per 1x10^6/ml increase
# 18. AMP.HIST: 1 == history of amputation, 0 == otherwise. HR referent == no prior amputation
# 19. CHF.HIST: 1 == history of CHF, 0 == otherwise. HR referent == no prior CHF
# 20. IHD.HIST: 1 == history of IHD, 0 == otherwise. HR referent == no prior IHD
# 21. STROKE.HIST: 1 == history of stroke, 0 == otherwise. HR referent == no prior stroke
# 22. ULCER.HIST: 1 == history of diabetic ulcer, 0 == otherwise. HR referent == no prior diabetic ulcer

# Microvascular complications: UKPDS paper ESM Table 5. Microvascular complications include blindness, diabetic ulcer amputation and renal failure.
# The risk factors (the ones that are not defined above for macrovascular complications) used to predict microvascular complications 
# are the following (notation used in the model, definition from UKPDS paper ESM Table 2):

# 23. eGFR60more: Same as eGFR. Continuous spline (knot at 60) and further DIVIDED by 10. HR per 10 ml/min/1.73m^2 increase if > 60
# 24. HAEM: haemoglobin g/dL measured continuously. HR per 1 g/dL increase
# 25. HEART.R: heart rate (beats per minute) determined from inspiration/expiration RR on ECG. Continuously measured and further DIVIDED by 10.
#              HR per 10 bpm increase.
# 26. BLIND.HIST: 1 == history of blindness, 0 == otherwise. HR referent == no prior blindness

# Risk of death: UKPDS paper ESM Table 6. Four equations depending on events and history. The risk factors (the ones that are not 
# defined above for macro/microvascular complications) are the following (notation used in the model, definition from UKPDS paper ESM Table 2):

# 27. YEAR: duration of diabetes in years measured continuously. HR per year increase in duration of diabetes. 
#           Note it increases + 1 after each simulated year.
# 28. BMI1: 1 == BMI < 18.5 m/kg^2, 0 == otherwise. HR referent 18.5 m/kg^2 <= BMI < 25 m/kg^2. It depends thus on the variable BMI.
# 29. BMI3: 1 == BMI >= 25 m/kg^2, 0 == otherwise. HR referent 18.5 m/kg^2 <= BMI < 25 m/kg^2. It depends thus on the variable BMI.
# 30. CURR.AGE: current age in years measured continuously. HR per year increase in current age. Note CURR.AGE = AGE.DIAG + YEAR at the 
#               beginning of the simulation. Also after each year in the simulation if YEAR is properly updated.
# 31. AMP1.EVENT: 1 == first amputation, 0 == otherwise. HR referent == no first amputation event
# 32. AMP2.EVENT: 1 == second amputation, 0 == otherwise. HR referent == no second amputation event
# 33. IHD.EVENT: 1 == IHD, 0 == otherwise. HR referent == no IHD event
# 34. MI.EVENT: 1 == MI, 0 == otherwise. HR referent == no MI event
# 35. MI.HIST: 1 == history of MI, 0 == otherwise. HR referent == no prior MI
# 36. RENAL.EVENT: 1 == renal failure event, 0 == otherwise. HR referent == no renal failure event
# 37. RENAL.HIST: 1 == history of renal failure, 0 == otherwise. HR referent == no prior renal failure
# 38. STROKE.EVENT: 1 == stroke, 0 == otherwise. HR referent == no stroke event
