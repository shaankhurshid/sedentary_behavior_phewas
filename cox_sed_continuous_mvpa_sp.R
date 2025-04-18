# Script to process phecode tables

# Depends
library(data.table)
library(stringr)
library(survival)
library(polspline)

# Index of file names
list <- paste0('~/phecode_tables_202401/',list.files('~/phecode_tables_202401'))

# Load exposure/covariate data
sed <- fread('cox_data_sed_121724.csv')
setkey(sed,sample_id)

# Format dates
for (j in (c('accel_date','phenotype_censor_date'))){set(sed,j=j,value=as.Date(sed[[j]],format='%Y-%m-%d'))}

# Init vars
out <- data.table(); n <- 1

# Looping cox model
for (i in list){
  # Merge
  phecode <- NULL; analysis_set <- NULL
  phecode <- read.table(i,sep='\t',header = TRUE); setDT(phecode)
  setkey(phecode,sample_id)
  analysis_set <- sed[phecode,nomatch=0]
  # Format variables
  analysis_set[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
  # Create analysis variables
  analysis_set[,time_to_event := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - accel_date)/365.25,as.numeric(phenotype_censor_date - accel_date)/365.25),
                                        as.numeric(censor_date - accel_date)/365.25)]
  # Remove prevalent disease or no follow-up
  analysis_set <- analysis_set[!is.na(time_to_event) & time_to_event > 0]
  # Define events and follow-up
  disease <- analysis_set$disease[1]
  n_events <- nrow(analysis_set[has_disease==1])
  fu_median <- quantile(analysis_set$time_to_event,0.50); fu_q1 <- quantile(analysis_set$time_to_event,0.25); fu_q3 <- quantile(analysis_set$time_to_event,0.75)
  # If less than 10 cases, abort
  if (n_events < 120){
    hr_sed <- NA; lower_sed <- NA; upper_sed <- NA; chisq_sed <- NA; p_sed <- NA
    result <- data.table(disease,n_events,fu_median,fu_q1,fu_q3,hr_sed,lower_sed,upper_sed,chisq_sed,p_sed)
    out <- rbind(out,result)
    print(paste0("Skipping phenotype ",analysis_set$disease[1]," since < 10 cases"))
    if (n %% 50 == 0){print(paste0("Just finished model ",n," out of ",length(list),"!"))}
    n <- n+1; next}
  # Fit cox model
  model <- coxph(Surv(time_to_event,has_disease) ~ sed_hours + accel_age + sex + race +
                   tob + etoh + tdi + employment_status + self_health + diet + qual_ea + pspline(mvpa_daily_total,df=0), data=analysis_set)
  
  hr_sed <- exp(model$coefficients[1]); lower_sed <- exp(confint(model)[1,1])
  upper_sed <- exp(confint(model)[1,2]); chisq_sed <- summary(model)$coefficients[1,4]
  p_sed <- summary(model)$coefficients[1,6]
  
  result <- data.table(disease,n_events,fu_median,fu_q1,fu_q3,hr_sed,lower_sed,upper_sed,chisq_sed,p_sed)
  out <- rbind(out,result)
  if (n %% 50 == 0){print(paste0("Just finished model ",n," out of ",length(list),"!"))}
  n <- n+1
}

# Save out
write.csv(out,file='cox_sed_continuous_mvpa_spline.csv',row.names=F)