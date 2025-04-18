# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)
library(sktools)
library(plyr)

# Load censor data
censor_data <- fread(file='censor_202401.csv')

# Load sedentary time data
sed <- fread(file='complete_time_sed_with_all_activities_090424.csv')

# Exposures
## Sed hours per day
sed[,':='(sed_hours = sed_daily_total/(7*60))]
## Sed quartiles
sed[,':='(sed_quartiles = factor(quantilize(sed_daily_total,4)))]
## Sed binary
sed[,':='(sed_binary = factor(ifelse(sed_quartiles==4,'High','Low')))]

# Load withdrawals
withdrawals <- fread(file='w7089_20241216.csv')

#### Fix censor data in censor file
# Load center categories
center <- fread(file='center0.csv')
center_lookup <- fread(file='enrollment_correspondences.csv')

# Add center value to dataset
setkey(censor_data,sample_id); setkey(center,sample_id)
censor_data[center,':='(center_code = i.value)]

setkey(censor_data,center_code); setkey(center_lookup,Code)
censor_data[center_lookup,':='(center_location = i.Region)]

# Now correct censor dates based on location
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(center_location=='England',phenotype_censor_date,
                                                         ifelse(center_location=='Scotland',pmin(phenotype_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                                pmin(phenotype_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]

# And set censor date to date of death for those who died
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(!is.na(death_date),pmin(death_date,phenotype_censor_date),phenotype_censor_date),origin='1970-01-01'))]

# Remove missing exposure data
setkey(sed,sample_id); setkey(censor_data,sample_id)
censor_data[sed,':='(activity_group = i.activity_group)]
censor_data <- censor_data[!is.na(activity_group)] 

# Remove withdrawals
censor_data <- censor_data[!(sample_id %in% withdrawals$V1)]

# Merges
setkey(sed,sample_id); setkey(censor_data,sample_id)
censor_data[sed,':='(sed_hours = i.sed_hours, sed_quartiles=i.sed_quartiles, sed_binary = i.sed_binary,
                     bmi = i.bmi,sbp = i.sbp, dbp = i.dbp,
                     accel_age = i.age_accel, sex = i.sex, race = i.race_category_adjust, tob = i.tob,
                     etoh = i.etoh_grams, tdi = i.tdi, employment_status = i.employment_status, 
                     self_health = i.self_health, diet = i.diet, qual_ea = i.qual_ea, accel_date = i.end_date,
                     mvpa_daily_total = i.mvpa_daily_total,prev_msk_accel = i.prev_msk_accel)]

# Scope columns
sed <- censor_data[,c('sample_id','sed_hours','sed_quartiles','sed_binary','accel_age','accel_date',
                      'sex','race','tob','etoh','tdi','employment_status',
                      'self_health','diet','qual_ea','phenotype_censor_date',
                      'bmi','sbp','dbp','mvpa_daily_total')]

# Write out 
write.csv(sed,file='cox_data_sed_121724.csv',row.names = F)

### ADD BLANKING PERIOD SENSITIVITY ANALYSIS
# Create blanked date (2 years after accelerometer)
sed[,blanked_date := accel_date + 365.25*2]

# Write out
write.csv(sed,file='cox_data_sed_blank_2y_121724.csv',row.names = F)

### BMI/SBP/DBP SENSITIVITY ANALYSIS
sed_anthro <- sed[!is.na(sbp) & !is.na(dbp) & !is.na(bmi)] # 89537 - 227 = 89310
write.csv(sed_anthro,file='cox_data_sed_anthro_121724.csv',row.names = F)

# Scope columns
sed_subgroups <- censor_data[,c('sample_id','sed_hours','sed_quartiles','sed_binary','accel_age','accel_date',
                                'sex','race','tob','etoh','tdi','employment_status',
                                'self_health','diet','qual_ea','phenotype_censor_date',
                                'bmi','sbp','dbp','mvpa_daily_total','prev_msk_accel')]

write.csv(sed_subgroups,file='cox_data_sed_subgroups_121724.csv',row.names = F)

