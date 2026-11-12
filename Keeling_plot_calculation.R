# Keeling calculations - BOKU
# Created by Lauren Gillespie July 2021 (post doc on the eLTER Plus project under the supervision of Dr. Eugenio Diaz-Pines)
# R version 4.0.4 (2021-02-15)

# User needs: 
# 1) to change the baseURL (line 23) to the URL of the folder where the N2O data subfolders are stored. 
#    In the example here, the link is to the BOKU Drive. If the BOKU Drive is used, the link needs to be generated through a 'share' link
#    otherwise the script cannot access the data.
# 2) indicate where the output files should be saved (i.e. 'export_files_to'; line 29)
#    I was unable to get it to upload direclty to the BOKU Drive. So, for the moment, it needs to be a location on a computer.
# 3) Once that is done, select the entire script with control 'a' and then press 'Run' or control enter.
# 4) After the script has run, the time it took for the script to run and where the files were exported to will appear in the console.

# The exported file is named on line 179 'Data_combined_file_name', the name is 'N2O_data_combined_cleaned_' and the start and end date of the data. 


#cleanning of workspace and time recording
rm(list = ls())
start_time <- Sys.time()

# import the data from BOKU drive
# URL for the folder where the N2O data subfolders are stored
baseURL = "https://drive.boku.ac.at/d/0dc1a30a1e484c44bd58/"

folders = c("files/?p=%2FGFT_output_N2O_ppb%2Fn2o_", "files/?p=%2FGFT_output_d15N%2Fn2o_", "files/?p=%2FGFT_output_d15Nalpha%2Fn2o_", 
            "files/?p=%2FGFT_output_d15Nbeta%2Fn2o_", "files/?p=%2FGFT_output_d18O%2Fn2o_")

# set the location where files should be exported to (i.e. the working directory)
export_files_to = "C:/Users/BOKU_/OneDrive/Desktop/Post_docBOKU/EXAFOR/Data/Keeling_test"

# create function to call data from BOKU drive
import_data = function(x) {
  chamberID = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16")
  suffixURL = ".txt&dl=1"
  url2fetch <- paste(baseURL, x, chamberID, suffixURL, sep='')
  data <- do.call("rbind",lapply(url2fetch, FUN=function(file){read.table(url(file), header=F, sep="\t")}))
  return(data)
}

# apply function
N2O_data = lapply(folders, import_data)

# rename the dataframes to be clearer
names(N2O_data) <- c("d15N_flux", "d15N_bulk", "d15N_alpha", "d15N_beta", "d18O")

# Primary links used to create above script
# https://www.r-bloggers.com/2015/03/getting-data-from-an-online-source/
# https://psychwire.wordpress.com/2011/06/05/testing-different-methods-for-merging-a-set-of-files-into-a-dataframe/

# add column heads
colnames <- c("timestamp",	"chamber",	"fluxLin",	"r2Lin",	"nmeas",	"ncgas",	"temp",	"conc1",	"conc2",
              "conc3",	"conc4",	"conc5",	"conc6",	"conc7",	"conc8",	"conc9",	"conc10",	"conc11",
              "conc12",	"conc13",	"conc14",	"m.dt1",	"m.dt2",	"m.dt3",	"m.dt4",	"m.dt5",	"m.dt6",
              "m.dt7",	"m.dt8",	"m.dt9",	"m.dt10",	"m.dt11",	"m.dt12",	"m.dt13",	"m.dt14", "inverseconc")

N2O_data = lapply(N2O_data, setNames, colnames)

# Remove rows that do not have 14 concentrations (i.e. 'all_vars') > 0, a zero potentially indicates a problem 
library(dplyr)
conc <- N2O_data[[1]] %>% names() %>% purrr::keep(~ stringr::str_detect(.,"^conc"))
N2O_data_zeros_removed <- lapply(N2O_data, function(x) filter_at(x, vars(all_of(conc)), all_vars(.!=0))) # https://www.r-bloggers.com/2019/10/filtering-a-data-frame-by-condition-on-multiple-columns/

# Convert the data from wide to long format to facilitate application of the linear model
N2O_data_long = lapply(N2O_data_zeros_removed, tidyr::gather, concentration, value, conc1:conc14, factor_key=TRUE)

cols_to_keep = c("chamber", "timestamp", "concentration", "value")
N2O_data_long = lapply(N2O_data_long, dplyr::select, all_of(cols_to_keep)) # remove superfluous columns to clean things up

# Add a column 'variable' identifying which variable (i.e. N2O gas flux or isotope)/the file name is associated with the data
# This is in preparation for when the dataframes are bound into a single dataframe (they are currently stored as a list)
N2O_data_long <- mapply(cbind, N2O_data_long, "variable"=names(N2O_data_long), SIMPLIFY=F)

# Time to bind the N2o ppb gas flux and all isotope (d15N 'total', d15Nalpha, d15beta, and d18O) data together
N2O_data_long_combined = dplyr::bind_rows(N2O_data_long)

# Set the timestamp column as a date.time object
N2O_data_long_combined$timestamp = as.POSIXct(N2O_data_long_combined$timestamp)

# spread the data from long to wide by 'variable' (aka N2O flux and associated isotopoes)
library(dplyr)
N2O_data_wide_combined = N2O_data_long_combined  %>% 
  tidyr::spread(key=variable, value=value)


# Create a column for the inverse of the N2O flux
# I originally tried to include this in the fit_by_timestamp function below to avoid having this extra line ("%s ~ (1/n2o_02)"), but it didn't work.
N2O_data_wide_combined$d15N_flux_inverse = ifelse(N2O_data_wide_combined$d15N_flux>0, 1/N2O_data_wide_combined$d15N_flux, NA) 

# Filter out the NA's from the columns: "d15N_bulk", "d15N_alpha", "d15N_beta", "d18O"
varlist = c( "d15N_bulk", "d15N_alpha", "d15N_beta", "d18O")
N2O_data_wide_combined_withoutNA = N2O_data_wide_combined %>% filter_at(vars(all_of(varlist)), all_vars(!is.na(.)))
#sapply(N2O_data_wide_combined_withoutNA, function(x) sum(is.na(x))) # to check the number of NAs in the dataframe

# combine 'chamber' and 'timestamp' into a single variable so that it can be used in the 'split' function below.
# this allows the data to be run by chamber and by timestamp. 
N2O_data_wide_combined_withoutNA$chamber_timestamp = paste(N2O_data_wide_combined_withoutNA$chamber, N2O_data_wide_combined_withoutNA$timestamp, sep = "_")


# Apply lm model
# uses code structure as suggested here: https://stackoverflow.com/questions/42445741/run-linear-models-by-group-over-list-of-variables-in-r
# want to group by timestamp and by chamber and apply the linear model across the four isotope columns (i.e. bulk, alpha, beta, and d18O)


# This function is to retrieve the linear model statistics
sumfun <- function(x) {
  model_stats = c(
    terms = as.character(summary(x)$terms), # lm structure
    slope = round(as.numeric(coef(x)[2]), 5), # slope
    intercept = round(as.numeric(summary(x)$coefficients[1,1]), 5), # intercept
    #t_value = summary(x)$coefficients[ 1, 3], # t-value
    R2 = round(as.numeric(summary(x)$r.squared), 5),
    adj_R2 = round(as.numeric(summary(x)$adj.r.squared), 5),
    SE = round(as.numeric(summary(x)$coefficients[1, 2]), 5), # standard error of the intercept
    #RMSE =sqrt(mean(resid(x) ^ 2, na.rm = TRUE)), # RMSE
    p_value = round(as.numeric(summary(x)$coefficients[1, 4]), 5) # p_value of the intercept
  )
  return(model_stats)
}

# Function to run the lm over the list of variables (i.e. varlist), varlist defined in the results = lapply() code immediately after function is created
fit_by_timestamp_by_chamber = function(by_timestamp_by_chamber, varlist) {
  models <- lapply(varlist, function(x) {
    # apply the modeling function to our list of air variables
    form = as.formula(sprintf("%s ~ d15N_flux_inverse", x))
    lm(form, data = by_timestamp_by_chamber, na.action = na.omit)  # linear model with air variable substituted
  })
  return(as.data.frame(t(sapply(models, sumfun))))
}

results = lapply(split(N2O_data_wide_combined_withoutNA, N2O_data_wide_combined_withoutNA$chamber_timestamp), FUN = fit_by_timestamp_by_chamber, varlist = c( "d15N_bulk", "d15N_alpha", "d15N_beta", "d18O"))

# add the timestamp as a column
results <- mapply(cbind, results, "timestamp"=names(results), SIMPLIFY=F)

# Time to bind the N2o ppb gas flux and all isotope (d15N 'total', d15Nalpha, d15beta, and d18O) data together
results_combined = results %>%
  dplyr::bind_rows() 

# separate the timestamp_chamber column into two columns
results_combined = results_combined %>% tidyr::separate(timestamp, c("chamber", "timestamp"), sep="_")

# change numeric values from 'chr' to 'num'
results_combined[, c("slope", "intercept", "R2", "adj_R2", "SE", "p_value")] <- lapply(results_combined[,c("slope", "intercept", "R2", "adj_R2", "SE", "p_value")], as.numeric) 


# spread the data from long to wide by 'gas' (aka N2O flux and associated isotopoes)
library(data.table)
results_combined_wide = dcast(setDT(results_combined), timestamp+chamber~terms2, value.var=c("slope", "intercept", "R2", "adj_R2", "SE", "p_value"))

# to get the site preference 'SP' column, the beta intercept should be subtracted from the alpha intercept (i.e. intercept_d15Nalpha_02 - intercept_d15Nbeta_02)
results_combined_wide$site_preference = results_combined_wide$intercept_d15N_alpha - results_combined_wide$intercept_d15N_beta


# add  the N2Olin flux ppb and r2linfrom the original dataframe (dset by timestamp and chamber) 
results_combined_wide$chamber = as.integer(results_combined_wide$chamber)
data.table::setDT(results_combined_wide)[N2O_data[[1]], c("fluxLin", "r2Lin") := .(fluxLin, r2Lin), on=c("chamber", "timestamp")]

# reorder columns to group the slope, r2, SE, etc. variables together by isotope
results_combined_wide = base::cbind(
  select(results_combined_wide, c("chamber", "timestamp", "fluxLin", "r2Lin", "site_preference")),
  select(results_combined_wide,contains("bulk")),
  select(results_combined_wide,contains("alpha")),
  select(results_combined_wide,contains("beta")),
  select(results_combined_wide,contains("d18O")))


# set working directory, where files should be saved
setwd(export_files_to)

# export the Keeling results to desired folder
# First create a file name so that there is the beginning timestamp date (i.e. min(timestamp)) and end timestamp date (i.e. mxa(timestamp))
Keeling_results_file_name = paste("Keeling_results_", as.Date(min(results_combined_wide$timestamp)), "_to_", as.Date(max(results_combined_wide$timestamp)),".csv", sep = "")
write.csv2(results_combined_wide, Keeling_results_file_name, row.names=FALSE)

# export the cleaned df
# don't need the chamber_timestamp variable now so best get rid of it first
N2O_data_wide_combined_withoutNA_without_chamber_timestamp<- dplyr::select(N2O_data_wide_combined_withoutNA, -c(chamber_timestamp))
# Create a file name so that there is the beginning timestamp date (i.e. min(timestamp)) and end timestamp date (i.e. mxa(timestamp))
Data_combined_file_name = paste("N2O_data_combined_cleaned_", as.Date(min(N2O_data_wide_combined_withoutNA_without_chamber_timestamp$timestamp)), "_to_", as.Date(max(N2O_data_wide_combined_withoutNA_without_chamber_timestamp$timestamp)),".csv", sep = "")
write.csv2(N2O_data_wide_combined_withoutNA_without_chamber_timestamp, Data_combined_file_name, row.names=FALSE)

# calculate the time it takes to run entire script
end_time <- Sys.time()
cat("Script run time: ")
end_time - start_time
cat("Files exported to:", getwd())
