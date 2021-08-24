
## Cole Oost
## Project Script 2021

# Estimating Long-Term Phenological Responses to Climate in the Ephemeral Flora of 
# Southern California Deserts ----

## ======================================================================================

# Loading the Climate Data ----

# Clear the global environment #
rm(list = ls())

# Set the working directory #
setwd("/Volumes/Personal_Files/climate_data")

# Load necessary packages #
library(dplyr)
library(tidyr)
library(ggplot2)
library(cruts)
library(rgdal)
library(sf)
library(sp)
library(units)
library(ncdf4)
library(raster)

# Load in the precipitation and temperature data sets. 
pre <- brick("cru_ts4.04.1901.2019.pre.dat.nc", varname = "pre")
tmp <- brick("cru_ts4.04.1901.2019.tmp.dat.nc", varname = "tmp")

# Crop the temperature and precipitation data sets to focus on the American Southwest.
southwest.area <- extent(-125, -105, 25, 45)
pre <- crop(pre, southwest.area)
tmp <- crop(tmp, southwest.area)

# Create a function that will more efficiently load in the future precipitation files.
load_future_pre <- function(rastlist){
  
  # Organize the list numerically by month.
  starter <- c(rastlist[1], rastlist[5:12])
  ender <- rastlist[2:4]
  rastlist <- c(starter, ender)
  
  # Import all of the raster files in that folder.
  allrasters <- lapply(rastlist, raster)
  
  # Turn the raster list into a stack of layers and reduce the spatial extent of the data.
  strawberry <- raster::stack(allrasters)
  southwest.area <- extent(-125, -105, 25, 45)
  captaincrunch <- crop(strawberry, southwest.area)
}

# The four priority GCM's for CA Fourth Climate Assessment are HadGEM2-ES (HE),
# CNRM-CM5 (CN), CanESM2 (Not Available), and MIROC5 (MC).

# Store all raster files for HE/precipitation/RCP 4.5 as a list, then use function.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/HE/he45pr50")
rastlist <- list.files(pattern = '.tif', all.files = TRUE, full.names = FALSE)
he45pr50 <- load_future_pre(rastlist)

# Store all raster files for CN/precipitation/RCP 4.5 as a list, then use function.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/CN/cn45pr50")
rastlist <- list.files(pattern = '.tif', all.files = TRUE, full.names = FALSE)
cn45pr50 <- load_future_pre(rastlist)

# Store all raster files for MC/precipitation/RCP 4.5 as a list, then use function.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/MC/mc45pr50")
rastlist <- list.files(pattern = '.tif', all.files = TRUE, full.names = FALSE)
mc45pr50 <- load_future_pre(rastlist)

# Store all raster files for HE/precipitation/RCP 8.5 as a list, then use function.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/HE/he85pr50")
rastlist <- list.files(pattern = '.tif', all.files = TRUE, full.names = FALSE)
he85pr50 <- load_future_pre(rastlist)

# Store all raster files for CN/precipitation/RCP 8.5 as a list, then use function.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/CN/cn85pr50")
rastlist <- list.files(pattern = '.tif', all.files = TRUE, full.names = FALSE)
cn85pr50 <- load_future_pre(rastlist)

# Store all raster files for MC/precipitation/RCP 8.5 as a list, then use function.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/MC/mc85pr50")
rastlist <- list.files(pattern = '.tif', all.files = TRUE, full.names = FALSE)
mc85pr50 <- load_future_pre(rastlist)

# Load future temperature data for HE/temperature/RCP 4.5.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/HE/he45bi50")
he45tmp50 <- raster("he45bi5011.tif")
he45tmp50 <- crop(he45tmp50, southwest.area)

# Load future temperature data for CN/temperature/RCP 4.5.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/CN/cn45bi50")
cn45tmp50 <- raster("cn45bi5011.tif")
cn45tmp50 <- crop(cn45tmp50, southwest.area)

# Load future temperature data for MC/temperature/RCP 4.5.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/MC/mc45bi50")
mc45tmp50 <- raster("mc45bi5011.tif")
mc45tmp50 <- crop(mc45tmp50, southwest.area)

# Load future temperature data for HE/temperature/RCP 8.5.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/HE/he85bi50")
he85tmp50 <- raster("he85bi5011.tif")
he85tmp50 <- crop(he85tmp50, southwest.area)

# Load future temperature data for CN/temperature/RCP 8.5.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/CN/cn85bi50")
cn85tmp50 <- raster("cn85bi5011.tif")
cn85tmp50 <- crop(cn85tmp50, southwest.area)

# Load future temperature data for MC/temperature/RCP 8.5.
setwd("/Volumes/Personal_Files/BAZONKA/WorldClim/MC/mc85bi50")
mc85tmp50 <- raster("mc85bi5011.tif")
mc85tmp50 <- crop(mc85tmp50, southwest.area)

## ======================================================================================

# Camissonia campestris Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Camissonia_campestris.csv")

# Load dplyr and tidyr.
library("dplyr")
library("tidyr")

# Create two blank vectors for later use.
human_observations_total <- c()
preserved_specimens_total <- c()

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
         basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
         basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Filter the data set to show only entries that were flowering.
vf_data <- filter(vf_data, VF == "YES")

# How many entries per decade?
my_summary <- count(vf_data, Decade, sort = TRUE) 
my_summary

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package 
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Camissonia campestris", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Camissonia_campestris_results <- results

## ======================================================================================

# Camissonia campestris Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
               mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Camissonia_campestris_results <- cbind(Camissonia_campestris_results, means_df)

## ======================================================================================

# Monoptilon bellioides Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Monoptilon_bellioides.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Monoptilon bellioides", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Monoptilon_bellioides_results <- results

## ======================================================================================

# Monoptilon bellioides Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Monoptilon_bellioides_results <- cbind(Monoptilon_bellioides_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(Camissonia_campestris_results, Monoptilon_bellioides_results)

## ======================================================================================

# Mohavea confertiflora Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Mohavea_confertiflora.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Mohavea confertiflora", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Mohavea_confertiflora_results <- results

## ======================================================================================

# Mohavea confertiflora Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Mohavea_confertiflora_results <- cbind(Mohavea_confertiflora_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Mohavea_confertiflora_results)

## ======================================================================================

# Atrichoseris platyphylla Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Atrichoseris_platyphylla.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Atrichoseris platyphylla", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Atrichoseris_platyphylla_results <- results

## ======================================================================================

# Atrichoseris platyphylla Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Atrichoseris_platyphylla_results <- cbind(Atrichoseris_platyphylla_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Atrichoseris_platyphylla_results)

## ======================================================================================

# Langloisia setosissima Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Langloisia_setosissima.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Langloisia setosissima", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Langloisia_setosissima_results <- results

## ======================================================================================

# Langloisia setosissima Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Langloisia_setosissima_results <- cbind(Langloisia_setosissima_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Langloisia_setosissima_results)

## ======================================================================================

# Chaenactis fremontii Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Chaenactis_fremontii.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Chaenactis fremontii", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Chaenactis_fremontii_results <- results

## ======================================================================================

# Chaenactis fremontii Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Chaenactis_fremontii_results <- cbind(Chaenactis_fremontii_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Chaenactis_fremontii_results)

## ======================================================================================

# Nama demissum Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Nama_demissum.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Nama demissum", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Nama_demissum_results <- results

## ======================================================================================

# Nama demissum Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Nama_demissum_results <- cbind(Nama_demissum_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Nama_demissum_results)

## ======================================================================================

# Chylismia brevipes Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Chylismia_brevipes.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Chylismia brevipes", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Chylismia_brevipes_results <- results

## ======================================================================================

# Chylismia brevipes Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Chylismia_brevipes_results <- cbind(Chylismia_brevipes_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Chylismia_brevipes_results)

## ======================================================================================

# Eremalche rotundifolia Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eremalche_rotundifolia.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eremalche rotundifolia", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eremalche_rotundifolia_results <- results

## ======================================================================================

# Eremalche rotundifolia Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eremalche_rotundifolia_results <- cbind(Eremalche_rotundifolia_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eremalche_rotundifolia_results)

## ======================================================================================

# Mentzelia involucrata Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Mentzelia_involucrata.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Mentzelia involucrata", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Mentzelia_involucrata_results <- results

## ======================================================================================

# Mentzelia involucrata Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Mentzelia_involucrata_results <- cbind(Mentzelia_involucrata_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Mentzelia_involucrata_results)

## ======================================================================================

# Calochortus kennedyi Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Calochortus_kennedyi.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Calochortus kennedyi", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Calochortus_kennedyi_results <- results

## ======================================================================================

# Calochortus kennedyi Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Calochortus_kennedyi_results <- cbind(Calochortus_kennedyi_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Calochortus_kennedyi_results)

## ======================================================================================

# Gilia cana Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Gilia_cana.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Gilia cana", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Gilia_cana_results <- results

## ======================================================================================

# Gilia cana Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Gilia_cana_results <- cbind(Gilia_cana_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Gilia_cana_results)

## ======================================================================================

# Eriophyllum lanosum Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eriophyllum_lanosum.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eriophyllum lanosum", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eriophyllum_lanosum_results <- results

## ======================================================================================

# Eriophyllum lanosum Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eriophyllum_lanosum_results <- cbind(Eriophyllum_lanosum_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eriophyllum_lanosum_results)

## ======================================================================================

# Lupinus odoratus Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Lupinus_odoratus.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Lupinus odoratus", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Lupinus_odoratus_results <- results

## ======================================================================================

# Lupinus odoratus Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Lupinus_odoratus_results <- cbind(Lupinus_odoratus_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Lupinus_odoratus_results)

## ======================================================================================

# Geraea canescens Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Geraea_canescens.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Geraea canescens", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Geraea_canescens_results <- results

## ======================================================================================

# Geraea canescens Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Geraea_canescens_results <- cbind(Geraea_canescens_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Geraea_canescens_results)
View(total_results)

## ======================================================================================

# Eschscholzia minutiflora Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eschscholzia_minutiflora.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eschscholzia minutiflora", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eschscholzia_minutiflora_results <- results

## ======================================================================================

# Eschscholzia minutiflora Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eschscholzia_minutiflora_results <- cbind(Eschscholzia_minutiflora_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eschscholzia_minutiflora_results)

## ======================================================================================

# Rafinesquia neomexicana Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Rafinesquia_neomexicana.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Rafinesquia neomexicana", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Rafinesquia_neomexicana_results <- results

## ======================================================================================

# Rafinesquia neomexicana Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Rafinesquia_neomexicana_results <- cbind(Rafinesquia_neomexicana_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Rafinesquia_neomexicana_results)
View(total_results)

## ======================================================================================

# Malacothrix glabrata Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Malacothrix_glabrata.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Malacothrix glabrata", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Malacothrix_glabrata_results <- results

## ======================================================================================

# Malacothrix glabrata Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Malacothrix_glabrata_results <- cbind(Malacothrix_glabrata_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Malacothrix_glabrata_results)

## ======================================================================================

# Loeseliastrum matthewsii Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Loeseliastrum_matthewsii.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Loeseliastrum matthewsii", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Loeseliastrum_matthewsii_results <- results

## ======================================================================================

# Loeseliastrum matthewsii Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Loeseliastrum_matthewsii_results <- cbind(Loeseliastrum_matthewsii_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Loeseliastrum_matthewsii_results)

## ======================================================================================

# Lupinus shockleyi Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Lupinus_shockleyi.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Lupinus shockleyi", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Lupinus_shockleyi_results <- results

## ======================================================================================

# Lupinus shockleyi Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Lupinus_shockleyi_results <- cbind(Lupinus_shockleyi_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Lupinus_shockleyi_results)

## ======================================================================================

# Eschscholzia glyptosperma Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eschscholzia_glyptosperma.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eschscholzia glyptosperma", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eschscholzia_glyptosperma_results <- results

## ======================================================================================

# Eschscholzia glyptosperma Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eschscholzia_glyptosperma_results <- cbind(Eschscholzia_glyptosperma_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eschscholzia_glyptosperma_results)

## ======================================================================================

# Eriogonum pusillum Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eriogonum_pusillum.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eriogonum pusillum", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eriogonum_pusillum_results <- results

## ======================================================================================

# Eriogonum pusillum Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eriogonum_pusillum_results <- cbind(Eriogonum_pusillum_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eriogonum_pusillum_results)

## ======================================================================================

# Senecio mohavensis Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Senecio_mohavensis.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Senecio mohavensis", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Senecio_mohavensis_results <- results

## ======================================================================================

# Senecio mohavensis Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Senecio_mohavensis_results <- cbind(Senecio_mohavensis_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Senecio_mohavensis_results)

## ======================================================================================

# Anisocoma acaulis Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Anisocoma_acaulis.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Anisocoma acaulis", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Anisocoma_acaulis_results <- results

## ======================================================================================

# Anisocoma acaulis Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Anisocoma_acaulis_results <- cbind(Anisocoma_acaulis_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Anisocoma_acaulis_results)

## ======================================================================================

# Eriophyllum ambiguum Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eriophyllum_ambiguum.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eriophyllum ambiguum", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eriophyllum_ambiguum_results <- results

## ======================================================================================

# Eriophyllum ambiguum Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eriophyllum_ambiguum_results <- cbind(Eriophyllum_ambiguum_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eriophyllum_ambiguum_results)

## ======================================================================================

# Eucrypta micrantha Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eucrypta_micrantha.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eucrypta micrantha", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eucrypta_micrantha_results <- results

## ======================================================================================

# Eucrypta micrantha Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eucrypta_micrantha_results <- cbind(Eucrypta_micrantha_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eucrypta_micrantha_results)

## ======================================================================================

# Gilia stellata Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Gilia_stellata.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Gilia stellata", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Gilia_stellata_results <- results

## ======================================================================================

# Gilia stellata Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Gilia_stellata_results <- cbind(Gilia_stellata_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Gilia_stellata_results)

## ======================================================================================

# Linanthus demissus Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Linanthus_demissus.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Linanthus demissus", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Linanthus_demissus_results <- results

## ======================================================================================

# Linanthus demissus Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Linanthus_demissus_results <- cbind(Linanthus_demissus_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Linanthus_demissus_results)

## ======================================================================================

# Boechera pulchra Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Boechera_pulchra.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Boechera pulchra", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Boechera_pulchra_results <- results

## ======================================================================================

# Boechera pulchra Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Boechera_pulchra_results <- cbind(Boechera_pulchra_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Boechera_pulchra_results)

## ======================================================================================

# Eriophyllum pringlei Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Eriophyllum_pringlei.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Eriophyllum_pringlei", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Eriophyllum_pringlei_results <- results

## ======================================================================================

# Eriophyllum pringlei Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Eriophyllum_pringlei_results <- cbind(Eriophyllum_pringlei_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Eriophyllum_pringlei_results)

## ======================================================================================

# Mohavea breviflora Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Mohavea_breviflora.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Mohavea breviflora", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Mohavea_breviflora_results <- results

## ======================================================================================

# Mohavea breviflora Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Mohavea_breviflora_results <- cbind(Mohavea_breviflora_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Mohavea_breviflora_results)

## ======================================================================================

# Delphinium parishii Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Delphinium_parishii.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Delphinium parishii", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Delphinium_parishii_results <- results

## ======================================================================================

# Delphinium parishii Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Delphinium_parishii_results <- cbind(Delphinium_parishii_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Delphinium_parishii_results)

## ======================================================================================

# Caulanthus cooperi Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Caulanthus_cooperi.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Caulanthus cooperi", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Caulanthus_cooperi_results <- results

## ======================================================================================

# Caulanthus cooperi Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Caulanthus_cooperi_results <- cbind(Caulanthus_cooperi_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Caulanthus_cooperi_results)

## ======================================================================================

# Chaenactis carphoclinia Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Chaenactis_carphoclinia.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Chaenactis carphoclinia", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Chaenactis_carphoclinia_results <- results

## ======================================================================================

# Chaenactis carphoclinia Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Chaenactis_carphoclinia_results <- cbind(Chaenactis_carphoclinia_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Chaenactis_carphoclinia_results)

## ======================================================================================

# Antirrhinum filipes Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Antirrhinum_filipes.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Antirrhinum filipes", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Antirrhinum_filipes_results <- results

## ======================================================================================

# Antirrhinum filipes Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Antirrhinum_filipes_results <- cbind(Antirrhinum_filipes_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Antirrhinum_filipes_results)

## ======================================================================================

# Hesperocallis undulata Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Hesperocallis_undulata.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Hesperocallis undulata", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Hesperocallis_undulata_results <- results

## ======================================================================================

# Hesperocallis undulata Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Hesperocallis_undulata_results <- cbind(Hesperocallis_undulata_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Hesperocallis_undulata_results)

## ======================================================================================

# Phacelia pedicellata Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Phacelia_pedicellata.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Phacelia pedicellata", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Phacelia_pedicellata_results <- results

## ======================================================================================

# Phacelia pedicellata Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Phacelia_pedicellata_results <- cbind(Phacelia_pedicellata_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Phacelia_pedicellata_results)

## ======================================================================================

# Phacelia rotundifolia Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Phacelia_rotundifolia.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Phacelia rotundifolia", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Phacelia_rotundifolia_results <- results

## ======================================================================================

# Phacelia rotundifolia Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Phacelia_rotundifolia_results <- cbind(Phacelia_rotundifolia_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Phacelia_rotundifolia_results)

## ======================================================================================

# Prenanthella exigua Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Prenanthella_exigua.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Prenanthella exigua", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Prenanthella_exigua_results <- results

## ======================================================================================

# Prenanthella exigua Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Prenanthella_exigua_results <- cbind(Prenanthella_exigua_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Prenanthella_exigua_results)

## ======================================================================================

# Syntrichopappus fremontii Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Syntrichopappus_fremontii.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Syntrichopappus fremontii", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Syntrichopappus_fremontii_results <- results

## ======================================================================================

# Syntrichopappus fremontii Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Syntrichopappus_fremontii_results <- cbind(Syntrichopappus_fremontii_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Syntrichopappus_fremontii_results)

## ======================================================================================

# Astragalus layneae Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Astragalus_layneae.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Astragalus layneae", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Astragalus_layneae_results <- results

## ======================================================================================

# Astragalus layneae Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Astragalus_layneae_results <- cbind(Astragalus_layneae_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Astragalus_layneae_results)

## ======================================================================================

# Calochortus flexuosus Flowering Analysis ----

# Reset the working directory and load the cleaned data with flowering information.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")
vf_data <- read.csv("Calochortus_flexuosus.csv")

# Collect all sampling dates for preserved specimens and human observations.
human_observations <- filter(vf_data, basisOfRecord == "HumanObservation" |
                               basisOfRecord == "Human Observation")
human_observations_total <- c(human_observations$DOY, human_observations_total)
preserved_specimens <- filter(vf_data, basisOfRecord == "PreservedSpecimen" |
                                basisOfRecord == "Preserved Specimen")
preserved_specimens_total <- c(preserved_specimens$DOY, preserved_specimens_total)

# Load dplyr and filter the data set to show only entries that were flowering.
library("dplyr")
vf_data <- filter(vf_data, VF == "YES")

# Create a summary table with coordinate information for the samples.
temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                          datum = geodeticDatum, gbifID = gbifID)

# Filter the data according to geodetic datum (samples with unknown or unspecified 
# coordinate systems will be assumed WGS84).
NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                    datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")

# Convert into sf objects and set the coordinate reference system accordingly.
NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)

# Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)

# Bind the fixed coordinates back together into a single data frame.
hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
hoopla <- arrange(hoopla, gbifID)

# Extract the re-projected coordinates as separate longitude and latitude information.
fixed.coords <- as.data.frame(st_coordinates(hoopla))
column_names <- c("lon", "lat")
colnames(fixed.coords) <- (column_names)

# Bind the fixed coordinates back to the original data frame.
vf_data <- arrange(vf_data, gbifID)
vf_data <- cbind(vf_data, fixed.coords)

# Filter the data into five different time periods for analysis.
p_1940_to_1959 <- filter(vf_data, Decade == 1940 | Decade == 1950)
p_1960_to_1979 <- filter(vf_data, Decade == 1960 | Decade == 1970)
p_1980_to_1999 <- filter(vf_data, Decade == 1980 | Decade == 1990)
p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)

# Load the "Phest" phenology estimator package #
library(phest)

# Create vectors to represent any negative flowering times (from December or earlier).
december_1 <- filter(p_1940_to_1959, DOY > 304)
negative_1 <- c(december_1$DOY)
negative_1 <- negative_1 - 365

december_2 <- filter(p_1960_to_1979, DOY > 304)
negative_2 <- c(december_2$DOY)
negative_2 <- negative_2 - 365

december_3 <- filter(p_1980_to_1999, DOY > 304)
negative_3 <- c(december_3$DOY)
negative_3 <- negative_3 - 365

december_4 <- filter(p_2000_to_2019, DOY > 304)
negative_4 <- c(december_4$DOY)
negative_4 <- negative_4 - 365

# Create a vector of DOY (flowering) for each time period.
DOY_vec_1 <- c(negative_1, p_1940_to_1959$DOY)
DOY_vec_2 <- c(negative_2, p_1960_to_1979$DOY)
DOY_vec_3 <- c(negative_3, p_1980_to_1999$DOY)
DOY_vec_4 <- c(negative_4, p_2000_to_2019$DOY)

# Sort the DOY vectors in ascending order.
DOY_vec_1 <- sort(DOY_vec_1, decreasing = FALSE)
DOY_vec_2 <- sort(DOY_vec_2, decreasing = FALSE)
DOY_vec_3 <- sort(DOY_vec_3, decreasing = FALSE)
DOY_vec_4 <- sort(DOY_vec_4, decreasing = FALSE)

# Create a vector that shows the number of samples used to calculate each estimate.
n_samples <- c()

if(length(DOY_vec_1) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_1))
}

if(length(DOY_vec_2) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_2))
}

if(length(DOY_vec_3) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_3))
}

if(length(DOY_vec_4) > 50){
  n_samples <- c(n_samples, 50)
}else{
  n_samples <- c(n_samples, length(DOY_vec_4))
}

# Create a data frame to store the results.
results <- data.frame()

# Estimates for the period of 1940 to 1959.
res_1 <- c(weib.limit(DOY_vec_1, k = 50, upper = FALSE, alpha = 0.05))
res_1 <- c(res_1, weib.limit.bootstrap(DOY_vec_1, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_1)

# Estimates for the period of 1960 to 1979.
res_2 <- c(weib.limit(DOY_vec_2, k = 50, upper = FALSE, alpha = 0.05))
res_2 <- c(res_2, weib.limit.bootstrap(DOY_vec_2, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_2)

# Estimates for the period of 1980 to 1999.
res_3 <- c(weib.limit(DOY_vec_3, k = 50, upper = FALSE, alpha = 0.05))
res_3 <- c(res_3, weib.limit.bootstrap(DOY_vec_3, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_3)

# Estimates for the period of 2000 to 2019.
res_4 <- c(weib.limit(DOY_vec_4, k = 50, upper = FALSE, alpha = 0.05))
res_4 <- c(res_4, weib.limit.bootstrap(DOY_vec_4, k = 50, n = 1000, max.iter = 10, upper = FALSE))
results <- rbind(results, res_4)

# Clean up the "results" data frame and create names for the columns and rows.
column_names <- c("Species", "Period", "N_Samples", "Estimated_DOY", "Lower_CI",
                  "Upper_CI", "Standard_Error")
row_names <- c("p_1940_to_1959", "p_1960_to_1979", "p_1980_to_1999", "p_2000_to_2019")
species_name <- c(rep("Calochortus flexuosus", times = 4))

results <- cbind(n_samples, results)
results <- cbind(row_names, results)
results <- results[,-c(6)]
results <- cbind(species_name, results)
colnames(results) <- (column_names)
Calochortus_flexuosus_results <- results

## ======================================================================================

# Calochortus flexuosus Historic Climate Analysis ----

# Create a summary table with coordinate information for the samples.
samples_1 <- summarise(p_1940_to_1959, lon = lon, lat = lat)
samples_2 <- summarise(p_1960_to_1979, lon = lon, lat = lat)
samples_3 <- summarise(p_1980_to_1999, lon = lon, lat = lat)
samples_4 <- summarise(p_2000_to_2019, lon = lon, lat = lat)

# Extract precipitation data from the raster brick as a new data frame for each time
# period.
pre.sites_1 <- data.frame(raster::extract(pre, samples_1, ncol = 2, na.rm = TRUE))
pre.sites_2 <- data.frame(raster::extract(pre, samples_2, ncol = 2, na.rm = TRUE))
pre.sites_3 <- data.frame(raster::extract(pre, samples_3, ncol = 2, na.rm = TRUE))
pre.sites_4 <- data.frame(raster::extract(pre, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(pre.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(pre.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
pre.sites_1 <- pre.sites_1[,469:708]
pre.sites_2 <- pre.sites_2[,709:948]
pre.sites_3 <- pre.sites_3[,949:1188]
pre.sites_4 <- pre.sites_4[,1189:1428]

# Create vectors with years represented in each time period.
years_1 <- 1940:1959
years_2 <- 1960:1979
years_3 <- 1980:1999
years_4 <- 2000:2019

# Calculate total.annual.precipitation for each site in each time period.
precip.year.total_1 <- as.data.frame(sapply(years_1, function(x) rowSums(pre.sites_1[, grep(x, names(pre.sites_1))])))
precip.year.total_2 <- as.data.frame(sapply(years_2, function(x) rowSums(pre.sites_2[, grep(x, names(pre.sites_2))])))
precip.year.total_3 <- as.data.frame(sapply(years_3, function(x) rowSums(pre.sites_3[, grep(x, names(pre.sites_3))])))
precip.year.total_4 <- as.data.frame(sapply(years_4, function(x) rowSums(pre.sites_4[, grep(x, names(pre.sites_4))])))

# Rename the columns in the new "total precipitation per year" data frames.
names(precip.year.total_1) <- years_1
names(precip.year.total_2) <- years_2
names(precip.year.total_3) <- years_3
names(precip.year.total_4) <- years_4

# Merge the resulting precipitation data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, precip.year.total_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, precip.year.total_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, precip.year.total_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, precip.year.total_4)

# Create a character vector containing the name of the year each sample was taken.
years_sampled_1 <- as.character(c(p_1940_to_1959$Year_2))
years_sampled_2 <- as.character(c(p_1960_to_1979$Year_2))
years_sampled_3 <- as.character(c(p_1980_to_1999$Year_2))
years_sampled_4 <- as.character(c(p_2000_to_2019$Year_2))

# Create a numeric vector that contains the number of samples taken.
number_sampled_1 <- c(1:length(p_1940_to_1959$Year_2))
number_sampled_2 <- c(1:length(p_1960_to_1979$Year_2))
number_sampled_3 <- c(1:length(p_1980_to_1999$Year_2))
number_sampled_4 <- c(1:length(p_2000_to_2019$Year_2))

# Use the two vectors to create a third vector that records the total rainfall at the
# location for each sample in the specific year that they were sampled. From that, find
# the average total rainfall per year across all sites.
p_j_d <- c()
for(i in number_sampled_1){
  p_j_d <- c(p_j_d, p_1940_to_1959[i, years_sampled_1[i]])
}
mean_annual_rainfall_1 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_2){
  p_j_d <- c(p_j_d, p_1960_to_1979[i, years_sampled_2[i]])
}
mean_annual_rainfall_2 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_3){
  p_j_d <- c(p_j_d, p_1980_to_1999[i, years_sampled_3[i]])
}
mean_annual_rainfall_3 <- mean(p_j_d, na.rm = TRUE)

p_j_d <- c()
for(i in number_sampled_4){
  p_j_d <- c(p_j_d, p_2000_to_2019[i, years_sampled_4[i]])
}
mean_annual_rainfall_4 <- mean(p_j_d, na.rm = TRUE)

# Collate the means into a single vector.
avg_annual_precip <- c(mean_annual_rainfall_1, mean_annual_rainfall_2, mean_annual_rainfall_3,
                       mean_annual_rainfall_4)

# Add the means to a new data frame.
means_df <- data.frame(avg_annual_precip)

# Extract temperature data from the raster brick as a new data frame for each time
# period.
tmp.sites_1 <- data.frame(raster::extract(tmp, samples_1, ncol = 2, na.rm = TRUE))
tmp.sites_2 <- data.frame(raster::extract(tmp, samples_2, ncol = 2, na.rm = TRUE))
tmp.sites_3 <- data.frame(raster::extract(tmp, samples_3, ncol = 2, na.rm = TRUE))
tmp.sites_4 <- data.frame(raster::extract(tmp, samples_4, ncol = 2, na.rm = TRUE))

# Create vectors to represent the years and months found in the data sets.
years <- 1901:2019
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
            "Sep", "Oct", "Nov", "Dec")

# Change the names of the columns in the data frame.
names(tmp.sites_1) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_2) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_3) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")
names(tmp.sites_4) <- paste(rep(years, each = 12), rep(months, times = 119), sep = "_")

# Select only the relevant columns for each time period.
tmp.sites_1 <- tmp.sites_1[,469:708]
tmp.sites_2 <- tmp.sites_2[,709:948]
tmp.sites_3 <- tmp.sites_3[,949:1188]
tmp.sites_4 <- tmp.sites_4[,1189:1428]

# Calculate mean monthly rainfall for each time period.
tmp.month.mean_1 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_1[, grep(x, names(tmp.sites_1))], na.rm = TRUE)))
tmp.month.mean_2 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_2[, grep(x, names(tmp.sites_2))], na.rm = TRUE)))
tmp.month.mean_3 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_3[, grep(x, names(tmp.sites_3))], na.rm = TRUE)))
tmp.month.mean_4 <- as.data.frame(sapply(months, function(x) rowMeans(tmp.sites_4[, grep(x, names(tmp.sites_4))], na.rm = TRUE)))

# Rename the columns in the new data frames using the "months" vector.
names(tmp.month.mean_1) <- months
names(tmp.month.mean_2) <- months
names(tmp.month.mean_3) <- months
names(tmp.month.mean_4) <- months

# Merge the resulting temperature data with the original data frames.
p_1940_to_1959 <- cbind(p_1940_to_1959, tmp.month.mean_1)
p_1960_to_1979 <- cbind(p_1960_to_1979, tmp.month.mean_2)
p_1980_to_1999 <- cbind(p_1980_to_1999, tmp.month.mean_3)
p_2000_to_2019 <- cbind(p_2000_to_2019, tmp.month.mean_4)

# Average temperature in December during each time period.
Dec_1 <- mean(p_1940_to_1959[,c("Dec")], na.rm = TRUE)
Dec_2 <- mean(p_1960_to_1979[,c("Dec")], na.rm = TRUE)
Dec_3 <- mean(p_1980_to_1999[,c("Dec")], na.rm = TRUE)
Dec_4 <- mean(p_2000_to_2019[,c("Dec")], na.rm = TRUE)
avg_tmp_Dec <- c(Dec_1, Dec_2, Dec_3, Dec_4)

# Average temperature in January during each time period.
Jan_1 <- mean(p_1940_to_1959[,c("Jan")], na.rm = TRUE)
Jan_2 <- mean(p_1960_to_1979[,c("Jan")], na.rm = TRUE)
Jan_3 <- mean(p_1980_to_1999[,c("Jan")], na.rm = TRUE)
Jan_4 <- mean(p_2000_to_2019[,c("Jan")], na.rm = TRUE)
avg_tmp_Jan <- c(Jan_1, Jan_2, Jan_3, Jan_4)

# Average temperature in February during each time period.
Feb_1 <- mean(p_1940_to_1959[,c("Feb")], na.rm = TRUE)
Feb_2 <- mean(p_1960_to_1979[,c("Feb")], na.rm = TRUE)
Feb_3 <- mean(p_1980_to_1999[,c("Feb")], na.rm = TRUE)
Feb_4 <- mean(p_2000_to_2019[,c("Feb")], na.rm = TRUE)
avg_tmp_Feb <- c(Feb_1, Feb_2, Feb_3, Feb_4)

# Average temperature in March during each time period.
Mar_1 <- mean(p_1940_to_1959[,c("Mar")], na.rm = TRUE)
Mar_2 <- mean(p_1960_to_1979[,c("Mar")], na.rm = TRUE)
Mar_3 <- mean(p_1980_to_1999[,c("Mar")], na.rm = TRUE)
Mar_4 <- mean(p_2000_to_2019[,c("Mar")], na.rm = TRUE)
avg_tmp_Mar <- c(Mar_1, Mar_2, Mar_3, Mar_4)

# Average temperature in April during each time period.
Apr_1 <- mean(p_1940_to_1959[,c("Apr")], na.rm = TRUE)
Apr_2 <- mean(p_1960_to_1979[,c("Apr")], na.rm = TRUE)
Apr_3 <- mean(p_1980_to_1999[,c("Apr")], na.rm = TRUE)
Apr_4 <- mean(p_2000_to_2019[,c("Apr")], na.rm = TRUE)
avg_tmp_Apr <- c(Apr_1, Apr_2, Apr_3, Apr_4)

# Average temperature in May during each time period.
May_1 <- mean(p_1940_to_1959[,c("May")], na.rm = TRUE)
May_2 <- mean(p_1960_to_1979[,c("May")], na.rm = TRUE)
May_3 <- mean(p_1980_to_1999[,c("May")], na.rm = TRUE)
May_4 <- mean(p_2000_to_2019[,c("May")], na.rm = TRUE)
avg_tmp_May <- c(May_1, May_2, May_3, May_4)

# Add the monthly temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_Mar, avg_tmp_Dec, avg_tmp_Jan, avg_tmp_Feb,
                  avg_tmp_Apr, avg_tmp_May)

# Calculate the average temperature in "winter" months (December - February).
avg_tmp_winter <- c()
for(i in 1:length(means_df)){
  jenny <- (means_df[i, 2] + means_df[i, 3] + means_df[i, 4])/3
  avg_tmp_winter <- c(avg_tmp_winter, jenny)
}
avg_tmp_winter <- avg_tmp_winter[1:4]

# Calculate the average temperature in "spring" months (Mar - May).
avg_tmp_spring <- c()
for(i in 1:length(means_df)){
  carly <- (means_df[i, 5] + means_df[i, 6] + means_df[i, 7])/3
  avg_tmp_spring <- c(avg_tmp_spring, carly)
}
avg_tmp_spring <- avg_tmp_spring[1:4]

# Add the seasonal temperature means to the "means_df" data frame.
means_df <- cbind(means_df, avg_tmp_winter, avg_tmp_spring)

# Combine the climate results with the flowering estimations.
Calochortus_flexuosus_results <- cbind(Calochortus_flexuosus_results, means_df)

# Combine the results for this species with the total results so far.
total_results <- rbind(total_results, Calochortus_flexuosus_results)

## ======================================================================================

# Save the Results as a .csv File ----
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_results")
write.csv(total_results, "final_42_species_results.csv")

## ======================================================================================

# Future Climate Analysis ----

# Load in the results file.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_results")
total_results <- read.csv("final_42_species_results.csv")
species_names <- sort(unique(total_results$Species))

# Set the working directory.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_files")

# Load in all flowering observation files as a list.
filelist <- list.files(pattern = '.csv', all.files = TRUE, full.names = FALSE)
list_of_csv <- lapply(filelist, read.csv)

# Make a blank data frame.
hopper <- data.frame()

# Create a loop.
for(i in list_of_csv){
  
  # Filter the data set to show only entries that were flowering in the 2000's/2010's.
  vf_data <- filter(i, VF == "YES")
  p_2000_to_2019 <- filter(vf_data, Decade == 2000 | Decade == 2010)
  
  # Create a summary table with coordinate information for the samples.
  temp_summary <- summarise(vf_data, lon = decimalLongitude, lat = decimalLatitude,
                            datum = geodeticDatum, gbifID = gbifID)
  
  # Filter the data according to geodetic datum (samples with unknown or unspecified 
  # coordinate systems will be assumed WGS84).
  NAD27_guys <- temp_summary%>%filter(datum == "NAD27" | datum == "NAD 27" | datum == "NAD 1927")
  NAD83_guys <- temp_summary%>%filter(datum == "NAD83" | datum == "NAD 83" | datum == "NAD 1983")
  WGS84_guys <- temp_summary%>%filter(datum != "NAD27" & datum != "NAD 27" & datum != "NAD 1927",
                                      datum != "NAD83" & datum != "NAD 83" & datum != "NAD 1983")
  
  # Convert into sf objects and set the coordinate reference system accordingly.
  NAD27_coords <- st_as_sf(NAD27_guys, coords = c("lon", "lat"), crs = 4267)
  NAD83_coords <- st_as_sf(NAD83_guys, coords = c("lon", "lat"), crs = 4269)
  WGS84_coords <- st_as_sf(WGS84_guys, coords = c("lon", "lat"), crs = 4326)
  
  # Re-project NAD27 and NAD83 coordinates into the WGS84 reference system.
  fixed_NAD27 <- st_transform(NAD27_coords, crs = 4326)
  fixed_NAD83 <- st_transform(NAD83_coords, crs = 4326)
  
  # Bind the fixed coordinates back together into a single data frame.
  hoopla <- rbind(WGS84_coords, fixed_NAD27, fixed_NAD83)
  hoopla <- arrange(hoopla, gbifID)
  
  # Extract the re-projected coordinates as separate longitude and latitude information.
  fixed.coords <- as.data.frame(st_coordinates(hoopla))
  column_names <- c("lon", "lat")
  colnames(fixed.coords) <- (column_names)
  
  # Create a summary table with coordinate information for the samples.
  samples_4 <- drop_na(fixed.coords)
  
  # Extract projected average annual precipitation for all three models under RCP 4.5.
  HE45pre <- data.frame(rowSums(raster::extract(he45pr50, samples_4, ncol = 2, na.rm = TRUE)))
  names(HE45pre) <- c("rainfall")
  HE45pre <- drop_na(HE45pre)
  mean_HE45pre <- mean(HE45pre$rainfall)
  
  CN45pre <- data.frame(rowSums(raster::extract(cn45pr50, samples_4, ncol = 2, na.rm = TRUE)))
  names(CN45pre) <- c("rainfall")
  CN45pre <- drop_na(CN45pre)
  mean_CN45pre <- mean(CN45pre$rainfall)
  
  MC45pre <- data.frame(rowSums(raster::extract(mc45pr50, samples_4, ncol = 2, na.rm = TRUE)))
  names(MC45pre) <- c("rainfall")
  MC45pre <- drop_na(MC45pre)
  mean_MC45pre <- mean(MC45pre$rainfall)
  
  mean_pre_45 <- mean(mean_HE45pre, mean_CN45pre, mean_MC45pre)
  
  # Extract projected average annual precipitation for all three models under RCP 8.5.
  HE85pre <- data.frame(rowSums(raster::extract(he85pr50, samples_4, ncol = 2, na.rm = TRUE)))
  names(HE85pre) <- c("rainfall")
  HE85pre <- drop_na(HE85pre)
  mean_HE85pre <- mean(HE85pre$rainfall)
  
  CN85pre <- data.frame(rowSums(raster::extract(cn85pr50, samples_4, ncol = 2, na.rm = TRUE)))
  names(CN85pre) <- c("rainfall")
  CN85pre <- drop_na(CN85pre)
  mean_CN85pre <- mean(CN85pre$rainfall)
  
  MC85pre <- data.frame(rowSums(raster::extract(mc85pr50, samples_4, ncol = 2, na.rm = TRUE)))
  names(MC85pre) <- c("rainfall")
  MC85pre <- drop_na(MC85pre)
  mean_MC85pre <- mean(MC85pre$rainfall)
  
  mean_pre_85 <- mean(mean_HE85pre, mean_CN85pre, mean_MC85pre)
  
  # Extract projected average winter temperature for all three models under RCP 4.5.
  HE45tmp <- data.frame(raster::extract(he45tmp50, samples_4, ncol = 2, na.rm = TRUE))
  names(HE45tmp) <- c("temperature")
  HE45tmp <- drop_na(HE45tmp)
  mean_HE45tmp <- mean(HE45tmp$temperature)
  
  CN45tmp <- data.frame(raster::extract(cn45tmp50, samples_4, ncol = 2, na.rm = TRUE))
  names(CN45tmp) <- c("temperature")
  CN45tmp <- drop_na(CN45tmp)
  mean_CN45tmp <- mean(CN45tmp$temperature)
  
  MC45tmp <- data.frame(raster::extract(mc45tmp50, samples_4, ncol = 2, na.rm = TRUE))
  names(MC45tmp) <- c("temperature")
  MC45tmp <- drop_na(MC45tmp)
  mean_MC45tmp <- mean(MC45tmp$temperature)
  
  mean_tmp_45 <- mean(mean_HE45tmp, mean_CN45tmp, mean_MC45tmp)/10
  
  # Extract projected average winter temperature for all three models under RCP 8.5.
  HE85tmp <- data.frame(raster::extract(he85tmp50, samples_4, ncol = 2, na.rm = TRUE))
  names(HE85tmp) <- c("temperature")
  HE85tmp <- drop_na(HE85tmp)
  mean_HE85tmp <- mean(HE85tmp$temperature)
  
  CN85tmp <- data.frame(raster::extract(cn85tmp50, samples_4, ncol = 2, na.rm = TRUE))
  names(CN85tmp) <- c("temperature")
  CN85tmp <- drop_na(CN85tmp)
  mean_CN85tmp <- mean(CN85tmp$temperature)
  
  MC85tmp <- data.frame(raster::extract(mc85tmp50, samples_4, ncol = 2, na.rm = TRUE))
  names(MC85tmp) <- c("temperature")
  MC85tmp <- drop_na(MC85tmp)
  mean_MC85tmp <- mean(MC85tmp$temperature)
  
  mean_tmp_85 <- mean(mean_HE85tmp, mean_CN85tmp, mean_MC85tmp)/10
  
  # Create a dataframe to store projectied infromation on future conditions.
  future_conditions <- data.frame(cbind(mean_pre_45, mean_pre_85, mean_tmp_45,mean_tmp_85))
  hopper <- rbind(hopper, future_conditions)
  
}

# Add time period and species names to the dataframe.
p_future <- c(rep("p_2041_to_2060", times = length(hopper$mean_pre_45)))
plant_future <- hopper[,c(1,3)]
plant_future <- cbind(species_names, p_future, plant_future)

# Isolate the environmental conditions from the earlier estimates
results_weather_only <- total_results[,c(2,3,9,16)]

# Combine the two data frames.
plant_future <- setNames(plant_future, names(results_weather_only))
past_and_future_weather <- rbind(results_weather_only, plant_future)
past_and_future_weather <- arrange(past_and_future_weather, Species)

# Write the results to a csv file.
setwd("~/Desktop")
write.csv(past_and_future_weather, "past_and_future_weather.csv")

## ======================================================================================

# Making Maps (Figures 1 and 2) ----

# Load necessary packages.
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(stringr)
library(sf)
library(RColorBrewer)
library(ggspatial)
library(cowplot)
library(magick)
library(dplyr)
library(tidyr)
library(cruts)
library(rgdal)
library(raster)
library(sp)
library(units)
library(ncdf4)
library(ggpubr)
library(gridExtra)

# Set working directory.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/deserts_sw")

# Load North American desert shapefiles and re-project coordinates to WGS84.
deserts_sw <- st_read("deserts_sw.shp")
deserts_sw <- st_transform(deserts_sw, crs = 4326)

# Get map data for the entire United States, California, and Mexico.
usa <- map_data("usa")
states <- map_data("state")
mexico <- map_data("world", region = "mexico")
california <- subset(states, region %in% c("california"))
world <- map_data("world")

# Get coordinates for major cities in the American Southwest.
major_us_cities <- filter(world.cities, country.etc == "USA")
southwest_cities <- filter(major_us_cities, name == "Las Vegas" |
                             name == "Phoenix" |
                             name == "Los Angeles")
ca_capital <- filter(major_us_cities, name == "Sacramento")

# Set formatting guidelines to get rid of the axes labels, titles, and ticks.
ditch_the_axes <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
)

# Set colors.
desert_colors <- brewer.pal(n = 6, name = "Set3")
desert_colors <- desert_colors[3:6]
some_tans <- brewer.pal(n = 9, name = "OrRd")
light_tan <- some_tans[1]

# Set names for the deserts and bind them to the polygons.
desert_names <- c("Arizona Sonoran", "Colorado Sonoran", "Great Basin", "Mojave")

# Make a map of the United States highlighting California.
small_us <- ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = "grey70", color = "white", size = 0.1) +
  geom_polygon(data = california, aes(x = long, y = lat, group = group),
               fill = "hotpink", color = "white", size = 0.1) +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = 0.2) +
  theme_void() +
  ditch_the_axes +
  theme(aspect.ratio = (5.8/10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"))
small_us

# Plot a map of North American deserts focused on Southern California.
total_map <- ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = light_tan, color = "grey60") + 
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group),
               fill = light_tan, color = "grey60") +
  geom_polygon(data = california, aes(x = long, y = lat, group = group),
               fill = "antiquewhite", color = NA) +
  geom_sf(data = deserts_sw, color = NA, aes(fill = NAME)) +
  scale_fill_manual(values = desert_colors,
                    labels = desert_names) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey60", size = .6) + 
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = NA,
               color = "gray20", size = .75) +
  coord_sf(xlim = c(-124, -111.9), ylim = c(32.2, 39.5), expand = FALSE) +
  geom_point(data = southwest_cities, aes(x = long, y = lat), size = 1.8) +
  geom_text(label = "★", size = 3.2, family = "HiraKakuPro-W3", 
            aes(x = -121.47, y = 38.57)) +
  annotation_scale(location = "bl", width_hint = 0.3, 
                   bar_cols = c("black", "aliceblue")) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.2, "in"),
                         pad_y = unit(0.4, "in"),
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  annotate(geom = "text", x = -123, y = 35.5, label = "Pacific\nOcean",
           fontface = "italic", color = "grey50", size = 4) +
  annotate(geom = "text", x = -119.47, y = 33.8, label = "Los Angeles",
           fontface = "bold", color = "grey10", size = 3) +
  annotate(geom = "text", x = -115.22, y = 36.5, label = "Las Vegas",
           fontface = "bold", color = "grey10", size = 3) +
  annotate(geom = "text", x = -112.85, y = 33.54, label = "Phoenix",
           fontface = "bold", color = "grey10", size = 3) +
  annotate(geom = "text", x = -120.8, y = 38.22, label = "Sacramento",
           fontface = "bold", color = "grey10", size = 3) +
  theme_bw() +
  ggtitle("Deserts of Southern California") +
  ditch_the_axes +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "aliceblue"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 10.5),
        plot.title = element_blank())

# Make the map of the United States an inset on the main deserts map.
desert_with_inset <- ggdraw() +
  draw_plot(total_map) +
  draw_plot(small_us, x = 0.675, y = 0.77, width = 0.21, height = 0.21)
desert_with_inset

# Set the working directory #
setwd("/Volumes/Personal_Files/new_climate_projections")

# Load in the projected future precipitation and temperature data sets.
fut_pre_4.5 <- brick("prdiff_30yavg_ens10_rcp45_2035-2064.LOCA_2016-04-02.16th.CA_NV.tif")
fut_pre_8.5 <- brick("prdiff_30yavg_ens10_rcp85_2035-2064.LOCA_2016-04-02.16th.CA_NV.tif")
fut_tmx_4.5 <- brick("tasmaxdiff_30yavg_ens10_rcp45_2035-2064.LOCA_2016-04-02.16th.CA_NV.tif")
fut_tmx_8.5 <- brick("tasmaxdiff_30yavg_ens10_rcp85_2035-2064.LOCA_2016-04-02.16th.CA_NV.tif")

# Convert the rasters into data frames for analysis.
fut_pre_4.5_df <- as.data.frame(fut_pre_4.5, xy = TRUE)
fut_pre_8.5_df <- as.data.frame(fut_pre_8.5, xy = TRUE)

# Try to break the continuous scale into a discrete one for visualization.
fut_pre_4.5_df <- fut_pre_4.5_df %>%
  mutate(fct_pre_change = cut(prdiff_30yavg_ens10_rcp45_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                              breaks = 5))
togogo <- ggplot() +
  geom_bar(data = fut_pre_4.5_df, aes(fct_pre_change))
unique(fut_pre_4.5_df$fct_pre_change)

fut_pre_8.5_df <- fut_pre_8.5_df %>%
  mutate(fct_pre_change = cut(prdiff_30yavg_ens10_rcp85_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                              breaks = 5))
togogo <- ggplot() +
  geom_bar(data = fut_pre_8.5_df, aes(fct_pre_change))
unique(fut_pre_8.5_df$fct_pre_change)

# Based on that trial, change the bins to an appropriate scale.
custom_bins <- c(-1, -0.01, 0.01, 1)
fut_pre_4.5_df <- fut_pre_4.5_df %>%
  mutate(fct_pre_change_2 = cut(prdiff_30yavg_ens10_rcp45_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                                breaks = custom_bins))
fut_pre_8.5_df <- fut_pre_8.5_df %>%
  mutate(fct_pre_change_2 = cut(prdiff_30yavg_ens10_rcp85_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                                breaks = custom_bins))

# Set the correct labels and fill colors for the raster data.
rain_names <- c("Decrease (mm/day)", "Near Zero Change", "Increase (mm/day)")
rain_colors <- brewer.pal(n = 4, name = "Set1")
rain_colors <- rain_colors[c(4, 2, 3)]

# Create a data frame with the coordinates of the outer left plot dimensions.
mask_box <- data.frame(x = c(-124, -124),
                       y = c(31.75801, 39.5))

# Create a data frame with the coastal boarder of California.
trying <- data.frame(x = california$long, y = california$lat)
trying_2 <- trying[73:449,]
trying_3 <- trying_2[rev(rownames(trying_2)),]
mask_box_2 <- rbind(mask_box, trying_3)

# Create a data frame with the coastal border of Mexico.
meji <- data.frame(x = mexico$long, y = mexico$lat)
meji_2 <- meji[866:872,]
meji_3 <- meji_2[rev(rownames(meji_2)),]
mask_box_3 <- rbind(mask_box_2, meji_3)
ender <- data.frame(x = c(-124), y = c(31.75801))

# Bind all the data frames together to create a polygon mask.
mask <- rbind(mask_box_3, ender)

# Plot two maps for projected precipitation change in Southern California.
med_rain_map <- ggplot() +
  geom_raster(data = fut_pre_4.5_df , aes(x = x, y = y, fill = fct_pre_change_2)) + 
  scale_fill_manual(na.value = NA,
                    values = rain_colors,
                    labels = rain_names,
                    na.translate = FALSE) +
  labs(fill = "Projected Average\nPrecipitation Rate (mm/day") +
  geom_polygon(data = mask,
               aes(x = x, y = y), color = "grey40", size = .8,
               fill = "white") +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) + 
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) +
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey10", size = 1.4) +
  coord_sf(xlim = c(-123.5, -113.375), ylim = c(31.8, 39.5), expand = FALSE) +
  theme_bw() +
  ggtitle("Medium Emissions Scenario (RCP 4.5)") +
  annotate(geom = "text", x = -122.7, y = 32.4, label = "A",
           fontface = "bold", color = "black", size = 6) +
  ditch_the_axes +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  size = 14),
        legend.title =  element_text(face = "bold", size = 14),
        legend.position = "none",
        plot.margin = unit(c(.5, 0, 0, .5), "cm"),
        aspect.ratio = (9/10))
med_rain_map

high_rain_map <- ggplot() +
  geom_raster(data = fut_pre_8.5_df , aes(x = x, y = y, fill = fct_pre_change_2)) + 
  scale_fill_manual(na.value = NA,
                    values = rain_colors,
                    labels = rain_names,
                    na.translate = FALSE) +
  labs(fill = "Projected Average\nPrecipitation Rate (mm/day") +
  geom_polygon(data = mask,
               aes(x = x, y = y), color = "grey40", size = .8,
               fill = "white") +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) + 
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) +
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey10", size = 1.4) +
  coord_sf(xlim = c(-123.5, -113.375), ylim = c(31.8, 39.5), expand = FALSE) +
  theme_bw() +
  ggtitle("High Emissions Scenario (RCP 8.5)") +
  annotate(geom = "text", x = -122.7, y = 32.4, label = "B",
           fontface = "bold", color = "black", size = 6) +
  ditch_the_axes +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  size = 14),
        legend.title =  element_text(face = "bold", size = 14),
        legend.position = "none",
        plot.margin = unit(c(.5, .5, 0, 0), "cm"),
        aspect.ratio = (9/10))
high_rain_map

# Plot a map of projected precipitation change in Southern California.
high_rain_legend <- ggplot() +
  geom_raster(data = fut_pre_8.5_df , aes(x = x, y = y, fill = fct_pre_change_2)) + 
  scale_fill_manual(na.value = NA,
                    values = rain_colors,
                    labels = rain_names,
                    na.translate = FALSE) +
  labs(fill = "Change in Average\nPrecipitation Rate") +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = .8) + 
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = .8) +
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey10", size = 1.3) +
  coord_sf(xlim = c(-124, -113.375), ylim = c(31.5625, 39.5), expand = FALSE) +
  theme_bw() +
  ggtitle("High Emissions Scenario (RCP 8.5)") +
  ditch_the_axes +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  size = 14),
        legend.title =  element_text(face = "bold", size = 14),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        aspect.ratio = (9/10))

# Extract the legend.
rain_leg <- get_legend(high_rain_legend)

# Convert to a ggplot object.
rain_legend <- as_ggplot(rain_leg)
rain_legend <- rain_legend + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Convert the rasters into data frames for analysis.
fut_tmx_4.5_df <- as.data.frame(fut_tmx_4.5, xy = TRUE)
fut_tmx_8.5_df <- as.data.frame(fut_tmx_8.5, xy = TRUE)

# Try to break the temperature scale into a discrete one for visualization.
fut_tmx_4.5_df <- fut_tmx_4.5_df %>%
  mutate(fct_tmx_change = cut(tasmaxdiff_30yavg_ens10_rcp45_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                              breaks = 5))
togogo <- ggplot() +
  geom_bar(data = fut_tmx_4.5_df, aes(fct_tmx_change))
unique(fut_tmx_4.5_df$fct_tmx_change)

fut_tmx_8.5_df <- fut_tmx_8.5_df %>%
  mutate(fct_tmx_change = cut(tasmaxdiff_30yavg_ens10_rcp85_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                              breaks = 5))
togogo <- ggplot() +
  geom_bar(data = fut_tmx_8.5_df, aes(fct_tmx_change))
unique(fut_tmx_8.5_df$fct_tmx_change)

# Based on that trial, change the bins to an appropriate scale.
custom_bins <- c(1, 1.75, 2.5, 3.25, 4, 4.75)
fut_tmx_4.5_df <- fut_tmx_4.5_df %>%
  mutate(fct_tmx_change_2 = cut(tasmaxdiff_30yavg_ens10_rcp45_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                                breaks = custom_bins))
fut_tmx_8.5_df <- fut_tmx_8.5_df %>%
  mutate(fct_tmx_change_2 = cut(tasmaxdiff_30yavg_ens10_rcp85_2035.2064.LOCA_2016.04.02.16th.CA_NV,
                                breaks = custom_bins))

# Set the correct labels and fill colors for the raster data.
tmx_names <- c("1.00 - 1.75 °C", "1.75 - 2.50 °C", "2.50 - 3.25 °C", "3.25 - 4.00 °C",
               "> 4.00 °C")
tmx_colors <- brewer.pal(n = 9, name = "YlOrRd")
tmx_colors <- tmx_colors[c(1, 3, 5, 7, 9)]

# Plot two maps for projected temperature change in Southern California.
med_tmx_map <- ggplot() +
  geom_raster(data = fut_tmx_4.5_df , aes(x = x, y = y, fill = fct_tmx_change_2)) + 
  scale_fill_manual(na.value = NA,
                    values = tmx_colors,
                    labels = tmx_names,
                    na.translate = FALSE) +
  labs(fill = "Increase in Average\nMaximum Temperature") +
  geom_polygon(data = mask,
               aes(x = x, y = y), color = "grey40", size = .8,
               fill = "white") +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) + 
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) +
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey10", size = 1.4) +
  coord_sf(xlim = c(-123.5, -113.375), ylim = c(31.8, 39.5), expand = FALSE) +
  theme_bw() +
  ggtitle("Medium Emissions Scenario (RCP 4.5)") +
  annotate(geom = "text", x = -122.7, y = 32.4, label = "C",
           fontface = "bold", color = "black", size = 6) +
  ditch_the_axes +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  size = 14),
        legend.title =  element_text(face = "bold", size = 14),
        legend.position = "none",
        plot.margin = unit(c(0, 0, .5, .5), "cm"),
        aspect.ratio = (9/10))
med_tmx_map

high_tmx_map <- ggplot() +
  geom_raster(data = fut_tmx_8.5_df , aes(x = x, y = y, fill = fct_tmx_change_2)) + 
  scale_fill_manual(na.value = NA,
                    values = tmx_colors,
                    labels = tmx_names,
                    na.translate = FALSE) +
  labs(fill = "Increase in Average\nMaximum Temperature") +
  geom_polygon(data = mask,
               aes(x = x, y = y), color = "grey40", size = .8,
               fill = "white") +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) + 
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = 1.2) +
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey10", size = 1.4) +
  coord_sf(xlim = c(-123.5, -113.375), ylim = c(31.8, 39.5), expand = FALSE) +
  theme_bw() +
  ggtitle("High Emissions Scenario (RCP 8.5)") +
  annotate(geom = "text", x = -122.7, y = 32.4, label = "D",
           fontface = "bold", color = "black", size = 6) +
  ditch_the_axes +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  size = 14),
        legend.title =  element_text(face = "bold", size = 14),
        legend.position = "none",
        plot.margin = unit(c(0, .5, .5, 0), "cm"),
        aspect.ratio = (9/10))
high_tmx_map

high_tmx_legend <- ggplot() +
  geom_raster(data = fut_tmx_8.5_df , aes(x = x, y = y, fill = fct_tmx_change_2)) + 
  scale_fill_manual(na.value = NA,
                    values = tmx_colors,
                    labels = tmx_names,
                    na.translate = FALSE) +
  labs(fill = "Increase in Average\nMaximum Temperature") +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = .8) + 
  geom_polygon(data = mexico, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40", size = .8) +
  geom_polygon(data = california, aes(x = long, y = lat, group = group), fill = NA,
               color = "grey10", size = 1.3) +
  coord_sf(xlim = c(-123.5, -113.375), ylim = c(31.5625, 39.5), expand = FALSE) +
  theme_bw() +
  ggtitle("High Emissions Scenario (RCP 8.5)") +
  ditch_the_axes +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5,
                                  size = 14),
        legend.title =  element_text(face = "bold", size = 14),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        aspect.ratio = (9/10))

# Extract the legend.
tmx_leg <- get_legend(high_tmx_legend)

# Convert to a ggplot object.
tmx_legend <- as_ggplot(tmx_leg)
tmx_legend <- tmx_legend + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Arrange all four plots in the device window at the same time.
all_maps <- ggarrange(med_rain_map, rain_legend, high_rain_map,
                      med_tmx_map, tmx_legend, high_tmx_map,
                      ncol = 3, nrow = 2, widths = c(1, .7, 1))
all_maps

## ======================================================================================

# Plotting Results by Time Period (Figure 4) ----

# Load in the results file.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_results")
total_results <- read.csv("final_42_species_results.csv")

# Load the plotting package.
library(ggplot2)
library(ggpubr)

# Create a simple linear regression.
time_periods <- c("1940_to_1959", "1960_to_1979",
                  "1980_to_1999", "2000_to_2019")
time_periods <- gsub("_", " ", time_periods)

# Plot the timelines.
time_plot <- ggplot(total_results, aes(x = Period, y = Estimated_DOY, color = Species,
                          group = Species)) +
  geom_point(size = 2) +
  geom_line(linetype = "dashed") +
  stat_summary(aes(y = Estimated_DOY, group = 1), fun = mean,
               colour = "black", geom = "line", size = 1, group = 1) +
  stat_summary(aes(y = Estimated_DOY, group = 1), fun = mean,
               colour = "black", geom = "point", size = 3, group = 1) +
  scale_x_discrete(name = "Time Period", 
                   labels = c("1940 to 1959", "1960 to 1979",
                              "1980 to 1999", "2000 to 2019")) +
  scale_y_continuous(name = "Estimated Day of First Flowering", 
                     breaks = c(-91.5, -60.5, -30.5, 0.5, 31.5, 59.5, 90.5, 120),
                     limits = c(-91.5, 120),
                     labels = c(" ", "OCT   ", "NOV   ", "DEC   ", 
                                "JAN   ", "FEB   ",
                                "MAR   ", "APR   ")) +
  guides(color = guide_legend(ncol = 3)) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.title = element_blank(),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid",
                                        colour = "lightgrey"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = -0.45, 
                                   hjust = 1.2),
        plot.margin = unit(c(.5, .5, .5, .5), "cm"),
        legend.position = "bottom",
        text = element_text(size = 13),
        legend.text = element_text(size = 11))
time_plot

## ======================================================================================

# Plotting Results by Rainfall and Temperature (Figure 5) ----

# Load in the results file.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_results")
total_results <- read.csv("final_42_species_results.csv")

# Load the necessary packages.
library(lme4)
library(MuMIn)
library(ggplot2)
library(ggpubr)

# Define weights for a weighted least squares regression based on CI width.
CI_weights <- 1 / sqrt(abs(total_results$Upper_CI - total_results$Lower_CI))

# Create a multiple weighted linear regression.
wt_lm <- lm(Estimated_DOY ~ avg_annual_precip, data = total_results, weights = CI_weights)
summary(wt_lm)

# Create a linear mixed effects model with Species and Period as random factors.
kablooey <- lmer(Estimated_DOY ~ avg_annual_precip + (1 | Species), data = total_results)
summary(kablooey)

# Make a plot showing the relationship between flowering and rainfall.
rain_plot <- ggplot(total_results, aes(x = avg_annual_precip, 
                                       y = Estimated_DOY, 
                                       color = Species)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = FALSE) +
  geom_abline(intercept = -24.82167, slope = 0.20351, size = 1, linetype = "dashed") +
  geom_abline(intercept = -43.38375, slope = 0.23604, size = 1) +
  scale_x_continuous(name = "Mean Annual Rainfall (mm)",
                     breaks = seq(150, 600, by = 50)) +
  scale_y_continuous(name = "Estimated Day of First Flowering", 
                     breaks = c(-91.5, -60.5, -30.5, 0.5, 31.5, 59.5, 90.5, 120),
                     limits = c(-91.5, 120),
                     labels = c(" ", "OCT  ", "NOV  ", "DEC  ", 
                                "JAN  ", "FEB  ",
                                "MAR  ", "APR  ")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid",
                                        colour = "lightgrey"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = -0.45, 
                                   hjust = 1.2),
        legend.position = "none",
        legend.title = element_blank(),
        text = element_text(size = 13))
rain_plot

# Define weights for a weighted least squares regression based on CI width.
CI_weights <- 1 / sqrt(abs(total_results$Upper_CI - total_results$Lower_CI))

# Create a multiple weighted linear regression.
wt_lm <- lm(Estimated_DOY ~ avg_tmp_winter, data = total_results, weights = CI_weights)
summary(wt_lm)

# Create a linear mixed effects model with Species and Period as random factors.
kablooey <- lmer(Estimated_DOY ~ avg_tmp_winter + (1 | Species), data = total_results)
summary(kablooey)

# Make a plot showing the relationship between flowering and winter temperature.
wint_plot <- ggplot(total_results, aes(x = avg_tmp_winter,
                                       y = Estimated_DOY, 
                                       color = Species)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = FALSE) +
  geom_abline(intercept = 155.44, slope = -13.00, size = 1, linetype = "dashed") +
  geom_abline(intercept = 170.917, slope = -15.686, size = 1) +
  scale_x_continuous(name = "Mean Winter Temperature (Celcius)",
                     breaks = seq(5, 12.5, by = 1)) +
  scale_y_continuous(name = "Estimated Day of First Flowering", 
                     breaks = c(-91.5, -60.5, -30.5, 0.5, 31.5, 59.5, 90.5, 120),
                     limits = c(-91.5, 120),
                     labels = c(" ", "OCT  ", "NOV  ", "DEC  ", 
                                "JAN  ", "FEB  ",
                                "MAR  ", "APR  ")) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid",
                                        colour = "lightgrey"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = -0.45, 
                                   hjust = 1.2),
        legend.position = "none",
        legend.title = element_blank(),
        text = element_text(size = 13))
wint_plot

# Arrange both plots in the device window simultaneously.
both_env_response <- ggarrange(rain_plot, wint_plot, nrow = 2, ncol = 1, labels = c("A", "B")) 
both_env_response

## ======================================================================================

# Making the Phylogeny/Timelines, and Testing for Phylogenetic Signal (Figure 6) ----

# Load necessary packages 
library(BIEN)
library(phytools)
library(phylogram)
library(ape)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(caper)
library(geiger)

# Load in the species names and traits file.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/species_for_vf")
species_traits <- read.csv("species_and_traits_DOI.csv")
species_traits <- filter(species_traits, Scientific.Name != "")

# Load in the results file.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_results")
total_results <- read.csv("final_42_species_results.csv")
total_results <- arrange(total_results, Species)

# Create a character vector with all the names of the study species.
all_species_names <- sort(unique(species_traits$Scientific.Name))
search_species_names <- sub(" ", "_", all_species_names)

# Download a complete phylogeny of all new world plant species.
complete_phylogeny <- BIEN_phylogeny_complete()

# Check to see which species aren't included in the phylogeny and delete them from
# the search vector.
for (i in search_species_names){
  if(sum(unique(complete_phylogeny$tip.label == i)) != 1){
    search_species_names <- search_species_names[search_species_names != i]
    print(paste(i, "was not found"))
  }
}

# Re-add species that are recorded under synonyms.
search_species_names <- c(search_species_names, "Nama_demissa", "Camissonia_brevipes",
                          "Eriophyllum_pringlei", "Amblyopappus_pusillus", 
                          "Chaenactis_stevioides", "Ericameria_ericoides")
search_species_names <- search_species_names[search_species_names != "Chaenactis_carphoclinia"]

# Trim the phylogeny to keep only the species of interest to this project.
trimmed_phylogeny <- keep.tip(phy = complete_phylogeny, tip = search_species_names)

# Fix the species labels to reflect recent taxonomic changes.
old.labels <- c(trimmed_phylogeny$tip.label)
fix.labels <- sub("_", " ", old.labels)
fix.labels <- replace(fix.labels, fix.labels == "Camissonia brevipes", "Chylismia brevipes")
fix.labels <- replace(fix.labels, fix.labels == "Amblyopappus pusillus", "Syntrichopappus fremontii")
fix.labels <- replace(fix.labels, fix.labels == "Chaenactis stevioides", "Chaenactis carphoclinia")
fix.labels <- replace(fix.labels, fix.labels == "Ericameria ericoides", "Monoptilon bellioides")
fix.labels <- replace(fix.labels, fix.labels == "Nama demissa", "Nama demissum")
new.labels <- paste(" ", fix.labels)
trimmed_phylogeny$tip.label <- new.labels

# Plot the trimmed phylogeny.
par(mar=c(0,0,0,0))
plot.phylo(trimmed_phylogeny, edge.width = 1.5, font = 3, underscore = FALSE,
           x.lim = NULL, adj = 0)

# Create a new column in the data frame where estimates without a confidence interval
# are denoted using the integer 1 and those with a confidence interval are 0.
CI_Fail <- as.integer(is.na(total_results$Lower_CI))
CI_Fail <- CI_Fail + 1
CI_Fail_2 <- c()
for(i in CI_Fail){
  if(i == 1){
    CI_Fail_2 <- c(CI_Fail_2, 19)
  }else{
    CI_Fail_2 <- c(CI_Fail_2, 17)
  }
}
CI_Fail <- CI_Fail_2

# Bind the CI_Fail column to the total_results to define point shape in the plot.
total_results <- total_results[2:16]
total_results <- cbind(CI_Fail, total_results)

# Add all the species names to a vector.
not_in_phylo <- c("Anisocoma acaulis", "Atrichoseris platyphylla", "Prenanthella exigua",
                  "Rafinesquia neomexicana")
not_in_phylo <- sort(not_in_phylo, decreasing = TRUE)
total_timeline_labels <- c(not_in_phylo, fix.labels)

# Bind a place marker to the original data frame so labels are sorted on the timeline
# according to their place in the phylogeny.
numberss <- c(1:42)
cookie <- data.frame(cbind(total_timeline_labels, numberss))
cookie <- rbind(cookie, cookie, cookie, cookie)
shake <- arrange(cookie, total_timeline_labels)
place_marker <- c(shake$numberss)
place_marker <- as.numeric(place_marker)
total_results <- arrange(total_results, Species)
total_results <- cbind(place_marker, total_results)

# Create labels for each time period.
time_periods <- c("1940_to_1959", "1960_to_1979",
                  "1980_to_1999", "2000_to_2019")
time_periods <- gsub("_", " ", time_periods)

# Create a plot that shows individual timelines for each species.
timeline_plot <- ggplot(total_results,
                        aes(x = Estimated_DOY,
                            y = place_marker,
                            color = Period,
                        )) +
  geom_point(size = 3.5, shape = CI_Fail) +
  scale_x_continuous(breaks = c(-91.5, -60.5, -30.5, 0.5, 31.5, 59.5, 90.5, 120),
                     limits = c(-91.5, 120),
                     labels = c(" ", "OCT      ", "NOV    ", "DEC    ", 
                                "JAN     ", "FEB    ",
                                "MAR    ", "APR    ")) +
  scale_y_continuous(breaks = c(1:42),
                     limits = c(1, 42)) +
  scale_color_hue(labels = time_periods) +
  ggtitle("Estimated Day of First Flowering") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(size = 0.55, 
                                          linetype = "solid",
                                          colour = "black"),
        panel.grid.major.x = element_line(size = 0.4, 
                                          linetype = "solid",
                                          colour = "lightgrey"),
        axis.text.x = element_text(hjust = 1),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom")
timeline_plot

# Create a reduced dataset with species and their most recent flowering times.
p_2000_to_2019_only <- filter(total_results, Period == "p_2000_to_2019")
p_2000_to_2019_only <- p_2000_to_2019_only[,c(2,5,9,16)]
p_2000_to_2019_only <- arrange(p_2000_to_2019_only, Species)
p_2000_to_2019_only$Species <- sub("_", " ", p_2000_to_2019_only$Species)

# Edit species names so they match the tip labels.
stormo <- c()
for(i in p_2000_to_2019_only$Species){
  stormo <- c(stormo, paste(" ", i))
}
p_2000_to_2019_only$Species <- stormo

# Check for mismatches between the data set and the phylogeny.
name.check(trimmed_phylogeny, p_2000_to_2019_only,
           data.names = p_2000_to_2019_only$Species)

# Remove species from the data set that aren't found in the phylogeny.
matches <- match(p_2000_to_2019_only$Species, trimmed_phylogeny$tip.label, nomatch = 0)
p_2000_to_2019_only <- subset(p_2000_to_2019_only, matches != 0)
length(p_2000_to_2019_only$Species)

# Make species names the row names instead of a column and remove period.
p_2000_to_2019_only <- data.frame(p_2000_to_2019_only[,-1], row.names = p_2000_to_2019_only[,1])
estimated_doys <- p_2000_to_2019_only[,1]
names(estimated_doys) <- rownames(p_2000_to_2019_only)

# Test for phylogenetic signal among recent flowering dates.
phylosig(trimmed_phylogeny, estimated_doys, method = "lambda", test = TRUE)

# Test for phylogenetic signal based on the response between flowering and rainfall.
data <- data.frame(species = names(estimated_doys), resp = estimated_doys,
                   trait = p_2000_to_2019_only$avg_annual_precip)
c.data <- comparative.data(trimmed_phylogeny, data, species)
rain.pgls <- pgls(resp ~ trait, data = c.data, lambda = "ML")
summary(rain.pgls)

# Test for phylogenetic signal based on the response between flowering and temperature
data <- data.frame(species = names(estimated_doys), resp = estimated_doys,
                   trait = p_2000_to_2019_only$avg_tmp_winter)
c.data <- comparative.data(trimmed_phylogeny, data, species)
tmp.pgls <- pgls(resp ~ trait, data = c.data, lambda = "ML")
summary(tmp.pgls)

# Plot the likelihood profiles for the estimates of lambda.
par(mar = c(10, 5, 4, 5))
lambda.prof <- pgls.profile(rain.pgls, 'lambda')
plot(lambda.prof)

par(mar = c(10, 5, 4, 5))
lambda.prof <- pgls.profile(tmp.pgls, 'lambda')
plot(lambda.prof)

## ======================================================================================

# Making Model Results Tables (Table 1) ----

# Load necessary packages.
library(dplyr)
library(tidyr)()
library(lme4)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(sjPlot)

# Load in the results file.
setwd("~/Desktop/Research Project/Data and Scripts/new_mojave_data/csv_results")
total_results <- read.csv("final_42_species_results.csv")
total_results <- total_results[2:length(total_results)]
total_results <- arrange(total_results, Species)
past_and_future_weather <- read.csv("past_and_future_weather.csv")
fut_weather <- filter(past_and_future_weather, Period == "p_2041_to_2060")

# Make an empty data frame.
space_jame <- data.frame()

# Create a loop to get predictions.
for(i in unique(total_results$Species)){
  
  # Filter data by species.
  solo_uno <- filter(total_results, Species == i)
  
  # Make three models: one for individual species and two for the community (weighted 
  # and un-weighted).
  one_species_model <- lm(Estimated_DOY ~ avg_annual_precip + avg_tmp_winter,
                          data = solo_uno)
  
  # Define weights for a weighted least squares regression based on CI width.
  CI_weights <- 1 / sqrt(abs(total_results$Upper_CI - total_results$Lower_CI))
  
  # Create a multiple weighted linear regression.
  wt_lm <- lm(Estimated_DOY ~ avg_annual_precip + avg_tmp_winter,
              data = total_results, weights = CI_weights)
  
  # Create a linear mixed effects model with Species and Period as random factors.
  kablooey <- lmer(Estimated_DOY ~ avg_tmp_winter + avg_annual_precip +
                     (1 | Species), data = total_results)
  
  # Filter data by species.
  fut_only_one <- filter(fut_weather, Species == i)
  new_weather <- fut_only_one[,c(2,4,5)]
  
  # Get predictions for all three models.
  answer_one_species <- predict(one_species_model, newdata = new_weather)
  answer_community_weighted <- predict(wt_lm, newdata = new_weather)
  answer_mixed_effects <- predict(kablooey, newdata = new_weather)
  
  # Bind results together.
  answer <- cbind(answer_one_species, answer_community_weighted, answer_mixed_effects)
  
  # Bind combined results to a data frame.
  space_jame <- rbind(space_jame, answer)
}

# Filter the total results retrieve the most recent time period.
p_2000_to_2019_df <- filter(total_results, Period == "p_2000_to_2019")

# Turn estimated DOY's for the fourth historical period into a vector.
fourth_period_onsets <- p_2000_to_2019_df$Estimated_DOY

# Average across all three models to get future predictions.
average_predictions <- rowMeans(space_jame)

# Bind the data frame together, rename columns, then gather into period and DOYs.
final_girl <- cbind(fourth_period_onsets, space_jame, average_predictions)
column_names <- c("p_2000_to_2019", "p_2041_to_2060_1", "p_2041_to_2060_2",
                  "p_2041_to_2060_3", "p_2041_to_2060_4")
colnames(final_girl) <- (column_names)
past_and_predics <- gather(final_girl, "Period", "Estimated_DOYs", 1:length(final_girl))

# Export mixed effects model results as a table.
tab_model(kablooey, show.se = TRUE, show.stat = TRUE, show.p = FALSE,
          pred.labels = c("Intercept", "Mean Winter Temperature",
                          "Mean Annual Rainfall"),
          dv.labels = c("Mixed Effects Model"),
          string.pred = "Fixed Effects",
          string.ci = "95% CI",
          string.stat = "t-value",
          string.se = "Standard Error")

# Export weighted model results as a table.
tab_model(wt_lm, show.se = TRUE,
          pred.labels = c("Intercept", "Mean Winter Temperature",
                          "Mean Annual Rainfall"),
          dv.labels = c("Weighted Model"),
          string.pred = "Predictors",
          string.ci = "95% CI",
          string.p = "p-value",
          string.se = "Std. Error")

## ======================================================================================

# Comparing Future Predictions and Recent Observations (Figure 7)----

# Set colors.
red_colors <- brewer.pal(n = 9, name = "Reds")
red_colors <- c(red_colors[6], red_colors[6], red_colors[6],
                red_colors[6], red_colors[2])
green_colors <- brewer.pal(n = 9, name = "Greens")
green_colors <- c(green_colors[7], green_colors[3], green_colors[3],
                  green_colors[3], green_colors[3])
blue_colors <- brewer.pal(n = 9, name = "Blues")
blue_colors <- c(blue_colors[7], blue_colors[7], blue_colors[7],
                 blue_colors[7], blue_colors[2])

# Create labels for all five time periods.
time_periods <- c("1940_to_1959", "1960_to_1979",
                  "1980_to_1999", "2000_to_2019",
                  "2041 to 2060")
time_periods <- gsub("_", " ", time_periods)

# Plot past and future annual rainfall in this system.
past_and_future_rains_plot <- ggplot(past_and_future_weather,
                                     aes(x = Period, y = avg_annual_precip, fill = Period)) +
  geom_boxplot() +
  scale_x_discrete(name = "Time Period", 
                   labels = time_periods) +
  scale_fill_manual(values = blue_colors) +
  scale_y_continuous(name = "Average Annual Rainfall (mm)") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm"),
        legend.position = "none",
        text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1))
past_and_future_rains_plot

# Plot past and future winter temperature in this system.
past_and_future_temps_plot <- ggplot(past_and_future_weather,
                                     aes(x = Period, y = avg_tmp_winter, fill = Period)) +
  geom_boxplot() +
  scale_x_discrete(name = "Time Period", 
                   labels = time_periods) +
  scale_fill_manual(values = red_colors) +
  scale_y_continuous(name = "Average Winter Temp (°C)") +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = unit(c(.5, .5, .5, .5), "cm"),
        legend.position = "none",
        text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1))
past_and_future_temps_plot

# Plot the model predictions for future flowering compared to contemporary observed values.
model_predictions_plot <- ggplot(past_and_predics,
                                 aes(x = Period, y = Estimated_DOYs, fill = Period)) +
  geom_boxplot() +
  scale_x_discrete(name = "Time Period", 
                   labels = c("2000 to 2019", "Ordinary Least Squares",
                              "Weighted Least Squares", "Mixed Effects Model",
                              "Three Model Averages"),
                   guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values = green_colors) +
  scale_y_continuous(name = "Estimated Day of First Flowering", 
                     breaks = c(-91.5, -60.5, -30.5, 0.5, 31.5, 59.5, 90.5, 120),
                     limits = c(-91.5, 120),
                     labels = c(" ", "OCT ", "NOV ", "DEC ", 
                                "JAN ", "FEB ",
                                "MAR ", "APR ")) +
  ggtitle("Recent Observations vs. Future Predictions") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = "solid",
                                        colour = "lightgrey"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = -0.45, 
                                   hjust = 1.2),
        plot.margin = unit(c(.5, .5, .5, .5), "cm"),
        legend.position = "none",
        text = element_text(size = 13))
model_predictions_plot

# Arrange all three boxplots in the device window at the same time.
all_boxplots <- ggarrange(model_predictions_plot,
                          ggarrange(past_and_future_rains_plot, past_and_future_temps_plot,
                                    ncol = 2, labels = c("B", "C")),
                          nrow = 2, labels = "A") 
all_boxplots

# Calculate average weather conditions in every time period.
past_and_future_weather %>% group_by(Period) %>% summarise(mean_pre = mean(avg_annual_precip),
                                                           mean_tmp = mean(avg_tmp_winter))
# Get means for future predictions.
mean(final_girl$p_2041_to_2060_1)
mean(final_girl$p_2041_to_2060_2)
mean(final_girl$p_2041_to_2060_3)
mean(final_girl$p_2041_to_2060_4)

## ======================================================================================
