################################################################################
#
# All necessary packages
#
################################################################################

list.of.packages <- c("tidyverse", "lme4", "parameters", "RHRV", "spiro","remotes")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

install.packages("remotes")
remotes::install_github("wviechtb/esmpack@v0.1-17")

lapply(c(list.of.packages, "esmpack"), require, character.only = TRUE)

################################################################################
#
# Mixed linear models
#
################################################################################

# load data
dat <- read.table("viechtbauer2022_data_esmda_example.txt", header=TRUE, sep="\t", na.strings="", as.is=TRUE)

# or

dat <- read.table(file.choose(), header=TRUE, sep="\t", na.strings="", as.is=TRUE)

# number of subjects and number of rows of data
nsub(dat$id)
nrow(dat)

# summary statistics for variable
summary(dat$mood_cheerf)

# frequency table for variable (to check for out-of-range values)
table(dat$mood_cheerf)



# Intercept only model of lonely mood with random intercept in subjects

m = lmer(mood_lonely ~ 1 +  (1 |id), dat)

# Intraclass Correlation Coefficient 

CovMat <- as.data.frame(VarCorr(m))

CovMat$vcov[1]/(CovMat$vcov[1] + CovMat$vcov[2])



# Model with interactions

m = lmer(mood_relaxed ~ sex*mood_cheerf + (1 |id), dat)

model_parameters(m)


# Modeling between effect

dat <- cbind(dat, demean(dat,"mood_down",group = "id", )) # Creating within and between mood down variables

m = lmer(mood_lonely ~ mood_down_within + mood_down_between + (1 |id), dat)

model_parameters(m,effects = "fixed")


################################################################################
#
# Heart Rate Variability
#
################################################################################

hrv.data  = CreateHRVData() # Create object to store data and analyses

hrv.data = LoadBeatAscii(hrv.data, "beat_ascii.txt") # Load data

hrv.data = SetVerbose(hrv.data,TRUE) # More messages from RHRV package

hrv.data = BuildNIHR(hrv.data) # Build raw Heart Rate time series

PlotNIHR(hrv.data, main = "niHR") # Plot raw Heart Rate time series


hrv.data = FilterNIHR(hrv.data) # Replacing improbable data points

PlotNIHR(hrv.data, main = "niHR") # Plot again

hrv.data = CreateTimeAnalysis(hrv.data,size=600,interval = 7.8125) # Just like the name of the function


# Spline interpolation of HR with 4 data points per second
hrv.data = InterpolateNIHR(hrv.data, freqhr = 4, method = "spline") 

PlotHR(hrv.data, main = "HR") # plot interpolated HR series


hrv.data <- CreateFreqAnalysis(hrv.data)

# Calculate Power Spectrum Density
hrv.data <- CalculatePSD(hrv.data, indexFreqAnalysis = 1,
                         method = "ar", doPlot = F)
#Plot it
PlotPSD(hrv.data, indexFreqAnalysis = 1)


# Calculate Spectrogram - How Power Spectrum vary with time
hrv.data = 
  CalculatePowerBand(hrv.data , indexFreqAnalysis = 1,
                     size = 120, shift = 30, type = "fourier",
                     ULFmin = 0, ULFmax = 0.03, VLFmin = 0.03, VLFmax = 0.05,
                     LFmin = 0.05, LFmax = 0.15, HFmin = 0.15, HFmax = 0.4 )
# Plot it
spectrogram <- PlotSpectrogram(HRVData = hrv.data,
                               size = 120, shift = 30,
                               scale = "logaritmic",
                               freqRange = c(0.04, 0.4))


# Indicate different episodes (conditions) in data. 
hrv.data  <- AddEpisodes(hrv.data,
                         InitTimes = c(0, 2000),
                         Durations = c(2000, 2000),
                         Tags = c("Hypertension", "After Medicine"),
                         Values = c(0, 0))

# Plot data with episodes info
PlotHR(hrv.data, main = "HR",
       Tags = "all",Indexes = c(1, 2))

#Extract Spectrogram from episodes
dHypertension <- SplitPowerBandByEpisodes(hrv.data,Tag = "Hypertension")
dAfter_Medicine <- SplitPowerBandByEpisodes(hrv.data,Tag = "After Medicine")

# Calculate LF/HF ratio for all time windows
Hypertension_LFHF <- dHypertension$InEpisodes$LF / dHypertension$InEpisodes$HF
After_Medicine_LFHF <- dAfter_Medicine$InEpisodes$LF / dAfter_Medicine$InEpisodes$HF

# Plot results
boxplot(Hypertension_LFHF, After_Medicine_LFHF, outline = FALSE,
        col = c("blue", "red"),
        names = c("Hypertension", "After_Medicine"),
        main = "Hypertension vs After_Medicine",
        ylab = "LF/HF ratio")



################################################################################
#
# Electrodermal activity
#
################################################################################

# Load data
EDA = read.delim(file.choose(), sep = " ", header = F)

# Plot raw data
plot(EDA_exp$V1,EDA_exp$V2, type = "l",
     xlab="Time (seconds)", ylab=expression(mu*S))

# Applying 3Hz low pass filter
filt_eda = bw_filter(EDA_exp$V2, n = 4, W = 0.006, zero_lag = TRUE)

# Plot clean data
plot(EDA_exp$V1,filt_eda, type = "l",xlab="Time (seconds)", ylab=expression(mu*S))


# Applying 0.05Hz low pass filter
tonic = bw_filter(filt_eda, n = 2, W = 1e-04, zero_lag = TRUE)

# Plot Tonic component
plot(EDA_exp$V1,tonic, type = "l",xlab="Time (seconds)", ylab="Tonic")


# Calculate Phasic component - we subtract Tonic component from cleaned signal
phasic = filt_eda - tonic

# Plot Phasic component
plot(EDA_exp$V1,phasic, type = "l",xlab="Time (seconds)", ylab="Phasic")

