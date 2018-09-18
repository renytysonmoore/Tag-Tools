# Remora/DSG-ST spectrum level analysis
# Provides plot of raw data to select region to analyze
# Data are analyzed with a rectangular window and FFT averaging

# Revised by R Tyson Moore, 2018 Aug 28, to be a function that can read in multiple files to analyze

st_spectrum <- function(path,filename){

library(tuneR)
library(stats)

#path = "/w/loggerhead/R/remora/data"
#filename = "1610903610.160611051721.wav"

# import file as data
rmswavfile <- file.path(path, filename)
gain = 2 # high gain = 33 dB
hydrophoneCal = -180 # spherical piezo = -211 dBV/uPa

wavObj <- readWave(rmswavfile)
# wavObj = sine (100000) # test sine wave
#clipObj <- extractWave(wavObj) # allow user to select range
#data = clipObj@left
data = wavObj@left
data = data - mean(data) #remove DC offset
data = data / 32768.0 #scale to +/-1
rms = 20 * log10(sd(data)) - gain - hydrophoneCal

# Use sample rate to get 1 Hz bin
fftpts = floor(wavObj@samp.rate / 100)
nchunks = floor(length(data) / fftpts)

# loop to average FFTs
SUMDATA = numeric(fftpts)
for (i in 1:nchunks){
  datachunk = data[(fftpts * (i - 1) + 1) : (fftpts * i)]
  DATA = numeric(fftpts)
  DATA = 20 * log10(abs(fft(datachunk)) / (fftpts/2)) - gain - hydrophoneCal
  SUMDATA = SUMDATA + DATA
}
MEANDATA = SUMDATA / nchunks

#plot(SUMDATA, type="l")
#plot(MEANDATA[1:fftpts/2], type="l", xlab="Frequency (Hz)", ylab="Spectrum Level (dB re uPa^2/100Hz)")
z <- list(MEANDATA,rms)
return(z)
}
