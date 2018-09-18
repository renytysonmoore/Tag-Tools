#convert time to seconds of the day
# i.e., 10:05:05 is 10*60*60 + 05*60 + 05

time2sec <- function(hr,min,sec){
  x <- hr*60*60 + min*60 + sec
  return(x)
}