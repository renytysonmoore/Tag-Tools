# Function to get pitch, roll and yaw (heading) from AMX data
# April 24, 2018, Reny Tyson Moore, rtysonmoore@mote.org
# Modified from otload.R 

prh <- function(A,M){
  ##----------- Calculate pitch, roll, yaw----------------------##
  radPerDeg = 0.0174532925
  
  # roll
  phi = atan2(A$accelY, A$accelZ) #original
  #phi - atan2(A$accelY, A$accelX)
  #phi - atan2(A$accelZ, A$accelX)
  sinAngle = sin(phi)
  cosAngle = cos(phi)
  
  # de-rotate by roll angle
  Bfy = (M$magY * cosAngle) - (M$magZ * sinAngle)
  Bz = (M$magY * sinAngle) + (M$magZ * cosAngle)
  
  Gz = A$accelY * sinAngle + A$accelZ * cosAngle
  
  # theta = pitch angle (-90 to 90 degrees)
  theta = atan(-A$accelX / Gz)
  sinAngle = sin(theta)
  cosAngle = cos(theta)
  
  # de-rotate by pitch angle theta
  Bfx = (M$magX * cosAngle) + (Bz * sinAngle)
  Bfz = (-M$magX * sinAngle) + (Bz * cosAngle)
  
  # Psi = yaw = heading
  psi = atan2(-Bfy, Bfx)
  
  pitch = theta / radPerDeg
  roll = phi / radPerDeg
  head = 180 + (psi / radPerDeg)
  
  prh.data <- as.data.frame(cbind(pitch,roll,head))
  return(prh.data)
  
}
