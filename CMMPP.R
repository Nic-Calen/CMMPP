
#Initial var. parameters 
NumDevices <-  500

NumIterations <-  60

seedNum <- 3

set.seed(seedNum)

#Variable declaration 
distx <- NULL
disty <- NULL
dis <-  NULL
master <-  matrix(0, nrow = 1, ncol = NumIterations)
State <-  array(0, dim = c(  NumDevices , NumIterations, 2))
ActualState <- matrix(0 , nrow = NumDevices , ncol = NumIterations)

NumArrivals <- matrix(0 , nrow = NumDevices , ncol = NumIterations)
PacketSize <- ActualState <- matrix(0 , nrow = NumDevices , ncol = NumIteration)
Arrival <-  NULL
#State <-matrix(0 , nrow = NumDevices , ncol = NumIterations)

masterDevice <-  matrix(0 , nrow = NumDevices , ncol = NumIterations)
#Set Transition Matricies 
Pc <- matrix((c(0,1, 1,0)) , nrow = 2, ncol = 2)
Pu <- matrix((c(1,0,1,0)), nrow = 2 , ncol = 2)

#Two-state model is assumed here 
#State 1 and 2 arrival rates 

StateInt <- matrix(0.5, nrow = 2)

#Functions

Multicast_EDF <-  function( input) {
  norm1 = rnorm(1:100, mean = input-0.03/(input+0.03), sd = 0.005) 
  norm2 = rnorm(1:90 , mean = input+0.1/(input+0.05), sd = 0.01)
  norm3 = rnorm(1:8 , mean = input+0.05/(input+0.05), sd = 0.005)
  norm4 = matrix(c(norm1 , norm2 , norm3), nrow = 1)
  return(norm4)
}

Unicast_EDF <- function(input){
  norm1 = rnorm(1:100, mean = input-0.04/(input+0.04), sd = 0.005) 
  norm2 = rnorm(1:25 , mean = input-0.1/(input+0.1), sd = 0.01)
  norm3 = rnorm(1:8 , mean = input-0.05/(input+0.05), sd = 0.005)
  norm4 = matrix(c(norm1 , norm2 , norm3), nrow = 1)
  return(norm4)
}


#Need Unicast and multicast here 
Unicast <- Unicast_EDF(0.31)
Multicast <- Multicast_EDF(0.3)
  
starttime <-  as.double(proc.time()[3])

#Generating distance parameter - constant over time 

for( device in 1:NumDevices){
  #Generate X and Y values for spatial correlation
  distx[device] = runif(1)
  disty[device] = runif(1)
  
  #Distance from epicenter 
  dis[device] <-  sqrt(distx[device]^2 + disty[device]^2)
}


for ( iteration in 1:NumIterations){ #i
  #Master background process to couple MTC devices to MMPP models 
  
  ArrivalRate <-  matrix(c(0.0005, 1/iteration), nrow =1 , ncol =2)
  
  master[iteration] = sample(Multicast, 1)/max(Multicast)
  
  for( device in 1:NumDevices) #j
  {
    masterDevice[device, iteration] = dis[device]*master[iteration]
    
    #Calculating state probabilities for Markov process 
    #Left stochastic matrices are generated here (column elements sum
    # up to 1)
    Pn =  masterDevice[device , iteration]*Pc + (1 - masterDevice[device, iteration])*Pu
    
    if( iteration == 1)
    {
      
      State[device,iteration, ] = Pn %*% StateInt 
      
      
      # Temporarily storing a 2D matrix for the current device at the
      # current time, and then re-shape it for 2D matrix
      # multiplication (data frame it??)
      
      tmpState =  State[device, iteration,  ] 
      tmpState = array(tmpState , c(2,1))
      
    }
    
    if (iteration > 1)
    {
    
      State[device, iteration,]   <-  Pn%*%tmpState
      
      tmpState =  State[device, iteration, ]
      tmpState <-  array(tmpState, c(2,1))
    }
    
    # Assign states to devices based on most probable outcome
    if (State[device, iteration, 1] >= State[device, iteration, 2])
      Arrival <-  ArrivalRate[1]
    
    if(State[device, iteration, 1] < State[device, iteration, 2])
      Arrival <-  ArrivalRate[2]
    
    # Poisson process modulates the arrival rate + additional
    # processing
    
    ActualState[device, iteration] <- rpois(1, Arrival)
    NumArrivals[device,iteration] = ActualState[device,iteration]*NumIterations;
    PacketSize[device,iteration] = ActualState[device,iteration]*NumDevices;
  }
}

endtime <- as.double(proc.time()[3]) - starttime



timeblock <-  endtime/NumIterations

#Plot processing 

MapX <- distx*1000
MapY <-  disty*1000

timeWindow <-  PacketSize[1:NumDevices , 1:10]


#Finding which devices are operating in regular and alarm mode out of the
#devices that are silent, and returns their spatial coordinates

dat <- which(timeWindow == NumDevices, arr.ind =TRUE)

findRegularD <- dat[,1]
findRegularT <- dat[,2]

regularOP <- matrix(0, nrow=1 , ncol = length(findRegularT))
regularOPX <- matrix(0, nrow=1 , ncol = length(findRegularT))
regularOPY <- matrix(0, nrow=1 , ncol = length(findRegularT))

for ( a in 1:length(findRegularT))
{
  regularOP[a] <- timeWindow[findRegularD[a], findRegularT[a]]
  regularOPX[a] <- MapX[findRegularD[a]]
  regularOPY[a] <- MapY[findRegularD[a]]
}

dat <- which(timeWindow > NumDevices, arr.ind =TRUE)

findAlarmD <- dat[,1]
findAlarmT <- dat[,2]


alarmOP <- matrix(0, nrow=1 , ncol = length(findAlarmT))
alarmOPX <- matrix(0, nrow=1 , ncol = length(findAlarmT))
alarmOPY <- matrix(0, nrow=1 , ncol = length(findAlarmT))


for ( b in 1:length(findAlarmT))
{
  alarmOP[b] <- timeWindow[findAlarmD[b], findAlarmT[b]]
  alarmOPX[b] <- MapX[findAlarmD[b]]
  alarmOPY[b] <- MapY[findAlarmD[b]]
}

InterArrival <- matrix(Inf, 1000,25 )

for(device in 1:NumDevices)
{
   
  tmpIndex <- which(ActualState[device,] > 0)
  
  for( d in 2:(1+length(tmpIndex)))
     InterArrival[device, d-1] <- tmpIndex[d]-tmpIndex[d-1]

}

tmp <- which(InterArrival != Inf )
InterArrivalTimes <-  NULL

for(c in 1:(length(tmp)))
{
  InterArrivalTimes[c] <- InterArrival[tmp[c]]
  
}



On_Duration <- rep(list(0),NumDevices)
off_Duration <- rep(list(0),NumDevices)
Burst_Duration <- rep(list(0),NumDevices)


  

for( device2 in 1:NumDevices)
  {
    tmpI <- which(ActualState[device2,] > 0)
     
    for( l in 1:length(tmpI))
    { 
      
  
      
      if(length(tmpI) == 0  ){
        
        On_Duration[[device2]][l] <- 0
        off_Duration[[device2]][l] <- 0
        Burst_Duration[[device2]][l] <- 0 
      }
        
      else{
        On_Duration[[device2]][l] <- timeblock
        off_Duration[[device2]][l] <- (tmpI[l]-1)*timeblock
        Burst_Duration[[device2]][l] <- list(c(On_Duration[[device2]][l],off_Duration[[device2]][l]))
      
      }
      
          
    }   
       
  }



 
plot(density(InterArrivalTimes))
par(2)
plot(MapX,MapY, xlab = "Location X" , ylab = "Location Y" , col = "red", cex = .2 , xlim = c(0,1000), ylim = c(0,1100))
par(new = TRUE)
plot(regularOPX,regularOPY,xlab = "Location X" , ylab = "Location Y", col ="blue", xlim = c(0, 1000),ylim = c(0,1100))
par(new = TRUE)
plot(alarmOPX,alarmOPY,xlab = "Location X" , ylab = "Location Y", cex = 2, xlim = c(0, 1000),ylim = c(0,1100))
