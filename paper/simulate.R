
# Set MDP parameters for using the basic parameters in the functions
param<-setParam()

# read the optimal polcy based on basic weights for trafficability, workability, and completion criteria
policy <- read.csv2("polices/based_weight/policyMDP.csv", stringsAsFactors = F)

# define the coeficants for the probability of a rainy days in three scenarios
coefPre1<-0.2
coefPre2<-0.4
coefPre3<-0.6

# define the coeficants for the average daily temperature in three scenarios
coefTem1<-16
coefTem2<-14
coefTem3<-12

# generate a set of random numbers to specify wet and dry days
set.seed(10000)
rndValues<-runif(n = param$tMax, min = 0, max = 1)
#length(rndValues[rndValues<0.25])

# initial soil water content
iniTrueWat<-30
givenWatInfo<-FALSE

# Simulate the weather data and find optimal actions in three scenarios.
dataOptimal1<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre1,coefTem1,rndValues)$temData, precData = simWeather(param,coefPre1,coefTem1,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )
dataOptimal2<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre2,coefTem2,rndValues)$temData, precData = simWeather(param,coefPre2,coefTem2,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )
dataOptimal3<-optimalSearch(param = param, policy = policy, temData = simWeather(param,coefPre3,coefTem3,rndValues)$temData, precData = simWeather(param,coefPre3,coefTem3,rndValues)$precData, iniTrueWat = iniTrueWat, givenWatInfo = givenWatInfo )

# Store the results in three csv files

datS1 <- dataOptimal1
datS2 <- dataOptimal2
datS3 <- dataOptimal3
datS1$scenario<-1
datS2$scenario<-2
datS3$scenario<-3

#write on the csv files

# datS1$optAction = as.character(datS1$optAction)
# write.csv2(datS1, file ="results_data_paper/datS1.csv",row.names=FALSE)
#
# datS2$optAction = as.character(datS2$optAction)
# write.csv2(datS2, file ="results_data_paper/datS2.csv",row.names=FALSE)
#
# datS3$optAction = as.character(datS3$optAction)
# write.csv2(datS3, file ="results_data_paper/datS3.csv",row.names=FALSE)


