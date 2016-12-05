# Simulation results used in the paper

This folder contains simulated data of three scenarios with different weather information. The probability of a rainy day and the average daily temperature for Scenarios 1, 2, and 3 are assumed to be 0.2, 0.4,
0.6 and 12, 14, 16, respectively.


Each data frame includes information regarding weather data, sensor data, results of the Guassian SSM and the optimal actions for tillage operations:   

  t: Day number in September.
  MW: Observed soil-water content in the field measured by moisture sensors  
  MP: Mean of posterior distribution for latent variable (forecast error of rainfall-runoff model) in the Gaussian SSM
  SP: Variance of posterior distribution for latent variable (forecast error of rainfall-runoff model) in the Gaussian SSM
  Tem: Temperature data related to average daily temperature in each day. 
  Pre: Precipitation data related to total daily precipitation in each day.
  iMW: Index of MW
  iMP: Index of MP
  iSP: Index of SP
  iTem: Index of Tem
  iPre: Index of Pre
  operation: Tillage operation number (1:ploughing, 2:disk harrowing (harrowing-1), 3:spike tooth harrowing (harrowing-2), 4:planting)
  dayLeft: Number of days remained to finish tillage operation "operation".
  optAction: Optimal action for tillage operations (do.:perform tillage operation, pos.:postpone tillage operation, doF:Forcing to perform              tillage operation due to latest finish times) 
  weight: The optimal value function. 
  scenario: Scenario number. 
  
