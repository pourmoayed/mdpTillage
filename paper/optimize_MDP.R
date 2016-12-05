
# Find optimal policy based on the basic weights for trafficability, workability and completion criteria
##################################################################################################################
# Set MDP parameters:
param<-setParam(weightCompletion=1,weightWorkable=1,weightTraffic=1)

# Solve the MDP model and restore the results in a csv file named "policyMDP" in the root directory.
SolveMDPModel(param)

# Read the optimal policy of the MDP
policy <-read.csv2("../policyMDP.csv", stringsAsFactors = F)

# Store the policy in the related directory
write.csv2(policy, file ="polices/based_weight/policyMDP.csv",row.names=FALSE)
##################################################################################################################



# Find optimal policy based on the differentweights for trafficability, workability and completion criteria (Group 1 in the paper)
##################################################################################################################
# Set MDP parameters:
param<-setParam(weightCompletion=1,weightWorkable=0.2,weightTraffic=0.8)

# Solve the MDP model and restore the results in a csv file named "policyMDP" in the root directory.
SolveMDPModel(param)

# Read the optimal policy of the MDP
policy <-read.csv2("../policyMDP.csv", stringsAsFactors = F)

# Store the policy in the related directory
write.csv2(policy, file ="polices/weight_high_work/policyMDP.csv",row.names=FALSE)
##################################################################################################################



# Find optimal policy based on the differentweights for trafficability, workability and completion criteria (Group 2 in the paper)
##################################################################################################################
# Set MDP parameters:
param<-setParam(weightCompletion=1,weightWorkable=0.8,weightTraffic=0.2)

# Solve the MDP model and restore the results in a csv file named "policyMDP" in the root directory.
SolveMDPModel(param)

# Read the optimal policy of the MDP
policy <-read.csv2("../policyMDP.csv", stringsAsFactors = F)

# Store the policy in the related directory
write.csv2(policy, file ="polices/weight_high_traf/policyMDP.csv",row.names=FALSE)
##################################################################################################################

