
# Results of the paper "Optimal scheduling of tillage operations using finite-horizon Markov decision process and Gaussian state space models"


###########################################################################################################################

## R script file for getting results used in the paper.
## Do the following
##
## 1) Find the optimal policy
## 2) Create data for the 3 scenarios with different weather information (stored in 3 csv files)
## 3) Plot the results (stored in pdf files)

# remember to set the working dir to ./paper/
library(mdpTillage)
set.seed(10000)
unzip(zipfile="polices.zip",exdir="polices")

useSimPaper <- TRUE   # use the already simulated data?
if (!useSimPaper) {
  message("Use a new simulation.")
  source("optimize_MDP.R")   # find optimal policy
  source("simulate.R")        # simulate the 3 pens (use optimal policy to identify feed-mix and number of pigs)
} else {
  message("Use the data generated in the results_data_paper folder.")
  datS1 <- read.csv2("results_data_paper/datS1.csv")
  datS2 <- read.csv2("results_data_paper/datS2.csv")
  datS3 <- read.csv2("results_data_paper/datS3.csv")
}
source("plot.R", echo = TRUE)    # plot results as tex files
tools::texi2pdf(file = "ScenariosPaper_plot.tex", clean = T)
tools::texi2pdf(file = "OptimalPaper_plot.tex", clean = T)
tools::texi2pdf(file = "ComparePolicyPaper_plot.tex", clean = T)













