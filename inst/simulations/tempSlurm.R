require(dplyr)
require(WalkingDead2)
require(slurmR)

NN <- seq(50, 200, length=3) # Sample size options
RR <- seq(0.25, 0.75, length=3) # Infection Rate options
FF <- c(TRUE, FALSE) # Flight or fight
NNRRFF <- expand.grid(NN,RR,FF) # Combinations of NN and RR
NNRRFFrep <- NNRRFF %>% slice(rep(1:n(), each=500))

simwd2 <- function(j){
  NN <- NNRRFFrep[j,1]
  RR <- NNRRFFrep[j,2]
  FF <- NNRRFFrep[j,3]
  testRun <- wd2(n = NN, # Population size
                 m = round(0.15*NN), # Initial Infected Population size
                 k = round(0.10*NN), # # Initial Immunized Population size
                 fixedAllocation = FALSE, # Allocation type: randomized, if FALSE
                 infect = RR, # infection rate
                 infect_radius = .075, # radius of infection
                 nIterations = 50, # Number of runs
                 flight = FF # Infected will flee from immunized if TRUE
  )
  keepThis <- snapshot(testRun)
  return(keepThis)
}

require(parallel)
cl <- parallel::makeCluster(16)
clusterExport(cl,c("NNRRFFrep", "simwd2"))#,"wd2","simwd2","snapshot"))
clusterEvalQ(cl, library(WalkingDead2))
temp <- do.call(rbind.data.frame,parLapply(cl,1:dim(NNRRFFrep)[1],simwd2))
stopCluster(cl)
saveRDS(temp, "SimWD2Out.rds")
#
# cl <- makeSlurmCluster(16, partition = "notchpeak-shared-short",
#                        account = "notchpeak-shared-short")
#
# clusterExport(cl,c("NNRRFFrep", "simwd2"))#,"wd2","simwd2","snapshot"))
# clusterEvalQ(cl, library(WalkingDead2))
# temp <- do.call(rbind.data.frame,parLapply(cl,1:16,simwd2))
# stopCluster(cl)
# saveRDS(temp, "SimWD2Out.rds")
#
# njobs <- 10L
# temp <- do.call(rbind.data.frame,
#                 Slurm_sapply(1:dim(NNRRFF)[1],
#                              simwd2,
#                              njobs    = njobs,
#                              plan     = "collect",
#                              tmp_path = "/scratch/general/nfs1/u1372213", # This is where all temp files will be exported
#                              sbatch_opt = list(account = "notchpeak-shared-short",
#                                                 partition = "notchpeak-shared-short"),
#
#                              )
#                 )
#
# saveRDS(temp, "SimWD2Out.rds")

