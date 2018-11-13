# In order for slurm to run, we need to set the wd to be home. Perhaps there's a
# way to improve this, but should be OK for now.
setwd("~")

# Calling rslurm
ans <- rslurm::slurm_call(
  function(i) source(paste0(Sys.getenv("RPROJECT"), "/simulations/dgp.r")),
  params=list(i=1),
  jobname="aphylo-test2"
)

