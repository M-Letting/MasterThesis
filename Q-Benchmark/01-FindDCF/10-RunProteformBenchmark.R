# Initialize time
start_time <- Sys.time()

# Set to not show warnings
options(warn = -1)

setwd(here::here("Q-Benchmark"))
# Run data processing
source("01-FindDCF/01-DataProcessing.R")

# Run COPF
source("01-FindDCF/02-RunCOPF.R")

# Run PeCorA
source("01-FindDCF/03-RunPeCorA.R")

setwd(here::here())

# Run ProteoForge
system("python3 Q-Benchmark/01-FindDCF/04-RunProteoForge.py")

setwd(here::here("Q-Benchmark"))

# Run DCF
source("01-FindDCF/05-RunDCF.R")

# Run DCF Baseline
source("01-FindDCF/06-RunDCFBaselineNormalize.R")

# Run BenchmarkAnalysis + plotting
source("01-FindDCF/07-BenchmarkAnalysis.R")
