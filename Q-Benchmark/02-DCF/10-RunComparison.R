# Set working directory to project root
setwd(here::here())

# Source scripts for creating datasets
source("Q-Benchmark/02-DCF/00-MakeProteinComplexData.R")
setwd(here::here())
source("Q-Benchmark/02-DCF/00-MakeProteoMakerData.R")
setwd(here::here())

# Source scripts for data processing
source("Q-Benchmark/02-DCF/01-DataProcessing.R")
setwd(here::here())

# Source scripts for running methods
source("Q-Benchmark/02-DCF/02-RunPecora.R")
setwd(here::here())
system("python3 Q-Benchmark/02-DCF/03-RunProteoForge.py")
setwd(here::here())
source("Q-Benchmark/02-DCF/04-RunDCF.R")
setwd(here::here())
source("Q-Benchmark/02-DCF/04-RunDCFBaseline.R")
setwd(here::here())
source("Q-Benchmark/02-DCF/05-RunCOPF.R")
setwd(here::here())
source("Q-Benchmark/02-DCF/04-RunRPC.R")
setwd(here::here())
# source("Q-Benchmark/02-DCF/04-RunRPC2-3.R")
# setwd(here::here())

# Source scripts for evaluation
source("Q-Benchmark/02-DCF/06-CompareDisocveryMethods.R")
setwd(here::here())
source("Q-Benchmark/02-DCF/07-CompareMapping.R")
setwd(here::here())
source("Q-Benchmark/02-DCF/08-PlotCompareDiscoveryMethods.R")
