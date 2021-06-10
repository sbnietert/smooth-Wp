# imports ----
library(pdist)
library(doParallel)
library(WassersteinGoF)
library(ggplot2)
library(reshape2)
library(parallel)
library(purrr)
library(pracma)

# data records ----
records = list()

# metadata for run ----
d <- 5 # 1,3
num_samples <- floor(logspace(log10(30),log10(1000),9))[1:7]
num_runs <- array(10, dim=length(num_samples)) #floor(num_samples[length(num_samples)] * min_num_runs / num_samples)
sigmas <- c(0, 0.02, 0.05)
scale <- 100
MC_its <- 20000
metadata = list(dim=d, num_samples=num_samples, num_runs=num_runs, sigmas=sigmas, MC_its=MC_its, scale=scale, onesample=TRUE)

# create labeled data matrix ----
sigma_names = map(sigmas, function(s) paste("sigma = ", toString(s)))
num_samples_names = map(num_samples, function(n) paste("num samples = ", toString(n)))
run_names = map(1:num_runs[1], function(r) paste("run # ", toString(r)))

data <- array(0L, dim=c(length(sigmas), length(num_samples), num_runs[1]), dimnames=list(sigma_names,num_samples_names,run_names))

# compute distances ---
for(i in 1:length(num_samples)) {
  n <- num_samples[i]
  print(n)
  N <- MC_its
  for(j in 1:num_runs[i]) {
    xs_base <- matrix(scale * runif(n*d),ncol=d)
    xs_base2 <- xs_base[sample(nrow(xs_base),size=N,replace=TRUE),]
    xt_base <- matrix(scale * runif(N*d),ncol=d)
    noise_s <- matrix(rnorm(N*d),ncol=d)
    noise_t <- matrix(rnorm(N*d),ncol=d)
    for(k in 1:length(sigmas)) {
      sigma <- sigmas[k]
      xs <- xs_base2 + scale * sigma*noise_s
      xt <- xt_base + scale * sigma*noise_t
      data[k,i,j] = sqrt(max(WassersteinDist(xs,xt,p=2),0))
    }
  }
}

# store data in records ----
timestamp = format(Sys.time(), "%a %b %d %H:%M:%S %Y")
metadata$time = timestamp
records[[length(records)+1]] <- list(metadata = metadata, data = data)
save(records, file=paste(timestamp,".RData"))
