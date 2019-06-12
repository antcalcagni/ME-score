# Set environment ---------------------------------------------------------
rm(list=ls()); graphics.off()
source("ext_functions.R")

# Finney data -------------------------------------------------------------
data("vaso",package = "robustbase")
glm(formula = vaso$Y~log(vaso$Volume)+log(vaso$Rate),family = binomial)

X = cbind(1,log(vaso$Volume),log(vaso$Rate))

M=7
z_a = seq(from=-5,to=5,length.out = M) #define supports for a
z_b = seq(from=-10,to=10,length.out = M) #define supports for b
K = ncol(X)-1;N = nrow(X); KK = ncol(X)
p0 = rep(1/M,M+(M*K)) #starting points 
ones_M = matrix(1,M,1); ones_N = matrix(1,N,1); ones_K = matrix(1,K,1)

res_me = ME(X = X[,-1],y = vaso$Y)
print(res_me)










