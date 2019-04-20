# Set environment ---------------------------------------------------------
rm(list=ls()); graphics.off()
load("simulated_data.rda")
source("ext_functions.R")
library(parallel)
#cl = makeCluster(spec=detectCores(),type = "FORK") # type of cluster

# Define global working variables ------------------------------------------------
B = dim(gendata[[1]]$Y)[2] #number of generated datasets
M = 7 #number of support points for ME

# Run analysis using parallel ---------------------------------------------
data_out = list()
data_out_long = data.frame()

z_a = seq(from=-10,to=10,length.out = M) #define supports for a
z_b = seq(from=-10,to=10,length.out = M) #define supports for b
B=10
for(j in 1:length(gendata)){ #loop over length of design
  cat(paste("\n Design cell n.: ",j),sep="")
  
  # working variables for ME
  K = ncol(gendata[[j]]$X)-1;N = nrow(gendata[[j]]$X); KK = ncol(gendata[[j]]$X)
  p0 = rep(1/M,M+(M*K)) #starting points 
  ones_M = matrix(1,M,1); ones_N = matrix(1,N,1); ones_K = matrix(1,K,1)
  
  cat("\n > Running algorithms")
  par_out = mclapply(X = 1:B, FUN = function(k){
    list(nrs=NRs(X = gendata[[j]]$X,y = gendata[[j]]$Y[,k]),
    me=ME(X = gendata[[j]]$X[,-1],y = gendata[[j]]$Y[,k]))
    },mc.cores = detectCores(),mc.preschedule = TRUE)
  
  cat("\n > Saving current data")
  Y = list(
    NR_out = matrix(unlist(lapply(1:B,function(k)par_out[[k]]$nr[1,])),nrow = B,byrow = TRUE),
    NRF_out = matrix(unlist(lapply(1:B,function(k)par_out[[k]]$nr[2,])),nrow = B,byrow = TRUE),
    ME_out = matrix(unlist(lapply(1:B,function(k)par_out[[k]]$me$res_me)),nrow = B,byrow = TRUE),
    ME_out_extra = list(p_a=matrix(unlist(lapply(1:B,function(k)par_out[[k]]$me$p_a)),nrow=B,byrow=TRUE),
                        P_b=array(unlist(lapply(1:B,function(k)par_out[[k]]$me$P_b)),c(B,M,K)))
  )
  save(Y,file=paste("design_",j,".rda",sep=""))
  data_out[[j]] = Y
  
  data_out_long = rbind(data_out_long,
                        data.frame(sample_id=rep(1:B,KK),
                                   xval=as.vector(Y$NR_out[,1:KK]),
                                   conv1=rep(Y$NR_out[,KK+1],KK),
                                   conv2=NA,
                                   par=rep(paste("b",(seq(1:KK)-1),sep=""),each=B),
                                   separ=rep(gendata[[j]]$iid_separation[1:B],KK),
                                   cell_design=j,
                                   algorithm="NR"),
                        data.frame(sample_id=rep(1:B,KK),
                                   xval=as.vector(Y$NRF_out[,1:KK]),
                                   conv1=rep(Y$NRF_out[,KK+1],KK),
                                   conv2=NA,
                                   par=rep(paste("b",seq(1:KK),sep=""),each=B),
                                   separ=rep(gendata[[j]]$iid_separation[1:B],KK),
                                   cell_design=j,
                                   algorithm="NRF"),
                        data.frame(sample_id=rep(1:B,KK),
                                   xval=as.vector(Y$ME_out[,1:KK]),
                                   conv1=rep(Y$ME_out[,KK+1],KK),
                                   conv2=rep(Y$ME_out[,KK+2],KK),
                                   par=rep(paste("b",seq(1:KK),sep=""),each=B),
                                   separ=rep(gendata[[j]]$iid_separation[1:B],KK),
                                   cell_design=j,
                                   algorithm="ME"))
                           
}
save(data_out,file="design_complete.rda")
save(data_out_long,file="design_complete_long.rda")





