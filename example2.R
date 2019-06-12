# Set working environment -------------------------------------------------
rm(list=ls());graphics.off()
lapply(c("stats4","alabama"), require, character.only = TRUE)

# Generate data -----------------------------------------------------------
# y = rpois(n=16,lambda=6.8)
y = c(5,7,7,4,4,8,15,7,7,4,7,3,8,5,4,7)
N=length(y)


# ML estimates ------------------------------------------------------------
ll_poisson = function(m){
  return(-sum(dpois(y,lambda=m,log = TRUE)))}

res_ML = stats4::mle(minuslogl = ll_poisson, start=list(m=0.5), method="L-BFGS-B", lower = c(0),control=list("trace" = 2))
lambda_ML = as.numeric(res_ML@coef)



# ME  ---------------------------------------------------------------------
M = 5
z_lambda = seq(from=0,to=max(y),length.out = M) #define supports for lambda

p0 = rep(1/M,M) #starting points 

eval_f = function(x){return(sum(x*log(x+1e-9)))} #obj function
sym_grad_f = deriv(expr=~x*log(x),namevec = "x",function.arg = TRUE); eval_grad_f = function(x,z_lambda,M,N,y){return(sym_grad_f(x))} # symbolic obj gradient
Hbnds = function(x){return(x)} #positive bounds for p

eval_g = function(x){ #Constraints of the problem
  p_lambda = x
  lambda = sum(z_lambda*p_lambda)
  
  ceq = c(
    1-sum(p_lambda),     #normalization constraint 
    -N + sum(y)/lambda #score function
  )
  return(ceq)
}

res = alabama::constrOptim.nl(par = p0,fn = eval_f,gr = eval_grad_f,heq = eval_g,hin=Hbnds,control.outer = list(trace=TRUE,kkt2=TRUE))
res$convergence
p_lambda = res$par
lambda_ME_1 = sum(z_lambda*p_lambda)











