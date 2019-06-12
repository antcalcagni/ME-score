# Set working environment -------------------------------------------------
rm(list=ls());graphics.off()
lapply(c("stats4","alabama"), require, character.only = TRUE)

# Generate data -----------------------------------------------------------
#y = rgamma(17,shape = 1.5,rate = 5.5)
y = c(0.09,0.35,0.98,0.20,0.44,0.13,0.25,0.48,0.09,0.45,0.03,0.06,0.18,0.26,0.79,0.36,0.26)
N=length(y)


# ML estimates ------------------------------------------------------------
ll_gamma = function(a,b){
  return(-((a-1)*sum(log(y)) - (b*sum(y)) + (N*a)*log(b) - N*log(gamma(a))))}

res_ML = stats4::mle(minuslogl = ll_gamma, start=list(a=0.5,b=0.5), method="L-BFGS-B", lower = c(1e-2,1e-2),control=list("trace" = 2))
parx_ML = as.numeric(res_ML@coef)



# ME  ---------------------------------------------------------------------
m = log(mean(y)) - mean(log(y))
a0 = 1/(2*m)
b0 = a0/mean(y)

M = 5
z_a = seq(from=1e-4,to=a0+3,length.out = M)
z_b = seq(from=1e-4,to=b0+3,length.out = M)

p0 = rep(1/M,M*2) 

eval_f = function(x){return(sum(x*log(x+1e-9)))} #obj function
sym_grad_f = deriv(expr=~x*log(x),namevec = "x",function.arg = TRUE); eval_grad_f = function(x){return(sym_grad_f(x))} # symbolic obj gradient
Hbnds = function(x){return(x)} #positive bounds for p

eval_g = function(x){ #Constraints of the problem
  p_a = x[1:M]; p_b = x[(1+M):(2*M)]
  a = sum(z_a*p_a);
  b = sum(z_b*p_b)
  
  ceq = c(
    1-sum(p_a),
    1-sum(p_b),
    -sum(y)+((N*a)/b),
    sum(log(y)) + N*log(b) - N*digamma(a)
  )
  return(ceq)
}

res = alabama::constrOptim.nl(par = p0,fn = eval_f,gr = eval_grad_f,heq = eval_g,hin=Hbnds,control.outer = list(trace=TRUE,kkt2=TRUE))
res$convergence
p_a = res$par[1:M]; p_b = res$par[(1+M):(2*M)]
parx_ME = c(sum(z_a*p_a), sum(z_b*p_b))

