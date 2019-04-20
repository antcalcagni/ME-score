eval_f = function(x,X,y){return(sum(x*log(x+1e-4)))} #obj function
sym_grad_f = deriv(expr=~x*log(x),namevec = "x",function.arg = TRUE); eval_grad_f = function(x,X,y){return(sym_grad_f(x))} # symbolic obj gradient
Hin = function(x,X,y){return(x)}

eval_g = function(x,X,y){ #Constraints of the problem
  p_a = x[1:M]; P_b = matrix(x[(M+1):(M+K*M)],K,M)
  ax = as.numeric(z_a%*%p_a)
  bx = kronecker(diag(K),t(z_b))%*%matrixcalc::vec(t(P_b))
  p = 1/(1+exp(-(ax+X%*%bx)))
  
  ceq = c(
    1-p_a%*%ones_M,     #normalization constraint 
    ones_K-P_b%*%ones_M,  #normalization constraint 
    t(ones_N)%*%(y-p), #score eq. for a
    t(X)%*%(y-p) #score eq. for b
  )
  return(ceq)
}

ME = function(X,y){
  res_me = matrix(NA,1,K+3)
  p_a = rep(NA,M)
  P_b = matrix(NA,K,M)
  
  tryCatch(
    expr = {
      res = alabama::constrOptim.nl(par = p0,fn = eval_f,gr = eval_grad_f,hin = Hin,heq = eval_g,control.outer = list(trace=FALSE,kkt2=TRUE),X=X,y=y)    
      res$convergence
      res$K #the smallest, the better
      p_a = res$par[1:M]; P_b = matrix(res$par[(M+1):(M+K*M)],K,M)
      a_me = z_a%*%p_a
      b_me = kronecker(diag(K),t(z_b))%*%matrixcalc::vec(t(P_b))
      res_me[1,] = c(a_me,b_me,res$convergence,res$K)
    },
    error = function(e){print(e)},
    #warning = function(w){},
    finally = {}
  )
  
  return(list(res_me=res_me,
              p_a=p_a,
              P_b=P_b))}

NRs = function(X,y){
  res_glm = matrix(NA,2,ncol(X)+1)
  
  tryCatch(
    expr = {
      fit_glm = glm(formula = y~-1+X,family = binomial)
      res_glm[1,] = c(fit_glm$coef,fit_glm$converged)
      fit_glmf = logistf::logistf(formula = y~-1+X)
      res_glm[2,] = c(fit_glmf$coef,NA)
    },
    error = function(e){print(e)},
    #warning = function(w){},
    finally = {}
  )
  return(res_glm)
}




