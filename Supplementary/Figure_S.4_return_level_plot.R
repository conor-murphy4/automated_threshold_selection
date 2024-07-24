#Cases 1-3 , p exceedance prob, par (sig, xi), u true
case123_true_quant <- function(p,par,u){
  x_p <- (par[1]/par[2])*((6*p/5)^(-par[2])-1) + u
  return(x_p)
}

#Case 4, p exceedance prob, p_1 P(<1), par (sig, xi), u true
case4_true_quant <- function(p, par, p_1, u){
  x_p <- ((par[1]+par[2]*(u))/par[2])*((p/(1-p_1))^(-par[2])-1) + u
  return(x_p)  
}

dev.new(width=9.17, height=6,noRStudioGD = TRUE)
par(mfrow=c(1,1),bg='transparent')

p_seq <- seq(0,0.2, by=0.000001)
true_pars <- matrix(c(0.5,0.1, 0.5,0.1, 0.5, -0.05, 0.5,0.1), nrow=4, ncol=2, byrow=T)
colours <- c('red', 'blue', 'grey', 'yellow')
line_types <- c('solid', 'dashed', 'dotted', 'dotdash')
y_p_case4 <- case4_true_quant(p=p_seq, par=true_pars[4,], p_1 = 0.7206618, u=1.0)
plot((1/-log(1-p_seq)), y_p_case4, log='x', lwd=2, lty=line_types[4], type='l', xlab="Return period", ylab = "Return value")
for(i in 1:3){
  par = true_pars[i,]
  print(par)
  y_p = case123_true_quant(p=p_seq, par=par, u=1.0)
  lines((1/-log(1-p_seq)),y_p, lwd=2, lty=line_types[i])
}

legend(x=5, y=8, legend = c("Case 1/2", "Case 3", "Case 4"), lty=c("solid", "dotted", "dotdash"), lwd=2, cex=1.5)

