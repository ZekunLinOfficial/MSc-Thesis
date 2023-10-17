#plot individually for theta in multivariate case

#plot

variables_in_theta_for_plot <- rep(0,iteration)
for (l in 1:iteration){
  variables_in_theta_for_plot[l] <- theta[,,l][1,7] 
  
}
plot(variables_in_theta_for_plot, xlab= "iteration",ylab= "theta",ylim = c(0,1) ,type = "l")

#histgram

variables_in_theta_for_plot <- rep(0,iteration)
for (l in 1:iteration){
  variables_in_theta_for_plot[l] <- theta[,,l][1,7] 
  
}
hist(variables_in_theta_for_plot, breaks = 30, main = 'Histogram of MCMC Samples', xlab = 'Parameter Values', ylab = 'Density')

#plot individually in multivariate case

plot(w[,1],type = "l",xlab = "iteration",ylab = "omega_1",ylim = c(0,1),cex.lab = 1.5)
plot(w[,2],type = "l",xlab = "iteration",ylab = "omega_2",ylim = c(0,1),cex.lab = 1.5)
plot(theta[,1],type = "l",xlab = "iteration",ylab = "theta_1",ylim = c(0,1),cex.lab = 1.5)
plot(theta[,2],type = "l",xlab = "iteration",ylab = "theta_2",ylim = c(0,1),cex.lab = 1.5)

hist(w[,1], breaks = 30, main = '', xlab = 'omega_1', ylab = 'Density',xlim = c(0,1))
hist(w[,2], breaks = 30, main = '', xlab = 'omega_2', ylab = 'Density',xlim = c(0,1))
hist(theta[,1], breaks = 30, main = '', xlab = 'theta_1', ylab = 'Density',xlim = c(0,1))
hist(theta[,2], breaks = 30, main = '', xlab = 'theta_2', ylab = 'Density',xlim = c(0,1))

