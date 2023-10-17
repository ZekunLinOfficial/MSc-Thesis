#specify the grade of data collected in gibbs for mixture of multivariate Bernoulli

grade <- 3

#plots will be stored in pdf files under current working directory

#plot for all parameters

pdf(paste0("non_iden_grade_",grade,"_plot_risk_averse_multivariate_Bernoulli_",iteration,"_simulations.pdf"),width = 8, height =6)
for (j in 1:k){
  plot(w[,j],type = "l",xlab = "iteration",ylab = paste0("mixture weight ",j),ylim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
}
for (j in 1:k) {
  for (i in 1:q){
    plot(theta[j,i,],type = "l",xlab = "iteration",ylab = paste0("joint probability ",paste0(hmatrix[i,], collapse = "")," of component ",j),ylim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
  }
}
dev.off()

#histgram for all variables

pdf(paste0("non_iden_grade_",grade,"_histgram_risk_averse_multivariate_Bernoulli_",iteration,"_simulations.pdf"),width = 8, height =6)
for (j in 1:k){
  hist(w[,j],breaks = 30, main = '',xlab = paste0("mixture weight ",j),ylab = "density",xlim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
}
for (j in 1:k) {
  for (i in 1:q){
    hist(theta[j,i,],breaks = 30, main = '',xlab = paste0("joint probability ",paste0(hmatrix[i,], collapse = "")," of component ",j),ylab = "density",ylim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
  }
}
dev.off()



#imposing an identifiability constraint, ordering by the joint probability of the permutation given by the first row of hmatrix (all elements take 0) into ascending order

for (l in 1:iteration){
  sorted_indices <- order(theta[,,l][,1])
  for (each_permutation in 1:q)
    theta[,,l][,each_permutation] <- theta[,,l][,each_permutation][sorted_indices]
  w[l, ] <- w[l, sorted_indices]
}



#plot for all variables again

pdf(paste0("iden_grade_",grade,"_plot_risk_averse_multivariate_Bernoulli_",iteration,"_simulations.pdf"),width = 8, height =6)
for (j in 1:k){
  plot(w[,j],type = "l",xlab = "iteration",ylab = paste0("mixture weight ",j),ylim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
}
for (j in 1:k) {
  for (i in 1:q){
    plot(theta[j,i,],type = "l",xlab = "iteration",ylab = paste0("joint probability ",paste0(hmatrix[i,], collapse = "")," of component ",j),ylim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
  }
}
dev.off()


#histgram for all variables again

pdf(paste0("iden_grade_",grade,"_histgram_risk_averse_multivariate_Bernoulli_",iteration,"_simulations.pdf"),width = 8, height =6)
for (j in 1:k){
  hist(w[,j],breaks = 30, main = '',xlab = paste0("mixture weight ",j),ylab = "density",xlim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
}
for (j in 1:k) {
  for (i in 1:q){
    hist(theta[j,i,],breaks = 30, main = '',xlab = paste0("joint probability ",paste0(hmatrix[i,], collapse = "")," of component ",j),ylab = "density",ylim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
  }
}
dev.off()
