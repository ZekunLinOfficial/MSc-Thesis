#have to put all collected Rdata file under the directory
  
  #the list of grades could be returned by the original data, but I just typed it here to make sure all files are there 
  
  list_of_grades <- c(3,5,7,8,10,11,13)
  
  for (grade in list_of_grades){
    
    all_objects <- ls()
    
    # keep the list_of_grades
    objects_to_keep <- c("list_of_grades")
    
    # everything else to be removed before the end of each loop
    objects_to_remove <- setdiff(all_objects, objects_to_keep)
    
    load(paste0("C:\\Users\\zekun\\OneDrive\\Desktop\\simulation_data\\multivariate\\risk_averse\\grade",grade,"_risk_averse_multivariate_Bernoulli_20000_simulations.RData"))
    
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
        hist(theta[j,i,],breaks = 30, main = '',xlab = paste0("joint probability ",paste0(hmatrix[i,], collapse = "")," of component ",j),ylab = "density",xlim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
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
        hist(theta[j,i,],breaks = 30, main = '',xlab = paste0("joint probability ",paste0(hmatrix[i,], collapse = "")," of component ",j),ylab = "density",xlim = c(0,1),cex.lab = 1.45, cex.axis = 1.5)
      }
    }
    dev.off()
    
    # Remove everything else
    rm(list = objects_to_remove)
  }