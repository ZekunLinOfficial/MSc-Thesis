#have to put all collected Rdata file under the directory

#load any of the Rdata file, it should have a data.frame named "fulldata" which contain the list of grades

load("C:\\Users\\zekun\\OneDrive\\Desktop\\simulation_data\\multivariate\\risk_averse\\by_grade\\grade3_risk_averse_multivariate_Bernoulli_20000_simulations.RData")

list_of_grades <- unique(fulldata$grade)

row_names <- c()

#add k variable names for mixture weight
for (j in 1:k){
  row_names <- append(row_names,paste0("w_",j))  
}

for (j in 1:k){
  for (i in 1:q){
    row_names <- append(row_names,paste0("theta_",j,"_",paste0(as.character(hmatrix[i,]),collapse = ""))) 
  }
}

table_of_statistics <- data.frame(
  mean = rep(0,k+k*q),
  variance = rep(0,k+k*q),
  confidence_low = rep(0,k+k*q),
  confidence_up = rep(0,k+k*q)
)

rownames(table_of_statistics) <- row_names

for (grade in list_of_grades){
  
  load(paste0("C:\\Users\\zekun\\OneDrive\\Desktop\\simulation_data\\multivariate\\risk_averse\\by_grade\\grade",grade,"_risk_averse_multivariate_Bernoulli_20000_simulations.RData"))
  
  #imposing an identifiability constraint, ordering by the Componential Bernoulli success probabilities into ascending order
  
  for (l in 1:iteration){
    sorted_indices <- order(theta[,,l][,1])
    for (each_permutation in 1:q)
      theta[,,l][,each_permutation] <- theta[,,l][,each_permutation][sorted_indices]
    w[l, ] <- w[l, sorted_indices]
  }

  for (i in 1:length(row_names)){
    if (substr(row_names[i], 1, 1)=="w"){
      the_index <- as.numeric(substr(row_names[i], 3, 3))
      table_of_statistics[i,1] <- round(mean(w[,the_index]),4)
      table_of_statistics[i,2] <- round(var(w[,the_index]),4)
      table_of_statistics[i,3] <- round(t.test(w[,the_index])$conf.int[1],4)
      table_of_statistics[i,4] <- round(t.test(w[,the_index])$conf.int[2],4)
    
    }
    
    else if(substr(row_names[i], 1, 1)=="t"){
      
      #need to be careful here, I'm actually a little bit of lazy, if the index is larger than 9 this won't work
      the_index <- as.numeric(substr(row_names[i], 7, 7))
      permutation_index <- which(hmatrix$Var1==as.numeric(substr(row_names[i], 9, 9)) & hmatrix$Var2==as.numeric(substr(row_names[i], 10, 10)) & hmatrix$Var3==as.numeric(substr(row_names[i], 11, 11)))
      table_of_statistics[i,1] <- round(mean(theta[the_index,permutation_index,]),4)
      table_of_statistics[i,2] <- round(var(theta[the_index,permutation_index,]),4)
      table_of_statistics[i,3] <- round(t.test(theta[the_index,permutation_index,])$conf.int[1],4)
      table_of_statistics[i,4] <- round(t.test(theta[the_index,permutation_index,])$conf.int[2],4)
    }
  }
  
  
  save(table_of_statistics,file=paste0("statistics_for_grade",grade,"_multivariate_Bernoulli.Rdata"))

}

