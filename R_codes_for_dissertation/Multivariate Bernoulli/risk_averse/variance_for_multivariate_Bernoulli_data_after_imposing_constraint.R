#have to put all collected Rdata file under the directory
library(xtable)
#load any of the Rdata file, it should have a data.frame named "fulldata" which contain the list of grades

load("C:\\Users\\zekun\\OneDrive\\Desktop\\simulation_data\\multivariate\\risk_averse\\by_grade\\grade3_risk_averse_multivariate_Bernoulli_20000_simulations.RData")

list_of_grades <- unique(fulldata$grade)

row_names <- c()

#add k variable names for mixture weight
for (j in 1:k){
  row_names <- append(row_names,paste0("omega_",j))  
}

for (j in 1:k){
  for (i in 1:q){
    row_names <- append(row_names,paste0("theta_",j,"_",paste0(as.character(hmatrix[i,]),collapse = ""))) 
  }
}

table_of_variance <- data.frame(
  grade3 = rep(0,k+k*q),
  grade5 = rep(0,k+k*q),
  grade7 = rep(0,k+k*q),
  grade8 = rep(0,k+k*q),
  grade10 = rep(0,k+k*q),
  grade11 = rep(0,k+k*q),
  grade13 = rep(0,k+k*q)
)

rownames(table_of_variance) <- row_names

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
    if (substr(row_names[i], 1, 1)=="o"){
      the_index <- as.numeric(substr(row_names[i], 7, 7))
      table_of_variance[[paste0("grade",grade)]][i] <- round(var(w[,the_index]),4)
      
    }
    
    else if(substr(row_names[i], 1, 1)=="t"){
      
      #need to be careful here, I'm actually a little bit of lazy, if the index is larger than 9 this won't work
      the_index <- as.numeric(substr(row_names[i], 7, 7))
      permutation_index <- which(hmatrix$Var1==as.numeric(substr(row_names[i], 9, 9)) & hmatrix$Var2==as.numeric(substr(row_names[i], 10, 10)) & hmatrix$Var3==as.numeric(substr(row_names[i], 11, 11)))
      table_of_variance[[paste0("grade",grade)]][i] <- round(var(theta[the_index,permutation_index,]),4)
    }
  }
}

all_objects <- ls()
objects_to_keep <- c("table_of_variance")
objects_to_remove <- setdiff(all_objects, objects_to_keep)
rm(list = objects_to_remove)

save(table_of_variance,file="variance_for_all_grades_multivariate_Bernoulli.Rdata")

latex_code <- xtable(table_of_variance, caption = "variance for parameters from 20000 iterations using Algorithm.1 when fitting
a mixture of two multivariate Bernoulli distributions to different grade", digits = 4)

