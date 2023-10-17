#have to put all collected Rdata file under the directory

#load any of the Rdata file, it should have a data.frame named "fulldata" which contain the list of grades

load("C:\\Users\\zekun\\OneDrive\\Desktop\\simulation_data\\sequential\\risk_averse\\by_grade\\grade3_risk_averse_sequential_Bernoulli_20000_simulations.RData")

list_of_grades <- unique(fulldata$grade)

row_names <- c()

#add k variable names for mixture weight
for (j in 1:k){
  row_names <- append(row_names,paste0("w_",j))  
}

#add k variable names for Bernoulli odd
for (j in 1:k){
  row_names <- append(row_names,paste0("theta_",j))  
}

table_of_statistics <- data.frame(
  mean = rep(0,2*k),
  variance = rep(0,2*k),
  confidence_low = rep(0,2*k),
  confidence_up = rep(0,2*k)
)

rownames(table_of_statistics) <- row_names

all_objects <- ls()
objects_to_keep <- c("list_of_grades","row_names","table_of_statistics")
objects_to_remove <- setdiff(all_objects, objects_to_keep)
rm(list = objects_to_remove)

for (grade in list_of_grades){
  
  all_objects <- ls()
  
  # keep the list_of_grades
  objects_to_keep <- c("list_of_grades","row_names","table_of_statistics")
  
  # everything else to be removed before the end of each loop
  objects_to_remove <- setdiff(all_objects, objects_to_keep)
  
  load(paste0("C:\\Users\\zekun\\OneDrive\\Desktop\\simulation_data\\sequential\\risk_averse\\by_grade\\grade",grade,"_risk_averse_sequential_Bernoulli_20000_simulations.RData"))
  
  #imposing an identifiability constraint, ordering by the Componential Bernoulli success probabilities into ascending order
  
  for (l in 1:iteration){
    sorted_indices <- order(theta[l, ])
    theta[l, ] <- theta[l, sorted_indices]
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
      the_index <- as.numeric(substr(row_names[i], 7, nchar(row_names[i])))
      table_of_statistics[i,1] <- round(mean(theta[,the_index]),4)
      table_of_statistics[i,2] <- round(var(theta[,the_index]),4)
      table_of_statistics[i,3] <- round(t.test(theta[,the_index])$conf.int[1],4)
      table_of_statistics[i,4] <- round(t.test(theta[,the_index])$conf.int[2],4)
    }
  }
  
  
  save(table_of_statistics,file=paste0("statistics_for_grade",grade,"_sequential_Bernoulli.Rdata"))
}

