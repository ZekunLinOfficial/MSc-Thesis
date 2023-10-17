load("C:\\Users\\zekun\\OneDrive\\Desktop\\Dissertation on Bayesian\\Data set.Rdata")

grades <- unique(data_for_presentation$grade)

num_observation <- nrow(data_for_presentation)  

for (i in grades){
  assign(paste0("grade",as.character(i),"_count"),0)
  assign(paste0("grade",as.character(i),"_0_count"),0)
  assign(paste0("grade",as.character(i),"_1_count"),0)
  assign(paste0("grade",as.character(i),"_2_count"),0)
  assign(paste0("grade",as.character(i),"_3_count"),0)
  
}

for (i in 1:num_observation){
  grade <- data_for_presentation[i,"grade"]
  grade <- as.character(grade)
  grade_count_character <- paste0("grade",grade,"_count")
  count <- get(grade_count_character)
  assign(grade_count_character,count+1)
  
  
  task2 <- data_for_presentation[i,"task_5"]
  task3 <-data_for_presentation[i,"task_6"]
  task4 <- data_for_presentation[i,"task_7"]
  tasksum <- task2+task3+task4
  for (j in grades){
    if (grade == j){
      if (tasksum == 0){
        count_0 <- paste0("grade",grade,"_0_count")
        assign(count_0,get(count_0) + 1)
      }else if (tasksum == 1){
        count_1 <- paste0("grade",grade,"_1_count")
        assign(count_1,get(count_1) + 1)
      }else if (tasksum == 2){
        count_2 <- paste0("grade",grade,"_2_count")
        assign(count_2,get(count_2) + 1)
      }else{
        count_3 <- paste0("grade",grade,"_3_count")
        assign(count_3,get(count_3) + 1)
      }
    }
  }
}


par(mfrow=c(3,3))

for (i in grades){
  zero <- get(paste0("grade",as.character(i),"_0_count"))
  one <- get(paste0("grade",as.character(i),"_1_count"))
  two <- get(paste0("grade",as.character(i),"_2_count"))
  three <- get(paste0("grade",as.character(i),"_3_count"))
  
  barplot(c(zero,one,two,three),
          names.arg= c("0","1","2","3"),
          main = paste("grade",as.character(i)))
  
}

