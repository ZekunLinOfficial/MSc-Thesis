library(DirichletReg)

#load data

load("C:\\Users\\zekun\\OneDrive\\Desktop\\Dissertation on Bayesian\\Data set.Rdata")

fulldata <- data_for_presentation[c("task_1","task_2","task_3","task_4","task_5","task_6","task_7","subject_id","wave","grade")]

school1_data <- fulldata[fulldata$wave==1,]
school2_data <- fulldata[fulldata$wave==2,]
school3_data <- fulldata[fulldata$wave==3,]

grade3_data <- fulldata[fulldata$grade==3,]
grade5_data <- fulldata[fulldata$grade==5,]
grade7_data <- fulldata[fulldata$grade==7,]
grade8_data <- fulldata[fulldata$grade==8,]
grade10_data <- fulldata[fulldata$grade==10,]
grade11_data <- fulldata[fulldata$grade==11,]
grade13_data <- fulldata[fulldata$grade==13,]

#choose the data of different grades for simulation

mydata <- grade8_data

#index for tasks of interest

interest <- 2:4

#run the following to reset the seed since we load a saved data

rm(.Random.seed, envir=globalenv())


#define inverse sampling

discrete_inverse_sampling <- function( prob ) {
  U  <- runif(1)
  if(U <= prob[1]){
    return(1)
  }
  for(state in 2:length(prob)) {
    if(sum(prob[1:(state-1)]) < U && U <= sum(prob[1:state]) ) {
      return(state)
    }
  }
}


#define an identity function

idenfunction = function(x,y){
  if (x==y){
    return(1)
  }
  else {
    return(0)
  }
} 

#iteration of simulation

iteration <- 5000

#number of observation

n <- nrow(mydata)

#dimension of observations

m <- 3

#number of components

k <- 2

#componential proportions

w <- matrix(0,iteration,k)

#initial setting of proportions

initial_value_for_w <- rep(1/k,k) 

for (i in 1:k){
  w[1,i] <- initial_value_for_w[i]           
}

#latent variables

z <- matrix(0,iteration,n)

#initial setting of membership

initial_value_for_z <- sample(c(1:k),n,replace = TRUE)

for (i in 1:n){
  z[1,i] <- initial_value_for_z[i]
}

#common success odds for sequential univariate Bernoulli Distributions

theta <- matrix(0,iteration,k)

#initial success odds

initial_value_for_theta <- rep(0.5,k)

for (i in 1:k){
  theta[1,i] <- initial_value_for_theta[i]           
}

#hyper-parameters for the priors

#prior for Dirichlet of componential proportions  
delta <- rep(2,k)

#prior for Beta of common success odds, and there is k components

alpha <- rep(0,k)
beta <- rep(0,k)

alpha <- rep(1,k)
beta <- rep(1,k)

#define a pmf for sequential Bernoulli with common odds

sequntialBernoulli = function(y,m,p){
  pmf<-(p^(sum(y)))*((1-p)^(m-sum(y)))
  return(pmf)
}



#gibbs

for (l in 2:iteration) {
  
  #sampling the proportions w
  
  #number of z_i equal to component j
  
  n_z_i <- rep(0,k)
  
  for (j in 1:k){
    
    for (i in 1:n){
      if (z[l-1,i]==j) {
        n_z_i[j] <- n_z_i[j]+1
      }
    }
  }
  updated_w <- rdirichlet(1,delta + n_z_i)
  
  for (j in 1:k){
    w[l,j] <- updated_w[j]
  }
  
  #sampling the latent variables z_i
  
  #define a vector store the probabilities such that Z_i=1 ... Z_i=k
  
  
  for (i in 1:n) {
    p_z_i <- rep(0,k)
    pi_z_i_equal_j_given_rest <- rep(0,k)
    for (j in 1:k) {
      
      #set the indicator for the tasks here
      
      p_z_i[j] <- sequntialBernoulli(as.numeric(mydata[i,][interest]),m,theta[l-1,j]) * w[l,j]
      
    }
    
    for (j in 1:k){
      pi_z_i_equal_j_given_rest[j] <- p_z_i[j]/sum(p_z_i)
    }
    
    z[l,i] <- discrete_inverse_sampling(pi_z_i_equal_j_given_rest)
  }
  
  #sampling the parameters theta_j
  

  
  for (j in 1:k) {
    
    sums_over_n_first_term <- rep(0,n)
    sums_over_n_second_term <- rep(0,n)
    
    for (i in 1:n){
        sums_over_n_first_term[i] <- idenfunction(z[l,i],j) * sum(as.numeric(mydata[i,][interest]))
        sums_over_n_second_term[i] <- idenfunction(z[l,i],j) * (m-sum(as.numeric(mydata[i,][interest])))
    }    
    
    alpha_prime <- alpha[j] + sum(sums_over_n_first_term)
    beta_prime <- beta[j] + sum(sums_over_n_second_term)
    
    theta[l,j] <- rbeta(1,alpha_prime,beta_prime)
    
  }
  
  print(l)
}


#identifiability constraint

for (l in 1:iteration){
  sorted_indices <- order(theta[l, ])
  sorted_row <- theta[l, sorted_indices]
  theta[l, ] <- sorted_row
}


