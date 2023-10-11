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

mydata <- fulldata

#index for tasks of interest

interest <- 5:7


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

iteration <- 20000

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

for (j in 1:k){
  theta[1,j] <- initial_value_for_theta[j]           
}

#hyper-parameters for the priors

#prior for Dirichlet of componential proportions  
delta <- rep(2,k)

#prior for Beta of common success odds, and there is k components

alpha <- rep(0,k)
beta <- rep(0,k)

alpha <- rep(1,k)
beta <- rep(1,k)

#define a set store most recent state of parameters

parameter_state <- list(omega=initial_value_for_w)

for (i in 1:n){
  
  latent_variable <- paste0("zeta_", i)
  parameter_state[[latent_variable]] <- z[1,i]
  
}


for (j in 1:k){
  
  success_odd <- paste0("theta_",j)
  parameter_state[[success_odd]] <- theta[1,][j]
  
}


#define a pmf for sequential Bernoulli with common odds

sequntialBernoulli = function(y,m,p){
  pmf<-(p^(sum(y)))*((1-p)^(m-sum(y)))
  return(pmf)
}

#gibbs sampling

for (l in 2:iteration) {
  
  #first generate a random permutation of sampling for each iteration
  
  sampling_permutation <- sample(names(parameter_state))
  number_of_parameters <- length(sampling_permutation)
  
  for (each_parameter in 1:number_of_parameters){
    
    parameter <- sampling_permutation[each_parameter]
    
    #sampling the proportions w
    
    if (substr(parameter,1,1)=="o"){
      
      #number of z_i equal to component j
      
      n_z_i <- rep(0,k)
      
      for (j in 1:k){
        
        for (i in 1:n){
          if (parameter_state[[paste0("zeta_",i)]]==j) {
            n_z_i[j] <- n_z_i[j]+1
          }
        }
      }
      
      updated_w <- rdirichlet(1,delta + n_z_i)
      
      for (j in 1:k){
        w[l,j] <- updated_w[j]
      }
      
      parameter_state$omega <- updated_w
      
    }
    
    #sampling the latent variables z_i
    
    
    else if(substr(parameter,1,1)=="z"){
      
      #return value i for indicating which latent z_i
      value_i <- as.numeric(substr(parameter,6,nchar(parameter)))
      
      p_z_i <- rep(0,k)
      pi_z_i_equal_j_given_rest <- rep(0,k)
      
      for (j in 1:k) {
        
        #set the indicator for the tasks here
        
        p_z_i[j] <- sequntialBernoulli(as.numeric(mydata[value_i,][interest]),m,parameter_state[[paste0("theta_",j)]]) * parameter_state$omega[j]
      
      }
      
      for (j in 1:k){
        pi_z_i_equal_j_given_rest[j] <- p_z_i[j]/sum(p_z_i)
      }
      
      updated_z <- discrete_inverse_sampling(pi_z_i_equal_j_given_rest)
      
      z[l,value_i] <-updated_z
      
      parameter_state[[parameter]] <- updated_z
      
    }
    
    #sampling the joint probabilities theta_j
    
    else if(substr(parameter,1,1)=="t"){
      
      #return value j for indicating which componential theta_j
      
      value_j <- as.numeric(substr(parameter,7,nchar(parameter)))

      sums_over_n_first_term <- rep(0,n)
      sums_over_n_second_term <- rep(0,n)
      
      for (i in 1:n){
        sums_over_n_first_term[i] <- idenfunction(parameter_state[[paste0("zeta_",i)]],value_j) * sum(as.numeric(mydata[i,][interest]))
        sums_over_n_second_term[i] <- idenfunction(parameter_state[[paste0("zeta_",i)]],value_j) * (m-sum(as.numeric(mydata[i,][interest])))
      }    
      
      alpha_prime <- alpha[value_j] + sum(sums_over_n_first_term)
      beta_prime <- beta[value_j] + sum(sums_over_n_second_term)
      
      updated_theta <- rbeta(1,alpha_prime,beta_prime)
      
      theta[l,value_j]<-updated_theta
      
      parameter_state[[parameter]] <- updated_theta
      
    }
    
    else {
      print("something is wrong")
    }
    
  }
  print(l)
}


pdf("plots_all_grade_risk_averse_sequential_Bernoulli_20000_random_permutation_simulations.pdf",width = 8, height = 6)
plot(w[,1],type = "l",xlab = "iteration",ylab = "omega_1",ylim = c(0,1),cex.lab = 1.5)
plot(w[,2],type = "l",xlab = "iteration",ylab = "omega_2",ylim = c(0,1),cex.lab = 1.5)
plot(theta[,1],type = "l",xlab = "iteration",ylab = "theta_1",ylim = c(0,1),cex.lab = 1.5)
plot(theta[,2],type = "l",xlab = "iteration",ylab = "theta_2",ylim = c(0,1),cex.lab = 1.5)
dev.off()

hist(w[,2], breaks = 30, col = 'blue', main = 'Histogram of MCMC Samples', xlab = 'Parameter Values', ylab = 'Density')