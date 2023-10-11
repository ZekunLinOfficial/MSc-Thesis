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

mydata <- grade3_data

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

iteration <- 20000

#number of observation

n <- nrow(mydata)


#m is the dimension of observations

m <- 3

#generate the H_m matrix

set_of_generation<- rep(list(0:1), m)
hmatrix <- expand.grid(set_of_generation)

#number of permutations given by hmatrix

q <- nrow(hmatrix)

#initial value for the variables:

#number of components

k <- 2

#proportions

w <- matrix(0,iteration,k)

initial_value_for_w <- rep(1/k,k)

for (i in 1:k){
  w[1,i] <- initial_value_for_w[i]           
}

#latent variables

z <- matrix(0,iteration,n)

initial_value_for_z <- sample(c(1:k),n,replace = TRUE)

for (i in 1:n){
  z[1,i] <- initial_value_for_z[i]
}

#joint probabilities for multivariate Bernoulli

#k by q matrix, ith row gives the joint probabilities of ith component, so each row should add up to 1

theta <- array(0,dim = c(k,q,iteration))

for (i in 1:k) {
  for (j in 1:q){
    theta[,,1][i,j] <- 1/q 
  }
}

#hyper-parameters for the priors
delta <- rep(2,k)
cappa <- rep(2,q)

#define a density for multivariate Bernoulli

multiber = function(iteration,y,component){
  permutation_index <- which(hmatrix$Var1==as.numeric(y[1]) & hmatrix$Var2==as.numeric(y[2]) & hmatrix$Var3==as.numeric(y[3]))
  mass <- theta[,,iteration][component,permutation_index]
  return(mass)
}



#gibbs sampling
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
      
      p_z_i[j] <- multiber(l-1,mydata[i,][interest],j)*w[l,j]
    }
    for (j in 1:k){
      pi_z_i_equal_j_given_rest[j] <- p_z_i[j]/sum(p_z_i)
    }
    z[l,i] <- discrete_inverse_sampling(pi_z_i_equal_j_given_rest)
  }
  
  
  #sampling the parameters theta_j
  
  for (j in 1:k) {
    
    
    cappa_prime <- rep(0,q)
    
    for (each_permutation in 1:q){
      items_in_sum <- rep(0,n)
      
      for (i in 1:n){
        items_in_product <- rep(0,m)
        
        for (each_dimension in 1:m){ 
          items_in_product[each_dimension] <- idenfunction(mydata[i,each_dimension+1],hmatrix[each_permutation,each_dimension])
          
        }
        the_product <- prod(items_in_product)
        
        items_in_sum[i] <- idenfunction(z[l,i],j)*the_product
        
      }
      
      the_sum <- sum(items_in_sum)
      cappa_prime[each_permutation] <- cappa[each_permutation]+the_sum
      
    }
    updated_theta_j <- rdirichlet(1,cappa_prime)
    
    
    for (each_permutation in 1:q) {
      theta[,,l][j,each_permutation] <- updated_theta_j[each_permutation]
    }
    
  }
  print(l)
}

variables_in_theta_for_plot <- rep(0,iteration)

for (l in 1:iteration){
  variables_in_theta_for_plot[l] <- theta[,,l][2,1] 
  
}

plot(variables_in_theta_for_plot, xlab= "iteration",ylab= "theta" ,type = "l")

plot(w[,1],type = "l")
plot(w[,2],type = "l")


