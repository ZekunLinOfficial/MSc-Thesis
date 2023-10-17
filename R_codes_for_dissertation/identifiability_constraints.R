#identifiability constraint



#by ordering theta from  into ascending order, and use this permutation to relabel other parameters as well

for (l in 1:iteration){
  sorted_indices <- order(theta[l, ])
  theta[l, ] <- theta[l, sorted_indices]
  w[l, ] <- w[l, sorted_indices]
}


#by ordering theta_000 into ascending order, correspond to the joint probability of permutation 000, all risk-averse or all imprudent and use this permutation to relabel other parameters as well

for (l in 1:iteration){
  sorted_indices <- order(theta[,,l][,1])
  for (each_permutation in 1:q)
    theta[,,l][,each_permutation] <- theta[,,l][,each_permutation][sorted_indices]
  w[l, ] <- w[l, sorted_indices]
}


#by ordering omega into  ascending order, and use this permutation to relabel other parameters as well

for (l in 1:iteration){
  sorted_indices <- order(w[l, ])
  w[l, ] <- w[l, sorted_indices]
  theta[l, ] <- theta[l, sorted_indices]
}