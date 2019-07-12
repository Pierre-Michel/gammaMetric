# Script for functions calculating the rate of increase up to order n


# Function calculating the rate of increase in order 1
diff.1 = function(x, t)
{
  # x is an element of a list containing the samples
  x = unlist(x)
  l = length(x)
  
  res = (x[2:l] - x[1:(l-1)])/t[1:(l-1)]
  return(res)
}



# Recursive function calculating the rate of increase up to order n
diff.n = function(x, t, n)
{
  x = unlist(x)
  
  # If order = 0 we return the same sample
  if(n==0){res = x}
  
  # If n = 1 we apply diff.1
  else if(n==1){res = diff.1(x, t)}
  
  # If n > 1 we apply diff.n with n = n-1
  else
  {
    # Computation of order n-1
    x.df = diff.n(x, t, n-1) 
    l = length(x.df)
    
    # Computation for order n
    res = (x.df[2:l] - x.df[1:l-1])/t[2:l]
  }
  return(res)
}

# Function computing mean and standard deviation of an element of a list
calcul.ind = function(lst)
{
  x = unlist(lst)
  res = cbind(mean(x), sd(x))
  return(res)
}

# Function computing the mean/sd of the rate of increase up to ordre 'ordre"
create_ts_df = function(ts_list, ordre = 1)
{
  # Initialization
  res = NULL
  
  for(i in 0:ordre)
  {
    cat("Calcul de l'ordre",i,"\n")
    # Computation the rate of increase of order i
    ll = lapply(ts_list, FUN = function(x){diff.n(x, t = unlist(x), n = i)})
    
    # Computation of the mean/sd
    ll = lapply(ll, calcul.ind)
    ll_mat = do.call("rbind",ll)
    colnames(ll_mat) = paste(c("mn.","sd."), rep(i,2),"df",sep="")
    res = cbind(res, ll_mat)
  }  
  return(data.frame(res))
}
