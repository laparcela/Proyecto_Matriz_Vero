
## Function to project the population dynamics of each species and species type

## ARGUMENTS:
## X0     : Initial abundance vector, size 1 × S
## K      : Carrying capacity vector, size 1 × S. Carrying capacity values K = [k_i]
## R      : Intrinsic growth rate vector, size 1 × S. Growth rate values R = [r_i]
## A      : Trophic interaction matrix
## w      : Type of weight to use: if w = 1, then W is a vector of ones, size 1 × S
## Ti     : Initial time
## Tf     : Final time
## num    : Number of time steps between Ti and Tf
## bound  : Lower bound for abundance values (used to count surviving species).
##          Typically set to 10^(-6); see Irene’s work for reference

## OUTPUT:
## evol_cont  : Matrix of size S × Tf. Shows changes in species abundances over time
## abundancia  : Vector with the final abundances of each species
## Sp         : Number of species that survive in the final state (num + 1)

## Equation describing population dynamics
# dxi/dt = xi [(1-xi/Ki) (ri + ei*Sum(wji*aji*xj)) - Sum(wij*aij*xj) ]

source("./tipoTIBA_VZ.R")

solEuler_peso <- function(X0,K,R,A,w,Ti,Tf,num,bound){
  S <- dim(A)[1]# Number of species
  
  ## For the purpose of this study W = 1 
  W <- rep(1,S)
  
 ## To assign the efficiency of different species guilds:
## basal        : ei = 0
## intermediate : ei = 0.66  # reference from Kéfi
## top          : ei = 0.85                                                             
  
  ei <- rep(0.66,S)
  
  tipos <- tipoTIBA(A)  
  ei[tipos$tyTop] <- 0.85  
  ei[tipos$tyBasal] <- 0
  
  ## Size of the integration step
  h <- (Tf-Ti)/(num)
  
  evol_cont <- matrix(data=0, nrow=num+1,S) #                                   
  evol_cont[1,] = X0                       
                                                                        
  for (t in 2:(num+1)){                                                                       
    for (i in 1:S){                                                                       
      cap <- 1-(evol_cont[t-1,i]/K[i])                                                  
      pos <- sum(W[i]*A[,i]*evol_cont[t-1,])   
      #arreglo : wij * aij                                                             
      wa <- W * A[i,]# renglon i : los que se comen a i  
      neg <- sum(wa * evol_cont[t-1,])                                                 
      
      f1 <- evol_cont[t-1,i] * (cap * (R[i] + ei[i]*pos) - neg)                                                                                                          
      evol_cont[t,i] <- evol_cont[t-1,i] + h*f1
    }
    # print(paste("t = ", t))
    # print(evol_cont)
  }
  
  
  # For species richness                            
  Sp <- 0                                                                                   
  for (i in 1:S){

    if (!is.nan(evol_cont[num+1,i]) && evol_cont[num+1,i] >= bound && evol_cont[num+1,i] <= 1){
      Sp <- Sp+1
    }
  }
  
  #For final abundance
  abundancia <- evol_cont[num+1,]
  return(list(evol_cont = evol_cont, Sp = Sp, abundancia = abundancia))
}
         
