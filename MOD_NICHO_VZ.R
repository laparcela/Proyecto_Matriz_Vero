
# Function to generate interaction networks based on the niche model of Williams and Martínez (2000)

# Function arguments

# S : number of species (default = 10)

# C : connectance (default = 0.2)

#

# Function output

# A : trophic interaction matrix

mod_nicho <- function(S,C){
  alpha = 1
  beta = (alpha/(2*C))-alpha
  A = matrix(data=0, nrow=S, ncol=S)
  x = rbeta(1,alpha,beta) #alpha=1 and beta=alpha/2C - alpha 
  nicho=runif(S) 
  
  
  interv = rep(0, S)
for(i in 1:S){
  interv[i]= x * nicho[i]
}

  cent = rep(0,S)
 for(i in 1:S){
  ai = interv[i]/2
  bi = min(nicho[i], 1-interv[i]/2)
  cent[i] = ai + (bi-ai)*runif(1)
}
  
  for( i in 1:S){
    for(j in 1:S){
      if((cent[i] - interv[i]/2 <= nicho[j]) && (nicho[j] <= cent[i] + interv[i]/2)) {A[j,i]=1}
    }
  }

  # To remove canibalism
  for(i in 1:S){
    A[i,i] = 0
  }
  A
}


#PMatrix of trophic interactions
matriz_interacciones = mod_nicho(10,0.2)


#Save multiple networks in an array
redes <- function(S, C, n_redes){
  array_redes = array(NA,dim = c(S, S, n_redes))
  for(i in 1:n_redes){
    array_redes[,,i] = mod_nicho(S,C)
  }
  array_redes
}

mul_redes=redes(10, 0.2, 150)

library(GGally)
library(network)
library(sna)
library(ggplot2)

net = network(matriz_interacciones, directed = TRUE)
ggnet2(net)
ggnet2(matriz_interacciones)

#ggnet2 can also size nodes by calculating their in-degree, out-degree, or total (Freeman) degree, using the degree function of the sna package.
ggnet2(matriz_interacciones, size = "degree")

# remove any isolated nodes
x = ggnet2(matriz_interacciones, size = "degree", size.min = 1)


