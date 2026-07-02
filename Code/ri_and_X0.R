
# Function to assign the growth rates of species
# The sign of the growth rate is chosen according to the species type
# Positive for basal species
# Negative for all other species

# S: number of species. S is always equal to 10
# A: trophic interaction matrix

nicho_taza_crecimiento <- function(S,A){
  ri <- rep(NA,S)
  tipos <- tipoTIBA(A)
  
  ri[tipos$tyTop] <- -1 
  ri[tipos$tyInter] <- -1 
  ri[tipos$tyBasal] <- 1
  ri
}


r_aleatorizado <- function(ri) {
  raleatorio <- rep(NA, length(ri))
  for (i in 1:length(ri)) {
    if (ri[i] == 1) {
      raleatorio[i] <- runif(1, 0, 2.5) # r positivas entre 0 y 2.5
    } else {
      raleatorio[i] <- runif(1, -1, 0) # r negativas entre -1 y 0
    }
  }
  return(raleatorio)
}


#For all the networks. They must be in an array.
#n_redes = dim(new_array2)[3]
mat_ri <- function(S, n_redes,array_redes){
  matriz_signos = matrix(NA, ncol=S, nrow = n_redes)
  matriz_tasas = matrix(NA, ncol=S, nrow = n_redes)
  for(i in 1:n_redes){
    matriz_signos[i,] <- nicho_taza_crecimiento(10, array_redes[,,i])
    matriz_tasas[i,] <- r_aleatorizado2(matriz_signos[i,])
  }
  matriz_tasas
}

## Function to generate the initial abundances of species (X0)
abundancias0 <- function(S,n_redes){
  ab <- matrix(NA, nrow = n_redes,ncol= S)
  for(i in 1:n_redes){
    ab[i,] <- runif(10, 0, 1)
  }
  ab
}



