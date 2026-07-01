
# Function to count species in each trophic level 
filtro <- function(array_A, n_redes){
  n_tipo = matrix(NA, ncol=4,nrow = n_redes)
  colnames(n_tipo)=c("top","inter","basal","ais")
  for(i in 1:n_redes){
    tipos <- tipoTIBA(array_A[,,i])
    n_tipo[i,1] = tipos$numtipo[1]
    n_tipo[i,2] = tipos$numtipo[2]
    n_tipo[i,3] = tipos$numtipo[3]
    n_tipo[i,4] = tipos$numtipo[4]
  }
  n_tipo
}


aa=filtro(array_redes,n_redes = dim(array_redes)[3])

sp_ais=which(aa[,4] > 0) #Networks with isolated species

new_array = array_redes[,,-sp_ais] 

bb=filtro(new_array,n_redes = dim(new_array)[3]) 

which(bb[,4] > 0) #check for isolated species
which(bb[,3] >= 1) #networks with at least one basal species
no_bas=which(bb[,3] == 0)

redes_filtradas = new_array[,,-no_bas] #array of network witouth isolated species and at least one basal species
pru=filtro(redes_filtradas,n_redes = dim(redes_filtradas)[3])
which(pru[,4] > 0)
which(pru[,3] == 0)

