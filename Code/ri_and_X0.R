#Verónica Zepeda
#verozepeda@ciencias.unam.mx


#Función para asignar las tasas de crecimiento de las especies
#Primero se elige el signo de la tasa de crecimiento de acuerdo al tipo de especie
#Positivo para especies basales
#Negativo para todas las otras especies


#S : número de especies. S siempre es igual a 10
#A : matriz de interacciones tróficas

nicho_taza_crecimiento <- function(S,A){
  ri <- rep(NA,S)
  tipos <- tipoTIBA(A)
  
  #Asigna los signos para cada gremio de especies
  ri[tipos$tyTop] <- -1 
  ri[tipos$tyInter] <- -1 
  ri[tipos$tyBasal] <- 1
  ri
}

signos <- nicho_taza_crecimiento(10, A)
signos <- nicho_taza_crecimiento(10, redes_filtradas[,,1])

#Función que asigna la magnitud de la ri de manera aleatoria dadas ciertas restricciones. 

#r_aleatorizado <- function(ri) {
  #raleatorio <- rep(NA, length(ri))
  #for (i in 1:length(ri)) {
    #if (ri[i] == 1) {
      #en el código original se multiplica por ri pero según yo no tiene sentido 
      #raleatorio[i] <- ri[i] * runif(1, 0, 2.5) # r positivas entre 0 y 2.5 
    #} else {
      #raleatorio[i] <- ri[i] * runif(1, -1, 0) # r negativas entre -1 y 0
    #}
  #}
  #return(raleatorio)
#}

r_aleatorizado2 <- function(ri) {
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



#Para obtener las ri´s de cada especie de una matriz de interacciones tróficas
tasas <- r_aleatorizado2(signos)

#Para obtener las ri´s de cada especie para un array de matrices de interacciones tróficas
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

#Guarda las ri generadas en un objeto llamado ris

ris2 <- mat_ri(10,n_redes = dim(redes_filtradas)[3], redes_filtradas)
save(ris2,file="/Users/veronicazepeda/Documents/PosdoctLANCIS/Datos/ris_BUENO2.Rdata")


#Función para generar las abundancias iniciales de las especies (X0)
abundancias0 <- function(S,n_redes){
  ab <- matrix(NA, nrow = n_redes,ncol= S)
  for(i in 1:n_redes){
    ab[i,] <- runif(10, 0, 1)
  }
  ab
}

X0 = abundancias0(10,149)
save(X0,file="/Users/veronicazepeda/Documents/PosdoctLANCIS/Datos/X0_BUENO.Rdata")

dim(redes_filtradas)
dim(ris)
dim(X0)

