#Verónica Zepeda
#verozepeda@ciencias.unam.mx

## Función para proyectar la dinámica poblacional de cada especie y tipo de especie.
#Este código es una traducción y adaptación del código en MatLab reportado en CITA ARTICULO LETY

## ARGUMENTOS:                                                                                 
## X0     : Vector de abundancia inicial, tamaño 1 X S                                  
## K      : Vector de capacidad de carga, tamaño 1 X S. Valores de capacidad de karga K=[ki's.]                                        
## R      : Vector de tasa de crecimiento intrínseco, tamaño 1 X S. Valores de tasa de crecimento intrinseco R=[ri's].                             
## A      : Matriz de interacciones tróficas. 
## w      : tipo de peso que se usará: si w=1, entonces es un vector W = 1, de tamaño 1 X S 
## Ti     : tiempo inicial
## Tf     : tiempo final 
## num    : número de interacciones ¿? numero de pasos desde Ti hasta Tf. ¿a qué se refiere con el número de pasos?
## bound  : cota inferior para valor de abundancia (para contar las especies que sobreviven). Esta la tomaremos como 10^(-6), en general. revisar trabajo irene  

## SALIDA:
## evol_cont  : Matriz de S X Tf. Muestra los cambios en las abundancias de las especies a lo largo del tiempo. 
## abundancia : Vector con las abundancias finales de cada especie. 
## Sp         : número de especies que sobreviven en el estado final (num+1) 

##Ecuación que describe la dinámica poblacional
# dxi/dt = xi [(1-xi/Ki) (ri + ei*Sum(wji*aji*xj)) - Sum(wij*aij*xj) ]

#Cargar la función que clasifica a las especies según su tipo (top, intermedia, basales, aisladas)
source("/Users/veronicazepeda/Documents/PosdoctLANCIS/Código/tipoTIBA_VZ.R")

solEuler_peso <- function(X0,K,R,A,w,Ti,Tf,num,bound){
  S <- dim(A)[1]# num. de especies	
  
  ## Asignando pesos de las especies. Para este estudio W = 1 dado que el peso es igual para todas las especies. 
  W <- rep(1,S)
  
  ## Para asignar la eficiencia de los distintos gremios de especies:                                                                            
  ## basal        :  ei=0                                                                    
  ## intermedia   :  ei=0.66  #referencia de Kefi                                                              
  ## top          :  ei=0.85                                                                 
  
  ei <- rep(0.66,S)
  
  tipos <- tipoTIBA(A)  
  ei[tipos$tyTop] <- 0.85  
  ei[tipos$tyBasal] <- 0
  
  ## tamaño de paso ¿QUÉ ES EL TAMAÑO DE PASO? 
  h <- (Tf-Ti)/(num)
  
  evol_cont <- matrix(data=0, nrow=num+1,S) #Matriz que contiene los cambios en las abundancias para todas las especies a lo largo del tiempo. 
  
  #corresponde al valor inicial de las especies                                        
  evol_cont[1,] = X0                       
  
  ### iteraciones. Empieza en 2, en evol_cont cada renglón corresponde a un paso de tiempo                                                                          
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
  
  
  # contando las especies que sobreviven y no explotan. Cuantifica la riqueza.                             
  Sp <- 0                                                                                   
  for (i in 1:S){
    #Para que no se enoje R cuando explotan los valores
    if (!is.nan(evol_cont[num+1,i]) && evol_cont[num+1,i] >= bound && evol_cont[num+1,i] <= 1){
      Sp <- Sp+1
    }
  }
  
  ### Para obtener la abundancia final
  abundancia <- evol_cont[num+1,]
  return(list(evol_cont = evol_cont, Sp = Sp, abundancia = abundancia))
}
         
