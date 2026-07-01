#Verónica Zepeda
#verozepeda@ciencias.unam.mx

## Función que clasifica a las especies de una matríz de adyacencia de acuerdo a su tipo (Top, Intermedias, Basales, Aisladas)
#Este código es una traducción y adaptación del código en MatLab reportado en CITA ARTICULO LETY

## Argumentos:                                                           
##   A : matriz de interacciones

## SALIDA:                                                            
##   numtipo : arreglo de 1X4, da el numero de especies T,I,B y A     
##   tyTop   : arreglo de las especies que son Top                    
##   tyInter : arreglo de las especies que son Intermedias            
##   tyBasal : arreglo de las especies que son Basales                
##   tyAis   : arreglo de las especies que son Aisladas

library(sna) #para usar outdegree, indegree and degree

tipoTIBA <- function(A){
    S <- dim(A)[1] # numero de especies                                   
    numtipo <- rep(0,4)
    
    ## identifica la especie que corresponde a cada tipo y lo guarda en un vector                                             
    tyTop <- c()   # vector para especies superiores                        
    tyInter <- c() # vector para especies intermedias                      
    tyBasal <- c() # vector para especies basales                          
    tyAis <- c()   # vector para especies aisladas 
    
    ## calcular grado de entrada y salida de la red
    outd <- degree(A, cmode = "outdegree", gmode = "digraph")        
    intd <- degree(A, cmode = "indegree", gmode = "digraph")
    
    # contando tipos de especies y agregandolas a los vectores de cada tipo   
    for (j in 1:S){                                                      
        if (outd[j] == 0 && intd[j] > 0){# contando especies superior : T
            numtipo[1] <- numtipo[1] + 1                              
            tyTop <- c(tyTop, j)                                           
        } else if (outd[j] > 0 && intd[j] > 0){# especie intermedia: I"   
            numtipo[2] <- numtipo[2] + 1                             
            tyInter <- c(tyInter, j)                                       
        } else if (outd[j] > 0 && intd[j] == 0){# especie basal: B"       
            numtipo[3] <- numtipo[3] + 1                         
            tyBasal <- c(tyBasal, j)                                    
        } else if (outd[j] == 0 && intd[j] == 0){# especie aislada: A"    
            numtipo[4] <- numtipo[4] + 1                           
            tyAis <- c(tyAis, j)                                           
        }                                                                
    }
    return(list(numtipo = numtipo, tyTop = tyTop, tyInter = tyInter, tyBasal = tyBasal, tyAis = tyAis))
} 
