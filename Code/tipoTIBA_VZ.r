
## Function to classify species in an adjacency matrix according to their type (Top, Intermediate, Basal, Isolated)

## ARGUMENTS:
##   A : interaction (adjacency) matrix

## OUTPUT:
##   numtype : 1 × 4 array giving the number of species in each category (Top, Intermediate, Basal, Isolated)
##   tyTop   : array of species classified as Top
##   tyInter : array of species classified as Intermediate
##   tyBasal : array of species classified as Basal
##   tyIsol  : array of species classified as Isolated

library(sna) 

tipoTIBA <- function(A){
    S <- dim(A)[1] # number of species                                 
    numtipo <- rep(0,4)
                                               
    tyTop <- c()   # vector for top species                       
    tyInter <- c() # vector for intermediate species                     
    tyBasal <- c() # vector for basal species                         
    tyAis <- c()   # vector for isolated species
    
    outd <- degree(A, cmode = "outdegree", gmode = "digraph")        
    intd <- degree(A, cmode = "indegree", gmode = "digraph")
     
    for (j in 1:S){                                                      
        if (outd[j] == 0 && intd[j] > 0){
            numtipo[1] <- numtipo[1] + 1                              
            tyTop <- c(tyTop, j)                                           
        } else if (outd[j] > 0 && intd[j] > 0){  
            numtipo[2] <- numtipo[2] + 1                             
            tyInter <- c(tyInter, j)                                       
        } else if (outd[j] > 0 && intd[j] == 0){       
            numtipo[3] <- numtipo[3] + 1                         
            tyBasal <- c(tyBasal, j)                                    
        } else if (outd[j] == 0 && intd[j] == 0){    
            numtipo[4] <- numtipo[4] + 1                           
            tyAis <- c(tyAis, j)                                           
        }                                                                
    }
    return(list(numtipo = numtipo, tyTop = tyTop, tyInter = tyInter, tyBasal = tyBasal, tyAis = tyAis))
} 
