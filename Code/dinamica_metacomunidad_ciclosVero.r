library(abind)
library(deSolve)

#' Simulation of metacommunity dynamics in a landscape
#'
#' For an ecological community, simulates metacommunity dynamics
#' where at each time iteration a community dynamic occurs
#' along with multiple iterations of migration and mortality dynamics.
#' 
#' @param paisaje A matrix of character strings, where each string indicates the type of patch in the landscape
#' @param condiciones_iniciales A vector of length S with the initial population sizes
#' for each species (the order of these values must correspond to the same order of species in the
#' interaction matrix)
#' @param funcion_dinamica_comunitaria Function that describes the differential equation
#' of population dynamics. Required parameter for the
#' ode function (see [deSolve::ode()])
#' @param parametros_dinamica_comunitaria List of parameters received by
#' the population dynamics function. Parameters in the list
#' must be named in the same way as they are named in the function.
#' Required parameter for the ode function (see [deSolve::ode()]).
#' @param parches_con_dinamica_comunitaria A character vector with the
#' names of the patch types where community dynamics are executed.
#' @param tasas_migracion A data.frame with two columns. One column
#' must be named "tipo" (type) and contain the different patch types present
#' in the landscape. The other column must be named "tasa_dispersion" (dispersal rate) and
#' contain the dispersal rates corresponding to each patch type
#' @param tasas_mortalidad A data.frame with two columns. One column
#' must be named "tipo" and contain the different patch types present
#' in the landscape. The other column must be named "tasa_mortalidad" (mortality rate) and
#' contain the mortality rates corresponding to each patch type
#' @param tiempo_total Number of times the community dynamics
#' will be executed during the simulation
#' @param pasos_mig_mort Number of times the migration and mortality dynamics
#' will be executed for each iteration of the community dynamics
#'
#' @return The first dimension contains
#' ((total_time * (migration_mortality_steps + 1)) + 1) entries that
#' correspond to the total number of iterations of the community dynamics and
#' migration-mortality dynamics of the model, plus the initial conditions.
#' The second and third dimensions represent the position of a patch in the
#' landscape, and therefore have sizes nrow(landscape) and ncol(landscape), respectively.
#' The fourth dimension contains ncol(interaction_matrix), where each entry
#' represents the population size of a species. Thus, the entry (k, x, y, i)
#' of the resulting array represents the population size of species i
#' in the patch with coordinates (x, y) in the landscape during iteration k
#' of the model dynamics.

simular_dinamica_metacomunidad <- function(paisaje,
                                           metacomunidad,
                                           funcion_dinamica_comunitaria,
                                           parametros_dinamica_comunitaria,
                                           parches_con_dinamica_comunitaria,
                                           tasas_migracion,
                                           tasas_mortalidad,
                                           tiempo_total,
                                           pasos_mig = 5,
                                           pasos_mort = 1,
                                           tiempo_final_int = 1,
                                           metodo_integracion="euler",
                                           tasas_reproduccion){
    
    num_especies <- length(tasas_reproduccion)
    num_elementos_por_paso <- (pasos_mig + pasos_mort) + 1 #This last one is the step of the population dynamics.
    num_total_indices <- (tiempo_total * num_elementos_por_paso) + 1
    ## The main array where all data are stored is initialized.
    dinamica_metacomunidad <- array(0, dim = c(num_total_indices,
                                               nrow(paisaje),
                                               ncol(paisaje),
                                               num_especies))
    dinamica_metacomunidad[1, , , ] <- metacomunidad
    tiempos <- c(0)
    indices <- c(1)
    procesos <- c("init")
    
    for (t in (1:tiempo_total)) {
      ## T defines the indices in the first dimension of `dinamica_metacomunidad`
      ## where community dynamics results are stored.
      ## Intermediate indices correspond to migration–mortality steps.
      ## +1 ensures indexing starts at 2, since 1 stores initial conditions.
      
        tiempos <- append(tiempos, rep(t, num_elementos_por_paso))
        T <- (t * num_elementos_por_paso) - (num_elementos_por_paso - 1) + 1
        dinamica_metacomunidad[T, , , ] <-
            dinamica_comunitaria(
                paisaje,
                dinamica_metacomunidad[T-1, , , ],
                parches_con_dinamica_comunitaria,
                funcion_dinamica_comunitaria,
                parametros_dinamica_comunitaria,
                tiempo_final_int)
        procesos <- append(procesos, rep("DP",1))
        indices <- append(indices, T)
        indice <- tail(indices,1) + 1
        if (pasos_mig > 0){
            for (k in 1:pasos_mig){
                dinamica_metacomunidad[indice, , , ] =
                    migracion_metacomunidad(
                        paisaje,
                        dinamica_metacomunidad[indice-1, , , ],
                        tasas_migracion)
                indice <- indice + 1
                indices <- append(indices, indice)
                procesos <- append(procesos, "migracion")
               
            }
          dinamica_metacomunidad[indice, , , ] =
            mortalidad_metacomunidad(
              paisaje,
              dinamica_metacomunidad[indice-1, , , ],
              tasas_mortalidad)
          indice <- indice + 1
          indices <- append(indices, indice)
          procesos <- append(procesos, "mortalidad")
          
        } else {
            print("Los pasos_mig_mort son cero! No hay dinámica espacial")
        }
    }
    return(list(dinamica=dinamica_metacomunidad,
                tiempos=tiempos,
                indices=indices,
                procesos=procesos))
}

#' Fill a landscape with initial species populations
#'
#' Assigns initial population sizes of each species to every patch
#' in the landscape according to patch type.
#'
#' @param paisaje Character matrix indicating patch types.
#' @param num_especies Number of species in the community.
#' @param poblaciones_iniciales Named list of numeric vectors with
#' initial population sizes for each patch type. Names must match
#' the patch types in paisaje.
#'
#' @return 3D array representing the initial metacommunity:
#' [x, y, i] = population of species i at patch (x, y).
#' 
generar_metacomunidad_inicial <- function(paisaje,
                                          num_especies,
                                          poblaciones_iniciales){
    for (i in 1:length(poblaciones_iniciales)){
        if (length(poblaciones_iniciales[[i]]) != num_especies){
            stop ("No todas las poblaciones iniciales tienen la misma longitud")
        }
    }
    if (! all( unique( as.vector(paisaje) %in% names(poblaciones_iniciales) ) ) ){
        stop ("Falta definir poblaciones iniciales para algunos de parche")
    }
    
    metacomunidad_inicial <- array(0, dim = c(nrow(paisaje),
                                              ncol(paisaje),
                                              num_especies))
    for (x in 1:nrow(paisaje)){
        for (y in 1:ncol(paisaje)){
            tipo_parche <- paisaje[x,y]
            metacomunidad_inicial[x,y,] <- poblaciones_iniciales[[tipo_parche]]
        }
    }
    return(metacomunidad_inicial)
}

#' Run community dynamics
#'
#' Applies community dynamics to each patch in a landscape using
#' a specified function. Uses the ode function from the deSolve package.
#' See [deSolve::ode()].
#'
#' @param paisaje Character matrix indicating patch types.
#' @param metacomunidad 3D array representing the metacommunity at a given time.
#' The first and second dimensions correspond to patch position in the landscape,
#' and the third to species. Entry (x, y, i) represents the population
#' of species i at patch (x, y).
#' @param parches_con_dinamica_comunitaria Character vector with patch types
#' where community dynamics occur.
#' @param tiempos Numeric vector of time points for which population dynamics
#' are computed. Required by ode (see [deSolve::ode()]).
#' @param funcion Function describing the differential equation of population
#' dynamics. Required by ode (see [deSolve::ode()]).
#' @param parametros List of parameters passed to the population dynamics
#' function. Names must match those used in funcion. Required by ode
#' (see [deSolve::ode()]).
#'
#' @return 3D array representing the metacommunity after applying
#' community dynamics. Dimensions correspond to [x, y, i], where each
#' entry is the population of species i at patch (x, y).
dinamica_comunitaria <- function(paisaje,
                                 metacomunidad,
                                 parches_con_dinamica_comunitaria,
                                 funcion,
                                 parametros,
                                 tiempo_final_int
                                 ){
  metacomunidad_dinamica <- array(0, dim = dim(metacomunidad))
    for (x in 1:nrow(paisaje)){
        for (y in 1:ncol(paisaje)){
            if (paisaje[x,y] %in% parches_con_dinamica_comunitaria){
              metacomunidad_dinamica[x,y, ] <- as.matrix.data.frame(
                ode(y = metacomunidad[x,y, ],
                    times = seq(0,tiempo_final_int,by=0.01),
                    func = funcion,
                    parms = parametros, method = "euler"))[(tiempo_final_int/0.01) + 1,-1]#Para tomar el ultimo valor de la abundancia
            } else {
              metacomunidad_dinamica[x,y, ] = metacomunidad[x,y, ]
            }
        }
    }
  return(metacomunidad_dinamica)
}


#' Compute species migration in the landscape
#'
#' Computes migration dynamics across the landscape based on
#' species population sizes and migration rates.
#'
#' @param paisaje Character matrix indicating patch types.
#' @param metacomunidad 3D array representing the metacommunity at a given time.
#' The first and second dimensions correspond to patch position in the landscape,
#' and the third to species. Entry (x, y, i) represents the population
#' of species i at patch (x, y).
#' @param tasas_migracion A data.frame with two columns. One column
#' must be named "tipo" and contain the different patch types in the landscape.
#' The other column must be named "tasas_migracion" and contain the
#' dispersal rates corresponding to each patch type.
#'
#' @return A 3D array where entry (x, y, i) represents the population
#' size of species i at patch (x, y) after applying migration dynamics.
migracion_metacomunidad <- function(paisaje,
                                    metacomunidad,
                                    tasas_migracion){

tipos_en_paisaje <- unique(as.vector(paisaje))
    tipos_en_tasas_migracion <- tasas_migracion$tipo
    if (! all(tipos_en_paisaje %in% tipos_en_tasas_migracion)){
        stop ("No se ha definido la tasa de migración para todos los tipos de celda definidos en el paisaje")
    }
        
    
    metacomunidad_perdida = array(0, dim = dim(metacomunidad))
    metacomunidad_ganancia = array(0, dim = dim(metacomunidad))
    metacomunidad_migracion = array(0, dim = dim(metacomunidad))
    num_especies = dim(metacomunidad)[3]
    
    for (i in 1:num_especies){
        metapoblacion <- metacomunidad[, , i]
        metapoblacion_perdida = array(0, dim = dim(metapoblacion))
        metapoblacion_ganancia = array(0, dim = dim(metapoblacion))
        
        for (x in 1:nrow(paisaje)) {
            for (y in 1:ncol(paisaje)) {
                metapoblacion_perdida[x, y] <- metapoblacion[x, y] *
                    tasas_migracion$tasa_migracion[
                                        tasas_migracion$tipo == paisaje[x, y]]
            }
        }
     
        for (x in 1:nrow(paisaje)) {
            for (y in 1:ncol(paisaje)) {
               
                suma_perdida_vecinos <-
                    arreglo_enrroscado(metapoblacion_perdida, x+1, y-1) + # ar-iz
                    arreglo_enrroscado(metapoblacion_perdida, x+1, y  ) + # ar
                    arreglo_enrroscado(metapoblacion_perdida, x+1, y+1) + # ar-de
                    arreglo_enrroscado(metapoblacion_perdida, x  , y-1) + # iz
                    arreglo_enrroscado(metapoblacion_perdida, x  , y+1) + # de 
                    arreglo_enrroscado(metapoblacion_perdida, x-1, y-1) + # ab-iz
                    arreglo_enrroscado(metapoblacion_perdida, x-1, y  ) + # ab
                    arreglo_enrroscado(metapoblacion_perdida, x-1, y+1)   # ab-de
                metapoblacion_ganancia[x, y] <- suma_perdida_vecinos / 8 
            }
        }
   
        metacomunidad_perdida[, , i] <- metapoblacion_perdida
        metacomunidad_ganancia[, , i] <- metapoblacion_ganancia
        metacomunidad_migracion[, , i] <- metapoblacion + (metapoblacion_ganancia -
                                                           metapoblacion_perdida)
    }
    return(metacomunidad_migracion)
}

#' Get elements from an array using wrap-around indexing
#'
#' Helper function for migration calculations. 
#' @param B array
#' @param i row index i
#' @param j column index j
#'
#' @return The element (i, j) of array B, where indexing wraps
#' horizontally and vertically (e.g., in a 3 x 2 array, indices (0, -1)
#' return entry (3, 1)).
#' 
#' 
arreglo_enrroscado <- function(B, i, j){
    n_rows = dim(B)[1]
    n_cols = dim(B)[2]
    i = i %% n_rows
    j = j %% n_cols
    if ( i == 0 ){ i = n_rows }
    if ( j == 0 ){ j = n_cols }
    return ( B[i,j] )
}


#' Compute species mortality in the landscape
#'
#' Computes mortality dynamics across the landscape based on
#' species population sizes and mortality rates.
#' 
#' Arguments: 
#' @param paisaje 
#' @param metacomunidad 
#' @param tasas_migracion 
#' 
#' 
#' @return #' A 3D array where entry (x, y, i) represents the population
#' size of species i at patch (x, y) after applying mortality dynamics.
#' 
mortalidad_metacomunidad <- function(paisaje,
                                     metacomunidad,
                                     tasas_mortalidad){
    metacomunidad_mortalidad <- array(0, dim = dim(metacomunidad))
    for (x in 1:nrow(paisaje)){
        for (y in 1:ncol(paisaje)){
            tasa_muerte <-
                tasas_mortalidad$tasa_mortalidad[
                                     tasas_mortalidad$tipo == paisaje[x,y]]
            metacomunidad_mortalidad[x,y,] <-
                (1 - tasa_muerte) * metacomunidad[x,y, ]
        }
    }
    return(metacomunidad_mortalidad)
}

simular_dinamica_metacomunidad_ciclos <- function(paisaje,
                                                  condiciones_iniciales,
                                                  funcion_dinamica_comunitaria,
                                                  parametros_dinamica_comunitaria,
                                                  parches_con_dinamica_comunitaria,
                                                  secuencia_tasas_migracion,
                                                  secuencia_tasas_mortalidad,
                                                  secuencia_tiempo_total,
                                                  pasos_mig = 5,
                                                  pasos_mort = 1,
                                                  tiempo_final_int = 1,
                                                  metodo_integracion = "euler",tasas_reproduccion, metacomunidad){
                                                      #num_especies <- 10
    #metacomunidad <- generar_metacomunidad_inicial(paisaje,
    #num_especies,
    #                                               condiciones_iniciales)
    dinamica_metacomunidad_ciclos <- list()
    tiempos <- c()
    indices <- c()
    procesos <- c()
    ciclos <- c()
    ultimo_tiempo_ciclo <- 0
    ultimo_indice_ciclo <- 0
    
    for(i in 1:length(secuencia_tiempo_total)){
        tiempo_total = secuencia_tiempo_total[i]
        tasas_migracion = as.data.frame(secuencia_tasas_migracion[i])
        tasas_mortalidad = as.data.frame(secuencia_tasas_mortalidad[i])
        res_dinamica_metacomunidad = simular_dinamica_metacomunidad(
            paisaje,
            metacomunidad,
            funcion_dinamica_comunitaria,
            parametros_dinamica_comunitaria,
            parches_con_dinamica_comunitaria,
            tasas_migracion,
            tasas_mortalidad,
            tiempo_total,
            pasos_mig,
            pasos_mort,
            tiempo_final_int,
            metodo_integracion,
            tasas_reproduccion)
        dinamica_metacomunidad = res_dinamica_metacomunidad$dinamica
        tiempos = append(tiempos,
                         res_dinamica_metacomunidad$tiempos + ultimo_tiempo_ciclo )
        ultimo_tiempo_ciclo <- tail(tiempos, 1)
        indices = append(indices,
                         res_dinamica_metacomunidad$indices + ultimo_indice_ciclo )
        ultimo_indice_ciclo <- tail(indices, 1) + 1
        procesos <- append(procesos, res_dinamica_metacomunidad$procesos)
        ciclos <- append(ciclos,
                         rep(i,tiempo_total*
                               ((pasos_mig + pasos_mort) + 1)+1)) 
        ultima_entrada <- dim(dinamica_metacomunidad)[1]
        metacomunidad <- dinamica_metacomunidad[ultima_entrada, , , ]

        longitud <- length(dinamica_metacomunidad_ciclos)
        dinamica_metacomunidad_ciclos[[longitud + 1]] <- dinamica_metacomunidad
    }

    return( list(dinamica=abind(dinamica_metacomunidad_ciclos, along = 1),
                 tiempos=tiempos,
                 indices=indices,
                 procesos=procesos,
                 ciclo=ciclos))
}
