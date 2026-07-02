
#Metacommunity dynamics
dinamica_metacomunidad_ciclos_todo <- function(
  paisaje,
  matriz_abundancias,
  ecuacion_lety,
  parametros_ecuacion_lety,
  parches_con_dinamica_comunitaria,
  secuencia_tasas_migracion,
  secuencia_tasas_mortalidad,
  secuencia_tiempo_total,
  pasos_mig,
  pasos_mort,
  tiempo_final_int,
  metodo_integracion,
  array_redes,
  matriz_ris
){

  
  n_redes = dim(array_redes)[3] 
  num_especies <- dim(matriz_abundancias)[2]
  num_elementos_por_paso <- (pasos_mig + pasos_mort) + 1 
  num_total_indices <- (tiempo_total * num_elementos_por_paso) + 1
  array_dinamica_meta <- array(NA, dim = c((num_total_indices*(repeticiones*2)), nrow(paisaje), ncol(paisaje), num_especies, n_redes)) #Array para resultados dinámica metapoblacional.
  
  for(j in 1:n_redes){
    
    
    cat(paste0("j = ", j, "\n"))
    
    #Matrix of initial abundances
    pobs_iniciales <- list(p = matriz_abundancias[j,],
                                  t = rep(0,10),
                                  r = rep(0,10))
    
    
    #Matrix of species growth rates
    tasas_reproduccion = matriz_ris[j,]
    
    #Networks array
    matriz_interacciones = array_redes[,,j]
    
    #EFICIENCY VALUES
    e <- rep(0.66,S)
    e[tipoTIBA(matriz_interacciones)$tyTop] <- 0.85
    e[tipoTIBA(matriz_interacciones)$tyBasal] <- 0
    
    #list of parameters for differential equation of population dynamics
    parametros_ecuacion_lety <- list(k = k,
                                     r = tasas_reproduccion, 
                                     A = matriz_interacciones, 
                                     W = W,
                                     e = e)
    ##To create an initial metacommunity
    metacomunidad <- generar_metacomunidad_inicial(paisaje,num_especies,pobs_iniciales)
    
    
    ####Simulation of metacommunity dynamics
    
    guarda <-  simular_dinamica_metacomunidad_ciclos(
      paisaje,
      pobs_iniciales,
      ecuacion_lety,
      parametros_ecuacion_lety,
      parches_con_dinamica_comunitaria,
      secuencia_tasas_migracion,
      secuencia_tasas_mortalidad,
      secuencia_tiempo_total,
      pasos_mig,
      pasos_mort,
      tiempo_final_int,
      metodo_integracion,
      tasas_reproduccion,
      metacomunidad)
    
    array_dinamica_meta[,,,,j] <- guarda$dinamica
  
    }
  array_dinamica_meta
}

