#Code to run all the simulated experiments with different periods 
#of temporal variability and different initial conditions (dry or rainy)
#for all network configurations

#interactions networks
load("./array_redes_filtradas.Rdata")

#Growth rates of each species in the communities
load("./ris.Rdata")

##Initial abundances
load("./X0.Rdata")

matriz_abundancias = X0 #matrix of species abundances
matriz_ris = ris2 #matrix of population growth rates
array_redes = redes_filtradas


#Metacommunity functions
source("./dinamica_metacomunidad_ciclosVero.r")
source("./tipoTIBA_VZ.r")
source("./loop_metaCiclos.R")



###Required arguments:

## differential equation of population dynamics
ecuacion_lety <- function(t,x,parametros){
  with(as.list(c(x, parametros)), {
    dx <- x * ( (1 - (x / k)) * ( r + ( e * ((t(W * A)) %*% x))) - (( W * A ) %*% x))
    return(list(dx))
  })
}

S <- 10 #number of species
k <- rep(1,S) #Carrying capacity was always 1
W <- matrix(rep(1,S*S), nrow=S, ncol=S, byrow=TRUE)

#Simulated landscape
paisaje <- matrix(c("p",  "t",  "r",  "t",  "r",  "t",  "p",  "p",  "p",  "p",  
                    "r",  "p",  "r",  "p",  "r",  "t",  "r",  "t",  "r",  "p",  
                    "p",  "t",  "r",  "t",  "r",  "t",  "r",  "t",  "r",  "t",  
                    "r",  "p",  "p",  "t",  "p",  "p",  "r",  "t",  "r",  "t",  
                    "r",  "t",  "p",  "t",  "r",  "p",  "p",  "p",  "r",  "p",  
                    "r",  "t",  "r",  "t",  "r",  "t",  "p",  "p",  "r",  "p",  
                    "p",  "t",  "p",  "t",  "r",  "t",  "r",  "t",  "r",  "t",  
                    "r",  "p",  "r",  "t",  "r",  "p",  "r",  "t",  "r",  "t",  
                    "r",  "t",  "r",  "p",  "r",  "t",  "r",  "t",  "r",  "t",  
                    "r",  "t",  "r",  "p",  "p",  "t",  "r",  "p",  "r",  "p"), nrow=10, ncol=10,byrow = T)

#indices patches with community dynamics
parches_con_dinamica_comunitaria <- c("p")

tiempo_final_int <- 1 #units of time to evaluated the differential equation of population dynamics 
metodo_integracion = "euler" #integration method
pasos_mig=5 #number of migration steps
pasos_mort=1 #number of mortality steps

######################################
###Migration and mortality rates
#####################################


#Rainy conditions
tasas_migracion_A <- data.frame(tipo=c("p","t","r"),
                                tasa_migracion=c(0.2,0.5,0.9))
tasas_mortalidad_A <- data.frame(tipo=c("p","t","r"),
                                 tasa_mortalidad=c(0.01,0.25,0.70))

#Dry conditions
tasas_migracion_B <- data.frame(tipo=c("p","t","r"),
                                tasa_migracion=c(0.4,0.7,0.9))

tasas_mortalidad_B <- data.frame(tipo=c("p","t","r"),
                                 tasa_mortalidad=c(0.05,0.6,0.7))

#######################################
########SIMULATED EXPERIMENTS
#######################################

#########################
#####Rainy conditions
##########################
#low temporal variability starting with rainy conditions 

tiempo_total <- 50 
repeticiones <- 1
secuencia_tiempo_total <- rep(c(tiempo_total,
                                tiempo_total),
                              times = repeticiones)
secuencia_tipo_ciclo <- rep(c("A","B"), times = repeticiones)


secuencia_tasas_migracion <- rep(list(tasas_migracion_A,
                                      tasas_migracion_B),
                                 times = repeticiones)


secuencia_tasas_mortalidad <- rep(list(tasas_mortalidad_A,
                                       tasas_mortalidad_B),
                                  times = repeticiones)


resultado.ciclos_50_A <- simular_dinamica_metacomunidad_ciclos(
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
  metodo_integracion)

resultado.ciclos_50_A = dinamica_metacomunidad_ciclos_todo(paisaje,matriz_abundancias,ecuacion_lety,parametros_ecuacion_lety,parches_con_dinamica_comunitaria,secuencia_tasas_migracion,secuencia_tasas_mortalidad,secuencia_tiempo_total,pasos_mig,pasos_mort,tiempo_final_int,metodo_integracion,array_redes,matriz_ris)
save(resultado.ciclos_50_A,file="./resultado.ciclos_50_A.Rdata")


#Moderate temporal variability starting with rainy conditions
tiempo_total <- 20
repeticiones <- 3
secuencia_tiempo_total <- rep(c(tiempo_total,
                                tiempo_total),
                              times = repeticiones)
secuencia_tipo_ciclo <- rep(c("A","B"), times = repeticiones)


secuencia_tasas_migracion <- rep(list(tasas_migracion_A,
                                      tasas_migracion_B),
                                 times = repeticiones)


secuencia_tasas_mortalidad <- rep(list(tasas_mortalidad_A,
                                       tasas_mortalidad_B),
                                  times = repeticiones)


resultado.ciclos_20_A = dinamica_metacomunidad_ciclos_todo(paisaje,matriz_abundancias,ecuacion_lety,parametros_ecuacion_lety,parches_con_dinamica_comunitaria,secuencia_tasas_migracion,secuencia_tasas_mortalidad,secuencia_tiempo_total,pasos_mig,pasos_mort,tiempo_final_int,metodo_integracion,array_redes,matriz_ris)
save(resultado.ciclos_20_A,file="./resultado.ciclos_20_A.Rdata")



#High temporal variability starting with rainy conditions
tiempo_total <- 10
repeticiones <- 5
secuencia_tiempo_total <- rep(c(tiempo_total,
                                tiempo_total),
                              times = repeticiones)
secuencia_tipo_ciclo <- rep(c("A","B"), times = repeticiones)

secuencia_tasas_migracion <- rep(list(tasas_migracion_A,
                                      tasas_migracion_B),
                                 times = repeticiones)


secuencia_tasas_mortalidad <- rep(list(tasas_mortalidad_A,
                                       tasas_mortalidad_B),
                                  times = repeticiones)


resultado.ciclos_10_A = dinamica_metacomunidad_ciclos_todo(paisaje,matriz_abundancias,ecuacion_lety,parametros_ecuacion_lety,parches_con_dinamica_comunitaria,secuencia_tasas_migracion,secuencia_tasas_mortalidad,secuencia_tiempo_total,pasos_mig,pasos_mort,tiempo_final_int,metodo_integracion,array_redes,matriz_ris)
save(resultado.ciclos_10_A,file="./resultado.ciclos_10_A_rep.Rdata")



#########################
#####DRY CONDITIONS
##########################

#Low temporal variability starting with dry conditions
tiempo_total <- 50
repeticiones <- 1
secuencia_tiempo_total <- rep(c(tiempo_total,
                                tiempo_total),
                              times = repeticiones)
secuencia_tipo_ciclo <- rep(c("B","A"), times = repeticiones)

secuencia_tasas_migracion <- rep(list(tasas_migracion_B,
                                      tasas_migracion_A),
                                 times = repeticiones)


secuencia_tasas_mortalidad <- rep(list(tasas_mortalidad_B,
                                       tasas_mortalidad_A),
                                  times = repeticiones)


resultado.ciclos_50_B = dinamica_metacomunidad_ciclos_todo(paisaje,matriz_abundancias,ecuacion_lety,parametros_ecuacion_lety,parches_con_dinamica_comunitaria,secuencia_tasas_migracion,secuencia_tasas_mortalidad,secuencia_tiempo_total,pasos_mig,pasos_mort,iteraciones_integracion,metodo_integracion,array_redes,matriz_ris)
save(resultado.ciclos_50_B,file="./resultado.ciclos_50_B.Rdata")


#Moderate temporal variability starting with dry conditions
tiempo_total <- 20
repeticiones <- 3
secuencia_tiempo_total <- rep(c(tiempo_total,
                                tiempo_total),
                              times = repeticiones)
secuencia_tipo_ciclo <- rep(c("B","A"), times = repeticiones)

secuencia_tasas_migracion <- rep(list(tasas_migracion_B,
                                      tasas_migracion_A),
                                 times = repeticiones)


secuencia_tasas_mortalidad <- rep(list(tasas_mortalidad_B,
                                       tasas_mortalidad_A),
                                  times = repeticiones)


resultado.ciclos_20_B = dinamica_metacomunidad_ciclos_todo(paisaje,matriz_abundancias,ecuacion_lety,parametros_ecuacion_lety,parches_con_dinamica_comunitaria,secuencia_tasas_migracion,secuencia_tasas_mortalidad,secuencia_tiempo_total,pasos_mig,pasos_mort,tiempo_final_int,metodo_integracion,array_redes,matriz_ris)
save(resultado.ciclos_20_B,file="./resultado.ciclos_20_B.Rdata")



#High temporal variability starting with dry conditions
tiempo_total <- 10
repeticiones <- 5
secuencia_tiempo_total <- rep(c(tiempo_total,
                                tiempo_total),
                              times = repeticiones)
secuencia_tipo_ciclo <- rep(c("B","A"), times = repeticiones)

secuencia_tasas_migracion <- rep(list(tasas_migracion_B,
                                      tasas_migracion_A),
                                 times = repeticiones)


secuencia_tasas_mortalidad <- rep(list(tasas_mortalidad_B,
                                       tasas_mortalidad_A),
                                  times = repeticiones)


resultado.ciclos_10_B = dinamica_metacomunidad_ciclos_todo(paisaje,matriz_abundancias,ecuacion_lety,parametros_ecuacion_lety,parches_con_dinamica_comunitaria,secuencia_tasas_migracion,secuencia_tasas_mortalidad,secuencia_tiempo_total,pasos_mig,pasos_mort,tiempo_final_int,metodo_integracion,array_redes,matriz_ris)
save(resultado.ciclos_10_B,file="./resultado.ciclos_10_B_rep.Rdata")






