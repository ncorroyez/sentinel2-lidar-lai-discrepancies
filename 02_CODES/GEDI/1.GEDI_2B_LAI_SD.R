# Sylvie Durrieu -2025-10-30

#install.packages('rGEDI', repos = c('https://carlos-alberto-silva.r-universe.dev', 'https://cloud.r-project.org'))


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")

library("rhdf5")
library(rGEDI)
library("data.table")  # to read data faster than read.cvs
library("XML")
library("sf")
library("stringr")
#library("leaflet")

library(bit64)
library(lubridate)



#### fonctions 

## recherche du max d'un profil de pad
# Hmax_padz<-function(padz , dz=5)
# {
#   Hmax=0
#   if (length(padz)>0) 
#   { padz= padz[length(padz):1]
#   for (i in (1:length(padz)))
#   {
#     index=i
#     if (padz[i]>0) {break}
#   }     
#   Hmax= (length(padz)-index)*dz
#   } 
#   #print("youpi")
#   return(Hmax)
# }

infos_beam <- function(fileh5, stru2, num_beam) 
  {
    fileH5=fileh5
    nom_beam= names(stru2)[num_beam]
    
    infos_sub=NULL
  
    if (length(stru2[[nom_beam]]) > 1) {
        
      shot_number= H5Dread(eval(parse(text=paste0("fileH5&'/",nom_beam,"/geolocation/shot_number'"))), bit64conversion='bit64')
      elev_bin0= eval(parse(text=paste0("fileH5$'/",nom_beam,"/geolocation/elevation_bin0'")))
      elev_dem= eval(parse(text=paste0("fileH5$'/",nom_beam,"/geolocation/digital_elevation_model'")))
      diff_alt0 = elev_bin0- elev_dem 
  
      subset= which( eval(parse(text=paste0("fileH5$'/",nom_beam,"/algorithmrun_flag'"))) == 1         
                 & diff_alt0 <100 & diff_alt0 >0 
                 & eval(parse(text=paste0("fileH5$'/",nom_beam,"/geolocation/degrade_flag'")))==0)
  
      print(length(subset))
      if (length(subset)>0)   # >0 pb avec hmax si pas as.data.frame
        {
        infos_sub$shot_number <- as.character(shot_number[subset])  
        infos_sub$delta_time = eval(parse(text=paste0("fileH5$'/",nom_beam,"/geolocation/delta_time'")))[subset]
        infos_sub$lat =  eval(parse(text=paste0("fileH5$'/",nom_beam,"/geolocation/lat_lowestmode'")))[subset]      # [subset[1]:subset[length(subset)]] erreur[subset[1:length(subset)]]
        infos_sub$long = eval(parse(text=paste0("fileH5$'/",nom_beam,"/geolocation/lon_lowestmode'")))[subset]
        infos_sub$RH100 = eval(parse(text=paste0("fileH5$'/",nom_beam,"/rh100'")))[subset]/100
        infos_sub$pai= eval(parse(text=paste0("fileH5$'/",nom_beam,"/pai'")))[subset]
        infos_sub$cover= eval(parse(text=paste0("fileH5$'/",nom_beam,"/cover'")))[subset]
        #infos_sub$max_h_cutof= eval(parse(text=paste0("fileH5$'/",nom_beam,"/ancillary/maxheight_cuttoff'")))[subset[1]:subset[length(subset)]]
        #
        # for (i in 1:length(subset))
        #   {
        #   padz= eval(parse(text=paste0("fileH5$'/",nom_beam,"/pavd_z'")))[1:HT, i]
        #   infos_sub$Hmax[i] = Hmax_padz(padz, dz) 
        #   }
        
        #   if (length(subset) > 1)
        #       {
        #     infos_sub$Hmax = apply(eval(parse(text=paste0("fileH5$'/",nom_beam,"/pavd_z'")))[1:HT, subset] ,2,Hmax_padz, dz=dz)
        #       } else {
        #     infos_sub$Hmax = Hmax_padz(eval(parse(text=paste0("fileH5$'/",nom_beam,"/pavd_z'")))[1:HT, 1],dz)  
        #       }
        # ### pb de dimension dans Hmax si subset de taille 1 alors consid?r? comme vecteur
        ##"pb aussi de calcul pour des beam avec RH100 = 0 ou PAI non d?fini => mettre un seuil de RH100 pour le calcul, NULL 
        
      } 
      } else {print("subset vide") }
       
  return(infos_sub)
  
  }



infos_allbeams <- function(fileh5,stru2) 
{
  infos=NULL                  
  for (i in 1:(length(names(stru2))-1))  # -1 car dernier nom est "metadata" et non un beam
  { 
      infos_i <- infos_beam(fileh5, stru2,i)
      
    if (!is.null(infos_i)) {
      
    infos$shot_number = c(infos$shot_number, infos_i$shot_number)
    infos$delta_time = c(infos$delta_time, infos_i$delta_time)
    infos$lat=c(infos$lat, infos_i$lat)
    infos$long=c(infos$long, infos_i$long)
    infos$RH100= c(infos$RH100, infos_i$RH100)
    infos$pai= c(infos$pai, infos_i$pai)
    infos$cover= c(infos$cover, infos_i$cover)
    infos$Hmax=c(infos$Hmax, infos_i$Hmax)
    #infos$max_h_cutof = c(infos$max_h_cutof, infos_i$max_h_cutof)
    }
   
  }
  return(infos)
}



###### liste des fichiers 

file_name = file.choose()    # choisir un fichier pour recuperer le chemin et nom du fichier

nom_rep=dirname(path=file_name)  # recupere le nom du dossier ou se trouvent les fichiers
#liste_fic= dir(nom_rep, pattern = "*.h5")

liste_fic=paste0(nom_rep,"/",dir(nom_rep, pattern = "*.h5"))  # liste les fichiers a parcourir

# fileH5 <- H5Fopen(liste_fic[1])  # pour test
# 
# #h5dump(fileH5, recursive=TRUE, all=FALSE, load= FALSE)   # permet d'avoir la tt la hi?rarchie 
# stru_H5 =h5dump(fileH5, recursive=2, all=FALSE, load= FALSE)  
# stru_H5_1 =h5dump(fileH5, recursive=1, all=FALSE, load= FALSE) # permet d'avoir la  hi?rarchie jusqu'au 2?me niveau
# 
# stru2 <- h5dump(fileH5, recursive=2, all=FALSE, load= FALSE)  


infos_all=NULL
for (j in 1:length(liste_fic))
  #for (j in 15:15)  # test sur les  premiers fichiers telecharges
{
  
  file_name <- basename(liste_fic[j])
  orbit <- substr(unlist(stringr::str_split(file_name, "_"))[4], 2, 6)
  
  fileH5 <- tryCatch({
    # Tentative d'ouverture du fichier HDF5
    H5Fopen(liste_fic[j])
  }, error = function(e) {
    # Si une erreur se produit, afficher un message et passer au fichier suivant
    cat("Erreur lors de l'ouverture du fichier", liste_fic[j], ":", e$message, "\n")
    return(NULL)  # Retourne NULL pour signaler une erreur
  })
  
  # si l'ouverture a réussi
  if (!is.null(fileH5)) {
    print(paste("fic:",j, "open"))
    
    struc_H5 <- h5dump(fileH5, recursive=2, all=FALSE, load= FALSE) 
    infosfile <-infos_allbeams(fileH5,struc_H5)
    infosfile$orbit <-rep(orbit, times = length(infosfile$shot_number))
    
    infos_all$orbit <- c(infos_all$orbit, infosfile$orbit)
    infos_all$shot_number <- c(infos_all$shot_number, infosfile$shot_number)
    infos_all$delta_time <- c(infos_all$delta_time, infosfile$delta_time)
    infos_all$lat=c(infos_all$lat, infosfile$lat)
    infos_all$long=c(infos_all$long, infosfile$long)
    infos_all$RH100=c(infos_all$RH100, infosfile$RH100)
    infos_all$pai=c(infos_all$pai, infosfile$pai)
    infos_all$cover=c(infos_all$cover, infosfile$cover)
    infos_all$Hmax=c(infos_all$Hmax, infosfile$Hmax)
    #infos_all$max_h_cutof= c(infos_all$max_h_cutof, infosfile$max_h_cutof)
    h5closeAll()
    print(paste("fic:",j,"closed"))
    
  } else {
    # Si le fichier n'a pas pu être ouvert, passer au suivant
    cat("Passage au fichier suivant...\n")
  }
  
}

h5closeAll()

df <- as.data.frame(infos_all)

# master_time_epoch
# Number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z) and the GEDI epoch (2018-01-01T00:00:00Z). Add this value to delta_time parameters to compute full GPS seconds (relative to the GPS epoch) for each data point.

#gps_epoch <- as.POSIXct("1980-01-06 00:00:00", tz = "UTC")
gedi_epoch <- as.POSIXct("2018-01-01 00:00:00", tz="UTC")
#master_time_epoch <- (gedi_epoch - gps_epoch)*24*60*60

# Conversion brute
time <- gedi_epoch + df$delta_time

df$year <- format(time, "%Y") 
df$month <- format(time, "%m") 
df$day <- format(time, "%d") 
df$heure <- format(time, "%H:%M:%S") 

# Modifier le systeme de projection et sauver le dta frame 