# # # Prioritisation # # #

# Library
library(raster)
library(sp)
library(sf)
library(ape)
library(phytools)
library(plyr)
library(reshape)
library(dplyr)

# Load in sample data - https://doi.org/10.6084/m9.figshare.25736703
load(file = './EDGE_Zone_prioritisation_data/EDGE_Zone_mammal_sample.RData')

# dist = species distribution data in 96.5 x 96.5 km2 Mollweide equal area grid
# mam.data = 100 extinction-risk weighted phylogenies from Gumbs et al. 2023
# EDGE.data = median EDGE scores for mammals from Gumbs et al. 2023 
# mam.edge.2 = distribution of 1000 EDGE scores from Gumbs et al. 2023


#############################################################################

# # # Prioritisation on median EDGE scores # # #

#############################################################################


# Calculate median ePD loss from phylogenies
pd.tot<-list()
for(t in 1:100){
  pd <- sum(c(sum(drop.tip(mam.data[[4]][[t]],mam.data[[4]][[t]]$tip.label[which(!mam.data[[4]][[t]]$tip.label %in% mam.data[[5]])])$edge.length)))
  pd.tot[[length(pd.tot)+1]] <- pd
}
med.pd <- median(unlist(pd.tot)) 
pd.25 <- med.pd * 0.25 # 25% threshold


# Merge the EDGE scores with gridded data
cont <- merge(dist, EDGE.data, by = 'Species')

# Rasterise
rast <- cast(melt(cont[,c('X', 'Y', 'EDGE.score')], id=c("X","Y")), fun.aggregate=sum)

# Remove top cell
top.cells <- c()
spp.to.drop <- c()
perc<-c()
pdcap <- 0
top.cells <- rbind(top.cells, rast[which(rast$EDGE == max(rast$EDGE)),])
cont$Species <- as.character(cont$Species)
spp.to.drop <- c(spp.to.drop, cont$Species[cont$X == top.cells$X[NROW(top.cells)] & cont$Y == top.cells$Y[NROW(top.cells)]])
cont <- cont[!cont$Species %in% spp.to.drop,]

# Calculate PD percentage captured
pd.prop <- list()
for(t in 1:100){
  pd.i <- sum(c(sum(drop.tip(mam.data[[4]][[t]],mam.data[[4]][[t]]$tip.label[which(!mam.data[[4]][[t]]$tip.label %in% spp.to.drop)])$edge.length)))
  pd.prop[[length(pd.prop)+1]] <- pd.i
}
pdcap <- median(unlist(pd.prop)) 
first <- pdcap/med.pd; print(first)

while (pdcap < pd.25){ # Iteratively remove grid cells until the threshold is crossed
  
  rast <- cast(melt(cont[,c('X', 'Y', 'EDGE.score')], id=c("X","Y")), fun.aggregate=sum)
  
  # get the new top grid cell
  top.cells <- rbind(top.cells, rast[which(rast$EDGE == max(rast$EDGE)),])
  spp.to.drop <- c(spp.to.drop, cont$Species[cont$X == top.cells$X[NROW(top.cells)] & cont$Y == top.cells$Y[NROW(top.cells)]])
  # remove species in new top grid cell from dataset
  cont <- cont[!cont$Species %in% spp.to.drop,]
  
  # get proportion of total ePD loss captured
  pd.prop<-list()
  for(t in 1:100){
    pd.i <- sum(c(sum(drop.tip(mam.data[[4]][[t]],mam.data[[4]][[t]]$tip.label[which(!mam.data[[4]][[t]]$tip.label %in% spp.to.drop)])$edge.length)))
    pd.prop[[length(pd.prop)+1]] <- pd.i
  }
  pdcap <- median(unlist(pd.prop)) 
  perc[[length(perc)+1]] <- pdcap/med.pd
  
  print(paste(signif(pdcap/med.pd, 3 ),
              "captured for median EDGE Score set"))
  
  #save(top.cells, spp.to.drop, file = "median_EDGE_prioritisation.RData")
  print(unlist(perc))
}


top.cells$perc <- c(first, unlist(perc)) # data frame with all priority cells

save(top.cells, spp.to.drop, file = "./EDGE_Zone_prioritisation_data/median_EDGE_prioritisation.RData")


#############################################################################

# # # Prioritisation on the distribution of EDGE scores # # #

#############################################################################

# Calculate median ePD loss
pd.tot<-list()
for(t in 1:100){
  pd <- sum(c(sum(drop.tip(mam.data[[4]][[t]],mam.data[[4]][[t]]$tip.label[which(!mam.data[[4]][[t]]$tip.label %in% mam.data[[5]])])$edge.length)))
  pd.tot[[length(pd.tot)+1]] <- pd
}
med.pd <- median(unlist(pd.tot))
pd.25 <- med.pd * 0.25 # Threshold


# Selecting over 1000 possible values
for(i in 1:1000){
  print(i)
  tetr <- rbind(mam.edge.2[[i]]) # pulling set i of 1000 possible sets of EDGE scores for each clade
  cont <- merge(tetr, dist, by = 'Species') # join with distribution data
  contI <- cont # First iteration 
  
  
  # Initial raster
  rast <- cast(melt(contI[,c('X', 'Y', 'EDGE')], id=c("X","Y")), fun.aggregate=sum)
  
  # Remove top cell
  top.cells <- c()
  spp.to.drop <- c()
  perc<-c()
  pdcap <- 0
  top.cells <- rbind(top.cells, rast[which(rast$EDGE == max(rast$EDGE)),])
  contI$Species <- as.character(contI$Species)
  spp.to.drop <- c(spp.to.drop, contI$Species[contI$X == top.cells$X[NROW(top.cells)] & contI$Y == top.cells$Y[NROW(top.cells)]])
  contI <- contI[!contI$Species %in% spp.to.drop,]
  
  # Calculate first iteration PD percentage captured
  pd.prop <- list()
  for(t in 1:100){
    pd.i <- sum(c(sum(drop.tip(mam.data[[4]][[t]],mam.data[[4]][[t]]$tip.label[which(!mam.data[[4]][[t]]$tip.label %in% spp.to.drop)])$edge.length)))
    pd.prop[[length(pd.prop)+1]] <- pd.i
  }
  pdcap <- median(unlist(pd.prop)) #50451.66 MY
  first <- pdcap/med.pd; print(first)
  
  while (pdcap < pd.25){
    
    rast <- cast(melt(contI[,c('X', 'Y', 'EDGE')], id=c("X","Y")), fun.aggregate=sum)
    
    # New top grid cell
    nr <- nrow(rast[which(rast$EDGE == max(rast$EDGE)),])
    top.cells <- rbind(top.cells, rast[which(rast$EDGE == max(rast$EDGE)),])
    spp.to.drop <- c(spp.to.drop, contI$Species[contI$X == top.cells$X[NROW(top.cells)] & contI$Y == top.cells$Y[NROW(top.cells)]])
    # Remove species in new top grid cell from dataset
    contI <- contI[!contI$Species %in% spp.to.drop,]
    
    # Get proportion of total ePD loss captured
    pd.prop<-list()
    for(t in 1:100){
      pd.i <- sum(c(sum(drop.tip(mam.data[[4]][[t]],mam.data[[4]][[t]]$tip.label[which(!mam.data[[4]][[t]]$tip.label %in% spp.to.drop)])$edge.length)))
      pd.prop[[length(pd.prop)+1]] <- pd.i
    }
    pdcap <- median(unlist(pd.prop)) 
    
    for(n in 1:nr){
      perc[[length(perc)+1]] <- pdcap/med.pd
    }
    
    print(paste(signif(pdcap/med.pd, 3 ),
                "captured for EDGE Score set",i, sep= " "))
    
    save(top.cells, spp.to.drop, file = paste("./EDGE_Zone_prioritisation_data/EDGE_dist_iteration_",i,".RData",sep = ""))
    print(unlist(perc))
  }
  
  
  
  top.cells$perc <- c(first, unlist(perc))
  
  save(top.cells, spp.to.drop, file = paste("./EDGE_Zone_prioritisation_data/EDGE_dist_iteration_",i,".RData",sep = ""))
  
}


