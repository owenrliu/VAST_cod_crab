#### overlap functions to accompany Carroll et al. (2019) 'A review of methods for quantifying predator-prey overlap,' Global Ecology and Biogeography 
## please see accompanying paper for detailed descriptions of metrics and their ecological interpretations


## area overlap
## for binary data
## measures proportion of an area where two species co-occur
area_overlapfn <- function(prey, pred, area){
  total_area <- sum(area, na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/total_area
}

## range overlap
## for binary data
## measures the proportion of one species range where the other co-occurs
range_overlapfn<-function(prey, pred, area){
  area_prey <- sum(area[prey > 0], na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/area_prey
}

## local index of collocation
## estimates correlation of predator and prey densities
loc_collocfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(p_prey*p_pred, na.rm = T)/(sqrt(sum(p_prey^2, na.rm = T)*sum(p_pred^2, na.rm = T)))
}

## asymmetrical alpha
## measures pressure of predator on prey relative to underlying prey density
asymmalpha_overlapfn <-function(prey, pred){
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(p_pred*p_prey, na.rm = T)/sum(p_prey^2, na.rm = T)
}

## biomass-weighted overlap (scaled to max)
## measures amount of predator biomass interacting with prey relative to underlying prey biomass
biomass_overlapfn <- function(prey, pred) {
  sum((prey/max(prey, na.rm = T)) * (pred/max(pred, na.rm = T)), na.rm = T)/sum(prey/max(prey, na.rm = T), na.rm = T)
}

## Hurlbert's overlap
## measures interspecific encounter rate between predator and prey
hurlbert_overlapfn <- function(prey, pred, area) {
  area_occupied <- sum(area[pred > 0 | prey > 0], na.rm = T)
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum((p_pred*p_prey)/(area/area_occupied), na.rm = T)
}

## Schoener's D
## density or probability of occurrence data
## measures how equally predator and prey share available resources
schoeners_overlapfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  1 - 0.5 * (sum(abs(p_prey-p_pred), na.rm = T))
}

## Bhattacharyya's coefficient
## density or probability of occurrence data
## measures whether two species use space independently
bhatta_coeffn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(sqrt(p_prey*p_pred), na.rm = T)
}

## global index of collocation
## measures geographic distinctness by comparing centres of gravity and dispersion of sampled individuals
glob_collocfn <- function(prey_x, prey_y, prey, pred_x, pred_y, pred){
  prey_cgx <- sum(prey_x*prey, na.rm = T)/sum(prey, na.rm = T)
  prey_cgy <- sum(prey_y*prey, na.rm = T)/sum(prey, na.rm = T)
  prey_ix <- prey_x - prey_cgx
  prey_iy <- prey_y - prey_cgy
  prey_i <- sqrt(prey_ix^2 + prey_iy^2)
  prey_inert <- sum(prey * (prey_i^2), na.rm = T)/sum(prey, na.rm = T)
  pred_cgx <- sum(pred_x*pred, na.rm = T)/sum(pred, na.rm = T)
  pred_cgy <- sum(pred_y*pred, na.rm = T)/sum(pred, na.rm = T)
  pred_ix <- pred_x - pred_cgx
  pred_iy <- pred_y - pred_cgy
  pred_i <- sqrt(pred_ix^2 + pred_iy^2)
  pred_inert <- sum(pred * (pred_i^2), na.rm = T)/sum(pred, na.rm = T)
  GIC <- (((prey_cgx - pred_cgx)^2+(prey_cgy - pred_cgy)^2)/ (((prey_cgx-pred_cgx)^2+(prey_cgy-pred_cgy)^2)+prey_inert + pred_inert))
  if(!is.na(GIC))
    GIC <- 1-GIC
  else GIC <- 1
  GIC
}

## AB ratio
## measures predator production that can be attributed to spatial overlap with prey
AB_overlapfn <- function(prey, pred) { 
  mean((pred - mean(pred, na.rm = T)) * (prey - mean(prey, na.rm = T)), na.rm = T)/(mean(pred, na.rm = T) * mean(prey, na.rm = T)) 
}




