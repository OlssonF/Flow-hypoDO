#' @param wtr column or vector with a single water temperature profile
#' @param depths column or vector with the associated depths
#' @param inflowT column or object containing the inflow temp
#' @param use_density logical. Should the nominal intrusion depth be based on water density? (not temperature) 

intrusion.depth <- function(wtr, depths, inflowT, use_density = T) {
  
  inflowT <- unique(inflowT)
  
  if (length(inflowT) != 1) {
    stop('Only 1 inflow temperature per day')
  } else {
   
    if (is.na(inflowT) | sum(is.na(wtr)) > 3) {
      return(NA)
    } else {
      
      
      insituT <- tibble(depth = depths,
                        temp = wtr) |> 
        mutate(depth = as.numeric(str_remove(depth, 'wtr_')))
      
      if (use_density) {
        inflowT <- water.density(inflowT)
        
        insituT$temp <- water.density(insituT$temp)
      }
      
      inflowD <- approx(x = insituT$temp, 
                        y = insituT$depth,
                        # want the depth f the lake equivalent to the temp of the inflow
                        xout = inflowT,
                        # when the inflowT is outside of the range in the lake,   
                        # assign max and min values
                        yleft = 6, yright = 1)$y
      
      return(inflowD)
    } 
  }
  
  
  
}
