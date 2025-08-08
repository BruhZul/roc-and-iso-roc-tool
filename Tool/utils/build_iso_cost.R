#This function works exactly like the build_iso_curves() function (see the
#build_iso_curves.R file). The only difference is that this function uses two values:
#  cFN (the cost of False Negatives)
#  cFP (the cost of False Positives)
#
#The function generates the coordinates for a list of iso-MC curves (one for each
#value in the cost_range parameter). MC is the misclassification cost. That is, 
#MC = FP*cFP + FN*cFN.
build_iso_norm_cost_curve = function(rho, cost_range, cFN, cFP){
  result = c()
  
  lam = cFN/(cFN+cFP)
  for(cost in cost_range){
    inter_y0 = (cost - lam*rho)/((1-lam)*(1-rho))
    if(inter_y0 > 0){
      first_x = inter_y0
      first_y = 0
    }else{
      first_x = 0
      first_y = 1 - (cost/(lam*rho))
    }
    
    inter_y1 = cost/((1-lam)*(1-rho))
    if(inter_y1 <= 1){
      second_x = inter_y1
      second_y = 1
    }else{
      second_x = 1
      second_y = ((1-lam)*(1-rho) + lam*rho - cost)/(lam*rho)
    }
    
    x_vector = c(first_x, second_x)
    y_vector = c(first_y, second_y)
    
    points = list(x_vector, y_vector)
    names(points) = c("x","y")
    result = append(result, list(points))
  }
  
  return(result)
}
