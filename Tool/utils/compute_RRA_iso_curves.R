library(iRRA)

# Used for approximation
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x+1)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#used for comparisons
less_e_than = function(a, b, precision=1e-8){
  return(a < b || abs(a-b) <= precision)
}

greater_e_than = function(a, b, precision=1e-8){
  return(a > b || abs(a-b) <= precision)
}

is_equal = function(a, b, precision=1e-8){
  return(abs(a-b) <= precision)
}

compute_iso_value_from_rra = function(metric, rra, rho, precision=0.001){
  result = c()
  if(metric == "MCC"){
    result = iso_mcc_from_rra(rra, rho, precision)
  }
  else if(metric == "BA"){
    result = iso_ba_from_rra(rra, rho)
  }
  else if(metric == "Gmean"){
    result = iso_gmean_from_rra(rra, rho, precision)
  }
  else if(metric == "GM"){
    result = iso_gm_from_rra(rra, rho, precision)
  }
  else if(metric == "D2H"){
    result = iso_D2H_from_rra(rra, rho, precision) 
  }
  else if(metric == "precision"){
    result = iso_precision_from_rra(rra, rho) 
  }
  else if(metric == "TPR"){
    result = iso_tpr_from_rra(rra, rho) 
  }
  else if(metric == "F-score"){
    result = iso_fscore_from_rra(rra, rho) 
  }
  else if(metric == "NPV"){
    result = iso_npv_from_rra(rra, rho)
  }
  else if(metric == "TNR"){
    result = iso_tnr_from_rra(rra, rho) 
  }
  else if(metric == "FPR"){
    result = iso_fpr_from_rra(rra, rho) 
  }
  else if(metric == "NM"){
    result = iso_nm_from_rra(rra, rho) 
  }
  else if(metric == "J"){
    result = iso_j_from_rra(rra, rho) 
  }
  else if(metric == "Markedness"){
    result = iso_markedness_from_rra(rra, rho, precision)
  }
  else{
    stop("Value for parameter 'metric' is not valid")
  }
  
  return(result)
}

iso_mcc_from_rra = function(rra, rho, precision=0.001){
  if(rra == 0){
    value = 0 # This value gives the iso-curve that intersecates the point (rho,rho)
  }else if(rra==1){
    value = 1
  }else{
    found = FALSE
    upper_value = 1
    lower_value = 0
    while(!found){
      mcc = (upper_value + lower_value)/2
      curve = iso_MCC(mcc, rho)
      curve$x = c(0,curve$x,1)
      curve$y = c(0,curve$y,1)
      val = suppressWarnings(rra(curve$x, curve$y, 1, (1-rho)/rho, fallout = TRUE, recall = TRUE, plot = FALSE, print=FALSE)$rra)
      
      if(val == rra){
        value = mcc
        found = TRUE
      }else if(val < rra){
        lower_value = mcc
      }else{
        upper_value = mcc
      }
      if(upper_value-lower_value < precision){
        value = (upper_value + lower_value)/2
        found = TRUE
      }
    }
  }
  
  return(round(value, digits = decimalplaces(precision)))
}

iso_ba_from_rra = function(rra, rho){
  if(rra == 0){
    value = 0.5 # under 0.5 RRA is always 0. this returns the diagonal, as it is the highest value where the RRA of the iso-BA is 0
  }else{
    #case 1
    value = (1+suppressWarnings(sqrt(rra*(2*rho*(1-rho)))))/2
    if(!is.nan(value) && greater_e_than(rho, 2*value - 1) && less_e_than(rho, 2 - 2*value)){
      return(value)
    }
    
    #case 2
    value = (2*rra - 2*rho*rra + rho + 2)/4
    if(less_e_than(rho, 2*value - 1) && less_e_than(rho, 2 - 2*value)){
      return(value)
    }
    
    #case 3
    value = (2*rho*rra - rho + 3)/4
    if(greater_e_than(rho, 2*value - 1) && greater_e_than(rho, 2 - 2*value)){
      return(value)
    }
    
    #case 4
    value = 1 - (suppressWarnings(sqrt(2*rho*(1-rho)*(1-rra))))/2
    if(less_e_than(rho, 2*value - 1) && greater_e_than(rho, 2 - 2*value)){
      return(value)
    }
    
    #should never go here
    stop("Valid value for iso-BA not found. Make sure that the RRA and rho values are correct.")
  }
  
  return(value)
}

iso_gmean_from_rra = function(rra, rho, precision=0.001){
  if(rra == 0){
    value = sqrt((1-rho)*rho) # This value gives the iso-curve that intersecates the point (rho,rho)
  }else if(rra==1){
    value = 1
  }else{
    gmean_vals = seq(0, 1-precision, precision)
    for(i in 1:(length(gmean_vals)-1)){
      gmean_1 = gmean_vals[i]
      gmean_2 = gmean_vals[i+1]
      h_1 = gmean_1^2
      h_2 = gmean_2^2
      
      gmean_function_1 = function(x){
        h_1/(1-x)
      }
      
      gmean_function_2 = function(x){
        h_2/(1-x)
      }
      
      lower_limit = function(h){
        return(max(0, (rho-h)/rho))
      }
      
      if(h_1/(1-rho) > 1){
        val_1 = (integrate(gmean_function_1, lower = lower_limit(h_1), upper = 1-h_1)$value - (1-h_1-lower_limit(h_1))*rho + (rho-1+h_1)*(1-rho))/(rho*(1-rho))
      }else{
        val_1 = (integrate(gmean_function_1, lower = lower_limit(h_1), upper = rho)$value - (rho-lower_limit(h_1))*rho)/(rho*(1-rho))
      }
      
      if(h_2/(1-rho) > 1){
        val_2 = (integrate(gmean_function_2, lower = lower_limit(h_2), upper = 1-h_2)$value - (1-h_2-lower_limit(h_2))*rho + (rho-1+h_2)*(1-rho))/(rho*(1-rho))
      }else{
        val_2 = (integrate(gmean_function_2, lower = lower_limit(h_2), upper = rho)$value - (rho-lower_limit(h_2))*rho)/(rho*(1-rho))
      }
      
      if(val_1 == rra){
        value = gmean_1
        break
      }
      if(val_1 < rra && val_2 > rra){
        value = ifelse(abs(rra-val_1) > abs(rra-val_2), gmean_2, gmean_1)
        break
      }
    }
  }
  
  return(value)
}

iso_gm_from_rra = function(rra, rho, precision=0.001){
  if(rra == 0){
    value = 2*rho*(1-rho) # This value gives the iso-curve that intersecates the point (rho,rho)
  }else if(rra==1){
    value = 1
  }else{
    gm_vals = seq(0, 1-precision, precision)
    for(i in 1:(length(gm_vals)-1)){
      gm_1 = gm_vals[i]
      gm_2 = gm_vals[i+1]
      h_1 = gm_1/2
      h_2 = gm_2/2
      
      gm_function_1 = function(x){
        (h_1*(1-x))/(1-x-h_1)
      }
      
      gm_function_2 = function(x){
        (h_2*(1-x))/(1-x-h_2)
      }
      
      lower_limit = function(h){
        if((h-rho+rho*h)/(h-rho) > 1-h){
          return(0)
        }
        return(max(0, (h-rho+rho*h)/(h-rho)))
      }
      
      if((h_1*(1-rho))/(1-rho-h_1) > 1 || (h_1*(1-rho))/(1-rho-h_1) < 0){
        val_1 = (integrate(gm_function_1, lower = lower_limit(h_1), upper = (2*h_1-1)/(h_1-1))$value - ((2*h_1-1)/(h_1-1)-lower_limit(h_1))*rho + ((rho-(2*h_1-1)/(h_1-1)))*(1-rho))/(rho*(1-rho))
      }else{
        val_1 = (integrate(gm_function_1, lower = lower_limit(h_1), upper = rho)$value - (rho-lower_limit(h_1))*rho)/(rho*(1-rho))
      }
      
      if((h_2*(1-rho))/(1-rho-h_2) > 1 || (h_2*(1-rho))/(1-rho-h_2) < 0){
        val_2 = (integrate(gm_function_2, lower = lower_limit(h_2), upper = (2*h_2-1)/(h_2-1))$value - ((2*h_2-1)/(h_2-1)-lower_limit(h_2))*rho + ((rho-(2*h_2-1)/(h_2-1)))*(1-rho))/(rho*(1-rho))
      }else{
        val_2 = (integrate(gm_function_2, lower = lower_limit(h_2), upper = rho)$value - (rho-lower_limit(h_2))*rho)/(rho*(1-rho))
      }
      
      if(val_1 == rra){
        value = gm_1
        break
      }
      if(val_1 < rra && val_2 > rra){
        value = ifelse(abs(rra-val_1) > abs(rra-val_2), gm_2, gm_1)
        break
      }
    }
  }
  
  return(value)
}

iso_D2H_from_rra = function(rra, rho, precision=0.001){
  if(rra == 0){
    value = sqrt(((1-rho)^2 + (rho^2))/2)
  }else if(rra == 1){
    value = 0
  }else{
    d2h_vals = seq(precision, 1-precision, precision)
    for(i in 1:(length(d2h_vals)-1)){
      d2h_1 = d2h_vals[i]
      d2h_2 = d2h_vals[i+1]
      r_1 = d2h_1*sqrt(2)
      r_2 = d2h_2*sqrt(2)
      
      d2h_function_1 = function(x){
        1-sqrt(r_1^2 - x^2)
      }
      
      d2h_function_2 = function(x){
        1-sqrt(r_2^2 - x^2)
      }
      
      lower_limit = function(r){
        if(r <= 1-rho){
          return(0)
        }
        return(sqrt(r^2 - (1-rho)^2))
      }
      
      if(rho >= r_1){
        val_1 = (integrate(d2h_function_1, lower = lower_limit(r_1), upper = r_1)$value - (r_1-lower_limit(r_1))*rho + (rho-r_1)*(1-rho))/(rho*(1-rho))
      }else{
        val_1 = (integrate(d2h_function_1, lower = lower_limit(r_1), upper = rho)$value - (rho-lower_limit(r_1))*rho)/(rho*(1-rho))
      }
      
      if(rho >= r_2){
        val_2 = (integrate(d2h_function_2, lower = lower_limit(r_2), upper = r_2)$value - (r_2-lower_limit(r_2))*rho + (rho-r_2)*(1-rho))/(rho*(1-rho))
      }else{
        val_2 = (integrate(d2h_function_2, lower = lower_limit(r_2), upper = rho)$value - (rho-lower_limit(r_2))*rho)/(rho*(1-rho))
      }
      
      if(val_1 == rra){
        value = d2h_1
        break
      }
      if(val_1 > rra && val_2 < rra){
        value = ifelse(abs(rra-val_1) > abs(rra-val_2), d2h_2, d2h_1)
        break
      }
    }
  }
  
  return(value)
}

iso_precision_from_rra = function(rra, rho){
  k = (1-rho)/rho
  if(rra == 0){
    value = 1/(k+1)
  }else if(rra == 1){
    value = 1
  }else{
    #caso 1
    value = (1+rho)/(2*(1-rho)*(1-rra) + 1 + rho)
    if(greater_e_than(value, 1/(2-rho))){ 
      return(value)
    }
    
    #case 2
    A = rra*(2*(1-rho)-2*k)-(k+1)
    B = rra*(2*k - 2*(1-rho)) + 2
    C = -rho
    value = (-B-sqrt(B^2 - 4*A*C))/(2*A)
    if(less_e_than(value, 1/(2-rho))){
      return(value)
    }
    
    #should never go here
    stop("Valid value for iso-Precision not found. Make sure that the RRA and rho values are correct.")
  }
  
  return(value)
}

iso_tnr_from_rra = function(rra, rho){
  if(rra==0){
    value = 1-rho
  }else{
    value = rho*(rra-1) + 1
  }
  
  return(value)
}

iso_tpr_from_rra = function(rra, rho){
  if(rra==0){
    value = rho
  }else{
    value = rra*(1-rho)+rho
  }
  
  return(value)
}

iso_fscore_from_rra = function(rra, rho){
  k = (1-rho)/rho
  if(rra==0){
    value = rho
  }else if(rra==1){
    value = 1
  }else{
    #case 1
    A = 4 + rra*(2+2*(rho^2)-4*rho)
    B = -(8*rho + 4*rra*((1-rho)^2))
    C = 4*(rho^2)
    value = (-B+suppressWarnings(sqrt(B^2 - 4*A*C)))/(2*A)
    if(!is.nan(value) && less_e_than(value,2/(3-rho)) && less_e_than(value,(2*rho)/(1+rho))){
      return(value)
    }
    
    #case 2 No
    #case 3
    value = (4*(rra*(1-rho) + rho))/(3+rho+2*rra*(1-rho))
    if(less_e_than(value,2/(3-rho)) && greater_e_than(value,(2*rho)/(1+rho))){
      return(value)
    }
    
    #case 4
    A = 2*(2*rho - rho^2 - 3 + rra*((1-rho)^2))
    B = 4*(((1-rho)^2)*(1-rra) + 2)
    value = (-B + sqrt(B^2 + 16*A))/(2*A)
    if(!is.nan(value) && greater_e_than(value,2/(3-rho)) && greater_e_than(value,(2*rho)/(1+rho))){
      return(value)
    }
    
    #should never go here
    stop("Valid value for iso-Precision not found. Make sure that the RRA and rho values are correct.")
  }
  
  return(value)
}

iso_npv_from_rra = function(rra, rho){
  k = (1-rho)/rho
  if(rra==0){
    value = (rho + k -1)/k
  }else if(rra==1){
    value = 1
  }else{
    #case 1
    A = k^2 + 2*rra*((1-rho)^2)
    B = 2*((((rho-1)^3)/(rho^2)) - rra*((1-rho)^2))
    C = ((rho-1)^4)/(rho^2)
    value = (-B+suppressWarnings(sqrt(B^2 - 4*A*C)))/(2*A)
    if(!is.nan(value) && less_e_than(value, (-k)/(rho-1-k))){
      return(value)
    }
    
    #case 2
    value = (3*rho - rho^2 - 2)/(2*rho*rra*(1-rho) + rho^2 + rho - 2)
    if(greater_e_than(value, (-k)/(rho-1-k))){
      return(value)
    }
    
    #should never go here
    stop("Valid value for iso-NPV not found. Make sure that the RRA and rho values are correct.") 
  }
  
  return(value)
}

iso_fpr_from_rra = function(rra, rho){
  if(rra == 0){
    value = rho
  }else{
    value = rho*(1-rra)
  }
  
  return(value)
}

iso_nm_from_rra = function(rra, rho){
  k = (1-rho)/rho
  if(rra==0){
    value = 1-rho
  }else if(rra==1){
    value = 1
  }else{
    #caso 1
    A = 2*(2+(rho^2)*rra)
    B = -4*(2-2*rho+(rho^2)*rra)
    C = 4*((1-rho)^2)
    value = (-B+suppressWarnings(sqrt(B^2 - 4*A*C)))/(2*A)
    if(!is.nan(value) && less_e_than(value, (2*((1-rho)^2))/(rho^2-3*rho+2))){
      return(value)
    }
    
    #case 2
    value = (4*((1-rho)^2 + rho*rra*(1-rho)))/((rho-1)*(rho-4) + 2*rho*rra*(1-rho))
    if(less_e_than(value, (-2*(1-rho))/(rho^2+rho-2)) && 
       greater_e_than(value, (2*((1-rho)^2))/(rho^2-3*rho+2))){
      return(value)
    }
    
    #case 3
    A = 2*(rho^2-rho-2*k)+2*rho*rra*(1-rho)
    B = 4*(rho-rho^2+2*k)-4*rho*rra*(1-rho)
    C = -4*k
    value = (-B + suppressWarnings(sqrt(B^2 - 4*A*C)))/(2*A)
    if(!is.nan(value) && greater_e_than(value, (-2*(1-rho))/(rho^2+rho-2))){
      return(value)
    }
    
    #should never go here
    stop("Valid value for iso-NM not found. Make sure that the RRA and rho values are correct.") 
  }
  
  return(value)
}

iso_j_from_rra = function(rra, rho){
  if(rra==0){
    value = 0
  }else{
    #case 1
    value = sqrt(2*rho*(1-rho)*rra)
    if(greater_e_than(rho,value) && less_e_than(rho,1-value)){
      return(value)
    }
    
    #case 2
    value = (1-rho)*rra + (rho/2)
    if(less_e_than(rho,value) && less_e_than(rho,1-value)){
      return(value)
    }
    
    #case 3
    value = rho*rra + ((1-rho)/2)
    if(greater_e_than(rho,value) && greater_e_than(rho,1-value)){
      return(value)
    }
    
    #case 4
    value = 1 - sqrt(2*rho*(1-rho)*(1-rra))
    if(less_e_than(rho,value) && greater_e_than(rho,1-value)){
      return(value)
    }
    
    #should never go here
    stop("Valid value for iso-J not found. Make sure that the RRA and rho values are correct.")
  }
  
  return(value)
}

iso_markedness_from_rra = function(rra, rho, precision=0.001){
  if(rra == 0){
    value = 0 # This value gives the iso-curve that intersecates the point (rho,rho)
  }else if(rra==1){
    value = 1
  }else{
    found = FALSE
    upper_value = 1
    lower_value = 0
    while(!found){
      mark = (upper_value + lower_value)/2
      curve = iso_Markedness(mark, rho)
      curve$x = c(0,curve$x,1)
      curve$y = c(0,curve$y,1)
      val = suppressWarnings(rra(curve$x, curve$y, 1, (1-rho)/rho, fallout = TRUE, recall = TRUE, plot = FALSE, print=FALSE)$rra)
      #cat(lower_value,"-",upper_value,"->", mark,"->",val,"\n") #for debug
      
      if(val == rra){
        value = mark
        found = TRUE
      }else if(val < rra){
        lower_value = mark
      }else{
        upper_value = mark
      }
      if(upper_value-lower_value < precision){
        value = (upper_value + lower_value)/2
        found = TRUE
      }
    }
  }
  
  return(round(value, digits = decimalplaces(precision)))
}
