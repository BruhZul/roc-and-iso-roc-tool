# This file is a part of the https://github.com/BruhZul/roc-and-iso-roc-tool repository.
#
# The tool is free software and can be redistributed and/or modified under the 
# GNU General Public License v3.0.
#
# This tool is distributed but without any warranty or implied warranty of
# merchantability or fitness for a particular purpose. See the LICENSE file for
# more details, or refer to gnu.org/licenses for more details.



#Function that generates the family of iso-pm curves. iso-pm curves are ROC curves
#that have the same value for a given PM for all their points.
#For example, iso-MCC curve for MCC=0.4 is the ROC curve so that all points have
#a MCC value of 0.4.
#
#Parameters:
#  metric: The PM that we want to generate iso-PM curves from.
#  metric_values: a vector of PM values. Each value will generate one iso-PM curve
#    with a specific PM value.
#  rho: prevalence rho (proportion of positive elements in the test dataset).
#  precision: When the iso curve is not a line (like MCC), the curve is generated
#    as a collection of points. Precision indicates how many points to generate 
#    (the default value is one point for each x=0.01, which should be more than
#      enough).
#
#The result is a list of lists (each one is one iso-PM curve for each value in 
#metric_values). It is used like this:
#  res = build_iso_curves(......)
#  res[[1]]$x is the list of x coordinates (FPR) for the first curve (1)
#  res[[1]]$y is the list of y coordinates (TPR)
build_iso_curves = function(metric, metric_values, rho, precision = 0.01){
  result = c()
  
  for(i in metric_values){
    points = c()
    if(metric == "MCC"){
      points = iso_MCC(i, rho, precision)
    }
    else if(metric == "BA"){
      points = iso_BA(i)
    }
    else if(metric == "Gmean"){
      points = iso_Gmean(i, precision)
    }
    else if(metric == "GM"){
      points = iso_GM(i, precision)
    }
    else if(metric == "D2H"){
      points = iso_D2H(i, precision)
    }
    else if(metric == "precision"){
      points = iso_precision(i, rho)
    }
    else if(metric == "TPR"){
      points = iso_TPR(i)
    }
    else if(metric == "F-score"){
      points = iso_Fscore(i, rho)
    }
    else if(metric == "NPV"){
      points = iso_NPV(i, rho)
    }
    else if(metric == "TNR"){
      points = iso_TNR(i)
    }
    else if(metric == "FPR"){
      points = iso_FPR(i)
    }
    else if(metric == "NM"){
      points = iso_NM(i, rho)
    }
    else if(metric == "J"){
      points = iso_J(i)
    }
    else if(metric == "Markedness"){
      points = iso_Markedness(i, rho, precision)
    }
    else{
      stop("Value for parameter 'metric' is not valid")
    }
    
    result = append(result, list(points))
  }
  
  return(result)
}

iso_MCC = function(value, rho, precision=0.01){
  x_vector = c()
  y_vector = c()
  
  p = value^2
  k = (1-rho)/rho
  for(x in seq(0,1,precision)){
    if(value == 0){
      y = x
    }else if(value==1){
      y_vector = c(y_vector,1, 1)
      x_vector = c(x_vector, 0,0)
      break
    }else{
      y = ((-1*(2*(p*k-k)*x-p*(k+1)))/(2*(p+k)))+
        ((suppressWarnings(sqrt((2*(p*k-k)*x-p*(k+1))^2-4*(k+p*k^2)*(p+k)*x^2 + 
                                  4*(p+k)*(p*(k^2+k))*x)))/(2*(p+k))) 
    }
    if(is.nan(y)){
      break
    }
    if(y > 1){
      break
    }
    y_vector = c(y_vector, y)
    x_vector = c(x_vector, x)
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_BA = function(value){
  c = 2*value - 1
  first_x = ifelse(c<0, -c, 0)
  first_y = ifelse(c<0, 0, c)
  second_x = ifelse(1-c<1, 1-c, 1)
  second_y = ifelse(1-c<1, 1, 1+c)
  x_vector = c(first_x, second_x)
  y_vector = c(first_y, second_y)
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_Gmean = function(value, precision=0.01){
  y_vector = c()
  x_vector = c()
  
  h = value^2
  for(x in seq(0,1,precision)){
    if(x == 1){
      break
    }
    y = h/(1-x)
    if(y > 1){
      y_vector = c(y_vector, 1)
      x_vector = c(x_vector, 1-h)
      break
    }
    y_vector = c(y_vector, y)
    x_vector = c(x_vector, x)
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_GM = function(value, precision = 0.01){
  y_vector = c()
  x_vector = c()
  
  h = value/2
  for(x in seq(0,1,precision)){
    y = h*((1-x)/(1-x-h))
    if(is.nan(y)){
      break
    }
    if(y < 0){
      y_vector = c(y_vector, 1)
      x_vector = c(x_vector, (2*h-1)/(h-1))
      break
    }
    if(y > 1){
      y_vector = c(y_vector, 1)
      x_vector = c(x_vector, (2*h-1)/(h-1))
      break
    }
    y_vector = c(y_vector, y)
    x_vector = c(x_vector, x)
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_D2H = function(value, precision=0.01){
  y_vector = c()
  x_vector = c()
  
  r = value*sqrt(2)
  ctr = 0
  for(x in seq(0,r,precision)){
    y = 1-suppressWarnings(sqrt(r^2-x^2))
    if(is.nan(y)){
      break
    }
    if(y < 0){
      next
    }
    if(length(x_vector) == 1 && ctr > 1){
      x_vector = c(sqrt(r^2-1), x_vector)
      y_vector = c(0, y_vector) 
    }
    if(x > 1){
      y_vector = c(y_vector, 1-sqrt(r^2-1))
      x_vector = c(x_vector, 1)
      break
    }else{
      y_vector = c(y_vector, y)
      x_vector = c(x_vector, x)
    }
    
    if(x+precision >= r){
      y_vector = c(y_vector, 1)
      x_vector = c(x_vector, r)
      break
    }
    ctr = ctr + 1
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_precision = function(value, rho){
  k = (1-rho)/rho
  
  first_x = 0
  first_y = 0
  if(value == 0){
    second_x = 1
    second_y = 0
  }else if(value==1){
    second_x = 0
    second_y = 1
  }else{
    second_x = ifelse(value<rho, 1, (1-value)/(value*k))
    second_y = ifelse(value<rho, (value*k)/(1-value), 1)
  }
  x_vector = c(first_x, second_x)
  y_vector = c(first_y, second_y)
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_TPR = function(value){
  x_vector = c(0, 1)
  y_vector = c(value, value)
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_Fscore = function(value, rho){
  k = (1-rho)/rho
  c = value/(2-value)
  
  first_x = 0
  first_y = c
  if(value == 0){
    second_x = 1
    second_y = 0 
  }else{
    second_x = ifelse(c<rho, 1, (1-c)/(c*k))
    second_y = ifelse(c<rho, c*k+c, 1)
  }
  
  x_vector = c(first_x, second_x)
  y_vector = c(first_y, second_y)
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_NPV = function(value, rho){
  k = (1-rho)/rho
  
  if(value == 0){
    x_vector = c(1, 1)
    y_vector = c(0, 1)
  }else if(value == 1){
    x_vector = c(0, 1)
    y_vector = c(1, 1)
  }else{
    first_x = ifelse(value > k/(1+k), 0, 1 - (value/((1-value)*k)))
    first_y = ifelse(value > k/(1+k), 1 - (((1-value)*k)/value), 0)
    
    second_x = 1
    second_y = 1
    
    x_vector = c(first_x, second_x)
    y_vector = c(first_y, second_y) 
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_TNR = function(value){
  x_vector = c(1-value, 1-value)
  y_vector = c(0, 1)
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_FPR = function(value){
  x_vector = c(value, value)
  y_vector = c(0, 1)
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_NM = function(value, rho){
  k = (1-rho)/rho
  
  if(value == 0){
    x_vector = c(1, 1)
    y_vector = c(0, 1)
  }else if(value == 1){
    x_vector = c(1, 1)
    y_vector = c(1, 1)
  }else{
    first_x = ifelse(value >= (2*k)/(1+2*k), 0, ((2*(1-value))/(2-value))-(value/((2-value)*k)))
    first_y = ifelse(value >= (2*k)/(1+2*k), 1 - ((2*(1-value)*k)/value), 0)
    
    second_x = (2*(1-value))/(2-value)
    second_y = 1
    
    x_vector = c(first_x, second_x)
    y_vector = c(first_y, second_y) 
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_J = function(value){
  if(value>=0){
    x_vector = c(0, 1-value)
    y_vector = c(value, 1)
  }else{
    x_vector = c(-value, 1)
    y_vector = c(0, 1+value)
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}

iso_Markedness = function(value, rho, precision=0.01){
  x_vector = c()
  y_vector = c()

  k = (1-rho)/rho
  for(x in seq(0,1,precision)){
    if(value == 0){
      y = x
    }else if(value==1){
      y_vector = c(y_vector,1, 1)
      x_vector = c(x_vector, 0,0)
      break
    }else{
      b = 2*value*k*x - value*k + k - value
      c = value*(k^2)*(x^2) - (value*(k^2) + value*k + k)*x
      y =(-b + suppressWarnings(sqrt((b^2)-(4*value*c))))/(2*value)
    }
    if(is.nan(y)){
      break
    }
    if(y > 1){
      break
    }
    y_vector = c(y_vector, y)
    x_vector = c(x_vector, x)
  }
  
  points = list(x_vector, y_vector)
  names(points) = c("x","y")
  return(points)
}
