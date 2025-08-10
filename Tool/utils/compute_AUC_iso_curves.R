# This file is a part of the https://github.com/BruhZul/roc-and-iso-roc-tool repository.
#
# The tool is free software and can be redistributed and/or modified under the 
# GNU General Public License v3.0.
#
# This tool is distributed but without any warranty or implied warranty of
# merchantability or fitness for a particular purpose. See the LICENSE file for
# more details, or refer to gnu.org/licenses for more details.



#Function that, given the Area Under the Curve (AUC) of a ROC curve, computes
#the value of the iso-metric curve having the same AUC.
#For example, if AUC=0.5, the iso-MCC having the same AUC will have a value
#MCC=1.
#
#Parameters:
#  metric. The metric of the iso-metric curve to compute.
#  auc. The AUC value of the ROC curve.
#  rho. Prevalence rho (the proportion of positive elements in the test dataset)
#  precision. Some metrics are computed using simulation instead of mathematical
#    formulas. In these cases precision is used, and represents the number of 
#    digits when rounding the result.
#
#Examples:
#  value = compute_iso_value_from_auc("MCC", 0.75, 0.3) returns 0.303. We can 
#    then build the iso-MCC curve by using the build_iso_curves() function (see
#    the file build_iso_curves.R for more details). 
compute_iso_value_from_auc = function(metric, auc, rho, precision=0.001){
  result = c()
  if(metric == "MCC"){
    result = iso_mcc_from_auc(auc, rho, precision)
  }
  else if(metric == "BA"){
    result = iso_ba_from_auc(auc)
  }
  else if(metric == "Gmean"){
    result = iso_gmean_from_auc(auc,precision)
  }
  else if(metric == "GM"){
    result = iso_gm_from_auc(auc, precision)
  }
  else if(metric == "D2H"){
    result = iso_D2H_from_auc(auc, precision)
  }
  else if(metric == "precision"){
    result = iso_precision_from_auc(auc, rho)
  }
  else if(metric == "TPR"){
    result = iso_tpr_from_auc(auc)
  }
  else if(metric == "F-score"){
    result = iso_fscore_from_auc(auc, rho)
  }
  else if(metric == "NPV"){
    result = iso_npv_from_auc(auc, rho)
  }
  else if(metric == "TNR"){
    result = iso_tnr_from_auc(auc)
  }
  else if(metric == "FPR"){
    result = iso_fpr_from_auc(auc)
  }
  else if(metric == "NM"){
    result = iso_nm_from_auc(auc, rho)
  }
  else if(metric == "J"){
    result = iso_j_from_auc(auc)
  }
  else if(metric == "Markedness"){
    result = iso_markedness_from_auc(auc, rho, precision)
  }
  else{
    stop("Value for parameter 'metric' is not valid")
  }
  
  return(result)
}

iso_mcc_from_auc = function(auc, rho, precision=0.001){
  k = (1-rho)/rho
  
  if(auc==0.5){
    value = 0
  }else if(auc == 1){
    value = 1
  }else if(auc < 0.5){ #Per simulazione
    mcc_vals = seq(-1,0-precision,precision)
    
    for(i in 1:(length(mcc_vals)-1)){
      mcc_1 = mcc_vals[i]
      mcc_2 = mcc_vals[i+1]
      
      mcc_function_1 = function(x){
        a = (mcc_1^2)+k
        b = 2*((mcc_1^2)*k-k)*x - (mcc_1^2)*(k+1)
        c = (k+(mcc_1^2)*(k^2))*(x^2) - ((mcc_1^2)*((k^2)+k))*x
        (-b - sqrt((b^2)-4*a*c))/(2*a)
      }

      mcc_function_2 = function(x){
        a = (mcc_2^2)+k
        b = 2*((mcc_2^2)*k-k)*x - (mcc_2^2)*(k+1)
        c = (k+(mcc_2^2)*(k^2))*(x^2) - ((mcc_2^2)*((k^2)+k))*x
        (-b - sqrt((b^2)-4*a*c))/(2*a)
      }
      
      lower_limit = function(mcc){
        return(((mcc^2)*((k^2)+k))/(k+(mcc^2)*(k^2)))
      }
      
      val_1 = integrate(mcc_function_1, lower = lower_limit(mcc_1), upper = 1)$value
      val_2 = integrate(mcc_function_2, lower = lower_limit(mcc_2), upper = 1)$value
      
      if(val_1 == auc){
        value = mcc_1
        break
      }
      if(val_1 < auc && val_2 > auc){
        value = ifelse(abs(auc-val_1) > abs(auc-val_2), mcc_2, mcc_1)
        break
      }
    }
  }else{
    mcc_vals = seq(precision,1,precision)
    
    for(i in 1:(length(mcc_vals)-1)){
      mcc_1 = mcc_vals[i]
      mcc_2 = mcc_vals[i+1]
      
      mcc_function_1 = function(x){
        a = (mcc_1^2)+k
        b = 2*((mcc_1^2)*k-k)*x - (mcc_1^2)*(k+1)
        c = (k+(mcc_1^2)*(k^2))*(x^2) - ((mcc_1^2)*((k^2)+k))*x
        (-b + sqrt((b^2)-4*a*c))/(2*a)
      }
      
      mcc_function_2 = function(x){
        a = (mcc_2^2)+k
        b = 2*((mcc_2^2)*k-k)*x - (mcc_2^2)*(k+1)
        c = (k+(mcc_2^2)*(k^2))*(x^2) - ((mcc_2^2)*((k^2)+k))*x
        (-b + sqrt((b^2)-4*a*c))/(2*a)
      }
      
      upper_limit = function(mcc){
        p = mcc^2
        a = k + p*(k^2)
        b = 2*(p*k - k)
        c = p + k
        d = -p*((k^2) + k)
        e = -p*(k+1)
        return((-(b+d)-sqrt(((b+d)^2)-4*a*(c+e)))/(2*a))
      }
      
      val_1 = integrate(mcc_function_1, lower = 0, upper = upper_limit(mcc_1))$value + 1 - upper_limit(mcc_1)
      val_2 = integrate(mcc_function_2, lower = 0, upper = upper_limit(mcc_2))$value + 1 - upper_limit(mcc_2)
      
      if(val_1 == auc){
        value = mcc_1
        break
      }
      if(val_1 < auc && val_2 > auc){
        value = ifelse(abs(auc-val_1) > abs(auc-val_2), mcc_2, mcc_1)
        break
      }
    }
  }
  
  return(value)
}

iso_ba_from_auc = function(auc){
  if(auc <= 0.5){
    value = sqrt(auc/2)
  }else{
    value = 1 - sqrt((1-auc)/2)
  }
  
  return(value)
}

iso_gmean_from_auc = function(auc, precision=0.001){
  if(auc==0){
    value = 0
  }else if(auc==1){
    value = 1
  }else{
    # per simulazione
    gmean_vals = seq(precision, 1, precision)
    for(i in 1:(length(gmean_vals)-1)){
      gmean_1 = gmean_vals[i]
      gmean_2 = gmean_vals[i+1]
      h_1 = gmean_1^2
      h_2 = gmean_2^2
      
      val_1 = h_1*(1-log(h_1))
      val_2 = h_2*(1-log(h_2))
      
      if(val_1 == auc){
        value = gmean_1
        break
      }
      if(val_1 < auc && val_2 > auc){
        value = ifelse(abs(auc-val_1) > abs(auc-val_2), gmean_2, gmean_1)
        break
      }
    }
  }
  
  return(value)
}

iso_gm_from_auc = function(auc, precision=0.001){
  if(auc==0){
    value = 0
  }else if(auc==1){
    value = 1
  }else{
    gm_vals = seq(precision, 1, precision)
    for(i in 1:(length(gm_vals)-1)){
      gm_1 = gm_vals[i]
      gm_2 = gm_vals[i+1]
      h_1 = gm_1/2
      h_2 = gm_2/2
      
      val_1 = 2*h_1 + 2*(h_1^2)*(log(1-h_1)-log(h_1))
      val_2 = 2*h_2 + 2*(h_2^2)*(log(1-h_2)-log(h_2))
      
      if(val_1 == auc){
        value = gm_1
        break
      }
      if(val_1 < auc && val_2 > auc){
        value = ifelse(abs(auc-val_1) > abs(auc-val_2), gm_2, gm_1)
        break
      }
    }
  }
  
  return(value)
}

iso_D2H_from_auc = function(auc, precision = 0.001){
  if(auc < 1-(pi/4)){
    if(auc == 0){
      value = 1
    }else{
      d2h_vals = seq(sqrt(2)/2, 1, precision)
      for(i in 1:(length(d2h_vals)-1)){
        d2h_1 = d2h_vals[i]
        d2h_2 = d2h_vals[i+1]
        r_1 = d2h_1*sqrt(2)
        r_2 = d2h_2*sqrt(2)
        
        d2h_function_1 = function(x){
          1-sqrt(r_1^2-x^2)
        }
        
        d2h_function_2 = function(x){
          1-sqrt(r_2^2-x^2)
        }
        
        val_1 = integrate(d2h_function_1, lower = sqrt(r_1^2-1), upper = 1)$value
        val_2 = integrate(d2h_function_2, lower = sqrt(r_2^2-1), upper = 1)$value
        
        if(val_1 == auc){
          value = d2h_1
          break
        }
        if(val_1 > auc && val_2 < auc){
          value = ifelse(abs(auc-val_1) > abs(auc-val_2), d2h_2, d2h_1)
          break
        }
      }
    }
  }else{
    value = sqrt((2*(1-auc))/pi)
  }
  
  return(value)
}

iso_precision_from_auc = function(auc, rho){
  k = (1-rho)/rho
  
  if(auc == 1){
    value = 1
  }else if(auc<=0.5){
    value = (2*auc)/(k+2*auc)
  }else{
    value = 1/(2*k+1-2*k*auc)
  }
  
  return(value)
}

iso_tpr_from_auc = function(auc){
  return(auc)
}

iso_fscore_from_auc = function(auc, rho){
  k = (1-rho)/rho
  
  # of fscore is greater than 2/(2+k)
  a = k + 2 - k*auc
  b = -(2*k + 4 - 2*k*auc)
  c = 2
  value_1 = (-b - suppressWarnings(sqrt(b^2-4*a*c)))/(2*a)
  
  # if fscore is lower than 2/(2+k)
  value_2 = (4*auc)/(2+k+2*auc)
  
  if(is.nan(value_1) || value_1 < 2/(2+k)){
    value = value_2
  }else{
    value = value_1
  }
  
  
  return(value)
}

iso_npv_from_auc = function(auc, rho){
  k = (1-rho)/rho
  
  if(auc == 0){
    value = 0
  }else if(auc<=0.5){
    value = (2*k*auc)/(1+(2*k*auc))
  }else{
    value = k/(2+k-2*auc)
  }
  
  return(value)
}

iso_fpr_from_auc = function(auc){
  return(1-auc)
}

iso_tnr_from_auc = function(auc){
  return(auc)
}

iso_nm_from_auc = function(auc, rho){
  k = (1-rho)/rho
  
  if(auc == 0){
    value = 0
  }else{
    # if nm is greater than 2k/(2k+1) 
    a = -2-4*k+2*auc
    b = 2+4*k-2*auc
    c = -4*k
    value_1 = (-b + suppressWarnings(sqrt(b^2-a*c)))/(a)
    
    # if nm is lower than 2/(2+k)
    value_2 = (4*k*auc)/(2*k+1+2*k*auc)
    
    if(is.nan(value_1) || value_1 < (2*k)/(2*k+1)){
      value = value_2
    }else{
      value = value_1
    } 
  }
  
  return(value)
}

iso_j_from_auc = function(auc){
  if(auc <= 0.5){
    value = -1 + sqrt(2*auc)
  }else{
    value = 1 - sqrt(2*(1-auc))
  }
  
  return(value)
}

iso_markedness_from_auc = function(auc, rho, precision=0.001){ #da fare
  k = (1-rho)/rho
  
  if(auc==1){
    value = 1
  }else if(auc == 0){
    value = -1
  }else if(auc == 0.5){
    value = 0
  }else{
    if(auc < 0.5){
      mark_vals = seq(-1,0-precision,precision)
    }else{
      mark_vals = seq(precision,1,precision)
    }
    for(i in 1:(length(mark_vals)-1)){
      m_1 = mark_vals[i]
      m_2 = mark_vals[i+1]
      
      m_function_1 = function(x){
        a = m_1*(k^2)
        b = 2*m_1*k
        d = -(m_1*(k^2) + m_1*k +k)
        e = -(m_1*k - k + m_1)
        ((-(b*x+e))+sqrt((b*x+e)^2 - 4*m_1*(a*(x^2)+d*x)))/(2*m_1)
      }
      
      m_function_2 = function(x){
        a = m_2*(k^2)
        b = 2*m_2*k
        d = -(m_2*(k^2) + m_2*k +k)
        e = -(m_2*k - k + m_2)
        ((-(b*x+e))+sqrt((b*x+e)^2 - 4*m_2*(a*(x^2)+d*x)))/(2*m_2)
      }
      
      if(auc < 0.5){
        lower_limit = function(m){
          a = m*(k^2)
          d = (m*(k^2) + m*k +k)
          return(max(0, d/a))
        }
        
        val_1 = integrate(m_function_1, lower = lower_limit(m_1), upper = 1)$value
        val_2 = integrate(m_function_2, lower = lower_limit(m_2), upper = 1)$value
      }else{
        upper_limit = function(m){
          a = m*(k^2)
          b = 2*m*k
          d = -(m*(k^2) + m*k +k)
          e = -(m*k - k + m)
          return((-(b+d)-sqrt((b+d)^2-4*a*(m+e)))/(2*a))
        }
        
        val_1 = integrate(m_function_1, lower = 0, upper = upper_limit(m_1))$value + 1 - upper_limit(m_1)
        val_2 = integrate(m_function_2, lower = 0, upper = upper_limit(m_2))$value + 1 - upper_limit(m_2)
      }
      
      
      if(val_1 == auc){
        value = m_1
        break
      }
      if(val_1 < auc && val_2 > auc){
        value = ifelse(abs(auc-val_1) > abs(auc-val_2), m_2, m_1)
        break
      }
    }
  }
  
  return(value)
}








