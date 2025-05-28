## Inoue+14 IGM attenuation 
## https://arxiv.org/abs/1402.0677
## Lyman continuum and Lyman series absorption given distribution function of z and NHI
## Consider both Lyman alpha forest and damped lyman systems

## Load in Lyman series coefficients for Eqs, 21-22
## Data obtained from https://github.com/gbrammer/eazy-py/tree/master/eazy/data

tau_IGM_LSLAF = function(l_obs, z){
 
   ##Lyman alpha forest lines
  Inoue14_LAFcoef = data.frame(
    j = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40), 
    lambda = c(1215.67, 1025.72, 972.537, 949.743, 937.803, 930.748, 926.226, 923.15, 920.963, 919.352, 918.129, 917.181, 916.429, 915.824, 915.329, 914.919, 914.576, 914.286, 914.039, 913.826, 913.641, 913.48, 913.339, 913.215, 913.104, 913.006, 912.918, 912.839, 912.768, 912.703, 912.645, 912.592, 912.543, 912.499, 912.458, 912.42, 912.385, 912.353, 912.324), 
    AJ1 = c(0.0168976, 0.00469229, 0.00223898, 0.00131901, 0.000870656, 0.000617843, 0.000460924, 0.000356887, 0.000284278, 0.000231771, 0.000192348, 0.000162155, 0.000138498, 0.000119611, 0.000104314, 9.17397e-05, 8.12784e-05, 7.25069e-05, 6.50549e-05, 5.86816e-05, 5.31918e-05, 4.84261e-05, 4.4274e-05, 4.06311e-05, 3.73821e-05, 3.45377e-05, 3.19891e-05, 2.9711e-05, 2.76635e-05, 2.58178e-05, 2.41479e-05, 2.26347e-05, 2.12567e-05, 1.99967e-05, 1.88476e-05, 1.77928e-05, 1.68222e-05, 1.59286e-05, 1.50996e-05), 
    AJ2 = c(0.00235379, 0.000653625, 0.000311884, 0.000183735, 0.00012128, 8.6064e-05, 6.42055e-05, 4.97135e-05, 3.95992e-05, 3.22851e-05, 2.67936e-05, 2.25878e-05, 1.92925e-05, 1.66615e-05, 1.45306e-05, 1.27791e-05, 1.13219e-05, 1.01e-05, 9.06198e-06, 8.17421e-06, 7.40949e-06, 6.74563e-06, 6.16726e-06, 5.65981e-06, 5.20723e-06, 4.81102e-06, 4.45601e-06, 4.13867e-06, 3.85346e-06, 3.59636e-06, 3.36374e-06, 3.15296e-06, 2.961e-06, 2.78549e-06, 2.62543e-06, 2.4785e-06, 2.3433e-06, 2.21882e-06, 2.10334e-06), 
    AJ3 = c(0.000102611, 2.8494e-05, 1.35962e-05, 8.00974e-06, 5.28707e-06, 3.75186e-06, 2.79897e-06, 2.1672e-06, 1.72628e-06, 1.40743e-06, 1.16804e-06, 9.84689e-07, 8.41033e-07, 7.2634e-07, 6.33446e-07, 5.57091e-07, 4.93564e-07, 4.40299e-07, 3.95047e-07, 3.56345e-07, 3.23008e-07, 2.94068e-07, 2.68854e-07, 2.46733e-07, 2.27003e-07, 2.09731e-07, 1.94255e-07, 1.80421e-07, 1.67987e-07, 1.56779e-07, 1.46638e-07, 1.3745e-07, 1.29081e-07, 1.2143e-07, 1.14452e-07, 1.08047e-07, 1.02153e-07, 9.67268e-08, 9.16925e-08)
  )
    
  tau_LAF_mat = matrix(0, ncol = length(l_obs),nrow = dim(Inoue14_LAFcoef)[1])
  for(i in 1:dim(Inoue14_LAFcoef)[1]){
    
    LJ = Inoue14_LAFcoef$lambda[i]
    
    ttau = rep(0, length(l_obs))
    idx1 = l_obs < 2.2*LJ & l_obs > LJ & l_obs < LJ*(1+z)
    ttau[idx1] = Inoue14_LAFcoef$AJ1[i] * (l_obs[idx1] / LJ)^1.2
    idx2 = l_obs >= 2.2*LJ & l_obs < 5.7*LJ & l_obs > LJ & l_obs < LJ*(1+z)
    ttau[idx2] = Inoue14_LAFcoef$AJ2[i] * (l_obs[idx2] / LJ)^3.7
    idx3 = l_obs >= 5.7*LJ & l_obs > LJ & l_obs < LJ*(1+z)
    ttau[idx3] = Inoue14_LAFcoef$AJ3[i] * (l_obs[idx3] / LJ)^5.5
    
    tau_LAF_mat[i,] = ttau
  }
  tau_LAF = colSums(tau_LAF_mat)

  return(tau_LAF)
}
tau_IGM_LSDLA = function(l_obs, z){
  
  ##Damped Lyman alpha lines
  Inoue14_DLAcoef <- data.frame(
    j = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40), 
    lambda = c(1215.67, 1025.72, 972.537, 949.743, 937.803, 930.748, 926.226, 923.15, 920.963, 919.352, 918.129, 917.181, 916.429, 915.824, 915.329, 914.919, 914.576, 914.286, 914.039, 913.826, 913.641, 913.48, 913.339, 913.215, 913.104, 913.006, 912.918, 912.839, 912.768, 912.703, 912.645, 912.592, 912.543, 912.499, 912.458, 912.42, 912.385, 912.353, 912.324), 
    AJ1 = c(0.000161698, 0.000154539, 0.000149767, 0.000146031, 0.000142893, 0.000140159, 0.000137714, 0.000135495, 0.000133452, 0.000131561, 0.000129785, 0.000128117, 0.00012654, 0.000125041, 0.000123614, 0.000122248, 0.000120938, 0.000119681, 0.000118469, 0.000117298, 0.000116167, 0.000115071, 0.000114011, 0.000112983, 0.000111972, 0.000111002, 0.000110051, 0.000109125, 0.00010822, 0.000107337, 0.000106473, 0.000105629, 0.000104802, 0.000103991, 0.000103198, 0.00010242, 0.000101657, 0.000100908, 0.000100168), 
    AJ2 = c(5.38995e-05, 5.15129e-05, 4.99222e-05, 4.86769e-05, 4.76312e-05, 4.67196e-05, 4.59048e-05, 4.5165e-05, 4.44841e-05, 4.38536e-05, 4.32617e-05, 4.27056e-05, 4.21799e-05, 4.16804e-05, 4.12046e-05, 4.07494e-05, 4.03127e-05, 3.98938e-05, 3.94896e-05, 3.90995e-05, 3.87225e-05, 3.83572e-05, 3.80037e-05, 3.76609e-05, 3.73241e-05, 3.70005e-05, 3.66836e-05, 3.63749e-05, 3.60734e-05, 3.57789e-05, 3.54909e-05, 3.52096e-05, 3.4934e-05, 3.46636e-05, 3.43994e-05, 3.41402e-05, 3.38856e-05, 3.36359e-05, 3.33895e-05)
  ) 
  
  tau_DLA_mat = matrix(0, ncol = length(l_obs),nrow = dim(Inoue14_DLAcoef)[1])
  for(i in 1:dim(Inoue14_DLAcoef)[1]){
    
    LJ = Inoue14_DLAcoef$lambda[i]
    
    ttau = rep(0, length(l_obs))
    idx1 = l_obs < 3.0*LJ & l_obs > LJ & l_obs < LJ*(1+z)
    ttau[idx1] = Inoue14_DLAcoef$AJ1[i] * (l_obs[idx1] / LJ)^2.0
    idx2 = l_obs >= 3.0*LJ & l_obs > LJ & l_obs < LJ*(1+z)
    ttau[idx2] = Inoue14_DLAcoef$AJ2[i] * (l_obs[idx2] / LJ)^3.0
    tau_DLA_mat[i,] = ttau
  }
  tau_DLA = colSums(tau_DLA_mat)

  return(tau_DLA)
}
tau_IGM_LCLAF = function(l_obs, z){
  
  ## Lyman continuum from Lyman alpha forest
  
  LL = 911.75 ## hardcode Lyman limit
  zshift_LL = LL * (1 + z)
  
  zLAF1 = 1.2
  zLAF2 = 4.7
  
  ## z < 1.2
  tau_LAF1 = rep(0, length(l_obs))
  LAF1_idx_1 = l_obs < zshift_LL & l_obs > LL
  tau_LAF1[LAF1_idx_1] = 0.325 * ((l_obs[LAF1_idx_1]/LL)^1.2 - (1+z)^(-0.9) * (l_obs[LAF1_idx_1]/LL)^2.1)
  
  ## 1.2 <= z < 4.7
  tau_LAF2 = rep(0, length(l_obs))
  LAF2_idx_1 = l_obs < 2.2*LL & l_obs > LL
  tau_LAF2[LAF2_idx_1] = 2.55e-2 * (1+z)^1.6 * (l_obs[LAF2_idx_1]/LL)^2.1 + 0.325*(l_obs[LAF2_idx_1]/LL)^1.2 - 0.250*(l_obs[LAF2_idx_1]/LL)^2.1
  LAF2_idx_2 = l_obs >= 2.2*LL & l_obs < zshift_LL
  tau_LAF2[LAF2_idx_2] = 2.55e-2 * ((1+z)^1.6 * (l_obs[LAF2_idx_2]/LL)^2.1 - (l_obs[LAF2_idx_2]/LL)^3.7)
  
  ## z >= 4.7
  tau_LAF3 = rep(0, length(l_obs)) 
  LAF3_idx_1 = l_obs < 2.2*LL & l_obs > LL
  tau_LAF3[LAF3_idx_1] = 5.22e-4 * (1+z)^3.4 * (l_obs[LAF3_idx_1]/LL)^2.1 + 0.325*(l_obs[LAF3_idx_1]/LL)^1.2 - 3.14e-2*(l_obs[LAF3_idx_1]/LL)^2.1
  LAF3_idx_2 = l_obs >= 2.2*LL & l_obs < 5.7*LL
  tau_LAF3[LAF3_idx_2] = 5.22e-4 * (1+z)^3.4 * (l_obs[LAF3_idx_2]/LL)^2.1 + 0.218*(l_obs[LAF3_idx_2]/LL)^2.1 - 2.55e-2*(l_obs[LAF3_idx_2]/LL)^3.7
  LAF3_idx_3 = l_obs >= 5.7*LL & l_obs < zshift_LL
  tau_LAF3[LAF3_idx_3] = 5.22e-4 * ((1+z)^3.4 * (l_obs[LAF3_idx_3]/LL)^2.1 - (l_obs[LAF3_idx_3]/LL)^5.5)
  
  if(z < zLAF1){
    tau_LAF = tau_LAF1
  }else if(z >= zLAF1 & z < zLAF2){
    tau_LAF = tau_LAF2
  }else{
    tau_LAF = tau_LAF3
  }
  
  return(tau_LAF)
}
tau_IGM_LCDLA = function(l_obs, z){
  
  ## Lyman continuum from damped lyman systems
  
  LL = 911.75 ## hardcode Lyman limit
  zshift_LL = LL * (1 + z)
  
  zDLA = 2.0
  
  ## z < 2.0
  tau_DLA1 = rep(0, length(l_obs))
  DLA1_idx_1 = l_obs < zshift_LL & l_obs > LL
  tau_DLA1[DLA1_idx_1] = 0.211*(1+z)^2.0 - 7.66e-2 * (1+z)^2.3 * (l_obs[DLA1_idx_1]/LL)^-0.3 - 0.135*(l_obs[DLA1_idx_1]/LL)^2.0
  
  ## z >= 2.0
  tau_DLA2 = rep(0, length(l_obs))
  DLA2_idx_1 = l_obs < 3.0*LL & l_obs > LL
  tau_DLA2[DLA2_idx_1] = 0.634 + 4.7e-2*(1+z)^3 - 1.78e-2*(1+z)^3.3*(l_obs[DLA2_idx_1]/LL)^-0.3 - 0.135*(l_obs[DLA2_idx_1]/LL)^2.0 - 0.291*(l_obs[DLA2_idx_1]/LL)^-0.3
  DLA2_idx_2 = l_obs >= 3.0*LL & l_obs < zshift_LL & l_obs > LL
  tau_DLA2[DLA2_idx_2] = 4.7e-2*(1+z)^3 - 1.78e-2*(1+z)^3.3*(l_obs[DLA2_idx_2]/LL)^-0.3 - 2.92e-2*(l_obs[DLA2_idx_2]/LL)^3.0
  
  if(z < zDLA){
    tau_DLA = tau_DLA1
  }else{
    tau_DLA = tau_DLA2
  }
  
  return(tau_DLA)
}
tau_IGM_Tot = function(l_obs, z){
  tau_IGM_LSLAF(l_obs, z) +
    tau_IGM_LSDLA(l_obs, z) + 
    tau_IGM_LCLAF(l_obs, z) + 
    tau_IGM_LCDLA(l_obs, z)
}
Inoue14_IGM = function(l_obs, z){
  # pmin(1.0, exp(-1 * tIGM(l_obs, z)))
  exp(-1 * tau_IGM_Tot(l_obs, z))
}
