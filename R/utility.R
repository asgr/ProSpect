interp_param=function(x, params, log=FALSE, method='linear'){
  flag=rep(2,length(x))
  flag[x<min(params)]=1
  flag[x>max(params)]=3
  x[x<min(params)]=min(params)
  x[x>max(params)]=max(params)
  if(log){
    temp=interp1(log10(params), 1:length(params), xi=log10(x), method=method)
  }else{
    temp=interp1(params, 1:length(params), xi=x, method=method)
  }
  flag[temp%%1==0 & flag==2]=0
  ID_lo=floor(temp)
  ID_hi=ceiling(temp)
  ID_mode=ID_lo
  ID_mode[temp%%1>0.5]=ID_hi[temp%%1>0.5]
  return(data.frame(x=x, param_lo=params[floor(temp)], param_hi=params[ceiling(temp)], ID_lo=ID_lo, ID_hi=ID_hi, ID_mode=ID_mode, weight_lo=1-temp%%1, weight_hi=temp%%1, flag=flag))
}