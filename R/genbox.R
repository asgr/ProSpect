genbox=function(sfunc = function(t, total, gasfrac, ssfr){1}, time = 10, step = 0.05, ssfr = 0,
                alpha = 0.93, total = 1e+10, gasfrac = 1, starfrac = 1 - gasfrac,
                infunc = function(t,total, sin){0}, sin = 0, outfunc = function(t, total, gas2stars, alpha, sout) {0},
                sout = 0, Zsn = 0.13, Zgas = 0, Zstars = 0, dgas = 0, Zin = 0, Chi = 0.16,
                Chiin = 0,destroy = 0.01, yield){
  # metalicity enrichment model
  # 2011 Aaaron Robotham St Andrews / UWA
  # If equal Chi model then Chi=0.16. destroy should be below 1 (0.01 ish)
  output={}
  
  timesteps=seq(0,time,by=step)
  
  # run the model
  for(t in timesteps){
    gas2stars=sfunc(t,total,gasfrac,ssfr)
    infall=infunc(t,total,sin)*step*1e9
    outflow=outfunc(t,total,gas2stars,alpha,sout)*step*1e9
    gas2stars=gas2stars*step*1e9
    temp=.genstep(infall=infall, outflow=outflow, gas2stars=gas2stars, alpha=alpha, Zsn=Zsn, total=total, gasfrac=gasfrac, starfrac=starfrac, Zgas=Zgas, Zstars=Zstars, dgas=dgas, Zin=Zin, Chi=Chi, Chiin=Chiin, destroy=destroy, yield=yield)
    total=as.numeric(temp['total'])
    gasfrac=as.numeric(temp['gasfrac'])
    starfrac=as.numeric(temp['starfrac'])
    Zgas=as.numeric(temp['Zgas'])
    Zstars=as.numeric(temp['Zstars'])
    dgas=as.numeric(temp['dgas'])
    
    temp=c(SFR=gas2stars/(step*1e9),temp)
    
    if(isTRUE(starfrac<=1)){
      output=rbind(output,temp)
    }else{break}
  }
  output=cbind(time=timesteps[1:length(output[,1])], output)
  output=cbind(output, sumstarsform=cumsum(output[,'addstars']))
  rownames(output)=rep('',length(output[,1]))
  
  invisible(as.data.frame(output))
}



.genstep=function(gas2stars=0.01, alpha=0.93, Zsn=0.13, infall=0, outflow=0, total=1, gasfrac=1, starfrac=0, Zgas=0, Zstars=0, dgas=0, Zin=0, Chi=0.16, Chiin=0, destroy=0.01, yield){
  # metalicity enrichment model
  # 2012 Aaron Robotham St Andrews / UWA
  # If equal Chi model then Chi=0.16. destroy should be below 1 (0.01 ish)
  
  if(!missing(yield)){
    Zsn = Zgas + (yield*alpha)/(1-alpha)
  }
  
  # All "sum" variables specify total mass they contain.
  
  # initial
  
  sumgas = total*gasfrac				                        #total gas
  sumstars = total*starfrac			                        #total stars
  sumZgas = sumgas*Zgas				                          #mass of metals in gas
  sumZstars = sumstars*Zstars			                      #mass of metals in stars
  sumZall = sumZgas + sumZstars		                      #mass of all metals
  sumdgas = sumgas*dgas				                          #mass of dust in gas
  
  # inflow
  addmass = infall            													#rate of infall in proportional to total mass or constant										
  # add inflow to gas
  sumgas = sumgas + addmass													    #all infall material is gas
  sumZgas = sumZgas + Zin * addmass											#add the metal component of infall to metals in gas
  sumdgas = sumdgas + Chiin * addmass										#add the dust component of infall to dust in gas
  # form stars
  addstars = gas2stars    												      #convert a fraction of the gas to stars
  recycle = ( 1 - alpha ) * addstars										#instantaneously feedback some fraction
  dead = alpha * addstars														    #lock the rest up permanently
  # Woo, forming stars!
  sumgas = sumgas - addstars														#mass of gas decreases by the SF
  sumstars = sumstars + addstars												#mass of stars increases by the SF
  sumZgas = sumZgas - Zgas * addstars										#mass of metals in gas decreases by the SF times Zgas (this will then instantaneously feedback some fraction)
  sumdgas = sumdgas - dgas * addstars										#mass of dust in dust decreases by the SF times dgas
  sumZstars = sumZstars + Zgas * addstars								#mass of metals in stars increases by the SF times Zgas
  # And bang, lets blow some up.
  sumgas = sumgas + recycle													    #mass of gas increases by the recycled amount
  sumstars = sumstars - recycle												  #mass of stars decrease by the recycled amount
  sumZgas = sumZgas + Zsn * recycle											#mass of metals in gas increases by the recycled amount times Zsn
  sumdgas = sumdgas + Zsn * Chi * recycle								#mass of dust in gas increases by the recycled amount times Zsn times Chi
  sumZstars = sumZstars - Zgas * recycle								#mass of metals in stars decreases by the recycled amount times Zgas
  # outflow
  minusmass = outflow         													#rate of outflow in proportional to total mass or constant
  # subtract outflow from gas	
  sumgas = sumgas - minusmass													  #remove outflow gas from total gas
  sumZgasout = Zgas * minusmass                         #outflow metal mass
  sumZgas = sumZgas - 	sumZgasout            					#remove outflow metals from gas metals
  sumdgasout = dgas * minusmass                         #outflow dust mass
  sumdgas = sumdgas - 	sumdgasout              				#remove outflow dist from gas dust
  # Just for dust we have dust destruction:
  sumdgas = sumdgas - sumdgas * destroy * recycle       #dust is destroyed in prop to current SN rate and current dust mass
  
  # new numbers
  yield = (Zsn - Zgas) * (1 - alpha) / alpha						#yield = mass of medals ADDED to ISM (adds-alpha*adds)*(Zsn-Zgas) divided by mass lost from ISM (adds*alpha)
  sumZall = sumZgas + sumZstars												  #total mass in metals
  Zgas = sumZgas / sumgas														    #metal fraction in gas recalc
  Zstars = sumZstars / sumstars												  #metal fraction in stars recalc
  dgas = sumdgas / sumgas														    #dust fraction in gas recalc
  total = sumgas + sumstars													    #total material in the SF region
  starfrac = sumstars / total													  #total fraction of mass in stars
  gasfrac = sumgas / total													    #total fraction of mass in gas
  Zfrac = sumZall / total														    #metal fraction in the SF region 
  dustfrac = sumdgas / total													  #dust fraction in the SF region
  
  output=c(sumgas=sumgas, sumstars=sumstars, total=total, sumZgas=sumZgas, sumZstars=sumZstars, sumZall=sumZall, sumdgas=sumdgas, Zgas=Zgas, Zstars=Zstars, Zfrac=Zfrac, dgas=dgas, gasfrac=gasfrac, starfrac=starfrac, dustfrac=dustfrac, yield=yield, addstars=addstars, dead=dead, recycle=recycle, infall=infall, outflow=outflow)
  invisible(output)
}

