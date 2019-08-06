\name{genbox}
\alias{genbox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leaky Box Star Formation Model
}
\description{
A leaky box star formation model, that allows for farily concept inflow and outflow models. It tracks metals and dust in a self consistent manner. In the simplest mode it is simply the classic closed box model.
}
\usage{
genbox(sfunc = function(t, total, gasfrac, argsfr){1}, time = 10, step = 0.05, argsfr = 0,
alpha = 0.8, total = 1e+10, gasfrac = 1, starfrac = 1 - gasfrac, infunc = function(t,
total, argin){0}, argin = 0, outfunc = function(t, total, gas2stars, alpha, argout) {0},
argout = 0, Zsn = 0.12, Zgas = 0, Zstars = 0, dgas = 0, Zin = 0, Chi = 0.16, Chiin = 0,
destroy = 0.01, yield)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sfunc}{
Star formation functions, where the output has to be the star formation rate in Msol/yr. This function has to have arguments \option{t} for time in years, \option{total} (total current baryonic mass of the system in Msol), \option{gasfrac} (the gas fraction), and \option{argsfr} (the specific star formation rate, where what this means depends on the function- see Examples). It is up to the user which of these are used to help produce a SFR, e.g. the default is constant and 1 Msol/yr.
}
  \item{time}{
The time to run the box for (Gyr).
}
  \item{step}{
Step size to take between periods of discrete star formation (Gyr).
}
  \item{argsfr}{
The desired star formation rate modfying argument to be passed into \option{sfunc} as an argument. Only relevant if this is used in \option{sfunc}.
}
  \item{alpha}{
The fraction of mass locked up in a period of star formation. The default of 0.8 is close the fraction of mass in stars above 10 Msol for a Chabrier IMF, so is appropriate for use with \code{ProSpect} since the various SSPs use the Chabrier IMF.
}
  \item{total}{
Total initial baryonic mass of the system (Msol).
}
  \item{gasfrac}{
The gas fraction of the system .
}
  \item{starfrac}{
the stellar fraction of the system .
}
  \item{infunc}{
Inflow function, where the the output has to be the gas mass inflow in Msol/yr. This function has to have arguments \option{t} for time in years, \option{total} (total current baryonic mass of the system in Msol), and \option{argin} which can mean whatever you want it to mean. The default is no inflow.
}
  \item{argin}{
The main modifying input to \option{infunc}.
}
  \item{outfunc}{
Outflow function, where the the output has to be the gas mass inflow in Msol/yr. This function has to have arguments \option{t} for time in years, \option{total} (total current baryonic mass of the system in Msol),  \option{gas2stars} (the mass of gas just converted to stars), \option{alpha} (the fraction of mass locked up in a period of star formation), \option{argout} which can mean whatever you want it to mean. The default is no outflow.
}
  \item{argout}{
The main modifying input to \option{outfunc}.
}
  \item{Zsn}{
Metallcity of the gas ejected by a SN event (a value \option{Zsn}=0.12 will produce a mean yield near 0.03 since \option{yield} = \option{Zsn} x (1 - \option{alpha}) / \option{alpha} for low gas phase metallcities).
}
  \item{Zgas}{
Initial metallicity of the gas.
}
  \item{Zstars}{
Initial metallicity of the stars.
}
  \item{dgas}{
Initial dust fraction of the gas.
}
  \item{Zin}{
The metallicity of any inflowing gas. Only relevant if \option{infunc} is non-zero.
}
  \item{Chi}{
Fraction of dust in the metals formed.
}
  \item{Chiin}{
Fraction of dust in the metals inflowing. Only relevant if \option{infunc} is non-zero and \option{Zin}>0.
}
  \item{destroy}{
Dust mass destroyed, where in a given step it is calculated as \option{destroy} x dustmass x recyclemass. Here dustmass is the total mass in the form of dust, and recyclemass is the total gas gas ejected by a SN event. This means more dust gets destroyed if there is more dust present, and also if there are more SN events.
}
  \item{yield}{
The yield, this is defined as mass of medals ADDED to ISM (gas2stars-alpha*gas2stars)*(Zsn-Zgas) divided by mass lost from ISM (gas2stars*alpha), i.e. yield = (Zsn - Zgas) * (1 - alpha) / alpha. Strictly it should be derived, not defined. But to be consistent with "classic" closed box star formation, it can be set to constant value using the \option{yield} argument. In this case \option{Zsn} is defined to be consistent with the desired yield. Since \option{yield} = \option{Zsn} x (1 - \option{alpha}) / \option{alpha} for low gas phase metallcities, the default \code{ProSpect} values will return a mean yield of roughly 0.03.
}
}
\details{
This function allows for fairly complex star formation scenarios, especially in regards to inflow and outflow models, e.g. check the examples to see how outflow can be coupled to SN events.
}
\value{
A data.frame containing information at each time step specified by \option{step}:

\item{time}{Time (Gyr).}
\item{SFR}{Star formation rate (Msol/yr).}
\item{sumgas }{Mass of gas (Msol).}
\item{sumstars }{Mass of stars (Msol).}
\item{total}{Total baryonic mass (Msol), i.e. \option{total} = \option{sumgas} + \option{sumstars}.}
\item{sumZgas}{Mass of metals in gas (Msol).}
\item{sumZstars}{Mass of metals in stars (Msol).}
\item{sumZall}{Mass of all metals (Msol), i.e. \option{sumZall} = \option{sumZgas} + \option{sumZstar}.}
\item{Zgas}{Fraction of gas mass in metals.}
\item{Zstars}{Fraction of stellar mass in metals.}
\item{Zfrac}{Fraction of all mass in metals.}
\item{dgas}{Fraction of gas mass in dust.}
\item{gasfrac}{Fraction of mass in gas, i.e. \option{sumgas}/\option{total}.}
\item{starfrac}{Fraction of mass in stars, i.e. \option{sumstars}/\option{total}.}
\item{dustfrac}{Fraction of mass in dust.}
\item{yield}{Computed yield, where yield = (Zsn - Zgas) * (1 - alpha) / alpha.}
\item{addstars}{Stars formed in the current step (Msol).}
\item{dead}{Formed stars locked up in the current step (Msol).}
\item{recycle}{Formed stars recylced as gas in the current step (Msol).}
\item{infall}{Infall mass in the current step (Msol).}
\item{outfall}{Outfall mass in the current step (Msol).}
\item{sumstarsform}{Total star formed, i.e. \option{SFR} x \option{step} x 1e9 (Msol).}
}
\references{
Largely based on code I wrote for chemical evolution lectures given in St Andrews and UWA.
}
\author{
Aaron Robotham
}
\note{
The classic Zgas = -yield x ln(gasfrac) only works for this constant yield model, which does not make physical sense. Basically, do not use the \option{yield} argument except for comparison to other simple closed box models models. It should be obvious why a constant yield breaks- SN produce close to constant Zsn (Zgas = 0.1-0.2 in the recycled gas ejected by a SNII) regardless of the metallicity they form with (there is a very weak dependency on this), so as your gas enriches you efficiency of enrichment should drop. The constant yield model can produce Zgas much larger than Zsn, which really cannot make sense! People still insist on using it, but do not be one of them.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{IMF}}, ~~~
}
\examples{
#Remember, do not use fixed yield! But this recreates the classic closed box result:

testbox=genbox(yield=0.03)

plot(testbox[,c('gasfrac','Zgas')], type='l', lwd=3, xlab='Gas Fraction', ylab='Zgas')
curve(-0.03*log(x), add=TRUE, col='red', lwd=3, lty=3)

#In reality for an even remotely realistic model the yield must vary with Zgas:

testbox=genbox()
summary(testbox$yield)

plot(testbox[,c('gasfrac','Zgas')], type='l', lwd=3, xlab='Gas Fraction', ylab='Zgas')
curve(-mean(testbox$yield)*log(x), add=TRUE, col='red', lwd=3, lty=3)

#Here we will make a more complex inflow and outflow leaky box model:

#Infall proportional to the current mass of system:
InfallTot=function(t,total,argin){return=total*argin}
#Outflow coupled to SN events, where argout is in effect the mass loading:
OutflowSN=function(t,total,gas2stars,alpha,argout){gas2stars*(1-alpha)*argout}

testbox=genbox(infunc=InfallTot, argin=2e-10, outfunc=OutflowSN, argout=5)

plot(testbox[,c('gasfrac','Zgas')], type='l', lwd=3, xlab='Gas Fraction', ylab='Zgas')
curve(-mean(testbox$yield)*log(x), add=TRUE, col='red', lwd=3, lty=3)

plot(testbox[,c('time','Zgas')], type='l', lwd=3, xlab='Time / Gyr', ylab='Zgas')
lines(testbox[,c('time','Zstars')], lwd=3, col='red')

plot(testbox[,c('time','total')], type='l', lwd=3, xlab='Time / Gyr',
ylab='Total Mass / Msol')

#Some other models:

#Constant inflow over time:
InfallConst=function(t,total,argin){return=argin}

#Star formation is in proportion to the stellar mass (traditional definition of sSFR):
sfuncSSFRstar=function(t,total,gasfrac,argsfr){return=total*argsfr*(1-gasfrac)}

#Star formation is in proportion to the gas mass:
sfuncSSFRgas=function(t,total,gasfrac,argsfr){return=total*argsfr*gasfrac}

#Star formation is in proportion to the total mass of the system:
sfuncSSFRtot=function(t,total,gasfrac,argsfr){return=total*argsfr}

#This can make quite complex SFHs with different metal enrichment rates:

teststar=genbox(sfunc=sfuncSSFRstar, argsfr=1e-10, gasfrac=0.5, infunc=InfallConst,
argin=1, outfunc=OutflowSN, argout=10)
testgas=genbox(sfunc=sfuncSSFRgas, argsfr=1e-10, gasfrac=0.5, infunc=InfallConst, argin=1,
outfunc=OutflowSN, argout=10)
testtot=genbox(sfunc=sfuncSSFRtot, argsfr=1e-10, gasfrac=0.5, infunc=InfallConst, argin=1,
outfunc=OutflowSN, argout=10)

plot(teststar[,c('time','SFR')], type='l', xlab='Time / Gyr', ylab='SFR / Msol/yr',
ylim=c(0,1.5), col='red')
lines(testgas[,c('time','SFR')], col='blue') #This gives a constant (and lowest) SFR!
lines(testtot[,c('time','SFR')])

plot(teststar[,c('gasfrac','Zgas')], type='l', xlab='Gas Fraction', ylab='Zgas',
  xlim=c(0,0.5), ylim=c(0,0.02), col='red')
lines(testgas[,c('gasfrac','Zgas')], col='blue') #Fastest enrichment despite lowest SFR
lines(testtot[,c('gasfrac','Zgas')]) #Note the linear gasfrac to Zgas relationship

#Compare the stellar mass formed:

max(teststar$sumstars)/1e9
max(testgas$sumstars)/1e9
max(testtot$sumstars)/1e9

plot(teststar[c('starfrac','Zgas')], type='l', xlab='Star Fraction', ylab='Zgas',
xlim=c(0.5,1), ylim=c(0,0.05), col='red')
lines(testgas[c('starfrac','Zgas')], col='blue')
lines(testtot[c('starfrac','Zgas')])
curve(-mean(testtot$yield)*log(1-x), add=TRUE, lty=2) #the classic yield equation
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ stars }% use one of  RShowDoc("KEYWORDS")
