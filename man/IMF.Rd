\name{IMF}
\alias{IMF}
\alias{IMF_Chabrier}
\alias{IMF_Kroupa}
\alias{binlims}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chabrier and Kroupa Initial Mass Functions (IMFs)
}
\description{
Simple interface to differential forms of the Chabrier and Kroupa IMFs. This is useful when checking the fraction of mass returned for different IMFs etc.
}
\usage{
IMF_Chabrier(mass, alpha = 2.3, a = 0.08, b = 0.69, masslow = 0.01, massmax = 150,
  massform = 1, massmult = FALSE, rel.tol = .Machine$double.eps^0.25)

IMF_Kroupa(mass, alpha1 = 0.3, alpha2 = 1.3, alpha3 = 2.3, masslow = 0.01, mass1 = 0.08,
  mass2 = 0.5, massmax = 150, massform = 1, massmult = FALSE,
  rel.tol = .Machine$double.eps^0.25)

binlims(input, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mass}{
Numeric vector or data.frame/list; the target stellar mass to be evaluated (Msol). If a numeric vector then it will be the value of the IMF at the mass specfied. If a data.frame/list it must have components mass$lo and mass$hi (same length), and the output is the IMF integral between these mass limits (either numeric integral if \option{massmult} = FALSE, or mass integral if \option{massmult} = TRUE).
}
  \item{input}{
Numeric vector; the mid-point of the targets bins (in log or linear space depending on the \option{log} arguments).
}
  \item{alpha}{
Numeric scalar; the high mass slope (mass > 1 Msol) of the Chabrier IMF.
}
  \item{a}{
Numeric scalar; controls the gradient of the low mass slope close to mass = 1 Msol for the Chabrier IMF (lower is steeper).
}
  \item{b}{
Numeric scalar; controls the curvature of the low mass slope for the Chabrier IMF (lower is more curved).
}
  \item{masslow}{
Numeric scalar; the extreme low mass limit (Msol), where mass < masslow = 0. \option{massform} is then the mass formed between \option{masslow} and \option{massmax}.
}
  \item{massmax}{
Numeric scalar; the max mass limit (Msol), where mass > massmax = 0. \option{massform} is then the mass formed between \option{masslow} and \option{massmax}.
}
  \item{alpha1}{
Numeric scalar; the low mass slope (masslow < mass < mass1) of the Kroupa IMF.
}
  \item{alpha2}{
Numeric scalar; the mid mass slope (mass1 < mass < mass2) of the Kroupa IMF.
}
  \item{alpha3}{
Numeric scalar; the high mass slope (mass2 < mass < massmax) of the Kroupa IMF.
}
  \item{mass1}{
Numeric scalar; the low mass limit (Msol), see \option{alpha1}.
}
  \item{mass2}{
Numeric scalar; the mid mass limit (Msol), see \option{alpha2}.
}
  \item{massform}{
Numeric scalar; the mass formed between \option{masslow} and \option{massmax}, i.e. what the IMF should integrate to (Msol). Since the default is 1, this means the IMF returned is a true PDF and can be used as such.
}
  \item{massmult}{
Should the returns value be dN/dM (FALSE) or M.dN/dM (TRUE)? dN/dM is useful when evaluating the number of stars formed within limits, whilst M.dN/dM is used for evaluating the mass of stars formed within limits.
}
  \item{rel.tol}{
Numeric scaler; relative accuracy requested. See \code{\link{integrate}}.
}
  \item{log}{
Logical; should the bins be evenly split in log or linear space?
}
}
\details{
Mostly useful for exploring the fraction of mass in different stellar mass ranges when determining feedback etc. See Examples.
}
\value{
For the IMF functions, if \option{massmult} = FALSE then the return is dN/dM for the specified \option{mass}, if \option{massmult} = TRUE then the return is M.dN/dM.

For \code{binlims} the output is a list with elements 'lo' and 'hi'. This can be uses as an input for the IMF functions for the \option{mass} arguments.
}
\references{
Chabrier, 2003, PASP, 115, 763 \cr
Kroupa, 2001, MNRAS, 322, 231
}
\author{
Aaron Robotham
}

\examples{
curve(IMF_Chabrier, 0.01, 150, log='xy', xlab='Mass / Msol', ylab='PDF')
curve(IMF_Kroupa, 0.01, 150, add=TRUE, col='red')

#Check the fraction of mass that will go SNII (around 20\%) and rapidly enrich:

integrate(IMF_Chabrier, lower=10, upper=150, massmult=TRUE)$value

#If we form a cloud with 10^6 Msol, how many stars more massive than Msol will be formed?

integrate(IMF_Chabrier, lower=1, upper=150, massform=1e6, massmult=FALSE)$value

#And how much mass will there be in these stars?:

integrate(IMF_Chabrier, lower=1, upper=150, massform=1e6, massmult=TRUE)$value

#How many stars in total?

integrate(IMF_Chabrier, lower=0.01, upper=150, massform=1e6, massmult=FALSE)$value

#So numerically 6\% of the stars will be more massive than Msol, but more than half the
#mass will be in these stars! This shows why rare processes in a massive stars can matter
#so much.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ IMF }% use one of  RShowDoc("KEYWORDS")

