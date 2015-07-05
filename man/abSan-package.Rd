\name{abSan-package}
\alias{abSan-package}
\alias{abSan}
\docType{package}
\title{
Computes the Abreu-Sannikov repeated game algorithm
}
\description{
Code to solve repeated games using the method described in 
    Abreu & Sannikov's 2013 Theoretical Economics paper "An algorithm for 
    two-player games with perfect monitoring".
}
\details{
\tabular{ll}{
Package: \tab abSan\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-08-27\cr
License: \tab GPL-3\cr
}

Solves for the equilibrium set of values of a 2-played repeated game with monitoring. Games must be described as a list containing the actions, payoffs, and discount rate for the two players.  See the documentation for \code{model.initiate} for details about how to create model files, or one the example model files (such as \code{examples.PD}). 
The equilibrium set can be found under the default options using only the \code{abSan.eqm(modelName)} command.
Common examples available for immediate use are the Prisoner's Dilemma, Battle of the Sexes, and a Cournot duopoly.  Model files for these are stored in \code{examples.PD}, \code{examples.sexes}, and \code{examples.cournot}.  Abreu & Sannikov's arbitrary 9-action game on p.11 of their paper is also available as \code{examples.AS}.  For these examples, \code{abSan} is accurate to a minimum of 7 decimal places (and in most cases, more than 10).  A final example computes the solution to Fuchs & Lippi's independent monetary policy game in their 2006 ReStud paper.  The accuracy of this last example has not been independently verified. 

Comments and suggestions are gratefully received by the author.
}
\author{
Philip Barrett <pobarrett@uchicago.edu>
}
\references{
Abreu & Sannikov (2013) "An algorithm for two-player repeated games with perfect monitoring", Theoretical Economics. http://econtheory.org/ojs/index.php/te/article/viewForthcomingFile/1302/8114/1tion
Fuchs & Lippi (2006) "Monetary Union with Voluntary Participation", Review of Economic Studies
}
\keyword{ package }
\seealso{
\code{\link{abSan.eqm}}, \code{\link{examples.PD}}, \code{\link{examples.sexes}}, \code{\link{examples.cournot}}, \code{\link{examples.AS}}, \code{\link{examples.FL.union}}

Richard Kratzwer's Java implementation (rgsolve):  http://www.princeton.edu/~rkatzwer/rgsolve/index.html
}
\examples{
## Compute Abreu-Sannikov's example in the paper
sol <- abSan.eqm( modelName=examples.AS )
sol$status
  # Should be 1 for success
sol$vStar$mZ
  # Print the vertices of the outcome
example <- examples.AS()
  # Load the example model
abs( example$ans$mZ - sol$vStar$mZ )
  # Compute differences oto (pre-loaded) exact solution

## Compute the Cournot duopoly game
sol <- abSan.eqm( modelName=examples.cournot, modelOpts=list( 'iActs' = 15 ), charts=TRUE )
  # Plots charts for the equilibrium set and convergence in the present working directory
sol$time
  # Time taken


}
