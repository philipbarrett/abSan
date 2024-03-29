\name{examples.cournot}
\alias{examples.cournot}
\title{Cournot duopoly example}
\description{Example model definition function. Records a Cournot duopoly game and returns the a list of actions, payoffs, and a discount factor to define the game.  These can then be used by \code{model.initiate} to generate a model description used in the body of the AS algorithm.  This example also contains answers in the 'ans' list entry generated by Richard Kratzwer's rgsolve algorithm, when the default model options are passed.  These are used in accuracy tests.}
\usage{examples.cournot(opts)}
\arguments{
The list \code{opts} should be passed with the following entries:
  \item{iActs}{is the integer number of actions available to each player. Default is 15.}
}
\seealso{
  \code{\link{model.initiate}} \code{\link{examples.PD}}, \code{\link{examples.sexes}}, \code{\link{examples.AS}}, \code{\link{examples.FL.union}}
}
\examples{
model.defn <- examples.cournot()
model.defn$delta 
  # Print the discount factor
model <- model.initiate( examples.cournot )
  # Pre-process the model
sol <- abSan.eqm( model=model )
  # Solve
model <- model.initiate( examples.cournot, opts=list( 'iActs' = 80 ) )
sol <- abSan.eqm( model=model )
  # Alternative use
}