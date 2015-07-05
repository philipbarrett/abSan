\name{abSan.eqm}
\alias{abSan.eqm}
\title{Compute the set of equilibrium values}
\description{Calculates the set of values of all subgame-perfect equilibria for a repeated game.}
\usage{abSan.eqm(modelName, model, set, pun, tol=1e-12, charts=FALSE, maxIt=150, detail=FALSE,
    print.output=TRUE, save.solution=FALSE, par=FALSE, cluster, modelOpts)}

\arguments{
  \item{modelName}{the model definition function.}
  \item{model}{a post-processed model list. Usually created as the output of \code{model.initiate(modelName)} for some model. Useful when re-computing the same model with different solution options, as prevents re-defining the model each time.}
  \item{set}{an (m,2) matrix of vertices of the initial set of continuation values. Must wholly contain the equilibrium set. If missing or \code{NULL}, the set of stage payoffs is used.}
  \item{pun}{a length-2 vector of initial punishments. Must be less than the equilibrium punishment. If missing or \code{NULL}, the minmax payoff for the stage game is used.}
  \item{tol}{numeric convergence tolerance. The threshold Hausdorff difference between successive sets at which the algorithm terminates.}
  \item{charts}{Boolean flag for saving charts. Currently only equilibrium set and the sequence of convergent sets are created. These are saved as equilibrium.pdf, and convergence.pdf in the current working directory.}
  \item{maxIt}{The maximum number of iterations.}
  \item{detail}{Boolean flag for retaining information about the iterations}
  \item{par}{Boolean flag for using multicore execution.}
  \item{print.output}{Boolean flag for output display.}
  \item{save.solution}{Currently inactive.}
  \item{cluster}{number of nodes to use in cluster execution. Currently inactive.}
  \item{modelOpt}{options to pass to the model.}
}

\value{
  Returns a list of solution components:
  \item{status}{is \code{1} if solution converges to required tolerance.  Otherwise \code{0}.}
  \item{vStar}{a list containing a description of the equilibrium set. \code{vStar$mZ} is a 2-column matrix of vertices of the equilibrium set.  \code{vStar$mG} is a 2-column matrix of unit-magnitude normals to the boundary of the equilibrium set, and \code{vStar$vC} is a vector of intercepts to those normals, such that for every row \eqn{g} of \code{vStar$mG} and each entry \eqn{c} of \code{vStar$vC} the conditions \eqn{g.z\le c} are a necessary and sufficient condition for \eqn{z} to be in the equilibrium set. }
  \item{vBar}{the equilibrium punishment}
  \item{iterations}{The number of iterations to solution}
  \item{lSet}{returned only if \code{details=TRUE}. The list of sets at each iteration of the algorithm. }
  \item{lPun}{returned only if \code{details=TRUE}. The list of punishments at each iteration of the algorithm. }
  \item{time}{time to compute the equilibrium set.}
}
\seealso{
  \code{\link{model.initiate}}
}
\examples{
## Compute the Cournot duopoly game
sol <- abSan.eqm( modelName=examples.cournot, modelOpts=list( 'iActs' = 15 ) )
  # Benchmark solution
sol$time
  # Time taken
sol <- abSan.eqm( modelName=examples.cournot, modelOpts=list( 'iActs' = 15 ), print.output=FALSE )
  # Turn off output
sol$time
  # Time taken
sol <- abSan.eqm( modelName=examples.cournot, modelOpts=list( 'iActs' = 40 ) )
  # Using model options to compute a finer discretization
sol$time
  # Time taken
sol <- abSan.eqm( modelName=examples.cournot, modelOpts=list( 'iActs' = 40 ), par=TRUE )
  # Using multicore execution to speed up larger example
sol$time
  # Time taken
}
