\name{model.initiate}
\alias{model.initiate}
\title{Processes a model for use in the algorithm}
\description{Processes a model definition function in two ways.  First, it repackages the payoffs and actions for easier use in the rest of the code.  Second, it computes the other properties of a game, such as the marginal gain from deviating, which the Abreu-Sannikov code uses. These are then returned as a list for use in other functions, such as, \code{abSan.eqm}.}
\usage{model.initiate(defn.fn, opts)}
\arguments{
  \item{defn.fn}{the model definition function.  This must return a list including the following four entries: \code{iActions}, the 2-vector of the number actions for each player; \code{payoffs1}, the normal-form matrix of payoffs of the stage game for player 1; \code{payoffs2}, the normal-form matrix of payoffs of the stage game for player 1; and \code{delta} the discount factor. See, for instance, \code{examples.PD} in this package. }
  \item{opts}{a list of options to pass to the model definition function \code{defn.fn}.}
}
\value{
  Returns a list of components used in the Abreu-Sannikov algorithm:
  \item{iActions}{2-vector of the number actions for each player.}
  \item{mA}{a \eqn{(2,m)}-matrix, where \eqn{m} is the number of *joint* actions. Each row contains the action indices for the two players associated with that joint action.}
  \item{mF}{payoffs associated with the joint actions listed in mA.}
  \item{mH}{the gain from deviating from each joint action for the two players. Abreu-Sannikov's h(a) function.}
  \item{delta}{the discount factor.}
  \item{iJointActs}{the number of joint actions.}
  \item{minmax}{the minmax values of the stage game.}
}
\seealso{
  \code{\link{examples.PD}}, \code{\link{examples.sexes}}, \code{\link{examples.cournot}}, \code{\link{examples.AS}}, \code{\link{examples.FL.union}}
}
\examples{
model <- model.initiate( examples.FL.union, opts=list( 'iActs' = 80, piRange=c(-8,8) ) )
  # Initiate a model with a large number of actions
}