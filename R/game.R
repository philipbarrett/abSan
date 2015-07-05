#######################################################################################
# Functions to manipulate the equilibrium set after it has been computed              #
# Philip Barrett, Chicago                                                             #
# Created: 11aug2013                                                                  #
#######################################################################################


game.efficient <- function( set, iPlayer, val ){
# Returns the efficient point for the game with equilibrium set given by set,
# when player iPlayer has value val
  
  vNorm <- if ( iPlayer == 1 ) c( 1, 0 ) else c( 0, 1 )
        # Normal to the intersetcion line
  mExtreme <- sets.lineIntersect( set, vNorm, val )
        # Compute the intersections of the line with the set
  if ( is.null(mExtreme) )
    return( NULL )
          # Return NULL if no answer
  else
    return( apply( mExtreme, 2, max ) )
          # Return the one with the best outcome
}

game.vertex.cont <- function( vV, set, vPun, model ){
# Computes an IC action-continuation pair with value vVal
# Currently *only* works for pure-strategy equilibria
  
  delta <- model$delta
        # Useful rewrite
  for ( iAct in 1:(model$iJointActs) ){
    vW = 1 / delta * ( vV - ( 1 - delta ) * model$mF[ iAct, ] )
          # Continuation value associated with the action
    boContInSet <- sets.inside.pt( vW, set$mG, set$vC )
          # Flag for vW inside the set
    boIC <- ( delta * ( vW - vPun ) >= ( 1 - delta ) * model$mH[ iAct, ] )
          # Flag for incentive compatibility
    if( all( c( boContInSet, boIC ) ) )
      return( list( 'iAct'=iAct, 'stage'=model$mF[ iAct, ], 'cont'=vW ) )
  }
  stop('No IC action-continuation pair found for value ( ', vV[1], ', ', vV[2], ' )' )
}

game.vertex.support <- function( set, vPun, model, ... ){
# Reports the continuing vertex supported by each action of the game
  
  return( lapply( 1:(model$iJointActs), game.vertex.support.pt, set, vPun, model, ...  ) )
        # The vertices supported by each action
}

game.vertex.support.pt <- function( iAct, set, ... ){
# Reports the action and continuation supporting a given action

  mSupports <- abSan.rComponent( iAct, set, ... )
        # Compute the extremal points supported by this action
  boIdx <- apply( mSupports, 1, function(x) any( apply( set$mZ, 1, function(y) identical( x, y ) ) ) )
        # The indices of mSupports which are extremal points of set
  return( mSupports[ boIdx, ] )
}

game.path.vertex <- function( sol, modelName, initVal, iPds=20, ... ){
# Computes the sequence of actions, payoffs, and continuation values given an 
# equilibrium set and model and a distribution over the vertices of the
# equilibrium set
  
  # 1. Calculate the actions supporting each vertex
  model <- model.initiate( modelName )
        # Create the model object from the model name
  mVals <- matrix( initVal, nrow=1, ncol=2 )
  vActs <- NULL
  mPayoffs <- NULL
        # Initialize the sequence of values, actions, and payoffs
  
  #2. Loop over each period
  for ( iNN in 1:iPds ){
    thisDecomp <- game.vertex.cont( mVals[ iNN, ], sol$vStar, sol$vBar, model )
          # Decompose the payoff into an action and continuation
    mVals <- rbind( mVals, thisDecomp$cont )
    vActs <- c( vActs, thisDecomp$iAct )
    mPayoffs <- rbind( mPayoffs, thisDecomp$stage )
  }
  return( list( 'mVals'=mVals, 'vActs'=vActs, 'mPayoffs'=mPayoffs ) )
}

# game.path.dist <- function( vPt, set ){
# # Decomposes a continuation pair into a distribution over (at most four) extreme
# # points of the equilibrium set
#   
#   dist <- c( 1, 1 ) %*% vPt
#         # Projecting vPTs onto (1,1)
#   mBoundary <- sets.lineIntersect( set, c( 1, 1 ), dist )
#         # Computing the extremal associated with vPt
#   if( nrow( mBoundary ) == 1 ) mBoundary <- rbind( mBoundary, mBoundary )
#         # Special case of a vertex
#   prob <- ( vPt - mBoundary[ 2, ] ) / ( mBoundary[ 1, ] - mBoundary[ 2, ] )
#         # The randomization probability
#   if 
# }
