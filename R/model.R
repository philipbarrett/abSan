#######################################################################################
# Function to describe the model in the Abreu-Sannikov model                          #
# Philip Barrett, Chicago                                                             #
# Created: 30jul2013                                                                  #
#######################################################################################

model.initiate <- function( defn.fn, opts=NULL ){
# Initiates the model. Mostly rearranging elements of model.define into more
# useful but less intuitive formats
  
  # 0. Utilities
  maxDev <- function(  jAction, iPlayer, jPlayer ){
  # Computes the best deviating payoff for player i given player j plays action j
    max( mF[ ( mA[, jPlayer] == jAction ), iPlayer ] )
  }
  
  # 1. Body
  defn <- defn.fn( opts )
        # Model definition
  mA <- cbind( rep( 1:( defn$iActions[ 1 ] ), defn$iActions[ 2 ] ), 
               rep( 1:( defn$iActions[ 2 ] ), each = defn$iActions[ 1 ] ) )
        # Matrix of actions for the two players
  mF <- cbind( c( defn$payoffs1 ), c( defn$payoffs2 ) )
        # Matrix of stage payoffs
  devPay <- list( sapply( 1:( defn$iActions[ 1 ] ), maxDev, 1, 2 ), 
                  sapply( 1:( defn$iActions[ 2 ] ), maxDev, 2, 1 ) )
        # The max deviating payoff for each player across alternative actions.
        # NB: This is a list rather than a matrix because players may have
        # different numbers of actions
  mH <- cbind( devPay[[ 1 ]][ mA[ , 2 ] ], devPay[[ 2 ]][ mA[ , 1 ] ] ) - mF
        # The best improvement from deviating at a given action
  minmax <- c( min( apply( defn$payoffs1, 2, max ) ), min( apply( defn$payoffs2, 1, max ) ) )
        # The minmax payoffs of the two players
  
  return( list( 'iActions'=defn$iActions, 'mA'=mA, 'mF'=mF, 'mH'=mH, 'delta'=defn$delta, 
                'iJointActs' = prod( defn$iActions ), 'minmax'=minmax ) )
}