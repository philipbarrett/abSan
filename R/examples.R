#######################################################################################
# Example model descriptions for the Abreu-Sannikov method                            #
# Philip Barrett, Chicago                                                             #
# Created: 30jul2013                                                                  #
#######################################################################################

examples.AS <- function( opts=NULL ){
  # Records and returns the payoff matrices for the two players
  
  # 1. Modify these lines to define the game
  iActions <- c( 3, 3)
  # Number of actions for each player
  payoffs1 <- matrix( c( 16,  3,  0, #1,
                         21, 10, -1,
                         9,  5, -5 ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=T )
  payoffs2 <- matrix( c(  9, 13,  3,
                          1,  4,  0,
                          0, -4, -15 ),
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=T )
  # Payoffs for each player from the stage game
  delta <- 0.4
  # Discount factor
  
  mZ <- matrix( c( .43169727312472, 2.4, 6.17267888662042, 0, 10.15, 0, 20.125, 2.4, 
                   16, 9, 6.17267888662042, 12.0237911118091 ), ncol=2, byrow=T )
        # Taken from Abreu & Sannikov's paper
  ans <- sets.setZ( mZ )
        # The answer (from AS paper p.19)
  
  return( list( 'iActions'=iActions , 'payoffs1'=payoffs1, 'payoffs2'=payoffs2,
                'delta'=delta, 'ans'=ans ) )
  # Return these definitions
}

examples.PD <- function( opts=NULL ){
# An exceptionally simple prinsoner's dilemma model
  
  # 1. Modify these lines to define the game
  iActions <- c( 2, 2)
        # Number of actions for each player
  payoffs1 <- matrix( c( 4, 0,
                         6, 2 ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=T )
  payoffs2 <- matrix( c(  4, 6,
                          0, 2 ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=T )
        # Payoffs for each player from the stage game
  delta <- 0.8
        # Discount factor
  mZ <- matrix( c( 5, 2, 4, 4, 2, 5, 2, 2 ), ncol=2, byrow=T )
        # Taken from rgsolve
  ans <- sets.setZ( mZ )
        # The answer
  return( list( 'iActions'=iActions , 'payoffs1'=payoffs1, 'payoffs2'=payoffs2, 
                'delta'=delta, 'ans'=ans ) )
        # Return these definitions
}

examples.sexes <- function( opts=NULL ){
  # Battle of the sexes example
  
  # 1. Modify these lines to define the game
  iActions <- c( 2, 2)
        # Number of actions for each player
  payoffs1 <- matrix( c( 8, 3,
                         0, 5 ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=T )
  payoffs2 <- matrix( c(  5, 3,
                          0, 8 ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=T )
  # Payoffs for each player from the stage game
  delta <- 0.8
        # Discount factor
  mZ <- matrix( c( 8, 5, 5, 5, 5, 8 ), ncol=2, byrow=T )
        # Taken from rgsolve
  ans <- sets.setZ( mZ )
        # The answer
  return( list( 'iActions'=iActions , 'payoffs1'=payoffs1, 'payoffs2'=payoffs2, 
                'delta'=delta, 'ans'=ans ) )
        # Return these definitions
}

examples.cournot <- function( opts=NULL ){
  # Cournot duopoly example
  
  # 1. Modify these lines to define the game
  iActs <- if( is.null( opts$iActs ) ) 15 else opts$iActs
         # Discretization
  cost <- if( is.null( opts$cost ) ) c(.6, .6) else opts$cost
  maxOut <- 6
        # Maximum output
  iActions <- c( iActs, iActs)
        # Number of actions for each player
  make.payoff <- function( q, iPlayer ){
  # Computes the payoff for player iPlayer given output vector q
    price <- max( c( maxOut - sum(q), 0 ) )
    return( q[iPlayer] * ( price - cost[ iPlayer ] ) )
  }
  vActs <- seq( 0, 6, length.out=iActs )
        # Quantities for each player
  mActs <- cbind( rep( vActs, iActs ), rep( vActs, each=iActs ) )
        # List of joint strategies
  payoffs1 <- matrix( c( sapply( 1:iActs^2, function(i) make.payoff( mActs[ i, ], 1 ) ) ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=F )
  payoffs2 <- matrix( c( sapply( 1:iActs^2, function(i) make.payoff( mActs[ i, ], 2 ) ) ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=F )
        # Payoffs for each player from the stage game
  delta <- 0.8
        # Discount factor
  mZ <- matrix( c( 6.9665306122448980, 0.2865306122448980, 6.8767346938775520, 0.3967346938775508,
                   0.3967346938775508, 6.8767346938775520, 0.2865306122448980, 6.9665306122448980,
                   0.1910204081632653, 6.8974149659863940, 0.1175510204081632, 6.5334170591313450,
                   0.0000000000000000, 4.9390685504971220, 0.0000000000000000, 0.0000000000000000,
                   4.9390685504971220, 0.0000000000000000, 6.5334170591313450, 0.1175510204081632,
                   6.8974149659863940, 0.1910204081632653 ), ncol=2, byrow=T )
  ans <- sets.setZ( mZ )
        # Taken from rgsolve
  return( list( 'iActions'=iActions , 'payoffs1'=payoffs1, 'payoffs2'=payoffs2, 'delta'=delta, 'ans'=ans ) )
        # Return these definitions
}

<<<<<<< HEAD
examples.cournot.CES <- function( opts=NULL ){
# Cournot duopoly example with CES products
  
  # 1. Modify these lines to define the game
  iActs <- if( is.null( opts$iActs ) ) 15 else opts$iActs
    # Discretization
  cost <- if( is.null( opts$cost ) ) c(.6, .6) else opts$cost
  income <- 1
    # Household income
  share <- if( is.null( opts$share ) ) .5 else opts$share
  v.share <- c( share, 1 - share )
  eps <- if( is.null( opts$eps ) ) 1 else opts$eps
    # Household utililty parameters
  iActions <- c( iActs, iActs)
    # Number of actions for each player
  make.payoff <- function( q, iPlayer ){
  # Computes the payoff for player iPlayer given output vector q
    denom <-  ( v.share[2] ^ ( -1 / eps ) * q[1] ^ ( 1 - 1 / eps ) +
                 v.share[1] ^ ( -1 / eps ) * q[2] ^ ( 1 - 1 / eps ) )
    price <- income / denom * v.share[ -iPlayer ] * q[ - iPlayer ]
    return( q[iPlayer] * ( price - cost[ iPlayer ] ) )
  }
  vActs <- seq( .1, 6, length.out=iActs )
  # Quantities for each player
  mActs <- cbind( rep( vActs, iActs ), rep( vActs, each=iActs ) )
  # List of joint strategies
=======
examples.cournot.subn <- function( opts=NULL ){
  # Cournot duopoly example
  
  # 1. Modify these lines to define the game
  iActs <- if( is.null( opts$iActs ) ) 15 else opts$iActs
      # Discretization
  v.cost <- if( is.null( opts$v.cost ) ) c(.1, .1) else opts$v.cost
      # Marginal cost
  p.max <- if( is.null( opts$p.max ) ) 6 else opts$p.max
      # Maximum price
  p.min <- if( is.null( opts$p.min ) ) 0.1 else opts$p.min
  p.min <- max( c( v.cost, c( p.min, p.min ) ) )
      # Minimum price
  income <- if( is.null( opts$income ) ) 1 else opts$income
      # Consumers' nominal income
  v.share <- if( is.null( opts$v.share ) ) c(.5, .5) else opts$v.share
      # The vector of market shares
  elas <- if( is.null( opts$elas ) ) 10 else opts$elas
      # The elasticity of substitution
  iActions <- c( iActs, iActs)
      # Number of actions for each player
  make.payoff <- function( p, iPlayer ){
  # Computes the payoff for player iPlayer given price vector p
    p.level <- ( v.share %*% p ^ ( 1 - elas ) ) ^ ( 1 / ( 1 - elas ) )
        # The aggregate price level
    profit <- v.share[ iPlayer ] * ( p[ iPlayer ] - v.cost[ iPlayer ] ) * income / p.level *
      ( p[ iPlayer ] / p.level ) ^ ( - elas )
    return( profit )
  }
  vActs <- seq( p.min, p.max, length.out=iActs )
      # Quantities for each player
  mActs <- cbind( rep( vActs, iActs ), rep( vActs, each=iActs ) )
      # List of joint strategies
>>>>>>> ef46f19489b74d7f9a0702d745f86ef7f18bcffd
  payoffs1 <- matrix( c( sapply( 1:iActs^2, function(i) make.payoff( mActs[ i, ], 1 ) ) ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=F )
  payoffs2 <- matrix( c( sapply( 1:iActs^2, function(i) make.payoff( mActs[ i, ], 2 ) ) ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=F )
<<<<<<< HEAD
  # Payoffs for each player from the stage game
  delta <- 0.95
  # Discount factor
#   mZ <- matrix( c( 6.9665306122448980, 0.2865306122448980, 6.8767346938775520, 0.3967346938775508,
#                    0.3967346938775508, 6.8767346938775520, 0.2865306122448980, 6.9665306122448980,
#                    0.1910204081632653, 6.8974149659863940, 0.1175510204081632, 6.5334170591313450,
#                    0.0000000000000000, 4.9390685504971220, 0.0000000000000000, 0.0000000000000000,
#                    4.9390685504971220, 0.0000000000000000, 6.5334170591313450, 0.1175510204081632,
#                    6.8974149659863940, 0.1910204081632653 ), ncol=2, byrow=T )
#   ans <- sets.setZ( mZ )
  # Taken from rgsolve
  return( list( 'iActions'=iActions , 'payoffs1'=payoffs1, 'payoffs2'=payoffs2, 'delta'=delta ) )
#           , 'ans'=ans ) )
  # Return these definitions
=======
      # Payoffs for each player from the stage game
  delta <- if( is.null( opts$delta ) ) .9 else opts$delta
      # Discount factor
  return( list( 'iActions'=iActions , 'payoffs1'=payoffs1, 'payoffs2'=payoffs2, 'delta'=delta ) )
      # Return these definitions
>>>>>>> ef46f19489b74d7f9a0702d745f86ef7f18bcffd
}

examples.FL.union <- function( opts=NULL ){
# Exaple similar to Fuchs & Lippi's monetary policy game
  
  # 1. Modify these lines to define the game
  iActs <- if( is.null( opts$iActs ) ) 40 else opts$iActs
        # Discretization
  piRange <- if( is.null( opts$piRange ) ) c( -5, 5 ) else opts$piRange 
        # Maximum output
  iActions <- c( iActs, iActs)
        # Number of actions for each player
  alpha <- 2
  
  make.payoff <- function( pi, piOther ){
    # Computes the payoff for any when they play pi and the other player plays piOther
    return( - .5 * pi^2 + alpha * ( pi - piOther ) )
  }
  
  vActs <- seq( piRange[ 1 ], piRange[ 2 ], length.out=iActs )
        # Quantities for each player
  mActs <- cbind( rep( vActs, iActs ), rep( vActs, each=iActs ) )
        # List of joint strategies
  payoffs1 <- matrix( c( sapply( 1:iActs^2, function(i) make.payoff( mActs[ i, 1 ], mActs[i, 2 ] ) ) ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=F )
  payoffs2 <- matrix( c( sapply( 1:iActs^2, function(i) make.payoff( mActs[ i, 2 ], mActs[i, 1 ] ) ) ), 
                      nrow=iActions[ 1 ], ncol=iActions[ 2 ], byrow=F )
        # Payoffs for each player from the stage game
  delta <- if( is.null( opts$delta ) ) 0.1 else opts$delta
        # Discount factor
  return( list( 'iActions'=iActions , 'payoffs1'=payoffs1, 'payoffs2'=payoffs2, 'delta'=delta) )
#           , 'ans'=ans ) )
  # Return these definitions
  
}