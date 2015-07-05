#######################################################################################
# Functions to run the main the Abreu-Sannikov algorithm                              #
# Philip Barrett, Chicago                                                             #
# Created: 31jul2013                                                                  #
#######################################################################################

abSan.R <- function( set, vPun, model, print.output=TRUE, par=FALSE ){
# Computes Abreu & Sannikov's R(W,u) object. Inputs are: set, a set of 
# continuation values; vPun is a pair of punishment values for the two players; 
# and model is a model description created by model.initiate. For each action we
# compute C(a,W,u) and add to the appropriate payoffs in g(a). Finally, we take
# a union of all the set of all the actions.
    
  mSets <- NULL
        # Container for the sets of the form (for each a):
        #     ( 1 - delta ) * g(a) + delta * C( a, W, u )
  
  # 1. No parallelization
  if ( par == FALSE ){
    mSets <- abSan.rCore(set, vPun, model )
          # Insert into the list of sets
  }
  
  # 2. Multicore parallelism
  if ( par == TRUE ){
    lsum.mZ <- mclapply( 1:(model$iJointActs), abSan.rComponent, set, vPun, model, print.output )
          # Calculate the points in the sum
    mSets <- sets.bind( lsum.mZ )
          # Insert into the list of sets
  }
  
  # 3. Cluster execution
  if ( class( par ) == 'numeric' ){
    lsum.mZ <- parLapply(cl, 1:(model$iJointActs), abSan.rComponent, set, vPun, model, print.output )
          # Calculate the points in the sum
    mSets <- do.call( rbind, lsum.mZ )
          # Insret into the list of sets
  }
  return( sets.setZ( mSets ) )
        # Return the union over all sets
}

abSan.uUpdate <- function( set, vPun ){
# Computes the new punishment given the updated set and the previous 
  newPun <- sets.P( set )
        # Punishment from the new set
  return( apply( rbind( newPun, vPun ), 2, max ) )
}

abSan.eqm <- function( modelName=NULL, model=NULL, set=NULL, pun=NULL, tol=1e-12, charts=FALSE, 
                         maxIt=150, detail=FALSE, print.output=TRUE, save.solution=FALSE,
                         par=FALSE, cluster=NULL, modelOpts=NULL ){
  
  # 0. Set up
  ptm <- proc.time()
  if (is.null( model ) ) model <- model.initiate( modelName, modelOpts )
  library( parallel )
  if ( is.null( set ) ) set <- sets.setZ( model$mF )
  if ( is.null( pun ) ) pun <- model$minmax # sets.P( set ) #
  lSet <- list( set$mZ )
  lPun <- list( pun )
  
  cl <- cluster
  if( !is.null( cluster ) ) par <- 1
        # To make sure that parallel execution is used in subsequent code
  if ( ( is.null( cl ) ) & ( class( par ) == 'numeric' ) ){
    message('Iniitializing cluster')
    cl <<- makeCluster( par , type = 'PSOCK' )   
  }
  if ( ( !is.null( cluster ) ) | ( class( par ) == 'numeric' ) ){
    message('Sourcing to cluster')
    clusterCall( cluster , function() { source('src/abreuSannikov.R') ; 
                                  source('src/sets.R'); NULL } )
  }

  
  # 1. Main Body
  iCounter <- 0
  dist <- 2 * tol
        # Initialize the counter and distance
  while( ( dist > tol ) & ( iCounter < maxIt ) ){
    iCounter <- iCounter + 1
          # Increment counter
    set <- abSan.R( set, pun, model, print.output, par )
    pun <- abSan.uUpdate( set, pun )
    lSet[[ iCounter + 1 ]] <- set$mZ
    lPun[[ iCounter + 1 ]] <- pun
    
    dist <- sets.Hausdorff( lSet[[ iCounter ]], set$mZ )
    if ( print.output ) message('Iteration ', iCounter, 
                                ': Hausdorff distance = ', format( dist, digits=3 ) )
  }
  if ( dist < tol ) status <- 1 else status <- 0
  
  # 2. Charting
  if ( charts ) abSan.charts( set, lSet )
    # Call charting function here
  
  # 3. Save solution
  
  
  # 4. Return output
  if ( ( is.null(cluster==NULL) ) & ( class( par ) == 'numeric' ) ) stopCluster(cl)
  if ( detail )
    return( list( 'status'=status, 'vStar'=set, 'vBar'=pun, 'iterations'=iCounter, 
                  'lSet'=lSet, 'lPun'=lPun, 'time'= proc.time() - ptm ) )
  return( list( 'status'=status, 'vStar'=set, 'vBar'=pun, 'iterations'=iCounter,
                'time'= proc.time() - ptm ) )
}

abSan.charts <- function( set, lSet ){
  
  iSets <- length( lSet )
        # Number of converging sets
  
  # 1. The equilibrium set #
  pdf('equilibrium.pdf')
        # Initiate the output
  abSan.plotSet( set$mZ, col='black', lwd=2 )
        # Plot the equilibrium set
  dev.off()
        # Suppress further output
  
  pdf('convergence.pdf')
        # Initiate the output
  abSan.plotSet( lSet[[ 1 ]], col='blue', lwd=.5 )
        # Plot the ifrst iteration
  for ( iII in 2:iSets )
    abSan.addSet( lSet[[ iII ]], col='blue', lwd=.5 )
        # Plot the ifrst iteration
  abSan.addSet( set$mZ, col='black', lwd=2 )
  dev.off()
        # Suppress further output
  
}

abSan.addSet <- function( mZ, ... ){
# Adds a set to a chart
  lines( c( mZ[ , 1 ], mZ[ 1, 1 ] ), c( mZ[ , 2 ], mZ[ 1, 2 ] ), ... )
}

abSan.plotSet <- function( mZ, ... ){
# Plots set on a new chart
  plot( c( mZ[ , 1 ], mZ[ 1, 1 ] ), c( mZ[ , 2 ], mZ[ 1, 2 ] ), type='l', 
        xlab='Player 1 payoff', ylab='Payer 2 payoff', ... )
  
}
