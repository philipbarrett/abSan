#######################################################################################
# Test functions for abreuSannikov source file                                        #
# Philip Barrett, Chicago                                                             #
# Created: 08apr2013                                                                  #
#######################################################################################


##** 1. abSan.C **##

test.abSan <- function(){
# Test abSan.C and abSan.R
  
  # 1. Set up
  model <- model.initiate( examples.PD )
  set <- sets.setZ( model$mF )
  pun <- sets.P( set ) # model$minmax
        # Initializing the model, set, and punishment vector
  
  # 2. Expected answers for simple PD example - Calculated by hand
  lAns <- list( c( 4, 4 ),
               matrix( c( 0, 6, 5, .5,  5.75, .5 ), nrow=3, ncol=2, byrow=T ),
               matrix( c( .5, 5, .5, 5.75, 6, 0 ), nrow=3, ncol=2, byrow=T ),
               c( 2, 2 ) )
        # Values of C for all actions
  union.mZ <- ( 1 - model$delta ) * model$mF[ c( 1, 2, 2, 2, 3, 3, 3, 4 ), ] + 
    model$delta * do.call( 'rbind', lAns )
  ans.R <- sets.setZ( union.mZ )
        # Values of R for all actions
    
  # 3. Run tests
  for ( iII in 1:model$iJointActs )
    checkTrue( all( abs( abSan.C( iII, set, pun, model ) - lAns[[ iII ]] ) < 1e-14 ) )
  checkEquals( abSan.R( set, pun, model)$mZ, ans.R$mZ )
  checkEquals( abSan.R( set, pun, model)$mG, ans.R$mG )
  checkEquals( abSan.R( set, pun, model)$vC, ans.R$vC )
}