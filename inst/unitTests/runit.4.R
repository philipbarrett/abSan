#######################################################################################
# Test functions for game manipulation                                                #
# Philip Barrett, Chicago                                                             #
# Created: 09aug2013                                                                  #
#######################################################################################

# 1. Prisoner's dilemma
test.game <- function(){
  
  DEACTIVATED('Deactivating this test function')
  # 0. Set up a test equilibrium set
  example <- examples.PD()
        # The example model.  Answer included from other source (eg. rgsolve)
  sol <- abSan.eqm( examples.PD, print.output=TRUE, tol=1e-15, maxIt=200, set=example$ans )
        # Compute the equilibrium set using the known answer as a starting point
  
  # 1. Test the efficient point
  vEfficient <- game.efficient( sol$vStar, 1, 3 )
        # Compute the efficient point
  checkEquals( vEfficient, c( 3, 4.5 ) )
        # Check that that the answer is as expected
  vEfficient <- game.efficient( sol$vStar, 2, 3 )
        # Compute the efficient point
  checkEquals( vEfficient, c( 4.5, 3 ) )
        # And for player 2
  vEfficient <- game.efficient( sol$vStar, 1, 2 )
        #  An example with player 1's minimum payoff - line intersect is
        #  parallel here
  checkEquals( vEfficient, c( 2, 5 ) )
        # Check answer
  vEfficient <- game.efficient( sol$vStar, 1, 1 )
        #  Should be no solution - worse than min payoff
  checkEquals( vEfficient, NULL )
        # Check answer
  vEfficient <- game.efficient( sol$vStar, 2, 6 )
        #  Should be no solution - better than max payoff
  checkEquals( vEfficient, NULL )
        # Check answer
  
  # 2. game.vertex.cont
  model.PD <- model.initiate( examples.PD )
        # Initialize the prisoner's dilemma model
  best <- game.vertex.cont( c( 4, 4 ), sol$vStar, sol$vBar, model.PD )
        # Continuation vertex supported by (4,4)
  checkEquals( best$stage, c(4,4) )
  checkEquals( best$cont, c(4,4) )
        # Supports best outcome
  worst <- game.vertex.cont( c( 2, 2 ), sol$vStar, sol$vBar, model.PD )
        # Continuation vertex supported by (2,2)
  checkEquals( worst$stage, c(2,2) )
  checkEquals( worst$cont, c(2,2) )
        # Supports worst outcome
  extreme1 <- game.vertex.cont( c( 2, 5 ), sol$vStar, sol$vBar, model.PD )
        # Continuation vertex supported by (2,5)
  checkEquals( extreme1$stage, c(0,6) )
  checkEquals( extreme1$cont, c(2.5,4.75) )
        # Supports other extreme
  extreme2 <- game.vertex.cont( c( 5, 2 ), sol$vStar, sol$vBar, model.PD )
        # Continuation vertex supported by (2,5)
  checkEquals( extreme2$stage, c(6,0) )
  checkEquals( extreme2$cont, c(4.75,2.5) )
        # Supports other extreme
  
  # 3. game.path.vertex
  game.path.vertex( sol, examples.PD, c( 1, 0, 0, 0 ) )
  
}