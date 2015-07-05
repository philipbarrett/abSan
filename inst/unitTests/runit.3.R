#######################################################################################
# Test functions for example cases                                                    #
# Philip Barrett, Chicago                                                             #
# Created: 09aug2013                                                                  #
#######################################################################################

# 1. Prisoner's dilemma
test.PD <- function(){
  message('TESTING: Prisoner\'s dilemma')
  tol <- 1e-14
        # Error tolerance on test
  example <- examples.PD()
        # The example model.  Answer included from other source (eg. rgsolve)
  sol <- abSan.eqm( examples.PD, print.output=TRUE, tol=1e-15, maxIt=200, par=FALSE )
        # abSan solution
  checkTrue( all( abs( example$ans$mZ - sol$vStar$mZ ) < tol ) )
}

# 2. Battle of the sexes
test.sexes <- function(){
  message('TESTING: Battle of the sexes')
  tol <- 1e-14
        # Error tolerance on test
  example <- examples.sexes()
        # The example model.  Answer included from other source (eg. rgsolve)
  sol <- abSan.eqm( examples.sexes, print.output=TRUE, tol=1e-20, maxIt=200, par=FALSE )
        # abSan solution
  checkTrue( all( abs( example$ans$mZ - sol$vStar$mZ ) < tol ) )
}

# 3. Abreu & Sannikov's example
test.AS <- function(){
  message('TESTING: Abreu-Sannikov arbitrary example')
  tol <- 1e-07
        # Error tolerance on test
  example <- examples.AS()
        # The example model.  Answer included from other source (eg. rgsolve)
  sol <- abSan.eqm( examples.AS, print.output=TRUE, tol=1e-24, maxIt=100, par=FALSE )
        # abSan solution
  checkTrue( all( abs( example$ans$mZ - sol$vStar$mZ ) < tol ) )
}

# 4. Judd, Yeltekin & Conklin's Cournot example
test.cournot <- function(){
#   DEACTIVATED('Deactivating this test function')
  tol <- 1e-10
        # Error tolerance on test
  example <- examples.cournot()
        # The example model.  Answer included from other source (eg. rgsolve)
  sol <- abSan.eqm( modelName=examples.cournot, print.output=TRUE, tol=1e-12,
                      maxIt=200, detail=FALSE, par=TRUE, modelOpts=list( 'iActs' = 15 ) )
        # abSan solution
  checkTrue( max( abs( example$ans$mZ - sol$vStar$mZ ) ) < tol )
}
