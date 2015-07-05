#######################################################################################
# Test functions for sets source file                                                 #
# Philip Barrett, Chicago                                                             #
# Created: 30jul2013                                                                  #
#######################################################################################


##** 1. INITIALIZATION **##


##** 2. ASSIGNMENT **##

test.setVars <- function(){
# Test functions to set data inside a set description object
  
  # 0. Set up
  mZ <- matrix( c(  1, 0, 0, 0, 0, 1, 1, 1 ), ncol=2, byrow=TRUE )
        # New mZ
  vC <- c( 0, 0, 1, 1 )
        # New vC
  vClong <- c( vC, vC )
        # Long vC to test error handling
  mG <- matrix( c( 0, -1, -1, 0, 0, 1, 1, 0 ), ncol=2, byrow=TRUE )
        # Quarters of the unit circle, tilted 45 degrees
  
  # 1. Test mZ replacement
  testSet <- sets.setZ( mZ )
  checkTrue( all( testSet$mZ == mZ ) )
  checkEquals( testSet$mG, mG )
  checkEquals( testSet$vC, vC )
          # Check simple assgnment
  testSet <- sets.setZ( rbind( mZ, c( 1, 1 ) ) )
  checkTrue( all( testSet$mZ == mZ ) )
          # Check elimination of duplicates  
  testSet <- sets.setZ( rbind( mZ, c( .5, .5 ) ) )
  checkTrue( all( testSet$mZ == mZ ) )
          # Check elimination of interior points
  
  # 2. Test vC replacement
  testSet <- sets.setZ( mZ )
  newTest <- sets.setC( testSet, vC * 2 )
        # Replace vC
  checkEquals( newTest$vC, vC * 2 )
  checkEquals( newTest$mG, testSet$mG )
  checkEquals( newTest$mZ, 2 * testSet$mZ )
  checkException( sets.setC( testSet,  vClong ) )
        # Should fail, vClong is not the same size as the original vC

  # 3. Test mG/vC replacement
  testSet <- sets.setZ( mZ )
  newTest <- sets.setGC( mG, vC )
        # Replace mG and vC simultaneously
  checkEquals( testSet, newTest )
        # Check replacement-assigment works with mG and vC together
  mGnew <- matrix( c( 1, 1, 1, -1, -1, -1, -1, 1 ), ncol=2, byrow=TRUE )
  vCnew <- rep( 1 , 4 )
  mZnew <- matrix( c( 1, 0 , 0, -1, -1, 0, 0, 1 ), ncol=2, byrow=TRUE )
        # Now try with a new set
  newTest <- sets.setGC( mGnew, vCnew )
        # Replace mG and vC simultaneously
  checkTrue( max( abs( newTest$mZ - mZnew ) ) < 1e-15 )
        # Check the replacement of mZ
  checkException( sets.setGC( mG, vClong ) )
        # Should fail, vClong is not the same length as mG
  
}


##** 3. CONVERSIONS **##

test.conversions <- function(){
# Test conversion functions

  # 0. Set up
  mZ <- matrix( c(  1, 0, 0, 0, 0, 1, 1, 1 ), ncol=2, byrow=TRUE )
        # New mZ
  vC <- c( 0, 0, 1, 1 )
        # New vC
  vClong <- c( vC, vC )
        # Long vC to test error handling
  mG <- matrix( c( 0, -1, -1, 0, 0, 1, 1, 0 ), ncol=2, byrow=TRUE )
        # Quarters of the unit circle, tilted 45 degrees
  testSet <- sets.setZ( mZ )
  
  # 1. sets.ZGtoC
  checkEquals( sets.ZGtoC( mZ, mG ), vC )
        # Should generate vC identical to the one there already (as set is consistent)
  
  # 4. sets.GCtoZ
  mZtest <- sets.GCtoZ( mG, vC )
        # Compute mZ from mG, vC
  checkTrue( all( mZ == mZtest[ chull( mZtest ), ] ) )
        # Check the test case exactly
  checkTrue( sets.inside.pts( mZ, mG, vC) )
        # Check all points lie inside the G,c description
  vCtest <- sets.ZGtoC( mZ, mG )
  checkEquals( vCtest, vC )
        # Check inputs against the inverse conversion
  
  # 5. sets.ZtoGC
  lGC <- sets.ZtoGC( mZ )
        # The function to be tested
  checkEquals( lGC$mG, mG )
        # Check the simple test case
  lGC <- sets.ZtoGC( mZ )
  checkTrue( sets.inside.pts( mZ, lGC$mG, lGC$vC) )
        # Check all points lie inside the G,c description
  mZtest <- sets.GCtoZ( mG, vC )
  checkTrue( all( mZ == mZtest[ chull( mZtest ), ] ) )
        # Check inputs against the inverse conversion
  
}


##** 4. CHECKS **##

test.checks <- function(){
 
  # 1. Set up
  vC <- rep( 1, 4 )
  mG <- matrix( c( 1, 0, 0, 1, -1, 0, 0, -1 ), ncol=2, byrow=TRUE )
        # Test case is for points inside the unit circle
  
  # 2. Test true/false for single points inside the unit circle
  checkTrue( sets.inside.pt( c( 0, 0 ), mG, vC ) )
  checkTrue( !sets.inside.pt( c( 10, 0 ), mG, vC ) )
  checkTrue( sets.inside.pt( c( 0.9999999, 0 ), mG, vC ) )
  checkTrue( !sets.inside.pt( c( 1.0000001, 0 ), mG, vC ) )
        # Positive direction
  checkTrue( sets.inside.pt( c( -0.9999999, 0 ), mG, vC ) )
  checkTrue( !sets.inside.pt( c( -1.0000001, 0 ), mG, vC ) )
        # Negative direction
  
  # 2. Test true/false for collections of points inside the unit circle
  checkTrue( sets.inside.pts( matrix( c( 0, 0 , 1, 0), ncol=2, byrow=TRUE ), mG, vC ) )
  checkTrue( sets.inside.pts( matrix( c( 0, .99999 , 1, 0), ncol=2, byrow=TRUE ), mG, vC ) )
  checkTrue( !sets.inside.pts( matrix( c( 0, 1.01 , 1, 0), ncol=2, byrow=TRUE ), mG, vC ) )
  checkTrue( !sets.inside.pts( matrix( c( 0, 0 , 1.01, 0), ncol=2, byrow=TRUE ), mG, vC ) )
  checkTrue( !sets.inside.pts( matrix( c( 0, 1.01 , 1.01, 0), ncol=2, byrow=TRUE ), mG, vC ) )
        # Positive direction
  checkTrue( sets.inside.pts( matrix( c( 0, 0 , -1, 0), ncol=2, byrow=TRUE ), mG, vC ) )
  checkTrue( !sets.inside.pts( matrix( c( 0, 0 , -1.01, 0), ncol=2, byrow=TRUE ), mG, vC ) )
        # Negative direction
  checkTrue( sets.inside.pts( matrix( c( 0, 0 , -1, 0, .5, 0, 0, -.5 ), ncol=2, byrow=TRUE ), mG, vC ) )
        # More points
}


##** 5. OPERATIONS **##
test.operations <- function(){
  
  # 1. sets.slice.chull
  mZ <- matrix( c( 1, 0, 0, 1, -1, 0, 0, -1, 0, 0 ), ncol=2, byrow=TRUE )
  lCGans <- sets.ZtoGC( mZ[ chull( mZ ), ] )
      # New data
  testSet <- sets.setZ( mZ )
      # Define the test set
  newSet <- sets.chull( testSet )
      # Generate the chull set
  checkTrue( all( testSet$mZ == mZ[ chull( mZ ), ] ) )
  checkEquals( testSet$mG, lCGans$mG )
  checkEquals( testSet$vC, lCGans$vC )
  
  # 2. sets.eliminate
  mZ <- matrix( c( 1, 0, 0, 1, 0, 0, 0, 1e-24, 1, 1e-24 ), ncol=2, byrow=TRUE )
  checkEquals( sets.eliminate( mZ, 1e-24 ), mZ[ -5, ] )
  
  # 3. sets.unique
  mZ <- matrix( c( 1, 0, 0, 0, 0, 1e-24, 0, 1 ), ncol=2, byrow=TRUE )
  checkTrue( all( sets.unique( mZ, 1e-24 ) == mZ[ -3, ] ) )
        # Very simple test case
  mZ <- matrix( c( 1, 0, 0, 1, 0, 0, 0, 1e-24, 1, 1e-24 ), ncol=2, byrow=TRUE )
  checkTrue( all( sets.unique( mZ, 1e-24 ) == mZ[ c(2,3,5), ] ) )
        # Merely simple test case - Funny re-ordering of the rows as not unique
  mZ <- matrix( c( 1, 0, 1, 0, 1, 0, 1, 0 ), ncol=2, byrow=TRUE )
  checkTrue( all( sets.unique( mZ, 1e-24 ) == mZ[ 1, ] ) )
        # A case wehre we eliminate everything
  
  # 4. sets.trim
  mZ <- matrix( c( 1, 0, 0, 0, 0, 0, 0, 1,  1e-24, 1, .1, .1 ), ncol=2, byrow=TRUE )
  mZans <- matrix( c( 1, 0, 0, 0, 0, 1, .1, .1 ), ncol=2, byrow=TRUE )
        # New data for the slice
  testSet$mZ <- mZ
        # Replace the appropriate data
  newSet <- sets.trim( testSet )
         # Generate the trimmed slice
  checkTrue( all( newSet$mZ == mZans ) )
        # Compare to the test case - should eliminate unnecessary points, but retain interior one.
        # Do not test mG, vC here, as sets.slice.trim requires mZ ordered before application. This is
        # Tested in sets.slice.slean, which uses sets.slice.trim
  
  # 5. sets.clean
  mZ <- matrix( c( 1, 0, 0, 1, 0, 0, 1e-24, 1, 0, 0, .1, .1 ), ncol=2, byrow=TRUE )
  mZans <- matrix( c( 1, 0, 1e-24, 1, 0, 0 ), ncol=2, byrow=TRUE )
        # New data for the slice
  mZans <- mZans[ chull( mZans ), ]
        # Order the data
  testSet <- sets.setZ( mZ )
        # Replace the appropriate slice
  lCGans <- sets.ZtoGC( mZans )
        # New data for the conversion
  newSet <- sets.clean( testSet )
        # Generate the cleaned slice
  checkTrue( all( newSet$mZ == mZans ) )
  checkEquals( newSet$mG, lCGans$mG )
  checkEquals( newSet$vC, lCGans$vC )
        # Compare to the test case - should eliminate unnecessary points and order points clockwise
  checkTrue( sets.inside.pt( c( .1, .1) , newSet$mG, newSet$vC ) ) 
        # Make sure that mG,vC are correctly oriented
  
  # 6. sets.scale
  mZ <- matrix( c(  1, 0, 0, 0, 0, 1, 1, 1 ), ncol=2, byrow=TRUE )
  testSet <- sets.setZ( mZ )
  doubleSet <- sets.scale( testSet, 2 )
        # Create double-sized set
  checkEquals( doubleSet$mZ, testSet$mZ * 2 )
  checkEquals( doubleSet$mG, testSet$mG )
  checkEquals( doubleSet$vC, testSet$vC * 2 )
        # Compare to the test case  
  
  # 7. sets.add
  A.mZ <- matrix( c(  1, 0, 0, 0, 0, 1 ), ncol=2, byrow=TRUE )
  ASet <- sets.setZ( A.mZ )
  B.mZ <- matrix( c(  -1, 0, 0, 0, 0, -1 ), ncol=2, byrow=TRUE )
  BSet <- sets.setZ( B.mZ )
  ans.mZ <- matrix( c( 0, -1, -1, 0, -1, 1, 0, 1, 1, 0, 1, -1 ), ncol=2, byrow=TRUE )
  ansSet <- sets.setZ( ans.mZ )
        # Set up the two sets to add and create the answer
  testSet <- sets.add( ASet, BSet )
        # Create the set sum
  checkEquals( testSet$mZ, ansSet$mZ )
  checkEquals( testSet$mG, ansSet$mG )
  checkEquals( testSet$vC, ansSet$vC )
        # Check that the test is equal to the anticipated answer
  
  # 8. sets.union
  ans.mZ <- matrix( c( 0, -1, -1, 0, 0, 1, 1, 0 ), ncol=2, byrow=TRUE )
  ansSet <- sets.setZ( ans.mZ )
        # Set up the two sets to add and create the answer
  testSet <- sets.union( ASet, BSet )
        # Create the set sum
  checkEquals( testSet$mZ, ansSet$mZ )
  checkEquals( testSet$mG, ansSet$mG )
  checkEquals( testSet$vC, ansSet$vC )
        # Check that the test is equal to the anticipated answer
  C.mZ <- matrix( c(  0, 2, -1, 0, 1, 0 ), ncol=2, byrow=TRUE )
  CSet <- sets.setZ( C.mZ )
        # Create a third set to take the union of
  ans2.mZ <- matrix(  c( 0, -1, -1, 0, 0, 2, 1, 0 ), ncol=2, byrow=TRUE )
  ansSet2 <- sets.setZ( ans2.mZ )
        # Create the answer
  testSet2 <- sets.union( ASet, BSet, CSet )
        # Create the set sum
  checkEquals( testSet2$mZ, ansSet2$mZ )
  checkEquals( testSet2$mG, ansSet2$mG )
  checkEquals( testSet2$vC, ansSet2$vC )
        # Check that the test is equal to the anticipated answer
  
  # 9. sets.intersect
  intersect <- sets.lineIntersect( ASet, c( 0, 1 ), .5 )
  checkEquals( intersect, matrix( c( 0, .5, .5, .5 ), 2, 2, byrow=T ) )
  intersect <- sets.lineIntersect( ASet, c( 0, 1 ), -.5 )
  checkEquals( intersect, matrix( 0, nrow=0, ncol=2 ) )
  intersect <- sets.lineIntersect( ASet, c( 0, 1 ), 1 )
  checkEquals( intersect, matrix( c( 0, 1 ), 1, 2, byrow=T ) )
  intersect <- sets.lineIntersect( ASet, c( -1, 1 ), 0 )
  checkEquals( intersect, matrix( c( 0, 0, .5, .5 ), 2, 2, byrow=T ) )
  
  # 10. sets.l2norm
  checkEquals( sets.l2norm( c(0,0), c(1,0) ), 1 )
  checkEquals( sets.l2norm( c(0,0), c(1,1) ), sqrt( 2 ) )
  checkEquals( sets.l2norm( c(0,0), c(3,4) ), 5 )
  checkEquals( sets.l2norm( c(1,2), c(4,6) ), 5 )
  checkException( sets.l2norm( c(0,0,1), c(1,1) ) )
  
  # 11. sets.nearest.vert
  C.mZ <- matrix( c(  1, 0, 0, 0, 0, 1 ), ncol=2, byrow=TRUE )
  CSet <- sets.setZ( C.mZ )
  checkEquals( sets.nearest.vert( c( 1, 0 ), CSet$mZ ), 0 )
  checkEquals( sets.nearest.vert( c( 1, 1 ), CSet$mZ ), 1 )
  checkEquals( sets.nearest.vert( c( .5, .5 ), CSet$mZ ), 1 / sqrt( 2 ) )
  
  # 12. sets.nearest.line
  checkEquals( sets.nearest.line( c( 1, 1 ), CSet$mG, CSet$vC ), 1 / sqrt( 2 ) )
  checkEquals( sets.nearest.line( c( .5, -1 ), CSet$mG, CSet$vC ), 1 )
  checkEquals( sets.nearest.line( c( -1, 2 ), CSet$mG, CSet$vC ), Inf )
  
  # 12. sets.Hausdorff
  D.mZ <- matrix( c(  -1, 0, 0, 0, 0, -1 ), ncol=2, byrow=TRUE )
  DSet <- sets.setZ( D.mZ )
  checkEquals( sets.Hausdorff( CSet$mZ, DSet$mZ ), 1 )
        # Basic test case
  E.mZ <- matrix( c(  -1, 0.5, 0, 1, 0, 0 ), ncol=2, byrow=TRUE )
  ESet <- sets.setZ( E.mZ )
  checkEquals( sets.Hausdorff( ESet$mZ, CSet$mZ ), 1 )
        # Slightly harder test - the line distance matters
  checkEquals( sets.Hausdorff( CSet$mZ, ESet$mZ ), 1 )
        # The same
  J.mZ <- rbind( C.mZ, c( 0, 0.5 ) )
  JSet <- sets.setZ( J.mZ )
  checkEquals( sets.Hausdorff( CSet$mZ, JSet$mZ ), 0 )
        # Slightly harder test - the line distance matters
  
  # 13. sets.between
  A.mZ <- matrix( c(  1, 0, 0, 0, 0, 1 ), ncol=2, byrow=TRUE )
  ASet <- sets.setZ( A.mZ )
  XX <- sets.between( c(.5, .5), ASet )
  checkTrue( all( XX$mZ == matrix( c(1,0,0,1), nrow=2 ) ) )
  checkEquals( XX$prob, .5 )
  YY <- sets.between( c(.2, .8), ASet )
  checkTrue( all( YY$mZ == matrix( c(1,0,0,1), nrow=2 ) ) )
  checkEquals( YY$prob, .2 )
  ZZ <- sets.between( c(1, 0), ASet )
  checkTrue( all( ZZ$mZ == c(1,0) ) )
  checkEquals( ZZ$prob, 1 )
  checkException( sets.between( c(1, 1), ASet ) )
  
}


