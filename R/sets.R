#######################################################################################
# Functions to describe the continuation sets in the Abreu-Sannikov method            #
# Philip Barrett, Chicago                                                             #
# Created: 30jul2013                                                                  #
# Notes: Based on modelSets.R in the credible monetary-fiscal policy project          #
#######################################################################################


##** 1. INITIALIZATION **##


##** 2. ASSIGNMENT **##

sets.setZ <- function( mZnew ){
# Returns a set with vertices mZnew, cleaned and ordered. mZ is a complete set
# description, so all other set data are updated to be consistent with mZnew
  set <- list( 'mZ' = mZnew, vC = NULL, mG = NULL )
  return( sets.clean( set ) )
}

sets.setC <- function( set, vCnew ){
# Replaces the distances vC with vCnew in set. Replaces mZ and mQ as these may 
# not be consistent, but assumes mG is unchanged.  Requires set as an argument
# as a complete set description requires vC *and* mG
  if ( length( vCnew ) != length( set$vC ) ) 
    stop( 'Replacement vC must be same size as existing one' )
        # Error handling
  mG <- set$mG
        # Record the old mG
  mZ <- sets.GCtoZ( mG, vCnew )
        # Construct the associated vertices
  return( sets.setZ( mZ ) )
        # Clean and return the new set
}

sets.setGC <- function( mGnew, vCnew ){
# Replaces the normal/distance set description pair (mG,vC) with (mGnew,vCnew) 
# Replaces mZ and mQ as these may not be consistent.
  if ( nrow( mGnew ) != length( vCnew ) ) stop( 'mGnew and vCnew must be same size' )
        # Error handling
  mZ <- sets.GCtoZ( mGnew, vCnew )
        # Compute the new values of mZ
  return( sets.setZ( mZ ) )
        # Clean and return the new set
}


##** 3. CONVERSIONS **##

##** 4. CHECKS **##

sets.inside.pt.set <- function( vPt, set ){
# Wrapper for sets.inside.pt
  return( sets.inside.pt( vPt ), set$mG, set$vC )
}

##** 5. OPERATIONS **##

sets.clean <- function( set ){
# Orders and eliminates redundancy from a particular slice in a set
  out <- sets.trim( sets.chull( set ) )
}

sets.chull <- function( set ){
# Removes and reorders elements of a slice to match the clockwise-counted convex hull
  mZ <- set$mZ[ chull( set$mZ ), ]
        # Reduce the vertices to the convex hull
  lGC <- sets.ZtoGC( mZ )
        # Compute the associated normals and distances
  out <- c( list( 'mZ'=mZ ), lGC )
        # Combining the aswers into a slice
}

sets.trim <- function( set ){
# Eliminates near-redundancy from a particular slice in a set.  Assumes that the points 
# are ordered clockwise already
  mZ <- sets.unique( set$mZ, 1e-13 )
        # Retain the vertices whose co-ordinates differ by less than 1e-16 in both dimensions
  set <- c( list( 'mZ'=mZ ), sets.ZtoGC( mZ )  )
        # Update the mG, vC description
  mZ <- sets.shape.preserve( set, 1e-14)
        # Eliminate non-extremal points
  lGC <- sets.ZtoGC( mZ )
        # Compute the associated normals and distances
  out <- c( list( 'mZ'=mZ ), lGC )
        # Combining the aswers into a slice
}

sets.eliminate <- function( mZ, tol, par=FALSE ){
# Eliminates rows of mZ which are less than tol different from the first row
  iNN <- nrow( mZ )
        # Number of vertices
  mFirstRow <- matrix( mZ[1, ], nrow=iNN, ncol=2, byrow=TRUE )
        # Marix of just the first row
  if (par==FALSE) return( rbind( mZ[1, ], mZ[sapply( 1:iNN, function(i) sets.l2norm( mZ[1, ], mZ[i, ] ) >tol ), ] ) )
#   if (par==TRUE) return( rbind( mZ[1, ], 
#                                 mZ[ unlist( mclapply( 1:iNN, function(i) sets.l2norm( mZ[1, ], mZ[i, ] ) > tol ) ), ] ) )
}

sets.scale <- function( set, scale ){
# Multiplies all vertices of a set by a scale factor (can be <1).
  mZnew <- set$mZ * scale
  return( sets.setZ( mZnew ) )
        # Clean and return the new set
}

sets.add <- function( A, B ){
# Function to return the convex hull of the set-sum, i.e. the points v = a + b
# for all a in A and b in B
  extreme.mZ <- sets.Zadd( A$mZ, B$mZ )
        # The extreme points of the convex hull
  return( sets.setZ( extreme.mZ ) )
        # Clean and return the new set
}

sets.Zadd <- function( A.mZ, B.mZ ){
# Function to return the extreme points of the convex hull of the set-sum  
  nA <- nrow( A.mZ )
  nB <- nrow( B.mZ )
        # Number of rows of each matrix  
  repA <- matrix( rep( A.mZ, each=nB ), nB * nA, 2 )
  repB <- matrix( rep( B.mZ, nA ), nA * nB, 2, byrow=F )
        # Replicate the points
  return( repA + repB )
        # The set of all possible extreme points
}
  
sets.union <- function( ... ){
# Returns the convex hull of the union of two sets
  lSets <- list( ... )
        # List of sets to be united
  return( sets.union.list( lSets ) )
}

sets.union.list <- function( lSets ){
# Returns the union of a list of sets
  iSets <- length(lSets)
        # Number of sets
  lMatrices <- lapply( 1:iSets, function(i) lSets[[i]]$mZ )
  mZnew <- do.call("rbind", lMatrices)
        # Compute the extreme points
  return( sets.setZ( mZnew ) )
        # Clean and return the new set
}

sets.P <- function( set ){
# Calculates the minimum payoffs for both players
  return( apply( set$mZ, 2, min ) )
}

sets.between <- function( pt, set ){
# Computes the adjacent vertices of a point on the boundary
  
  vProj <- set$mG %*% pt - set$vC
        # Projection of pt onto normal vectors, less distances
  if( min( abs( vProj ) ) > 1e-10 ) stop('Vertex not on boundary')
        # Check point on some boundary
  
  idxVert <- apply( set$mZ, 1, function(x) identical( pt, x ) )
        # Check if pt is a vertex
  if (any(idxVert)) return( list( 'mZ'=set$mZ[ idxVert ] , 'prob'=1,
                                  'vDist' = apply( set$mZ, 1, function(x) 1 * identical( x, set$mZ[ idxVert ] ) )) )
  
  idxLine <- ( abs( vProj ) <= 1e-10 )
        # The line index which the point lies on
  mZ <- set$mZ[ abs( set$mZ %*% set$mG[ idxLine, ]  - set$vC[ idxLine ] ) < 1e-20, ]
        # The points on the line
  if( nrow( mZ ) != 2 ) stop('Multiple points along boundary')
        # Error check
  p <- ( pt - mZ[ 2, ] ) / ( mZ[ 1, ] - mZ[ 2, ] )
        # Solves pt = p * v_1 + ( 1 - p ) * v_2, for v_1, v_2 rows of mZ
  boCheckP <- ( abs( min(p) - max(p) ) < 1e-14 )
  if( !boCheckP ) stop('Randomization probabilities not correctly computed')
        # Another error check, both els of p should be the same
  vDist <- p[ 1 ] * apply( set$mZ, 1, function(x) identical( x, mZ[ 1, ] ) ) + 
    ( 1 - p[ 1 ] ) * apply( set$mZ, 1, function(x) identical( x, mZ[ 2, ] ) )
  return( list( 'mZ'=mZ, 'prob'=p[ 1 ], 'vDist'=vDist ) )  
}

##** 6. HAUSDORFF DISTANCE **##

# 0. Utilities

# 1. Body
sets.Hausdorff <- function( A, B ){
  # Hausdorff distance between sets defines by two matrices of points (i.e. the 
  # furthest nearest neighbour). This is the maxmin distance *from* A *to* B. The 
  # algorithm relies on the fact that the Hausdorff line must extend from a 
  # *vertex* of A, but can land at either a vertex or a side of B. So we first 
  # work out the distance from each vertex of A to the lines defining the eadge of
  # B, but reject the result if it is not in B. Then, if necessary, we compute the
  # distance to the nearest vertex.
  
  B <- B[ chull(B), ]
  # Ordering the Vertices
  lGC <- sets.ZtoGC( B )
  # Computing the normal-distance representation
  return( sets.Hausdorff.core( A, B, lGC$mG, lGC$vC ) )
  
}