/**************************************************************************************
* Alalgam of code from: sets.cpp, setConv.cpp, and abreuSannikov.cpp                  *
* Used for package definitions                                                        *
* Philip Barrett, Chicago                                                             *
* Created: 27aug2013                                                                  *
***************************************************************************************/

// [[Rcpp::depends(RcppEigen)]]
#include<Rcpp.h>
#include<RcppEigen.h>
using namespace Rcpp ;
using namespace Eigen ;

// Useful typdef conversions
typedef Map<MatrixXd> eMatd ;
typedef Map<VectorXd> eVecd ;
typedef Map<VectorXi> eVeci ;

double sets_l2norm_cpp_eigen( const VectorXd v1, const VectorXd v2 ){
// C++ translation of sets.l2norm

  // Error checking
  int iLen1 = v1.size() ;
  int iLen2 = v2.size() ;
        // The length of the two vectors
  if( iLen1 != iLen2 )
    ::Rf_error("Input vectors must be same length");
  return ( v1 - v2 ).norm() ;
}

//[[Rcpp::export('sets.l2norm')]]
double sets_l2norm_cpp( NumericVector v1, NumericVector v2 ){
// Wrapper for sets_l2norm_cpp_eigen
    const eVecd ev1( as<eVecd>( v1 ) ) ;
    const eVecd ev2( as<eVecd>( v2 ) ) ;
    return sets_l2norm_cpp_eigen( ev1, ev2 ) ;
}

double sets_nearest_vert_cpp_eigen( VectorXd pt, MatrixXd mB ){
// C++ translation of sets.nearest.vert
  
  int iB = mB.rows() ;
        // Number of rows in B
  VectorXd vDists( iB ) ;
        // Initialize the distance vector
  for( int iII=0; iII<iB; iII++ ){
    vDists[ iII ] = ( pt.transpose() - mB.row( iII ) ).norm() ;
          // Compute the distance to pt
  }
  return vDists.minCoeff() ;
}

// [[Rcpp::export('sets.nearest.vert')]]
double sets_nearest_vert_cpp( NumericVector pt, NumericMatrix B ){
// Wrapper for sets_nearest_vert_cpp_eigen
    const eVecd evPt( as<eVecd>( pt ) ) ;
    const eMatd emB( as<eMatd>( B ) ) ;
    return sets_nearest_vert_cpp_eigen( evPt, emB ) ;
}

bool sets_inside_pt_cpp_eigen( VectorXd vPt, MatrixXd mG, VectorXd vC ){
// C++ translation of sets.inside.pt
// Tests if the point vPt is inside the set described by normal/distance pair mG, vC
  
  //Error checking
  if( vPt.size() != 2 ) 
    ::Rf_error("vPt must have length 2") ;
  
  // Do the necessary matrix mutliplication
  VectorXd vProd = mG * vPt ;
  
  // Now check if all the points are within the bounds
  bool boInside = FALSE ;
  for( int iII=0; iII<vProd.rows(); iII++ ){
    boInside = ( vProd[ iII ] <= vC[ iII ] + 1e-11 ) ;
    if ( !boInside ) return( boInside ) ;
  }
  return( boInside ) ;
}

// [[Rcpp::export('sets.inside.pt')]]
bool sets_inside_pt_cpp( NumericVector vPt, NumericMatrix mG, NumericVector vC ){
// Wrapper for sets_inside_pt_cpp_eigen to allow R-style matrix I/O
  
  // Convert to Eigen types
  const eMatd emG( as<eMatd>( mG ) ) ;
  const eVecd evPt( as<eVecd>( vPt ) ) ;
  const eVecd evC( as<eVecd>( vC ) ) ;
  return sets_inside_pt_cpp_eigen( evPt, emG, evC ) ;
}

// [[Rcpp::export('sets.inside.pts')]]
bool sets_inside_pts_cpp( NumericMatrix mPts, NumericMatrix mG, NumericVector vC ){
// C++ translation of sets.inside.pt
// Tests if the points mPts are inside the set described by normal/distance pair mG, vC
  
  bool boInside ;
        // Instantiate boolean output
  // Convert to Eigen types
  const eMatd emG( as<eMatd>( mG ) ) ;
  const eMatd emPts( as<eMatd>( mPts ) ) ;      
  const eVecd evC( as<eVecd>( vC ) ) ;      
        
  for( int iII; iII<mPts.nrow(); iII++ ){
    boInside = sets_inside_pt_cpp_eigen( emPts.row( iII ), emG, evC ) ;
    if ( !boInside ) return( boInside ) ;
  }
  return( boInside ) ;
}

LogicalVector sets_inside_pts_idx_cpp_eigen( MatrixXd mPts, MatrixXd mG, VectorXd vC ){
// Returns the sets of indices of points inside the set described by mG, vC
  int iRows = mPts.rows() ;
        // Number of points
  LogicalVector boInside( iRows ) ;
        // Instantiate boolean output
  for( int iII; iII<iRows; iII++ ){
    boInside[ iII ] = sets_inside_pt_cpp_eigen( mPts.row( iII ), mG, vC ) ;
  }
  return( boInside ) ;
}

// [[Rcpp::export('sets.inside.pts.idx')]]
LogicalVector sets_inside_pts_idx_cpp( NumericMatrix mPts, NumericMatrix mG, NumericVector vC ){
// Wrapper for sets_inside_pts_idx_cpp_eigen to allow R-style matrix I/O
  // Convert to Eigen types
  const eMatd emG( as<eMatd>( mG ) ) ;
  const eMatd emPts( as<eMatd>( mPts ) ) ;
  const eVecd evC( as<eVecd>( vC ) ) ;
  return sets_inside_pts_idx_cpp_eigen( emPts, emG, evC ) ;
}

double sets_nearest_line_cpp_eigen( VectorXd pt, MatrixXd mG, VectorXd vC ){
// Eigen-format c++ rewrite of sets.nearest.line
// Computes distance of a point to the lines satisfying g.z = c, then picks
// the nearest which is inside B.  Relies on mG all being unit length
  int iNN = vC.size() ;
        // Number of lines
  const VectorXd vDiff = vC - mG * pt ;
        // The difference of projections of p onto the mG vectors from vC
  const VectorXd vDiffAbs( vDiff.cwiseAbs() ) ;
        // Absolute value of the distance
  const MatrixXd mNearPts = VectorXd::Ones( iNN ) * pt.transpose() + 
                              mG.cwiseProduct( vDiff * VectorXd::Ones( 2 ).transpose() ) ;
        // Compute the matrix of nearest points to p       
  LogicalVector boIdx = sets_inside_pts_idx_cpp_eigen( mNearPts, mG, vC ) ;
        // Vector of logical indices for mNearPts inside mG, vC set
  if ( is_false( any( boIdx ) ) ) return( INFINITY ) ;
        // Return outsized distance if none of mNearPts is inside the set
  VectorXd vDistInside( iNN ) ;
        // Initialize the vector distances of mNearPts entries inside the set  
  for( int iII=0; iII<iNN; iII++ )
    vDistInside[ iII ] = ( boIdx[ iII ] == 1 ? vDiffAbs[ iII ] : INFINITY ) ;
        // Compute the vector distances of mNearPts entries inside the set
  return vDistInside.minCoeff() ;
}

//[[Rcpp::export('sets.nearest.line')]]
double sets_nearest_line_cpp( NumericVector pt, NumericMatrix mG, NumericVector vC ){
// Wrapper for sets_nearest_line_cpp_eigen - allows conversion to R-style matrices
  const eMatd emG( as<eMatd>( mG ) ) ;
  const eVecd evPt( as<eVecd>( pt ) ) ;
  const eVecd evC( as<eVecd>( vC ) ) ;
  return sets_nearest_line_cpp_eigen( evPt, emG, evC) ;
}


//[[Rcpp::export('sets.Hausdorff.core')]]
double sets_Hausdorrf_core_cpp( NumericMatrix mA, NumericMatrix mZ, NumericMatrix mG, NumericVector vC ){
// Computes the Hausdorf distance between points in mA and the set defined by mZ-mG-vC

  // Eigen Conversion
  const eMatd emA( as<eMatd>( mA ) ) ; 
  const eMatd emZ( as<eMatd>( mZ ) ) ; 
  const eMatd emG( as<eMatd>( mG ) ) ; 
  const eVecd evC( as<eVecd>( vC ) ) ;

  int iRows = emA.rows() ;
        // The number of points in mA
  double dLine = 0 ;
  double dVert = 0 ;
  double dNewDist = 0;
  double dHausdorff = 0 ;
        // Initializing the distance variables
  
  // Loop over points in mA
  for( int iII=0; iII<iRows; iII++ ){
    dLine = sets_nearest_line_cpp_eigen( emA.row( iII ), emG, evC ) ;
          // Distance to nearest line
    dVert = sets_nearest_vert_cpp_eigen( emA.row( iII ), emZ ) ;
          // Distance to nearest vertex
    dNewDist = ( dLine < dVert ? dLine : dVert ) ;
          // The minimum distance
    dHausdorff = ( dNewDist > dHausdorff ? dNewDist : dHausdorff ) ;
          // Update Hausdorff with a minmax
  }

  return dHausdorff ;
}

MatrixXd sets_filter_rows( const MatrixXd mX, const VectorXi boV ){
// Returns rows of mX denoted by boV

  // Error checking
  if( mX.rows() != boV.size() )
    ::Rf_error("Dimension mismatch in row filter") ;
    
  const int iV = boV.sum() ;
        // Number of rows to be selected
  MatrixXd mOut( iV, mX.cols() ) ;
        // Initialize the output matrix
  int iCounter = 0 ;
        // Counter for rows of boV
  for( int iII=0; iII<mX.rows(); iII++ ){
    if ( boV[ iII ] == 1 )
      mOut.row( iCounter ) = mX.row( iII ) ;
    iCounter = iCounter + boV[ iII ] ;
  }
  return mOut ;
}

MatrixXd sets_filter_rows( const MatrixXd mX, const IntegerVector borV ){
  const Map<MatrixXi> boV( as<Map<MatrixXi> >(borV)) ;
  return sets_filter_rows( mX, boV ) ;
}

MatrixXd sets_filter_rows( const MatrixXd mX, const LogicalVector borV ){
  const IntegerVector boV( borV ) ;
  return sets_filter_rows( mX, boV ) ;
}

NumericMatrix sets_filter_rows( const NumericMatrix mX, const LogicalVector borV ){
  const eMatd emX( as<eMatd>( mX ) ) ;
  const IntegerVector boV( borV ) ;
  return wrap( sets_filter_rows( emX, boV ) );
}


MatrixXd sets_unique_cpp_eigen( const MatrixXd mZ, double tol ){
// Eliminates near-dulplicate rows from a clockwise-ordered matrix of 2D vertices
// Eigen-format I/O format rewrite of sets.unique
  
  const int iNN = mZ.rows() ;
        // Number of points
  MatrixXd mZneighbors( mZ ) ;
  mZneighbors.row( 0 ) = mZ.row( iNN - 1 ) ;
  mZneighbors.bottomRows( iNN - 1 ) = mZ.topRows( iNN - 1 ) ;
        // Compute the matrix of offset rows
  const MatrixXd mDiff = mZ - mZneighbors ;
  const VectorXd vDiff = mDiff.cwiseProduct( mDiff ).rowwise().sum().cwiseSqrt() ;
        // The l2norm difference between each point and the next

  VectorXi boValid = VectorXi::Zero( iNN ) ;
  
  for( int iII=0; iII<iNN; iII++ )
    if ( vDiff[ iII ] > tol )
      boValid[ iII ] = 1 ;
            // Count the points which are sufficiently different
  
  if( boValid.sum() == 0 )
    return( mZ.row( 0 ) ) ;
            // Special case when all entries are the same

return( sets_filter_rows( mZ, boValid ) ) ;
}

//[[Rcpp::export('sets.unique')]]
NumericMatrix sets_unique_cpp( const NumericMatrix mZ, double tol ){
// Wrapper for sets_unique_cpp_eigen - allows conversion to R-style matrices
  const eMatd emZ( as<eMatd>( mZ ) );
        // Convert to Eigen Matrix
  return( wrap( sets_unique_cpp_eigen( emZ, tol ) ) ) ;
}

//[[Rcpp::export('sets.lineIntersect')]]
SEXP sets_lineIntersect_cpp( const List set, const NumericVector dir, const double intersect ){
// c++ rewrite of cets.lineInterset
// Returns the intersection of a line dir.x = int with the boundary of set.  May
// return two points, one point, or NULL (if there is no intersection)

  // 0. Set up
  const eMatd mG( as<eMatd>( set["mG"] ) ) ; 
  const eVecd vC( as<eVecd>( set["vC"] ) ) ;
  const eVecd vDir( as<eVecd>( dir ) ) ;
        // Eigen conversions
  const int iNN = vC.size() ;
        // Number of normal vectors
  MatrixXd mIntersect = MatrixXd::Zero( iNN, 2 ) ;
        // Matrix of intersecting points
  Matrix2d mA ;
  mA << 0, 0, vDir.transpose() ; 
        // Initialize the matrix of gradients used throughout
  Vector2d vB( 0, intersect ) ;
        // Initialize the vector in Ax=b used throughout
  double dRelErr ;
        // The relative error of matrix inversion
  int iCounter = 0 ;
        // Counter to check how many rows have meaningful solutions
  Vector2d vNewSol ;
        // Placeholder for the new solution
  
  
  // 1. Solve for the intersections of the given line with all (extended) boundaries
  for ( int iII=0; iII<iNN; iII++ ){
    mA.row( 0 ) = mG.row( iII ) ;
    vB[ 0 ] = vC[ iII ] ;
          // Creating matrix and vector of the form A*x=b to solve for x
//    vNewSol = mA.householderQr().solve(vB) ;
          // Faster, less accurate
    vNewSol = mA.fullPivHouseholderQr().solve(vB) ;
          // Slower, very accurate (actually ends up with a faster overall algorithm)
    if( vB.norm()==0 )
      dRelErr = ( mA * vNewSol - vB ).norm() ;
    else
      dRelErr = ( mA * vNewSol - vB ).norm() / vB.norm() ;
          // Relative error of the solution - special case for vB=(0,0)
    if( dRelErr <= 1e-14){
      mIntersect.row( iCounter ) = vNewSol.transpose() ;
      iCounter++ ;
    } 
  }
  
  const MatrixXd mIntersectValid = mIntersect.topRows( iCounter ) ;
  
  // 2. Eliminate the points outside the convex hull
  const int iInt = mIntersectValid.rows() ;
        // Numer of total intersections
  LogicalVector boBoundary = sets_inside_pts_idx_cpp_eigen( mIntersectValid, mG, vC ) ;
        // The logical vector of boundary values
        
    // DEBUG LINES
//    Rcout << "mIntersectValid:\n" << mIntersectValid << std::endl << std::endl;
//    Rcout << "is_true( any(boBoundary) ) =" << is_true( any( boBoundary ) ) << std::endl << std::endl;
//    Rcout << "sets_filter_rows( mIntersectValid, boBoundary ), 1e-14 ):\n" << 
//          sets_filter_rows( mIntersectValid, boBoundary ) << std::endl << std::endl;
        
  if ( is_false( any( boBoundary ) ) ) {
    NumericMatrix mOut( 0, 2 ) ; 
    return mOut ; // R_NilValue ;
            // Return 0-row matrix if no intercept
  }
        
//  Rcout << "Past NULL return" << std::endl ;
  
  IntegerVector iBoundary( boBoundary ) ;
        // Convert logical to integer for filtering
  return( wrap( sets_unique_cpp_eigen( sets_filter_rows( mIntersectValid, iBoundary ), 1e-14 ) ) ) ;
        // Extract the boundary rows, remove non-unique points, and then wrap for output
}

MatrixXd sets_GCtoZ_cpp_eigen( const MatrixXd mG, const VectorXd vC ){
// Given a normal/distance description of a set, returns a set of points mZ which lie at
// the vertices.  This assumes that the vector of normals is ordered clockwise already.  

  int iNN = vC.size() ;
        // Number of normal vectors
  MatrixXd mZ( iNN, 2 ) ;
        // Initialize the matrix of intersections
  Matrix2d mA ;
  Vector2d vB ;
  double dRelErr = 0 ;
  VectorXd vNewSol ;
  double iCounter = 0 ;
  
  // 1. Solve for the interior points
  for ( int iII=0; iII < ( iNN - 1 ) ; iII++ ){
    mA = mG.block( iII, 0, 2, 2 ) ;
    vB = vC.segment( iII, 2 ) ;
          // Creating matrix and vector of the form A*x=b to solve for x
    vNewSol = ( mA.fullPivHouseholderQr().solve(vB) ) ;
          // Solving and adding to mZ
    if( vB.norm()==0 )
      dRelErr = ( mA * vNewSol - vB ).norm() ;
    else
      dRelErr = ( mA * vNewSol - vB ).norm() / vB.norm() ;
          // Relative error of the solution - special case for vB=(0,0)
    if( dRelErr <= 1e-14){
      mZ.row( iCounter ) = vNewSol.transpose() ;
      iCounter++ ;
    }
  }
  
  // 2. Solve for the end points
  mA << mG.row( iNN - 1 ), mG.row( 0 );
  vB << vC[ iNN - 1 ], vC[ 0 ] ;
        // Creating matrix and vector of the form A*x=b to solve for x
  vNewSol = ( mA.fullPivHouseholderQr().solve(vB) ) ;
  if( vB.norm()==0 )
    dRelErr = ( mA * vNewSol - vB ).norm() ;
  else
    dRelErr = ( mA * vNewSol - vB ).norm() / vB.norm() ;
        // Relative error of the solution - special case for vB=(0,0)
  if( dRelErr <= 1e-14){
    mZ.row( iCounter ) = vNewSol.transpose() ;
    iCounter++ ;
  }
  return( mZ.topRows( iCounter ) ) ;
}

//[[Rcpp::export('sets.GCtoZ')]]
NumericMatrix sets_GCtoZ_cpp( NumericMatrix mG, NumericVector vC ){
// Wrapper for sets_GCtoZ_cpp_eigen
  const eMatd emG( as<eMatd>( mG ) ) ;
  const eVecd evC( as<eVecd>( vC ) ) ;
  return wrap( sets_GCtoZ_cpp_eigen( emG, evC ) ) ;
}

VectorXd sets_ZGtoC_cpp_eigen( MatrixXd mZ, MatrixXd mG ){
// Computes the distances vC for the points mZ and search directions mG
  const MatrixXd mProj = mZ * mG.transpose() ;
  return mProj.colwise().maxCoeff() ;
        // Maximum of projections onto each of the search directions
}

//[[Rcpp::export('sets.ZGtoC')]]
NumericVector sets_ZGtoC_cpp( NumericMatrix mZ, NumericMatrix mG ){
// Wrapper for sets_ZGtoC_cpp_eigen
  const eMatd emZ( as<eMatd>( mZ ) ) ;
  const eMatd emG( as<eMatd>( mG ) ) ;
  return wrap( sets_ZGtoC_cpp_eigen( emZ, emG ) ) ;
}

MatrixXd sets_ZtoGC_cpp_eigen( MatrixXd mZ ){
// Given a set of vertices, computes normals and distances mG, vC.  Assumes that
// the vertices are ordered clockwise. Here, mG and vC are concatenated into one matrix

  const int iNN = mZ.rows() ;
        // Number of vertices
  Matrix2d mFlip ;
  mFlip << 0, 1, -1, 0 ;
        // Matrix to convert gradient to normal
  MatrixXd mZneighbors( iNN, 2 ) ;
  mZneighbors.topRows( iNN - 1 ) = mZ.bottomRows( iNN - 1 ) ;
  mZneighbors.row( iNN - 1 ) = mZ.row( 0 ) ;
        // Compute the matrix of offset rows
  const MatrixXd mDiff = mZneighbors - mZ ;
        // Difference between one element and the next
  MatrixXd mG_Raw  = mDiff * mFlip ;
        // Compute directions of mG
  MatrixXd mG( iNN, 2 ) ;
        // Initialize mG
  for( int iII; iII<iNN; iII++ )
    mG.row( iII ) = mG_Raw.row( iII ) / mG_Raw.row( iII ).norm() ;
  const VectorXd vC = sets_ZGtoC_cpp_eigen( mZ, mG ) ;
        // Compute vC
  MatrixXd mOut( iNN, 3 ) ;
  mOut << mG, vC ;
        // Concatenate garadient and intercept for output
  return mOut ;
}

//[[Rcpp::export('sets.ZtoGC')]]
List sets_ZtoGC_cpp( NumericVector mZ ){
// Wrapper for sets_ZtoGC_cpp_eigen. Convert output into an R-style list with
// element mG, vC.
  const eMatd emZ( as<eMatd>( mZ ) )  ;
  const MatrixXd mOut = sets_ZtoGC_cpp_eigen( emZ ) ;
        // The Eigen form of the function
  List lOut ;
        // The output list
  NumericMatrix mG = wrap( mOut.leftCols( 2 ) ) ;
        // Format the mG output
  lOut["mG"] = mG ;
        // Write it to the list
  const VectorXd vC = mOut.rightCols( 1 ) ;
        // Format the vC entry
  lOut["vC"] = wrap( vC ) ;
        // Write it to the list
  return lOut ;
}

MatrixXd sets_shape_preserve_cpp_eigen( MatrixXd mZ, MatrixXd mG, double dTol ){
// Remove vertices which make no difference to the shape, ie. non-extremal points
  
  const int iNN = mZ.rows() ;
        // Number of rows
  MatrixXd mGneighbors( iNN, 2 ) ;
  mGneighbors.row( 0 ) = mG.row( iNN - 1 ) ;
  mGneighbors.bottomRows( iNN - 1 ) = mG.topRows( iNN - 1 ) ;
        // The offset-by-one matrix for mG
  VectorXd vDiff = mG.cwiseProduct( mGneighbors ).rowwise().sum() - VectorXd::Ones( iNN ) ;
        // Projection of mG onto its neighbors
  VectorXi iLarge = (vDiff.array().abs() > dTol ).cast<int>(); ;
        // Integer matrix of comparisons
  return sets_filter_rows( mZ, iLarge ) ;
}

//[[Rcpp::export(sets.shape.preserve)]]
NumericMatrix sets_shape_preserve_cpp( List set, double dTol ){
// Wrpper for sets_shape_preserve_cpp_eigen
  const eMatd emZ( as<eMatd>( set["mZ"] ) ) ;
  const eMatd emG( as<eMatd>( set["mG"] ) ) ;
  return wrap( sets_shape_preserve_cpp_eigen( emZ, emG, dTol ) ) ;
}

//[[Rcpp::export('abSan.C')]]
NumericMatrix abSan_C_cpp( int iA, const List set, const NumericVector vPun, const List model, 
      const bool print_output = true, bool boCountFrom1 = true ){
/* Computes Abreu & Sannikov's C(a,W,u) object. Inputs are: iA, an action index;
 * set is a set of continuation values; vPun is a pair of punishment values 
 * for the two players; and model is a model description created by model.initiate. */
  
/* TODO: Change the I-O to eliminate passing the whole model. Only part it used. 
 * Rewrite with Eigen containers (maybe?)
 */
  
  /* 1. Set up */
  if( boCountFrom1 ) iA-- ;
        // Transform from 1-based counting to 0-based
  NumericMatrix mF = model["mF"] ;
  NumericMatrix mH = model["mH"] ;
  const NumericVector gA = mF.row( iA ) ;
        // The payoff of the current action
  const NumericVector hA = mH.row( iA ) ;
        // The best deviating payoffs for the players given that the the other
        // player plays the action defined in the action index.
  const double delta = as<double>( model["delta"] ) ;
        // Saves extra writing later
        
  /* 2. Return g(a) if incentive compatible, else compute the extreme points  
   * assuming that the IC constraint binds                                    */
  const int iIC = as<int>( all( delta * ( gA - vPun ) >= ( 1 - delta ) * hA ) ) ;
        // Boolean. True if incentive compatibility holds for both agents at a
        
  if ( iIC == 1 ){
    NumericMatrix mgA( 1, 2 ) ;
    mgA( 0, _ ) = gA ;
    return( mgA ) ;
          // The action is IC. Return it.
  }else{ 
    /* 3. If g(a) not IC, calculate the extreme points of set where IC binds  */
    const NumericVector vW = vPun + ( 1 - delta ) / delta * hA ;
          // The pair of values in W where the IC constraint binds
    const NumericMatrix mExtremePts0 = sets_lineIntersect_cpp( set, NumericVector::create(1.0, 0.0), vW[ 0 ] ) ;
    const NumericMatrix mExtremePts1 = sets_lineIntersect_cpp( set, NumericVector::create(0.0, 1.0), vW[ 1 ] ) ;
          // Compute the intersections of the lines w_1 = vW[1] and w_2 = vW[2] with set
    const LogicalVector mboExtremeIC0 = mExtremePts0( _, 1 ) >= vW[ 1 ] ;
    const LogicalVector mboExtremeIC1 = mExtremePts1( _, 0 ) >= vW[ 0 ] ;
          // Index which rows which meet both IC conditions
    NumericMatrix mExtremeIC0 = sets_filter_rows( mExtremePts0, mboExtremeIC0 ) ;
    NumericMatrix mExtremeIC1 = sets_filter_rows( mExtremePts1, mboExtremeIC1 ) ;
          // Strip out the non-incentive-compatible points
    const bool boICinSet = sets_inside_pt_cpp( vW, set["mG"], set["vC"] ) ;
          // See if the punishment point is inside the set
    const int iRows0 = mExtremeIC0.nrow() ;
    int iRows = iRows0 + mExtremeIC1.nrow() ;
          // Number of rows in the output matrix
    if( boICinSet ) iRows++ ;
          // Add an extra row if the punishment is in the set
    NumericMatrix mExtremeIC( iRows , 2 ) ;
          // Initialize the output matrix
    if( boICinSet ) iRows-- ;
          // Knock 1 off the counter if boICinSet holds true
    for ( int iII=0; iII<iRows; iII++ ){
      if( iII < iRows0 ){
        mExtremeIC( iII, _ ) = mExtremeIC0( iII, _ ) ;
      }
      else{
        mExtremeIC( iII, _ ) = mExtremeIC1( iII - iRows0, _ ) ;
      }
    }
          // Fill the output matrix by concatenating mExtremeIC0 and mExtremeIC1
    if( boICinSet )
      mExtremeIC( iRows, _ ) = vW ;
          // Add the IC point to the set of extreme points if contained in set
    if( iRows==0 && boICinSet==false ){
      NumericMatrix mOut( 0, 2 ) ;
      return( mOut ) ;
    }
          // If mExtreme is empty, return a 0-row empty matrix. Means
          // this can be handled as a matrix in future without errors
    NumericMatrix mOut = sets_unique_cpp( mExtremeIC, 1e-14 ) ;
    if( mOut.nrow() > 4 )
      throw(Rcpp::exception("abSan.C should not return more than 4 values","abreuSannikov.cpp", 97) ) ;
    return( mOut ) ;
  }
}

MatrixXd abSan_rComponent_cpp_eigen( int iA, const List set, const NumericVector vPun, const List model, 
      const bool print_output = true, bool boCountFrom1=true ){
// Compute the component of R for the given action
  if( boCountFrom1 ) iA-- ;
        // Transform from 1-based counting to 0-based
  const NumericMatrix cSet = abSan_C_cpp( iA, set, vPun, model, print_output, false ) ;
        // The outcome of C for the current action
  const int iPts = cSet.nrow() ;
        // Number of rows of cSet
  NumericMatrix mF = model["mF"] ;
        // Copy the stage payoffs
  const NumericVector vF = mF.row( iA ) ;
  eVecd evF( as<eVecd>( vF ) ) ;
  eMatd emCset( as<eMatd>( cSet ) ) ;
  const double delta = as<double>( model["delta"] ) ;
  MatrixXd mOut =  ( 1 - delta ) * VectorXd::Ones( iPts ) * evF.transpose() + delta * emCset ;
  return( mOut ) ;
        // Calculate the points in the sum
}

//[[Rcpp::export('abSan.rComponent')]]
NumericMatrix abSan_rComponent_cpp( int iA, const List set, const NumericVector vPun, const List model, 
      const bool print_output = true, bool boCountFrom1=true ){
// Wrapper for abSan_rComponent_cpp_eigen
  return( wrap( abSan_rComponent_cpp_eigen( iA, set, vPun, model, print_output, boCountFrom1 ) ) ) ;
        // Calculate the points in the sum
}

//[[Rcpp::export('abSan.rCore')]]
NumericMatrix abSan_rCore( const List set, const NumericVector vPun, const List model, 
      const bool print_output = true, bool boCountFrom1=true ){
// The core loop for computing Abreu & Sannikov's R function
  int iActs( as<int>( model["iJointActs"] ) ) ;
        // The number of joint actions
//  List lComponents( iActs ) ;
        // The list of individual components of R
  VectorXi vComponentRows( iActs) ;
  vComponentRows = VectorXi::Zero( iActs ) ;
        // The vector of the number of rows in each component of abSan.R
  MatrixXd mComponent ;
        // Initializing the matrix of returns from rComponent
  std::vector<MatrixXd> lComponents ;
        // Initializing the container which stores these returns
  for( int iII=0; iII<iActs; iII++ ){
    mComponent = abSan_rComponent_cpp_eigen( iII, set, vPun, model, print_output, false ) ;
          // Compute the component
    lComponents.push_back( mComponent ) ;
          // Record it
    vComponentRows[ iII ] = mComponent.rows() ;
          // Note the numer of rows
  }
  MatrixXd mOut( vComponentRows.sum(), 2 ) ;
        // Initializing the matrix of outputs
  int iRowsSoFar = 0 ;
        // Initializing the number of rows of mOut already filled
  for( int iII; iII<iActs; iII++ ){
    mOut.block( iRowsSoFar, 0, vComponentRows[ iII ], 2 ) = lComponents[ iII ] ;
          // Fill in the rows 
    iRowsSoFar = iRowsSoFar + vComponentRows[ iII ] ;
  }
  return wrap( mOut ) ;
}

//[[Rcpp::export('sets.bind')]]
NumericMatrix bind_par_cpp( const List lA ){
/* Takes a list of matrices and binds the rows together.  Equivalent to do.call(
 * "rbind, lA ") */
  const int iA = lA.length() ;
        // The number of matrices
  VectorXi iRows( iA ) ;
        // The numer of rows of each matrix
  std::vector<MatrixXd> elA ;
        // Initializing the container for the matrices in lA
  SEXP sexpThisA ;
  NumericMatrix mThisA ;
//  eMatd emThisA( as<eMatd>( mThisA ) ) ;
  for( int iII=0; iII<iA; iII++ ){
    sexpThisA = lA[ iII ] ;
          // Make sure the entry of lA is defined as an Eigen matrix
    mThisA = sexpThisA ;
    eMatd emThisA( as<eMatd>( mThisA ) );
    elA.push_back( emThisA ) ;
          // Record the matrix entry of lA
    iRows[ iII ] = emThisA.rows() ;
          // And the number of rows too
  }
  MatrixXd emA( iRows.sum(), 2) ;
        // Initialize the output matrix
  int iTotalRows = 0 ;
        // Counter for number of rows completed so far
  for( int iII=0; iII<iA; iII++ ){
    emA.block( iTotalRows, 0, iRows[ iII ], 2 ) = elA[ iII ] ;
          // Paste the Eigen-form matrix in elA to the appropriate rows of emA
    iTotalRows = iTotalRows + iRows[ iII ] ;
          // Update the counter for the rows
  }
  return wrap( emA ) ;
}