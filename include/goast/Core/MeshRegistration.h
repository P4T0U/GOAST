// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef __MESHREGISTRATION_H
#define __MESHREGISTRATION_H

#include <chrono>
#include <ctime>
#include <goast/Optimization.h>
#include "FEM.h"

//===========================================================================================================================
//  HELPER FUNCTIONS FOR RBM REGISTRATION
//===========================================================================================================================
template<typename MatType>
void getRotationMatrix( double alpha, double beta, double gamma, MatType& Q ) {
    Q.set(0,0, cos(beta) * cos(gamma) );
    Q.set(0,1, -sin(gamma) * cos(beta) );
    Q.set(0,2, sin(beta) );
    Q.set(1,0, sin(alpha) * sin(beta) * cos(gamma) + cos(alpha) * sin(gamma) );
    Q.set(1,1, cos(alpha) * cos(gamma) - sin(alpha) * sin(beta) * sin(gamma) );
    Q.set(1,2, -sin(alpha) * cos(beta) );
    Q.set(2,0, sin(alpha) * sin(gamma) - sin(beta) * cos(alpha) * cos(gamma) );
    Q.set(2,1, cos(alpha) * sin(beta) * sin(gamma) + sin(alpha) * cos(gamma) );
    Q.set(2,2, cos(alpha) * cos(beta) );
}
  
template<typename VecType>
void getAlphaDerivative( double alpha, double beta, double gamma, const VecType& y, VecType& deriv ) {
    deriv[0] = 0.;
    deriv[1] =  (cos(alpha) * sin(beta) * cos(gamma) - sin(alpha) * sin(gamma) ) * y[0] - ( sin(alpha) * cos(gamma) + cos(alpha) * sin(beta) * sin(gamma) ) * y[1] - ( cos(alpha) * cos(beta) ) * y[2];
    deriv[2] = (cos(alpha) * sin(gamma) + sin(beta) * sin(alpha) * cos(gamma) ) * y[0] - ( sin(alpha) * sin(beta) * sin(gamma) - cos(alpha) * cos(gamma) ) * y[1] - ( sin(alpha) * cos(beta) ) * y[2];
}
  
template<typename VecType>
void getBetaDerivative( double alpha, double beta, double gamma, const VecType& y, VecType& deriv ) {
    deriv[0] = -sin(beta) * cos(gamma) * y[0] + sin(gamma) * sin(beta) * y[1]  + cos(beta) * y[2];
    deriv[1] = sin(alpha) * cos(beta) * cos(gamma) * y[0] - sin(alpha) * cos(beta) * sin(gamma) * y[1] + sin(alpha) * sin(beta) * y[2];
    deriv[2] = - cos(beta) * cos(alpha) * cos(gamma) * y[0] + cos(alpha) * cos(beta) * sin(gamma) * y[1] - cos(alpha) * sin(beta) * y[2];
}
 
template<typename VecType>
void getGammaDerivative( double alpha, double beta, double gamma, const VecType& y, VecType& deriv ) {
    deriv[0] = -cos(beta) * sin(gamma) * y[0] - cos(gamma) * cos(beta) * y[1];
    deriv[1] = (cos(alpha) * cos(gamma) - sin(alpha) * sin(beta) * sin(gamma) ) * y[0] - ( cos(alpha) * sin(gamma) + sin(alpha) * sin(beta) * cos(gamma) ) * y[1];
    deriv[2] = (sin(alpha) * cos(gamma) + sin(beta) * cos(alpha) * sin(gamma) ) * y[0] + ( cos(alpha) * sin(beta) * cos(gamma) - sin(alpha) * sin(gamma) ) * y[1];
}

template<typename VecType>
void getFirstDerivativeOfRotationAngle( int i, double alpha, double beta, double gamma, const VecType& y, VecType& deriv ) {
    if( i == 0 ) getAlphaDerivative( alpha, beta, gamma, y, deriv );
    if( i == 1 ) getBetaDerivative( alpha, beta, gamma, y, deriv );
    if( i == 2 ) getGammaDerivative( alpha, beta, gamma, y, deriv );
  }

//===========================================================================================================================
// RIGID BODY MOTION REGISTRATION
//===========================================================================================================================

/**
 * \brief Least squares penalty on rigid body motions
 * \author Heeren
 * \tparam ConfiguratorType Container with datatypes
 *
 * For two corresponding meshes \f$X\f$ and \f$Y\f$ with positions \f$(x_i)_i\f$ and \f$(y_i)_i\f$ the energy is given by
 * \f[ E[Q] = \sum_i (x_i - G(y_i) )^2,\f]
 * where \f$G(x) = Qx + z\f$, \f$Q = Q[a,b,c]\f$ is a rotation matrix and \f$z = (d,e,f)\f$ is a translation.
 * The argument vector is given by \f$(a, b, c, d, e, f) \in \R^6\f$
 * In the current notation \f$X\f$ is the template mesh whereas \f$Y\f$ is the reference mesh.
 */
template<typename ConfiguratorType>
class LeastSquaresRBMEnergy
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {

  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const VectorType &_defGeom, _tempGeom;
  std::vector<int> _nodeIndices;
  int _seedPointIdx, _sizeOfNeighborhood;

public:
  LeastSquaresRBMEnergy( const VectorType &TempGeom, const VectorType &DefGeom ) : _defGeom( DefGeom ),
                                                                                   _tempGeom( TempGeom ) {
    for ( int i = 0; i < _defGeom.size() / 3; i++ )
      _nodeIndices.push_back( i );
  }

  LeastSquaresRBMEnergy( const VectorType &TempGeom, const VectorType &DefGeom, const std::vector<int> &nodeIndices )
          : _defGeom( DefGeom ), _tempGeom( TempGeom ), _nodeIndices( nodeIndices ) {}

  void apply( const VectorType &AnglesAndTranslation, RealType &Dest ) const {

    if ( AnglesAndTranslation.size() != 6 )
      throw BasicException( "LeastSquaresRBMEnergy::apply(): argument has wrong size!!!" );
    Dest = 0.;

    MatType Q;
    getRotationMatrix( AnglesAndTranslation[0], AnglesAndTranslation[1], AnglesAndTranslation[2], Q );
    VecType x, y, Gy;

    for ( int i = 0; i < _nodeIndices.size(); i++ ) {
      getXYZCoord<VectorType, VecType>( _tempGeom, x, _nodeIndices[i] );
      getXYZCoord<VectorType, VecType>( _defGeom, y, _nodeIndices[i] );
      Q.mult( y, Gy );
      for ( int j = 0; j < 3; j++ )
        x[j] -= Gy[j] + AnglesAndTranslation[3 + j];
      Dest += x.normSqr();
    }
  }

};

/**
 * \brief Gradient of LeastSquaresRBMEnergy
 * \author Heeren
 *
 * First derivative w.r.t. \f$(a,b,c,d,e,f)\f$ of \f$E[Q] = \sum_i (x_i - Gy_i)^2\f$,
 * where \f$G(x) = Q + z\f$, \f$Q = Q[a,b,c]\f$ is a rotation matrix and \f$z = (d,e,f)\f$ is a translation.
 */
template<typename ConfiguratorType>
class LeastSquaresRBMDerivative : public BaseOp< typename ConfiguratorType::VectorType > {
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  
  const VectorType& _defGeom, _tempGeom;
  std::vector<int> _nodeIndices;
  
public:
  LeastSquaresRBMDerivative( const VectorType& TempGeom, const VectorType& DefGeom ) : _defGeom( DefGeom), _tempGeom( TempGeom ){
    for( int i = 0; i < _defGeom.size() / 3; i++ )
      _nodeIndices.push_back( i );
  }
  
  LeastSquaresRBMDerivative( const VectorType& TempGeom, const VectorType& DefGeom, const std::vector<int>& nodeIndices ) : _defGeom( DefGeom), _tempGeom( TempGeom ), _nodeIndices( nodeIndices ){} 
  
  void apply( const VectorType& AnglesAndTranslation, VectorType& Dest ) const {
    
    if( AnglesAndTranslation.size() != 6 )
      throw BasicException ( "LeastSquaresRBMDerivative::apply(): argument has wrong size!!!" );
    if( Dest.size() != 6 )
        Dest.resize(6);
    Dest.setZero();
    
    MatType Q;
    getRotationMatrix( AnglesAndTranslation[0], AnglesAndTranslation[1], AnglesAndTranslation[2], Q );
    VecType x, y, z;    
    
    for( int i = 0; i < _nodeIndices.size(); i++ ){
      getXYZCoord<VectorType, VecType>( _tempGeom, x, _nodeIndices[i]);
      getXYZCoord<VectorType, VecType>( _defGeom, y, _nodeIndices[i]);    
      Q.mult(y, z);
      for( int j = 0; j < 3; j++ )
	x[j] -= z[j] + AnglesAndTranslation[3 + j];

      // derivative w.r.t rotation
      for( int j = 0; j < 3; j++ ){
        getFirstDerivativeOfRotationAngle( j, AnglesAndTranslation[0], AnglesAndTranslation[1], AnglesAndTranslation[2], y, z );
        Dest[j] -= 2. * dotProduct( x, z );  
      }
      
      // derivatives w.r.t. translation
      for( int j = 0; j < 3; j++ )
	Dest[3 + j] -= 2. * x[j];
    }
  }

};


/*
// Second derivative w.r.t. (a,b,c) of E[Q] = \sum_i (x_i - Qy_i)^2, where Q = Q[a,b,c] is a rotation matrix
//TODO complete!
template<typename ConfiguratorType>
class RotationAboutOriginHessian : public BaseOp< typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType > {
    
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  
  const VectorType& _defGeom, _tempGeom;
  std::vector<int> _nodeIndices;
  
public:
  RotationAboutOriginHessian( const VectorType& TempGeom, const VectorType& DefGeom ) : _defGeom( DefGeom), _tempGeom( TempGeom ){
    for( int i = 0; i < _defGeom.size() / 3; i++ )
      _nodeIndices.push_back( i );
  }
  
  RotationAboutOriginHessian( const VectorType& TempGeom, const VectorType& DefGeom, const std::vector<int>& nodeIndices ) : _defGeom( DefGeom), _tempGeom( TempGeom ), _nodeIndices( nodeIndices ){} 
  
  void apply( const VectorType& Angles, MatrixType& Dest ) const{
    if( Angles.size() != 3 )
      throw BasicException ( "RotationAboutOriginHessian::applyAdd: wrong number of entries!!!", __FILE__, __LINE__ );
    
    MatType Q;
    getRotationMatrix( Angles[0], Angles[1], Angles[2], Q );
    VecType x, y, z;    
    
    for( int i = 0; i < _nodeIndices.size(); i++ ){
      getXYZCoord<VectorType, VecType>( _tempGeom, x, _nodeIndices[i]);
      getXYZCoord<VectorType, VecType>( _defGeom, y, _nodeIndices[i]);       
      Q.mult(y, z);
      x -= z;
      
      std::vector< VecType > derivatives(3);      
      getAlphaDerivative( Angles[0], Angles[1], Angles[2], y, derivatives[0] );
      getBetaDerivative( Angles[0], Angles[1], Angles[2], y, derivatives[1] );
      getGammaDerivative( Angles[0], Angles[1], Angles[2], y, derivatives[2] );
      
      //TODO \sum_i 2 \partial_k (Qy_i) * \partial_j (Qy_i) - 2 (x_i - Qy_i) * \partial_{kj}^2 (Qy_i)
      for( int j = 0; j < 3; j ++){
	for( int k = 0; k < 3; k++ ){
	  RealType aux = derivatives[j] * derivatives[k];
	  Dest.add( j, k, 2.*aux);
	}
      }
    }
  }
  
  static void getAlphaAlphaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = 0.;
    deriv[1] =  -(sin(alpha) * sin(beta) * cos(gamma) + cos(alpha) * sin(gamma) ) * y[0] - ( cos(alpha) * cos(gamma) - sin(alpha) * sin(beta) * sin(gamma) ) * y[1] + ( sin(alpha) * cos(beta) ) * y[2];
    deriv[2] = -1. * (sin(alpha) * sin(gamma) - sin(beta) * cos(alpha) * cos(gamma) ) * y[0] - ( cos(alpha) * sin(beta) * sin(gamma) + sin(alpha) * cos(gamma) ) * y[1] - ( cos(alpha) * cos(beta) ) * y[2];
  }
  
  static void getAlphaBetaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = 0.;
    deriv[1] = cos(alpha) * cos(beta) * cos(gamma) * y[0] - cos(alpha) * cos(beta) * sin(gamma) * y[1] + cos(alpha) * sin(beta) * y[2];
    deriv[2] = cos(beta) * sin(alpha) * cos(gamma) * y[0] - sin(alpha) * cos(beta) * sin(gamma) * y[1] + sin(alpha) * sin(beta) * y[2];
  }
    
  static void getAlphaGammaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = 0.;
    deriv[1] = -1. * (cos(alpha) * sin(beta) * sin(gamma) + sin(alpha) * cos(gamma) ) * y[0] + ( sin(alpha) * sin(gamma) - cos(alpha) * sin(beta) * cos(gamma) ) * y[1];
    deriv[2] = (cos(alpha) * cos(gamma) - sin(beta) * sin(alpha) * sin(gamma) ) * y[0] - ( sin(alpha) * sin(beta) * cos(gamma) + cos(alpha) * sin(gamma) ) * y[1];
  }
  
  static void getBetaAlphaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = 0;
    deriv[1] = cos(alpha) * cos(beta) * cos(gamma) * y[0] - cos(alpha) * cos(beta) * sin(gamma) * y[1] + cos(alpha) * sin(beta) * y[2];
    deriv[2] = cos(beta) * sin(alpha) * cos(gamma) * y[0] - sin(alpha) * cos(beta) * sin(gamma) * y[1] + sin(alpha) * sin(beta) * y[2];
  }
  
  static void getBetaBetaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = -cos(beta) * cos(gamma) * y[0] + sin(gamma) * cos(beta) * y[1]  - sin(beta) * y[2];
    deriv[1] = -sin(alpha) * sin(beta) * cos(gamma) * y[0] + sin(alpha) * sin(beta) * sin(gamma) * y[1] + sin(alpha) * cos(beta) * y[2];
    deriv[2] =  sin(beta) * cos(alpha) * cos(gamma) * y[0] - cos(alpha) * sin(beta) * sin(gamma) * y[1] - cos(alpha) * cos(beta) * y[2];
  }  
  
  static void getBetaGammaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = sin(beta) * sin(gamma) * y[0] + cos(gamma) * sin(beta) * y[1];
    deriv[1] = -sin(alpha) * cos(beta) * sin(gamma) * y[0] - sin(alpha) * cos(beta) * cos(gamma) * y[1];
    deriv[2] = cos(beta) * cos(alpha) * sin(gamma) * y[0] + cos(alpha) * cos(beta) * cos(gamma) * y[1];
  }
  
  static void getGammaAlphaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = 0.;
    deriv[1] = -1. * (sin(alpha) * cos(gamma) + cos(alpha) * sin(beta) * sin(gamma) ) * y[0] + ( sin(alpha) * sin(gamma) - cos(alpha) * sin(beta) * cos(gamma) ) * y[1];
    deriv[2] = (cos(alpha) * cos(gamma) - sin(beta) * sin(alpha) * sin(gamma) ) * y[0] - ( sin(alpha) * sin(beta) * cos(gamma) + cos(alpha) * sin(gamma) ) * y[1];
  }
  
  static void getGammaBetaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = sin(beta) * sin(gamma) * y[0] + cos(gamma) * sin(beta) * y[1];
    deriv[1] = -sin(alpha) * cos(beta) * sin(gamma) * y[0] - sin(alpha) * cos(beta) * cos(gamma) * y[1];
    deriv[2] = cos(beta) * cos(alpha) * sin(gamma) * y[0] + cos(alpha) * cos(beta) * cos(gamma) * y[1];
  }  
  
  static void getGammaGammaDerivative( RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    deriv[0] = -cos(beta) * cos(gamma) * y[0] + sin(gamma) * cos(beta) * y[1];
    deriv[1] = -(cos(alpha) * sin(gamma) + sin(alpha) * sin(beta) * cos(gamma) ) * y[0] - ( cos(alpha) * cos(gamma) - sin(alpha) * sin(beta) * sin(gamma) ) * y[1];
    deriv[2] = -(sin(alpha) * sin(gamma) - sin(beta) * cos(alpha) * cos(gamma) ) * y[0] - ( cos(alpha) * sin(beta) * sin(gamma) + sin(alpha) * cos(gamma) ) * y[1];
  }
  
  static void getSecondDerivative( int i, int j, RealType alpha, RealType beta, RealType gamma, const Vec3<RealType>& y, VecType& deriv ) {
    // symmetry
    if( i > j )
      getSecondDerivative( j, i, alpha, beta, gamma, y, deriv );
      
    if( i == 0 ){
      if ( j == 0 ) getAlphaAlphaDerivative( alpha, beta, gamma, y, deriv );
      if ( j == 1 ) getAlphaBetaDerivative( alpha, beta, gamma, y, deriv );
      if ( j == 2 ) getAlphaGammaDerivative( alpha, beta, gamma, y, deriv );
    }
    
    if( i == 1 ){
      if ( j == 1 ) getBetaBetaDerivative( alpha, beta, gamma, y, deriv );
      if ( j == 2 ) getBetaGammaDerivative( alpha, beta, gamma, y, deriv );
    }
    
    if( (i == 2) && (j == 2) ) getGammaGammaDerivative( alpha, beta, gamma, y, deriv );

  }

};
*/

//===========================================================================================================================
// BARYCENTER REGISTRATION
//===========================================================================================================================

/**
 * \brief Register barycenter wrt. \f$l^2\f$ metric
 * \author Heeren
 *
 * For a given template mesh \f$T\f$ with n vertices \f$x_i\f$  the \f$l^2\f$-barycenter is given by \f[b[T] = 1/n \sum_i x_i.\f]
 * A free mesh \f$M\f$ is translated to \f$M'\f$ such that \f$b[M'] = b[T]\f$.
 */
template<typename ConfiguratorType>
class BarycenterRegistrationOp{ 
    
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  
  VecType _tempBarycenter;
  int _numNodes;
  
public:
  BarycenterRegistrationOp( const VectorType& TempGeom ) : _tempBarycenter(0.,0.,0.), _numNodes(TempGeom.size() / 3){ 
    for( int i = 0; i < _numNodes; i++ ){
      VecType node;
      getXYZCoord<VectorType, VecType>( TempGeom, node, i);
      _tempBarycenter += node;
    }
    _tempBarycenter *= 1. / _numNodes;
  }
  
  // argument is geometry of free mesh M, solution is translated mesh M' such that b[M'] = b[T], where b denotes the l^2-barycenter.
  void execute( VectorType& FreeGeom ) const {  
      
    if( FreeGeom.size() != 3 * _numNodes )
        throw BasicException("BarycenterRegistrationOp::execute() argument has wrong size!");
      
    // compute free barycenter
    VecType freeBarycenter(0.,0.,0.);
    for( int i = 0; i < _numNodes; i++ ){
      VecType node;
      getXYZCoord<VectorType, VecType>( FreeGeom, node, i);
      freeBarycenter += node;
    }
    freeBarycenter *= 1. / _numNodes;
    
    // compute shift
    freeBarycenter -= _tempBarycenter;
    
    // apply shift
    for( int i = 0; i < _numNodes; i++ ){
      VecType node;
      getXYZCoord<VectorType, VecType>( FreeGeom, node, i);
      node -= freeBarycenter;
      setXYZCoord<VectorType, VecType>( FreeGeom, node, i);
    }
  }
  
};

//! \brief Registration by minimizing the LeastSquaresRBMEnergy (see above) via gradient descent
//! \author Heeren
template<typename ConfiguratorType>
class LeastSquaresRBMRegistrationOp{ 
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  
  const VectorType& _tempGeom;
  std::vector<int> _nodeIndices;
  int _numOfIterations;
  RealType _eps;
  bool _quiet;
  std::vector<int> _fixedDOFs;
  
public:
  LeastSquaresRBMRegistrationOp( const VectorType& TempGeom, int numOfIterations = 1000, RealType Eps = 1e-5, bool Quiet = true ) : _tempGeom( TempGeom ), _nodeIndices( TempGeom.size()/3 ),_numOfIterations( numOfIterations ), _eps( Eps ), _quiet(Quiet){ 
    for( int i = 0; i < TempGeom.size() / 3; i++ )
      _nodeIndices[i] = i;
  }
  
  LeastSquaresRBMRegistrationOp( const VectorType& TempGeom, const std::vector<int>& nodeIndices, int numOfIterations = 1000, RealType Eps = 1e-5, bool Quiet = true ) : _tempGeom( TempGeom ), _nodeIndices(nodeIndices), _numOfIterations( numOfIterations ), _eps( Eps ), _quiet(Quiet){ }
  
  void setBoundaryMask( const std::vector<int>& fixedDOFs ) {
    _fixedDOFs.resize( fixedDOFs.size() );
    _fixedDOFs = fixedDOFs;
  }
  
  void execute( VectorType& FreeGeom ) const {  
      VectorType OptimalRBMdofs;
      execute( FreeGeom, OptimalRBMdofs );
  }
  
  // arguments are free mesh and the optimal rbm degrees of freedeom corresponding two three translations and three angles
  void execute( VectorType& FreeGeom, VectorType& OptimalRBMdofs ) const {  

    // optimize
    LeastSquaresRBMEnergy<ConfiguratorType> E( _tempGeom, FreeGeom, _nodeIndices );
    LeastSquaresRBMDerivative<ConfiguratorType> DE ( _tempGeom, FreeGeom, _nodeIndices );
    
    // solve with gradient descent
    OptimalRBMdofs.resize( 6 );
    VectorType initialGuess(6);
    initialGuess.setZero();
    GradientDescent<ConfiguratorType> gradDesc( E, DE, _numOfIterations, _eps, ARMIJO, _quiet );
    if( _fixedDOFs.size() > 0 )
        gradDesc.setBoundaryMask( _fixedDOFs );
    gradDesc.solve( initialGuess, OptimalRBMdofs );
    
    // apply optimal rbm
    MatType Q;
    getRotationMatrix( OptimalRBMdofs[0], OptimalRBMdofs[1], OptimalRBMdofs[2], Q );
    for( int j = 0; j < FreeGeom.size() / 3; j++ ){
      VecType node, aux;
      getXYZCoord<VectorType, VecType>( FreeGeom, node, j);
      Q.mult( node, aux );
      for( int k = 0; k < 3; k++ )
        aux[k] += OptimalRBMdofs[ 3 + k ];
      setXYZCoord<VectorType, VecType>( FreeGeom, aux, j);
    }
  }

};


//===========================================================================================================================
// MOMENTUM REGISTRATION
//===========================================================================================================================

/**
 * \brief Registration energy to measure difference in zeroth and first momentum (each in all three spatial components)
 * \author Heeren
 *
 * Let \f$ x = (x^0, x^1, x^2) \f$ the geometry of a fixed reference (resp. template) mesh with $n$ vertices, where \f$ x^i \in \R^n \f$.
 * Let \f$ y = (y^0, y^1, y^2) \f$ the geometry of the mesh that is supposed to be registered, where \f$ y^i \in \R^n \f$.
 * For \f$ z \in R^3 \f$, let \f$ r[z] := Qz + b \f$ a rigid motion, determined by the matrix \f$ Q \in SO(3) \f$ and translation vector \f$ b \in R^3 \f$.
 * Note that $Q$ and $b$ are described by three angles and three scalars, respectively, hence there are 6 dofs in total.
 * Let \f$M \in R^{n,n}\f$ be a mass matrix (of \f$x\f$) and $1 \in R^n$ a vector containing only ones.
 * The zeroth momentum wrt. \f$x\f$ and \f$y\f$ is given by \f$M(x^i - y^i)1\f$, \f$i = 0,1,2\f$,
 * and the first momentum by \f$M(x^i - y^i)x^{i+1} - M(x^{i+1} - y^{i+1})x^i\f$, \f$i = 0,1,2 (mod 3)\f$
 *
 * The argument of the energy are the six dofs describing the angles and translation.
 * The first three components of the energy's result represent the mismatch wrt. the zeroth momentum, i.e.  \f$M(x^i - r[y]^i)1\f$ for \f$i = 0,1,2\f$,
 * and the last three components represent the mismatch wrt. the first momentum, i.e.  \f$M(x^i - (r[y])^i)x^{i+1} - M(x^{i+1} - (r[y])^{i+1})x^i\f$ for \f$i = 0,1,2 (mod 3)\f$.
 */
template<typename ConfiguratorType>
class MomentumRegistrationEnergy : public BaseOp< typename ConfiguratorType::VectorType > {
 
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  
  int _numV;
  const VectorType& _tempGeom, _defGeom;
  VectorType _massOpAllOnes, _massOpIdentity;
  
public:
  MomentumRegistrationEnergy(  const MeshTopologySaver& Topology, const VectorType& TempGeom, const VectorType& DefGeom ) 
      :  _numV( Topology.getNumVertices() ), _tempGeom( TempGeom ), _defGeom( DefGeom ){ 
       MatrixType LumpedMassMatrix;  
       computeNodalAreas<ConfiguratorType>(Topology, TempGeom, LumpedMassMatrix, true );
       _massOpAllOnes = LumpedMassMatrix * VectorType::Constant( 3 * _numV, 1. );
       _massOpIdentity = LumpedMassMatrix * TempGeom;
      }
      
  MomentumRegistrationEnergy(  const  MatrixType& LumpedMassMatrix, const VectorType& TempGeom, const VectorType& DefGeom ) 
      :  _numV( TempGeom.size() / 3 ), _tempGeom( TempGeom ), _defGeom( DefGeom){ 
       _massOpAllOnes = LumpedMassMatrix * VectorType::Constant( 3 * _numV, 1. );
       _massOpIdentity = LumpedMassMatrix * TempGeom;
  }
 
  // Argument are rotation and translation dofs, Dest are mass momentums and angular momentums
  void apply( const VectorType& AnglesAndTranslation, VectorType& Dest ) const {
    
    if( AnglesAndTranslation.size() != 6 )
      throw BasicException ( "MomentumRegistrationEnergy::apply(): wrong number of arguments!!!" );
    if( Dest.size() != 6 )
      Dest.resize( 6 );
    Dest.setZero();
    
    MatType Q;
    getRotationMatrix( AnglesAndTranslation[0], AnglesAndTranslation[1], AnglesAndTranslation[2], Q );
    VecType x, trans( AnglesAndTranslation[3], AnglesAndTranslation[4], AnglesAndTranslation[5]);
    
    // run over all vertices
    for( int i = 0; i < _numV; i++ ){
      getXYZCoord<VectorType, VecType>( _defGeom, x, i);
      VecType Qx = Q * x;
      Qx += trans;
      for( int k = 0; k < 3; k++ ){
        // mass momentum
        Dest[k]   += _massOpAllOnes[i] * (_tempGeom[k*_numV + i] - Qx[k]);
        // angular momentum
        int nextIdx = (k+1)%3;
        Dest[3+k] += _massOpIdentity[k*_numV + i] * Qx[nextIdx] - _massOpIdentity[nextIdx*_numV + i] * Qx[k];
      }
    }
  }  

};

/**
 * \brief Derivative of MomentumRegistrationEnergy (see above)
 * \author Heeren
 *
 * Since the "energy" measures differences in zeroth and first momentum in all three spatial components,
 * the energy's result is a vector with six components.
 * Hence the derivative is a 6x6-matrix.
 */
template<typename ConfiguratorType>
class MomentumRegistrationDerivative : public BaseOp< typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType > {
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType; 
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  int _numV;
  const VectorType& _tempGeom, _defGeom;
  VectorType _massOpAllOnes, _massOpIdentity;
  
public:
  MomentumRegistrationDerivative(  const MeshTopologySaver& Topology, const VectorType& TempGeom, const VectorType& DefGeom ) 
      :  _numV( Topology.getNumVertices() ), _tempGeom( TempGeom ), _defGeom( DefGeom ){ 
       MatrixType LumpedMassMatrix;  
       computeNodalAreas<ConfiguratorType>(Topology, TempGeom, LumpedMassMatrix, true );
       _massOpAllOnes = LumpedMassMatrix * VectorType::Constant( 3 * _numV, 1. );
       _massOpIdentity = LumpedMassMatrix * TempGeom;
      }
      
  MomentumRegistrationDerivative(  const  MatrixType& LumpedMassMatrix, const VectorType& TempGeom, const VectorType& DefGeom ) 
      :  _numV( TempGeom.size() / 3 ), _tempGeom( TempGeom ), _defGeom( DefGeom){ 
       _massOpAllOnes = LumpedMassMatrix * VectorType::Constant( 3 * _numV, 1. );
       _massOpIdentity = LumpedMassMatrix * TempGeom;
  }
 
  void apply( const VectorType& AnglesAndTranslation, MatrixType& Dest ) const {
    
    if( AnglesAndTranslation.size() != 6 )
      throw BasicException ( "MomentumRegistrationDerivative::apply(): wrong number of arguments!!!");
    if( Dest.rows() != 6 || Dest.cols() != 6 )
        Dest.resize( 6, 6 );
    Dest.setZero();
    
    TripletListType triplets;
    pushTriplets( AnglesAndTranslation, triplets );
    Dest.setFromTriplets( triplets.cbegin(), triplets.cend() );

  }  
  
  void pushTriplets( const VectorType& AnglesAndTranslation, TripletListType& triplets ) const {
      
    if( AnglesAndTranslation.size() != 6 )
      throw BasicException ( "MomentumRegistrationDerivative::apply(): wrong number of arguments!!!");
      
    triplets.clear();
    triplets.reserve( 36 );      
          
    VecType x, DQx;
    
    // run over all vertices
    for( int i = 0; i < _numV; i++ ){
      getXYZCoord<VectorType, VecType>( _defGeom, x, i);
      
      for( int k = 0; k < 3; k++ ){
        getFirstDerivativeOfRotationAngle( k, AnglesAndTranslation[0], AnglesAndTranslation[1], AnglesAndTranslation[2], x, DQx );
        
        // mass momentum: deriv. wrt. translation
        triplets.push_back( TripletType( k, 3+k, -1. * _massOpAllOnes[i]) );
        
        // angular momentum: deriv. wrt. translation
        int nextIdx = (k+1)%3;
        int prevIdx = (k+2)%3;
        triplets.push_back( TripletType( 3+k,       3+k, -1. * _massOpIdentity[nextIdx*_numV + i] ) );
        triplets.push_back( TripletType( 3+prevIdx, 3+k,       _massOpIdentity[prevIdx*_numV + i] ) );
        
        for( int l = 0; l < 3; l++ ){            
          // mass momentum: deriv. wrt. angles
          triplets.push_back( TripletType( l, k, -1. * _massOpAllOnes[i] * DQx[l] ) );
          
          // angular momentum: deriv. wrt. angles
          nextIdx = (l+1)%3;
          triplets.push_back( TripletType( 3+l, k, _massOpIdentity[l*_numV +i] * DQx[nextIdx] - _massOpIdentity[nextIdx*_numV + i] * DQx[l] ) );          
        }
      }
    }
      
  }


};

/**
 * \brief Operator to register meshes such that their zeroth and first momentum coincide.
 * \author Heeren
 *
 * Done by finding root of MomentumRegistrationEnergy by means of Newton's method.
 * Note: The "energy" measures differences in zeroth and first momentum in all three spatial components,
 * hence the energy's result is a vector with six components.
 */
template<typename ConfiguratorType>
class MomentumRegistrationOperator{
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType; 
  
  const MeshTopologySaver& _topology;
  const VectorType& _tempGeom;
  MatrixType _lumpedMass;
  int _numOfIterations;
  RealType _eps;
  bool _quiet;
  
public:
  MomentumRegistrationOperator( const MeshTopologySaver& Topology, const VectorType TempGeom, int numOfIterations = 1000, RealType Eps = 1e-8, bool Quiet = true ) 
    : _topology(Topology), _tempGeom( TempGeom ), _numOfIterations( numOfIterations ), _eps( Eps ), _quiet(Quiet){ 
      computeLumpedMassMatrix<ConfiguratorType>(Topology, TempGeom, _lumpedMass, true );
  }  
  
  // arguments are free mesh and the optimal rbm degrees of freedeom corresponding two three translations and three angles
  void execute( VectorType& FreeGeom ) const {  
    VectorType OptimalRBMdofs;
    execute( FreeGeom, OptimalRBMdofs );
  }
    
  //  
  // arguments are free mesh and the optimal rbm degrees of freedeom corresponding two three translations and three angles
  void execute( VectorType& FreeGeom, VectorType& OptimalRBMdofs ) const {  

    // optimize
    MomentumRegistrationEnergy<ConfiguratorType> E( _lumpedMass, _tempGeom, FreeGeom );
    MomentumRegistrationDerivative<ConfiguratorType> DE( _lumpedMass, _tempGeom, FreeGeom );
    
    // solve with gradient descent
    OptimalRBMdofs.resize( 6 );
    VectorType initialGuess(6);
    initialGuess.setZero();
    
    OptimizationParameters<ConfiguratorType> optPars;
    optPars.setNewtonIterations( _numOfIterations );
    optPars.setEpsilon( _eps );
    if( !_quiet ) optPars.setQuietMode( SHOW_ALL );

    //std::cerr << "Start derivative tests" << std::endl;
    //VectorValuedDerivativeTester<DefaultConfigurator> ( E, DE, 1e-5 ).plotAllDirections ( initialGuess, "rbm_hessTest" );
    //return;
    
    NewtonMethod<ConfiguratorType>( E, DE, optPars ).solve( initialGuess, OptimalRBMdofs ); 
    
    // apply optimal rbm
    MatType Q;
    getRotationMatrix( OptimalRBMdofs[0], OptimalRBMdofs[1], OptimalRBMdofs[2], Q );
    for( int j = 0; j < _topology.getNumVertices(); j++ ){
      VecType node, aux;
      getXYZCoord<VectorType, VecType>( FreeGeom, node, j);
      Q.mult( node, aux );
      for( int k = 0; k < 3; k++ )
        aux[k] += OptimalRBMdofs[ 3 + k ];
      setXYZCoord<VectorType, VecType>( FreeGeom, aux, j);
    }
  }
   
    
};

#endif
