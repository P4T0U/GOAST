// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Full disscrete geodesic calculus for low-dimensional, embedded manifolds.
 * \author Effland, Heeren
 * 
 * \note This file has been copied and adapted from Alexander Effland's QuocMesh code.
 * 
 * Let \f$ \Omega \subset \R^2 \f$ a compact parameter domain, \f$ E: \Omega \to \R^3 \f$,
 * a parametrization of a regular surface. Here, we consider the sphere and the torus.
 * On \f$ \Omega \f$, we define \f$ W[x,y] = 1/2 \| E(x) - E(y) \|^2 \f$ and provide the 
 * basic ingredients for the discrete geodesic calculus for this finite dimensional manifold.* 
 * 
 * 
 */

#ifndef __FINITEDIMDISCRETEGEODCALC_H
#define __FINITEDIMDISCRETEGEODCALC_H

#include <iostream>

#include <goast/Core.h>
#include <goast/Optimization.h>

//========================================================================================================
//! \brief Trait to handle fixed points in variational setup (e.g. end points when computing geodesics)
//! \author Effland, Heeren
template < typename ConfiguratorType >
class FixedPointsTrait {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  VectorType _startPoint;
  VectorType _endPoint;
  bool _pointsSet;

public:
  FixedPointsTrait ( ) : _pointsSet ( false ) { }

  FixedPointsTrait ( const VectorType& StartPoint, const VectorType& FixedPoint )
  : _startPoint ( StartPoint ), _endPoint ( FixedPoint ), _pointsSet ( true ) { }

  const VectorType& getStartVector ( ) const {
    check ( );
    return _startPoint;
  }

  const VectorType& getEndVector ( ) const {
    check ( );
    return _endPoint;
  }

  void setVectors ( const VectorType& StartPoint, const VectorType& FixedPoint ) {
    _startPoint = StartPoint;
    _endPoint = FixedPoint;
    _pointsSet = true;
  }

  void unsetVectors ( ) {
    _pointsSet = false;
  }

private:
  void check ( ) const {
    if ( !_pointsSet )
      throw BasicException( "FixedPointsTrait::check(): Vector not set!" );
  }
};


//========================================================================================================
//! \brief Linear interpolation between start and end point
//! \author Heeren
template <typename VectorType>
void linspace( const VectorType& Start, const VectorType& End, VectorType& linPath, int K ) {
  int dim = Start.size();
  linPath.resize( (K-1) * dim );
  for( int k = 1; k < K; k++ )
      for( int i = 0; i < dim; i++ )
          linPath[(k-1)*dim+i] = ( (K-k) * Start[i] + k * End[i] ) / (1.*K);
}

//========================================================================================================
//! \brief Class to store fourth order tensor
//! \author Effland
template < typename RealType, int DimOfDomain, int DimOfRange >
class FourthTensor {
  RealType _data[DimOfRange][DimOfDomain][DimOfDomain][DimOfDomain];

public:
  FourthTensor ( ) {
    setZero ( );
  }

  void setZero ( ) {
    for ( int i = 0; i < DimOfRange; ++i )
      for ( int j = 0; j < DimOfDomain; ++j )
        for ( int k = 0; k < DimOfDomain; ++k )
          for ( int l = 0; l < DimOfDomain; ++l )
            _data[i][j][k][l] = 0.;
  }

  void set ( const int i, const int j, const int k, const int l, const RealType Value ) {
    _data[i][j][k][l] = Value;
  }

  RealType get ( const int i, const int j, const int k, const int l ) {
    return _data[i][j][k][l];
  }

  void setSymmetric ( const int i, const int numberOfYDerivatives, const RealType Value ) {
    unsigned short int positionOfY[3] = { 0, 0, 0 };
    for ( unsigned short int i = 0; i < numberOfYDerivatives; ++i )
      positionOfY[i] = 1;
    std::sort ( positionOfY, positionOfY + 3 );
    do {
      _data[i][positionOfY[0]][positionOfY[1]][positionOfY[2]] = Value;
    } while ( std::next_permutation ( positionOfY, positionOfY + 3 ) );
  }

  void writeData ( ) const {
    for ( int i = 0; i < DimOfRange; ++i )
      for ( int j = 0; j < DimOfDomain; ++j )
        for ( int k = 0; k < DimOfDomain; ++k )
          for ( int l = 0; l < DimOfDomain; ++l )
            std::cout << "(" << i << ", " << j << ", " << k << ", " << l << ") = " << _data[i][j][k][l] << std::endl;
  }
};


//========================================================================================================
//========================================================================================================
//========================================================================================================

//! \brief Interface for a generic surface parametrization \f$ E: \Omega \to \R^n \f$
//! \author Effland and Heeren
//! Here \f$ \Omega \subset \R^m \$ where m is dimension of the domain and n is the dimension of range.
//! Derived class as template argument (Barton-Nackman trick).
template < typename ConfiguratorType, unsigned short int DimOfDomain, unsigned short int DimOfRange, typename Imp >
class SurfaceParametrizationInterface {
public:
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::VectorType VectorType;       
   typedef typename ConfiguratorType::FullMatrixType MatrixType;

  void applyAdd ( const VectorType& Arg, VectorType& Dest ) const {
    if ( Arg.size ( ) != getDimOfDomain ( ) )
      throw BasicException ( "SurfaceParametrizationInterface::applyAdd(): arg size mismatch!" );
    if( Dest.size() != getDimOfRange ( ) )
        Dest.resize( getDimOfRange() );
    this->asImp().computeParametrization ( Arg, Dest );
  }

  //! compute parametrization \f$ E(x) \in \R^n \f$
  void apply ( const VectorType& Arg, VectorType& Dest ) const {
    if( Dest.size() != getDimOfRange ( ) )
        Dest.resize( getDimOfRange() );
    Dest.setZero ( );
    applyAdd ( Arg, Dest );
  }

  void applyAddInverse ( const VectorType& Arg, VectorType& Dest ) const {
    if ( Arg.size ( ) != getDimOfRange ( ) )
      throw BasicException ( "SurfaceParametrizationInterface::applyAddInverse(): arg size mismatch!" );
    if ( Dest.size ( ) != getDimOfDomain ( ) )
       Dest.resize( getDimOfDomain() );
    this->asImp().computeInverseParametrization ( Arg, Dest );

    VectorType test;
    apply( Dest, test );
    test -= Arg;
    if ( test.normSqr ( ) > 1e-8 ) {
      std::cerr << "Difference of both vectors = " << test.norm ( ) << "\n" << "\n";
      throw BasicException ("SurfaceParametrizationInterface::applyAddInverse(): applyAddInverse caused problem!" );
    }
  }

  //!
  void applyInverse ( const VectorType& Arg, VectorType& Dest ) const {
    Dest.setZero ( );
    applyAddInverse ( Arg, Dest );
  }

  //! \f$ DE(x) \in \R^{n,m} \f$
  void getDE ( const VectorType& Arg, MatrixType& Dest ) const {
    if( Arg.size ( ) != getDimOfDomain ( ) )
      throw BasicException ( "SurfaceParametrizationInterface::getDE(): arg size mismatch!" );
    if(  Dest.rows() != getDimOfRange ( ) || Dest.cols() != getDimOfDomain ( ) )
        Dest.resize( getDimOfRange ( ), getDimOfDomain ( ) );
    Dest.setZero ( );
    this->asImp().computeDE ( Arg, Dest );
  }

  //! metric \f$ g(x) = DE(x)^T DE(x) \in \R^{m,m} \f$
  void getMetricTensor ( const VectorType& Arg, MatrixType& Dest ) const {
    if( Arg.size ( ) != getDimOfDomain ( ) )
      throw BasicException ( "SurfaceParametrizationInterface::getMetricTensor(): Arg size mismatch!" );
    if( Dest.rows()!= getDimOfDomain ( ) || Dest.cols() != getDimOfDomain ( )  )
      Dest.resize( getDimOfDomain(), getDimOfDomain() );

    MatrixType DE ( getDimOfRange ( ), getDimOfDomain ( ) );
    getDE ( Arg, DE );
    Dest = DE.transpose() * DE;
  }

  //!
  void getD2E ( const VectorType& Arg, std::vector< MatrixType >& Dest ) const {
    if ( Arg.size () != getDimOfDomain( ) )
      throw BasicException ( "SurfaceParametrizationInterface::getD2E(): Arg size mismatch!" );
    if ( Dest.size () != getDimOfRange ( ) )
        Dest.resize( getDimOfRange() );
    for ( unsigned short int i = 0; i < Dest.size ( ); ++i ){
      Dest[i].resize( getDimOfDomain( ), getDimOfDomain( ) );
      Dest[i].setZero ( );
    }
    this->asImp().computeD2E ( Arg, Dest );
  }

  //!
  void getD3E ( const VectorType& Arg, FourthTensor < RealType, DimOfDomain, DimOfRange >& Dest ) const {
    if( Arg.size ( ) != getDimOfDomain ( ) )
      throw BasicException ( "SurfaceParametrizationInterface::getD3E(): Arg size mismatch!" );
    Dest.setZero ( );
    this->asImp().computeD3E ( Arg, Dest );
  }

  //! DE(Arg)^T * direction
  void applyAddDE ( const VectorType& Arg, const VectorType& Direction, VectorType& Dest ) const {
    if( Arg.size ( ) != getDimOfDomain ( ) || Direction.size ( ) != getDimOfRange ( ) )
      throw BasicException ( "SurfaceParametrizationInterface::applyAddDE(): Input size mismatch!" );
    if( Dest.size ( ) != getDimOfDomain ( ) ){
      Dest.resize( getDimOfDomain ( ) );
      Dest.setZero();
    }
    
    MatrixType derivative;
    getDE ( Arg, derivative );
    Dest += derivative.transpose() * Direction;
  }

  // vector-valued result DE(Arg1)^T * direction
  void applyDE ( const VectorType& Arg, const VectorType& Direction, VectorType& Dest ) const {
    Dest.resize( getDimOfDomain ( ) );
    Dest.setZero();
    applyAddDE ( Arg, Direction, Dest );
  }

  // matrix-valued result DE(Arg1)^T * DE(Arg2)
  void applyDE ( const VectorType& Arg1, const VectorType& Arg2, MatrixType& Dest ) const {
    if ( Arg1.size ( ) != getDimOfDomain ( ) )
        throw BasicException ( "SurfaceParametrizationInterface::applyDE(): First arg size mismatch!" );
    if ( Arg2.size ( ) != getDimOfDomain ( ) )
        throw BasicException ( "SurfaceParametrizationInterface::applyDE(): Second arg size mismatch!" );
    if( Dest.rows() != getDimOfDomain ( ) || Dest.cols() != getDimOfDomain ( ) )
        Dest.resize( getDimOfDomain ( ), getDimOfDomain ( ) );      

    MatrixType DArg1, DArg2;
    getDE ( Arg1, DArg1 );
    MatrixType DArg1T = DArg1.transpose();
    getDE ( Arg2, DArg2 );
    Dest = DArg1.transpose() * DArg2;
  }

  void applyAddDE ( const VectorType& Arg1, const VectorType& Arg2, MatrixType& Dest ) const {
    MatrixType temp;
    applyDE ( Arg1, Arg2, temp );
    Dest += temp;
  }

  void applyD2E ( const VectorType& Position, VectorType& Arg, MatrixType& Dest ) const {
    if( Position.size() !=  getDimOfDomain( ) )
        throw BasicException("SurfaceParametrizationInterface::applyD2E(): Position has wrong size!");
    if( Arg.size() !=  getDimOfRange( ) )
        throw BasicException("SurfaceParametrizationInterface::applyD2E(): Arg has wrong size!");
    if( Dest.rows() != getDimOfDomain ( ) || Dest.cols() != getDimOfDomain ( ) )
        Dest.resize( getDimOfDomain ( ), getDimOfDomain ( ) );
    Dest.setZero ( );
    
    std::vector< MatrixType > tensor( getDimOfRange ( ) );
    getD2E ( Position, tensor );    
    for ( unsigned short int i = 0; i < getDimOfDomain ( ); ++i )
      for( unsigned short int j = 0; j < getDimOfDomain ( ); ++j )
        for ( unsigned short int k = 0; k < getDimOfRange ( ); ++k )
          Dest( i, j ) += tensor[k]( i, j ) * Arg[k];
  }

  unsigned short int getDimOfDomain ( ) const {
    return this->asImp().DimOfDomain;
  }

  unsigned short int getDimOfRange ( ) const {
    return this->asImp().DimOfRange;
  }
  
/*
  // R(X,Y)Z
  template < typename DiscreteVectorFieldTypeX, typename DiscreteVectorFieldTypeY, typename DiscreteVectorFieldTypeZ >
  void applyRiemannCurvatureTensor ( const VectorType& BasePoint, const DiscreteVectorFieldTypeX& XField, const DiscreteVectorFieldTypeY& YField, const DiscreteVectorFieldTypeZ& ZField, VectorType& Dest ) const {
    if ( BasePoint.size ( ) != getDimOfDomain ( ) || Dest.size ( ) != getDimOfDomain ( ) )
      throw BasicException ( "Size mismatch!" );
    Dest.setZero ( );
    aol::Vector < RealType > X ( getDimOfDomain ( ) );
    XField.apply ( BasePoint, X );
    aol::Vector < RealType > Y ( getDimOfDomain ( ) );
    YField.apply ( BasePoint, Y );
    aol::Vector < RealType > Z ( getDimOfDomain ( ) );
    ZField.apply ( BasePoint, Z );
    for ( unsigned short int i = 0; i < getDimOfDomain ( ); ++i ) {
      for ( unsigned short int j = 0; j < getDimOfDomain ( ); ++j ) {
        for ( unsigned short int k = 0; k < getDimOfDomain ( ); ++k ) {
          for ( unsigned short int l = 0; l < getDimOfDomain ( ); ++l ) {
            RealType Rijkl = getDChristoffelTensor ( i, k, l, j, BasePoint ) - getDChristoffelTensor ( i, j, l, k, BasePoint );
            for ( unsigned short int s = 0; s < getDimOfDomain ( ); ++s ) {
              Rijkl += getChristoffelTensor ( j, s, l, BasePoint ) * getChristoffelTensor ( i, k, s, BasePoint ) - getChristoffelTensor ( k, s, l, BasePoint ) * getChristoffelTensor ( i, j, s, BasePoint );
            }
            Dest[l] += Rijkl * Z[i] * X[j] * Y[k];
          }
        }
      }
    }
  }

  // \nabla_X ( Y )
  template < typename DiscreteVectorFieldTypeX, typename DiscreteVectorFieldTypeY >
  void applyCovariantDerivative ( const VectorType& BasePoint, const DiscreteVectorFieldTypeX& XField, const DiscreteVectorFieldTypeY& YField, VectorType& Dest ) const {
    Dest.setZero ( );
    MatrixType jacobian ( getDimOfDomain ( ), getDimOfDomain ( ) );
    YField.getJacobianMatrix ( BasePoint, jacobian );
    aol::Vector < RealType > Y ( getDimOfDomain ( ) );
    YField.apply ( BasePoint, Y );
    aol::Vector < RealType > X ( getDimOfDomain ( ) );
    XField.apply ( BasePoint, X );
    for ( unsigned short int k = 0; k < getDimOfDomain ( ); ++k ) {
      for ( unsigned short int i = 0; i < getDimOfDomain ( ); ++i ) {
        Dest[k] += ( X[i] * jacobian.get ( k, i ) );
        for ( unsigned short int j = 0; j < getDimOfDomain ( ); ++j ) {
          Dest[k] += X[i] * Y[j] * getChristoffelTensor ( i, j, k, BasePoint );
        }
      }
    }
  }

  // \nabla_Outer \nabla_Inner ( Y )
  // it is assumed that Inner corresponds to a constant vector field!
  template < typename DiscreteVectorFieldTypeOuter, typename DiscreteVectorFieldTypeInner, typename DiscreteVectorFieldTypeY >
  void applySecondOrderCovariantDerivative ( const VectorType& BasePoint, const DiscreteVectorFieldTypeOuter& OuterField, const DiscreteVectorFieldTypeInner& InnerField, const DiscreteVectorFieldTypeY& YField, VectorType& Dest ) const {
    Dest.setZero ( );
    aol::RandomAccessContainer < MatrixType > hessian ( getDimOfDomain ( ), getDimOfDomain ( ), getDimOfDomain ( ) );
    YField.getHessianMatrix ( BasePoint, hessian );
    MatrixType jacobian ( getDimOfDomain ( ), getDimOfDomain ( ) );
    YField.getJacobianMatrix ( BasePoint, jacobian );
    aol::Vector < RealType > currentValueOfVectorField ( getDimOfDomain ( ) );
    YField.apply ( BasePoint, currentValueOfVectorField );
    aol::Vector < RealType > outer ( getDimOfDomain ( ) );
    OuterField.apply ( BasePoint, outer );
    aol::Vector < RealType > inner ( getDimOfDomain ( ) );
    InnerField.apply ( BasePoint, inner );
    MatrixType jacobianInnerField ( getDimOfDomain ( ), getDimOfDomain ( ) );
    InnerField.getJacobianMatrix ( BasePoint, jacobianInnerField );

    for ( unsigned short int n = 0; n < getDimOfDomain ( ); ++n ) {
      for ( unsigned short int l = 0; l < getDimOfDomain ( ); ++l ) {
        for ( unsigned short int i = 0; i < getDimOfDomain ( ); ++i ) {
          Dest[n] += outer[l] * inner[i] * hessian[l].get ( n, i );
          for ( unsigned short int m = 0; m < getDimOfDomain ( ); ++m ) {
            Dest[n] += outer[l] * inner[i] * ( jacobian.get ( m, i ) * getChristoffelTensor ( l, m, n, BasePoint ) + jacobian.get ( m, l ) * getChristoffelTensor ( i, m, n, BasePoint ) + currentValueOfVectorField[m] * getDChristoffelTensor ( i, m, n, l, BasePoint ) ) + outer[l] * jacobianInnerField.get ( i, l ) * currentValueOfVectorField[m] * getChristoffelTensor ( i, m, n, BasePoint );
            for ( unsigned short int j = 0; j < getDimOfDomain ( ); ++j ) {
              Dest[n] += outer[l] * inner[i] * currentValueOfVectorField[j] * getChristoffelTensor ( i, j, m, BasePoint ) * getChristoffelTensor ( l, m, n, BasePoint );
            }
          }
        }
      }
    }
  }

  template < typename ScalarProductType, typename DiscreteVectorFieldTypeX, typename DiscreteVectorFieldTypeY >
  RealType applySectionalCurvature ( const VectorType& BasePoint, const DiscreteVectorFieldTypeX& XField, const DiscreteVectorFieldTypeY& YField, const ScalarProductType& g ) const {
    aol::Vector < RealType > dest ( getDimOfDomain( ) );
    this->applyRiemannCurvatureTensor ( BasePoint, XField, YField, YField, dest );
    aol::Vector < RealType > X ( getDimOfDomain ( ) );
    XField.apply ( BasePoint, X );
    aol::Vector < RealType > Y ( getDimOfDomain ( ) );
    YField.apply ( BasePoint, Y );
    return g.apply ( BasePoint, dest, X ) / ( g.apply ( BasePoint, X, X ) * g.apply ( BasePoint, Y, Y ) - aol::Sqr ( g.apply ( BasePoint, X, Y ) ) );
  }


  RealType getChristoffelTensor ( const unsigned short i, const unsigned short j, const unsigned short k, const aol::Vector < RealType >& BasePoint ) const {
    if ( i >= getDimOfDomain ( ) || j >= getDimOfDomain ( ) || k >= getDimOfDomain ( ) )
      throw BasicException ( "Size mismatch!" );
    return this->asImp().ChristoffelTensor ( i, j, k, BasePoint );
  }

  RealType getDChristoffelTensor ( const unsigned short i, const unsigned short j, const unsigned short k, const unsigned short l, const aol::Vector < RealType >& BasePoint ) const {
    if ( i >= getDimOfDomain ( ) || j >= getDimOfDomain ( ) || k >= getDimOfDomain ( ) || l >= getDimOfDomain ( ) )
      throw BasicException ( "Size mismatch!" );
    return this->asImp().DChristoffelTensor ( i, j, k, l, BasePoint );
  }
*/

protected:
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//========================================================================================================
//! \brief Standard torus parametrization with inner and outer radius
//! \author Effland and Heeren
template < typename ConfiguratorType >
class TorusParametrization : public SurfaceParametrizationInterface < ConfiguratorType, 2, 3, TorusParametrization <ConfiguratorType> > {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;       
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
   
  const RealType _smallRadius;
  const RealType _greatRadius;

public:
  static const unsigned short int DimOfDomain = 2;
  static const unsigned short int DimOfRange = 3;

  TorusParametrization ( const RealType SmallRadius, const RealType GreatRadius )
  : SurfaceParametrizationInterface < ConfiguratorType, 2, 3, TorusParametrization <ConfiguratorType> > ( ), _smallRadius ( SmallRadius ), _greatRadius ( GreatRadius ) {
    if ( _smallRadius > _greatRadius )
      throw BasicException( "TorusParametrization::TorusParametrization() Wrong parameter choice!" );
  }

//protected:
  // Torus coordinates to Euclidean coordinates
  void computeParametrization ( const VectorType& Arg, VectorType& Dest ) const {
    if( Arg.size() != DimOfDomain )
        throw BasicException( "TorusParametrization::computeParametrization(): wrong size in arg!" );
    if( Dest.size() != DimOfRange )
        Dest.resize(DimOfRange);
    const RealType coeff = _greatRadius + _smallRadius * std::cos( Arg[1] );
    Dest[0] += coeff * std::cos ( Arg[0] );
    Dest[1] += coeff * std::sin ( Arg[0] );
    Dest[2] += _smallRadius * std::sin( Arg[1] );
  }
  
  RealType getGaussCurvature ( const VectorType& Arg ) const {
    if( Arg.size() != DimOfDomain )
        throw BasicException( "TorusParametrization::getGaussCurvature(): wrong size in arg!" );
    const RealType coeff = _greatRadius + _smallRadius * std::cos( Arg[1] );
    return std::cos(Arg[1]) / (_smallRadius * coeff);
  }

  void computeDE ( const VectorType& Arg, FullMatrixType& Derivative ) const {
    if( Arg.size() != DimOfDomain )
        throw BasicException( "TorusParametrization::computeDE(): wrong size in arg!" );
    if( Derivative.rows() != DimOfRange || Derivative.cols() != DimOfDomain )
      Derivative.resize( DimOfRange, DimOfDomain );
    
    const RealType coeff1 = _greatRadius + _smallRadius * std::cos ( Arg[1] );
    const RealType coeff2 = _smallRadius * ( -std::sin ( Arg[1] ) );
    Derivative( 0, 0 ) = coeff1 * ( -std::sin ( Arg[0] ) );
    Derivative( 1, 0 ) = coeff1 * ( std::cos ( Arg[0] ) );
    Derivative( 2, 0 ) = 0.;
    Derivative( 0, 1 ) = coeff2 * std::cos ( Arg[0] );
    Derivative( 1, 1 ) = coeff2 * std::sin ( Arg[0] );
    Derivative( 2, 1 ) = _smallRadius * std::cos ( Arg[1] );
  }

  void computeD2E ( const VectorType& Arg, std::vector<FullMatrixType>& Dest ) const {
    //std::cerr << "computeD2E" << std::endl;
    if( Arg.size() != DimOfDomain )
        throw BasicException( "TorusParametrization::computeD2E(): wrong size in arg!" );
    if( Dest.size() != DimOfRange )
      Dest.resize( DimOfRange );

    Dest[0].resize( DimOfDomain, DimOfDomain );
    Dest[0]( 0, 0 ) = ( _greatRadius + _smallRadius * std::cos ( Arg[1] ) ) * ( -std::cos ( Arg[0] ) ) ;
    Dest[0]( 0, 1 ) = _smallRadius * std::sin ( Arg[1] ) * std::sin ( Arg[0] );
    Dest[0]( 1, 0 ) = Dest[0]( 0, 1 );
    Dest[0]( 1 , 1 ) = ( -_smallRadius) * std::cos ( Arg[1])  * std::cos ( Arg[0] ) ;

    Dest[1].resize( DimOfDomain, DimOfDomain );
    Dest[1]( 0, 0 ) = ( _greatRadius + _smallRadius * std::cos ( Arg[1] ) ) * ( -std::sin ( Arg[0] ) );
    Dest[1]( 0, 1 ) = ( -_smallRadius ) * std::sin ( Arg[1] ) * std::cos ( Arg[0] );
    Dest[1]( 1, 0 ) = Dest[1]( 0, 1 );
    Dest[1]( 1, 1 ) = ( -_smallRadius ) * std::cos ( Arg[1] ) * std::sin ( Arg[0] ) ;

    Dest[2].resize( DimOfDomain, DimOfDomain );
    Dest[2].setZero ( );
    Dest[2]( 1, 1 ) = ( -_smallRadius ) * std::sin ( Arg[1] );
    //std::cerr << "computeD2E done" << std::endl;
  }

  void computeD3E ( const VectorType& Arg, FourthTensor < RealType, DimOfDomain, DimOfRange >& Dest ) const {
    Dest.setZero ( );
    const RealType sinu = std::sin ( Arg[0] );
    const RealType cosu = std::cos ( Arg[0] );
    const RealType sinv = std::sin ( Arg[1] );
    const RealType cosv = std::cos ( Arg[1] );
    const RealType R = _greatRadius;
    const RealType r = _smallRadius;

    Dest.setSymmetric ( 0, 0, ( R + r * cosv ) * sinu );
    Dest.setSymmetric ( 0, 1, r * cosu * sinv );
    Dest.setSymmetric ( 0, 2, r * cosv * sinu );
    Dest.setSymmetric ( 0, 3, r * cosu * sinv );

    Dest.setSymmetric ( 1, 0, - cosu * ( R + r * cosv ) );
    Dest.setSymmetric ( 1, 1, r * sinu * sinv );
    Dest.setSymmetric ( 1, 2, -r * cosu * cosv );
    Dest.setSymmetric ( 1, 3, r * sinu * sinv );

    Dest.setSymmetric ( 2, 0, 0. );
    Dest.setSymmetric ( 2, 1, 0. );
    Dest.setSymmetric ( 2, 2, 0. );
    Dest.setSymmetric ( 2, 3, -r * cosv );
  }

  // Euclidean coordinates to Torus coordinates
  void computeInverseParametrization ( const VectorType& Arg, VectorType& Dest ) const {
    if( Arg.size() != DimOfRange )
      throw BasicException( "TorusParametrization::computeInverseParametrization(): wrong size in arg!" );
      
    VectorType DestCopy( DimOfDomain );
    DestCopy[0] = std::atan2 ( Arg[1], Arg[0] );
    DestCopy[1] = std::acos ( ( Arg[0] / std::cos( DestCopy[0] ) - _greatRadius ) / _smallRadius );

    VectorType temp;
    apply ( DestCopy, temp );
    temp -= Arg;
    const RealType normUnchanged = temp.norm();

    VectorType DestCopyChanged ( DestCopy );
    DestCopyChanged[1] *= -1.0;
    apply ( DestCopyChanged, temp );
    temp -= Arg;
    const RealType normChanged = temp.norm();
    if ( normUnchanged < normChanged )
      Dest += DestCopy;
    else
      Dest += DestCopyChanged;
  }

/*
  RealType ChristoffelTensor ( const short i, const short j, const short k, const VectorType& BasePoint ) const {
    if ( i == 0 && j == 1 && k == 0 )
      return - ( _smallRadius * sin ( BasePoint[1] ) ) / ( _greatRadius + ( _smallRadius * cos ( BasePoint[1] ) ) );
    if ( i == 1 && j == 0 && k == 0 )
      return - ( _smallRadius * sin ( BasePoint[1] ) ) / ( _greatRadius + ( _smallRadius * cos ( BasePoint[1] ) ) );
    if ( i == 0 && j == 0 && k == 1 )
      return ( sin ( BasePoint[1] ) * ( _greatRadius + ( _smallRadius * cos ( BasePoint[1] ) ) ) ) / ( _smallRadius );
    return 0.;
  }

  RealType DChristoffelTensor ( const short i, const short j, const short k, const short l, const VectorType& BasePoint ) const {
    if ( i == 0 && j == 1 && k == 0 && l == 1 )
      return - ( Square( _smallRadius * sin ( BasePoint[1] ) ) ) / ( Square( _greatRadius + ( _smallRadius * cos ( BasePoint[1] ) ) ) )
          - ( _smallRadius * cos ( BasePoint[1] ) ) / ( _greatRadius + ( _smallRadius * cos ( BasePoint[1] ) ) );
    if ( i == 1 && j == 0 && k == 0 && l == 1)
      return - ( Square( _smallRadius * sin ( BasePoint[1] ) ) ) / ( Square( _greatRadius + ( _smallRadius * cos ( BasePoint[1] ) ) ) )
          - ( _smallRadius * cos ( BasePoint[1] ) ) / ( _greatRadius + ( _smallRadius * cos ( BasePoint[1] ) ) );
    if ( i == 0 && j == 0 && k == 1 && l == 1)
      return ( _greatRadius * cos ( BasePoint[1] ) + _smallRadius * Square( cos ( BasePoint[1] ) ) - _smallRadius * Square( sin ( BasePoint[1] ) ) ) / ( _smallRadius );
    return 0.;
  }
*/
  void plotSurface ( std::ofstream& Out ) const {
    Out << "r = " << _smallRadius << ";" << "\n";
    Out << "R = " << _greatRadius << ";" << "\n";
    Out << "x = ( R + r * cos ( v ) ) .* cos ( u );" << "\n";
    Out << "y = ( R + r * cos ( v ) ) .* sin ( u );" << "\n";
    Out << "z = r * sin ( v );" << "\n";
  }

protected:
    inline RealType Square( const RealType& val ) const {
        return val*val;
    }

};

//! \brief Standard sphere parametrization
//! \author Effland and Heeren
template < typename ConfiguratorType>
class SphereParametrization : public SurfaceParametrizationInterface < ConfiguratorType, 2, 3, SphereParametrization <ConfiguratorType> > {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;       
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  
  const RealType _radius;

public:
  static const unsigned short int DimOfDomain = 2;
  static const unsigned short int DimOfRange = 3;

  SphereParametrization ( const RealType Radius )
  : SurfaceParametrizationInterface <ConfiguratorType, 2, 3, SphereParametrization<ConfiguratorType> > ( ), _radius ( Radius ) { }

  void computeParametrization ( const VectorType& Arg, VectorType& Dest ) const {
    if( Arg.size() != DimOfDomain )
        throw BasicException( "SphereParametrization::computeParametrization(): wrong size in arg!" );
    if( Dest.size() != DimOfRange )
        Dest.resize(DimOfRange);
    Dest[0] = _radius * sin ( Arg[0] ) * cos ( Arg[1] );
    Dest[1] = _radius * sin ( Arg[0] ) * sin ( Arg[1] );
    Dest[2] = _radius * cos ( Arg[0] );
  }

  void computeInverseParametrization ( const VectorType& Arg, VectorType& Dest ) const {
    if( Arg.size() != DimOfRange )
        throw BasicException( "SphereParametrization::computeInverseParametrization(): wrong size in arg!" );
    if( Dest.size() != DimOfDomain )
        Dest.resize(DimOfDomain);
    Dest[1] = atan ( Arg[1] / Arg[0] );
    Dest[0] = acos ( Arg[2] / _radius );
  }
  
  RealType getGaussCurvature ( const VectorType& Arg ) const {
    if( Arg.size() != DimOfDomain )
        throw BasicException( "SphereParametrization::getGaussCurvature(): wrong size in arg!" );
    return 1./(_radius*_radius);
  }

  void plotSurface ( std::ofstream& Out ) const {
    Out << "R = " << _radius << ";" << "\n";
    Out << "x = R * sin ( u ) .* cos ( v );" << "\n";
    Out << "y = R * sin ( u ) .* sin ( v );" << "\n";
    Out << "z = R * cos ( u );" << "\n";
  }

  void computeDE ( const VectorType& Arg, FullMatrixType& Derivative ) const {
    const RealType sinu = sin ( Arg[0] );
    const RealType cosu = cos ( Arg[0] );
    const RealType sinv = sin ( Arg[1] );
    const RealType cosv = cos ( Arg[1] );
    Derivative.resize(DimOfRange, DimOfDomain);
    Derivative( 0, 0 ) = _radius * cosu * cosv;
    Derivative( 1, 0 ) = _radius * cosu * sinv;
    Derivative( 2, 0 ) = - _radius * sinu;
    Derivative( 0, 1 ) = - _radius * sinu * sinv;
    Derivative( 1, 1 ) = _radius * sinu * cosv;
    Derivative( 2, 1 ) = 0.;
  }

  void computeD2E ( const VectorType& Arg, std::vector <FullMatrixType>& Dest ) const {
    if( Arg.size() != DimOfDomain )
        throw BasicException( "SphereParametrization::computeD2E(): wrong size in arg!" );
    if( Dest.size() != DimOfRange )
        Dest.resize(DimOfRange);
    
    const RealType sinu = sin ( Arg[0] );
    const RealType cosu = cos ( Arg[0] );
    const RealType sinv = sin ( Arg[1] );
    const RealType cosv = cos ( Arg[1] );
    for ( unsigned short int i = 0; i < DimOfRange; ++i )
      Dest[i].resize(DimOfDomain, DimOfDomain );
 
    Dest[0]( 0, 0 ) = - _radius * sinu * cosv;
    Dest[0]( 0, 1 ) = - _radius * cosu * sinv;
    Dest[0]( 1, 0 ) =  Dest[0]( 0, 1 );
    Dest[0]( 1, 1 ) = - _radius * sinu * cosv;

    Dest[1]( 0, 0 ) = - _radius * sinu * sinv;
    Dest[1]( 0, 1 ) = _radius * cosu * cosv;
    Dest[1]( 1, 0 ) = Dest[1]( 0, 1 );
    Dest[1]( 1, 1 ) = - _radius * sinu * sinv;

    Dest[2].setZero( );
    Dest[2]( 0, 0 ) = - _radius * cosu;
  }

  void getD3E ( const VectorType& Arg, FourthTensor< RealType, DimOfDomain, DimOfRange >& Dest ) const {
    Dest.setZero ( );
    const RealType sinu = sin ( Arg[0] );
    const RealType cosu = cos ( Arg[0] );
    const RealType sinv = sin ( Arg[1] );
    const RealType cosv = cos ( Arg[1] );
    const RealType r = _radius;

    Dest.setSymmetric ( 0, 0, - r * cosu * cosv );
    Dest.setSymmetric ( 0, 1, r * sinu * sinv );
    Dest.setSymmetric ( 0, 2, - r * cosu * cosv );
    Dest.setSymmetric ( 0, 3, r * sinu * sinv );

    Dest.setSymmetric ( 1, 0, - r * cosu * sinv );
    Dest.setSymmetric ( 1, 1, - r * sinu * cosv );
    Dest.setSymmetric ( 1, 2, - r * cosu * sinv );
    Dest.setSymmetric ( 1, 3, - r * sinu * cosv );

    Dest.setSymmetric ( 2, 0, r * sinu );
    Dest.setSymmetric ( 2, 1, 0. );
    Dest.setSymmetric ( 2, 2, 0. );
    Dest.setSymmetric ( 2, 3, 0. );
  }
/*
  RealType ChristoffelTensor ( const short i, const short j, const short k, const aol::Vector < RealType >& BasePoint ) const {
    if ( i == 1 && j == 1 && k == 0 )
      return - sin ( BasePoint[0] ) * cos ( BasePoint[0] );
    if ( i == 1 && j == 0 && k == 1 )
      return aol::ZOTrait<RealType>::one / tan ( BasePoint[0] ) ;
    if ( i == 0 && j == 1 && k == 1 )
      return aol::ZOTrait<RealType>::one / tan ( BasePoint[0] );
    return aol::ZTrait<RealType>::zero;
  }

  RealType DChristoffelTensor ( const short i, const short j, const short k, const short l, const aol::Vector < RealType >& BasePoint ) const {
    if ( i == 1 && j == 1 && k == 0 && l == 0 )
      return aol::Sqr ( sin ( BasePoint[0] ) ) - aol::Sqr ( cos ( BasePoint[0] ) );
    if ( i == 1 && j == 0 && k == 1 && l == 0)
      return - aol::ZOTrait<RealType>::one / aol::Sqr ( sin ( BasePoint[0] ) );
    if ( i == 0 && j == 1 && k == 1 && l == 0 )
      return - aol::ZOTrait<RealType>::one / aol::Sqr ( sin ( BasePoint[0] ) );
    return aol::ZTrait<RealType>::zero;
  }
*/
};


//========================================================================================================
//========================================================================================================
//========================================================================================================
//! \brief Base class for generic energy \f$ W \f$
//! \author Effland and Heeren
template < typename SurfaceParametrization >
class PairwiseEnergyBase {
protected:
  const SurfaceParametrization& _surface;
  typedef SurfaceParametrization SurfaceType;

public:
  PairwiseEnergyBase ( const SurfaceParametrization& Surface )
: _surface ( Surface) { }

  const SurfaceType& getSurfaceParam ( ) const {
    return _surface;
  }
};


//========================================================================================================
//! \brief Simple energy \f$ W(x,y)= 0.5 * \|E(x) - E(y)\|^2 \f$, where E is surface parametrization
//! \author Effland and Heeren
template < typename ConfiguratorType, typename SurfaceParametrization >
class PairwiseEnergy : public PairwiseEnergyBase < SurfaceParametrization > {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType;
  
public:
  PairwiseEnergy ( const SurfaceParametrization& Surface ) : PairwiseEnergyBase < SurfaceParametrization > ( Surface ) { }

  typedef SurfaceParametrization SurfaceType;

  void applyAdd ( const VectorType& Arg1, const VectorType& Arg2, RealType & Dest ) const {
    if( Arg1.size ( ) != Arg2.size ( ) )
      throw BasicException( "PairwiseEnergy::applyAdd(): Size mismatch!" );
    VectorType arg1Y, arg2Y;
    this->getSurfaceParam( ).apply ( Arg1, arg1Y );
    this->getSurfaceParam( ).apply ( Arg2, arg2Y );
    arg1Y -= arg2Y;
    Dest += 0.5 * arg1Y.squaredNorm();
  }

  void apply ( const VectorType& Arg1, const VectorType& Arg2, RealType & Dest ) const {
    Dest = 0.;
    applyAdd ( Arg1, Arg2, Dest );
  }

  // D_1 W(x,y) = DE(x)^T * (E(x) - E(y))
  void applyAddDFirstArg ( const VectorType& Arg1, const VectorType& Arg2, VectorType& Dest ) const {
    if ( Arg1.size ( ) != Arg2.size ( ) )
      throw BasicException( "PairwiseEnergy::applyAddDFirstArg(): Size mismatch!" );
    VectorType arg1Y, arg2Y;
    this->getSurfaceParam( ).apply ( Arg1, arg1Y );
    this->getSurfaceParam( ).apply ( Arg2, arg2Y );
    arg1Y -= arg2Y;
    this->getSurfaceParam( ).applyAddDE ( Arg1, arg1Y, Dest );
  }

  void applyDFirstArg ( const VectorType& Arg1, const VectorType& Arg2, VectorType& Dest ) const {
    Dest.setZero ( );
    applyAddDFirstArg ( Arg1, Arg2, Dest );
  }

  // D_2 W(x,y) = DE(y)^T * (E(y) - E(x))
  void applyAddDSecondArg ( const VectorType& Arg1, const VectorType& Arg2, VectorType& Dest ) const {
    if ( Arg1.size ( ) != Arg2.size ( ) )
      throw BasicException( "PairwiseEnergy::applyAddDSecondArg(): Size mismatch!" );
    VectorType arg1Y, arg2Y;
    this->getSurfaceParam( ).apply ( Arg1, arg1Y );
    this->getSurfaceParam( ).apply ( Arg2, arg2Y );
    arg2Y -= arg1Y;
    this->getSurfaceParam( ).applyAddDE ( Arg2, arg2Y, Dest );
  }

  void applyDSecondArg ( const VectorType& Arg1, const VectorType& Arg2, VectorType& Dest ) const {
    Dest.setZero ( );
    applyAddDSecondArg ( Arg1, Arg2, Dest );
  }

  // D_{11} W(x,y) = DE(x)^T * DE(x) + D^2E(x) * (E(x) - E(y))
  void D2XX ( const VectorType& Arg1, const VectorType& Arg2, MatrixType& Dest ) const {
//std::cerr << "D2XX" << std::endl;
    boundCheckHessian ( Arg1, Arg2, Dest );
    Dest.setZero ();

    VectorType arg1Y, arg2Y;
    this->getSurfaceParam( ).apply ( Arg1, arg1Y );
    this->getSurfaceParam( ).apply ( Arg2, arg2Y );
    arg1Y -= arg2Y;
//std::cerr << "applyD2E" << std::endl;
    MatrixType temp;
    this->getSurfaceParam( ).applyD2E ( Arg1, arg1Y, temp );
//std::cerr << "applyDE" << std::endl;
    this->getSurfaceParam( ).applyDE ( Arg1, Arg1, Dest );
//std::cerr << "done" << std::endl;
    Dest += temp;
//std::cerr << "leave" << std::endl;
  }

  // D_{12} W(x,y) = DE(x)^T * DE(y) 
  void D2XY ( const VectorType& Arg1, const VectorType& Arg2, MatrixType& Dest ) const {
    boundCheckHessian ( Arg1, Arg2, Dest );
    Dest.setZero ();
    this->getSurfaceParam( ).applyDE ( Arg1, Arg2, Dest );
    Dest *= -1.;
  }
  
  // D_{22} W(x,y) = DE(y)^T * DE(y) + D^2E(y) * (E(y) - E(x))
  void D2YY ( const VectorType& Arg1, const VectorType& Arg2, MatrixType& Dest ) const {
    boundCheckHessian ( Arg1, Arg2, Dest );
    Dest.setZero ();

    VectorType arg1Y, arg2Y;
    this->getSurfaceParam( ).apply ( Arg1, arg1Y );
    this->getSurfaceParam( ).apply ( Arg2, arg2Y );
    arg2Y -= arg1Y;

    MatrixType temp;
    this->getSurfaceParam( ).applyD2E ( Arg2, arg2Y, temp );
    this->getSurfaceParam( ).applyDE ( Arg2, Arg2, Dest );
    Dest += temp;
  }
  
  unsigned short int getDimOfDomain() const {
    return this->getSurfaceParam( ).getDimOfDomain();
  }
  
private:
  void boundCheckHessian ( const VectorType& Arg1, const VectorType& Arg2, MatrixType& Dest ) const {
    if ( Arg1.size ( ) != getDimOfDomain() )
        throw BasicException( "PairwiseEnergy::boundCheckHessian (): First argument has wrong size!");
    if ( Arg1.size ( ) != Arg2.size ( ) )
      throw BasicException( "PairwiseEnergy::boundCheckHessian (): Second argument has wrong size!");
    if( Dest.rows ( ) != this->getSurfaceParam( ).getDimOfDomain ( ) || Dest.cols ( ) != this->getSurfaceParam( ).getDimOfDomain ( ) )
        Dest.resize( this->getSurfaceParam( ).getDimOfDomain ( ), this->getSurfaceParam( ).getDimOfDomain ( ) );
  }

};


//========================================================================================================
//========================================================================================================
//========================================================================================================
//! \brief Discrete path energy \f$ E[x_0, ..., x_K] = K \sum_{k=1}^K W[x_{k-1},x_k] \f$
//! \author Effland and Heeren
//! Here x_0 and x_K are given as fixed point traits in the constructor.
template < typename ConfiguratorType, typename PairwiseEnergyType >
class PathEnergy : public BaseOp < typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType;    
    
  const PairwiseEnergyType& _pairwiseEnergy;
  const FixedPointsTrait<ConfiguratorType>& _fixedPointsTrait;
  int _K, _dimOfDomain;

public:
  PathEnergy( const PairwiseEnergyType& PairwiseEnergy, const FixedPointsTrait<ConfiguratorType>& FixedPointsTrait, int K )
: _pairwiseEnergy ( PairwiseEnergy ), _fixedPointsTrait ( FixedPointsTrait ), _K(K), _dimOfDomain( PairwiseEnergy.getSurfaceParam( ).getDimOfDomain() ) { }

  void apply( const VectorType& Arg, RealType & Dest ) const {
    if( Arg.size() != (_K-1) * _dimOfDomain )
      throw BasicException("PathEnergy::apply(): arg has wrong size!");
      
    // bring into more convenient form
    std::vector< Eigen::Ref<const VectorType> > Points;
    Points.reserve(_K-1);
    for( int k = 0; k < _K-1; k++ )
       Points.push_back( Arg.segment(k*_dimOfDomain, _dimOfDomain) );

    Dest = 0.;
    _pairwiseEnergy.applyAdd ( _fixedPointsTrait.getStartVector ( ), Points[0], Dest );
    for ( int i = 0; i < _K-2; ++i )
      _pairwiseEnergy.applyAdd ( Points[i], Points[i+1], Dest );
    _pairwiseEnergy.applyAdd ( Points[_K-2], _fixedPointsTrait.getEndVector ( ), Dest );
    Dest *= _K;
  }
};

//! \brief Gradient of Discrete path energy
//! \author Effland and Heeren
template < typename ConfiguratorType, typename PairwiseEnergyType >
class DPathEnergy : public BaseOp < typename ConfiguratorType::VectorType > {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType;  
  
  const PairwiseEnergyType& _pairwiseEnergy;
  const FixedPointsTrait <ConfiguratorType>& _fixedPointsTrait;
  int _K, _dimOfDomain;
  
public:
  DPathEnergy( const PairwiseEnergyType& PairwiseEnergy, const FixedPointsTrait<ConfiguratorType>& FixedPointsTrait, int K  )
: _pairwiseEnergy ( PairwiseEnergy ), _fixedPointsTrait ( FixedPointsTrait ), _K(K), _dimOfDomain( PairwiseEnergy.getSurfaceParam( ).getDimOfDomain() ) { }

  void apply ( const VectorType& Arg, VectorType& Dest ) const {
    if( Arg.size() != (_K-1) * _dimOfDomain )
      throw BasicException("DPathEnergy::apply(): arg has wrong size!");
    if( Arg.size() != Dest.size() )
        Dest.resize( Arg.size() );
    Dest.setZero();
    VectorType temp;
    
    // bring into more convenient form
    std::vector< Eigen::Ref<const VectorType> > Points;
    Points.reserve(_K-1);
    for( int k = 0; k < _K-1; k++ )
       Points.push_back( Arg.segment(k*_dimOfDomain, _dimOfDomain) );
    
    // compute path energy gradient
    _pairwiseEnergy.applyDSecondArg ( _fixedPointsTrait.getStartVector (), Points[0], temp );
    for( int i = 0; i < _dimOfDomain; i++ )
        Dest[i] += _K * temp[i];
    _pairwiseEnergy.applyDFirstArg( Points[_K-2], _fixedPointsTrait.getEndVector ( ), temp );
    for( int i = 0; i < _dimOfDomain; i++ )
        Dest[(_K-2)*_dimOfDomain+i] += _K * temp[i];

    // inner points
    for( int k = 0; k < _K-1; k++ ){     
          temp.setZero();
          if( k < _K-2 )
            _pairwiseEnergy.applyAddDFirstArg( Points[k], Points[k+1], temp );
          if( k > 0 )
            _pairwiseEnergy.applyAddDSecondArg ( Points[k-1], Points[k], temp );
          for( int i = 0; i < _dimOfDomain; i++ )
            Dest[k*_dimOfDomain+i] += _K * temp[i];
    }

  }
};

//! \brief Hessian of Discrete path energy
//! \author Effland and Heeren
template< typename ConfiguratorType, typename PairwiseEnergyType >
class D2PathEnergy : public BaseOp < typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType > {
    
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType; 
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType; 
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  const PairwiseEnergyType& _pairwiseEnergy;
  const FixedPointsTrait<ConfiguratorType>& _fixedPointsTrait;
  int _K, _dimOfDomain;
  
public:
  D2PathEnergy( const PairwiseEnergyType& PairwiseEnergy, const FixedPointsTrait<ConfiguratorType>& FixedPointsTrait, int K )
: _pairwiseEnergy ( PairwiseEnergy ), _fixedPointsTrait ( FixedPointsTrait ), _K(K), _dimOfDomain( PairwiseEnergy.getSurfaceParam( ).getDimOfDomain() ) { }

  void apply( const VectorType& Arg, SparseMatrixType& Dest ) const {
    if( Arg.size() != (_K-1) * _dimOfDomain )
      throw BasicException("D2PathEnergy::apply(): arg has wrong size!");
    if( Dest.rows() != Arg.size() || Dest.cols() != Arg.size() )
        Dest.resize( Arg.size(), Arg.size() );
    Dest.setZero();
    
    // bring into more convenient form
    std::vector< Eigen::Ref<const VectorType> > Points;
    Points.reserve(_K-1);
    for( int k = 0; k < _K-1; k++ )
       Points.push_back( Arg.segment(k*_dimOfDomain, _dimOfDomain) );
    
    MatrixType temp1, temp2;    
    TripletListType tripletList;      
    tripletList.reserve( 4 * _K * _dimOfDomain * _dimOfDomain );

    // 
    for ( int k = 0; k < _K-1; ++k ) {
        
      // diagonal blocks  
      if ( k != _K-2 )
        _pairwiseEnergy.D2XX ( Points[k], Points[k+1], temp1 );
      else
        _pairwiseEnergy.D2XX ( Points[_K-2], _fixedPointsTrait.getEndVector ( ), temp1 );
      
      if ( k != 0 )
        _pairwiseEnergy.D2YY ( Points[k-1], Points[k], temp2 );
      else
        _pairwiseEnergy.D2YY ( _fixedPointsTrait.getStartVector ( ), Points[0], temp2 );
      
      for( int i = 0; i < _dimOfDomain; i++ )
          for( int j = 0; j < _dimOfDomain; j++ )
              tripletList.push_back( TripletType( k*_dimOfDomain+i, k*_dimOfDomain+j, _K * ( temp1(i,j) + temp2(i,j) ) ) ); 
          
      // off diagonal blocks
      if ( k != 0 ) {
        _pairwiseEnergy.D2XY ( Points[k], Points[k-1], temp1 );
        for( int i = 0; i < _dimOfDomain; i++ )
          for( int j = 0; j < _dimOfDomain; j++ )
              tripletList.push_back( TripletType( k*_dimOfDomain+i, (k-1)*_dimOfDomain+j, _K * temp1(i,j) ) );
      }

      if ( k != _K-2 ) {
        _pairwiseEnergy.D2XY ( Points[k], Points[k+1], temp1 );
        for( int i = 0; i < _dimOfDomain; i++ )
          for( int j = 0; j < _dimOfDomain; j++ )
              tripletList.push_back( TripletType( k*_dimOfDomain+i, (k+1)*_dimOfDomain+j, _K * temp1(i,j) ) );
      }
    }
    
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
};



//========================================================================================================
//========================================================================================================
//========================================================================================================
//! \brief Class to compute discrete K-geodesics by minimization of discrete path energy
//! \author Effland and Heeren
//! The integer K is set in the constructor whereas the fixed end points are passed in the method execute().
template < typename ConfiguratorType, typename PairwiseEnergyType >
class DiscreteGeodesic {  
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType; 
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;

  const PairwiseEnergyType& _pairwiseEnergy;
  int _K;
  const bool _quiet;
  const int _maxIterations;
  const RealType _stopEpsilon;

public:
  explicit DiscreteGeodesic ( const PairwiseEnergyType& PairwiseEnergy, int K, const int MaxIterations, const RealType StopEpsilon, bool Quiet = true )
  : _pairwiseEnergy ( PairwiseEnergy ), _K(K), _quiet( Quiet ), _maxIterations ( MaxIterations ), _stopEpsilon ( StopEpsilon ) { }

  void execute( const VectorType& Arg, VectorType& Dest, const FixedPointsTrait<ConfiguratorType>& FixedPoints ) const {

    PathEnergy < ConfiguratorType, PairwiseEnergyType > E ( _pairwiseEnergy, FixedPoints, _K );
    DPathEnergy < ConfiguratorType, PairwiseEnergyType > DE ( _pairwiseEnergy, FixedPoints, _K );
    D2PathEnergy < ConfiguratorType, PairwiseEnergyType > D2E ( _pairwiseEnergy, FixedPoints, _K );

    RealType energy;
    VectorType gradient;
    SparseMatrixType Hessian;
    
    QUIET_MODE solverQuietMode = _quiet ? SUPERQUIET : SHOW_ALL;
    
    // initial
    if(!_quiet){
      E.apply( Arg, energy );
      std::cerr << "Initial energy = " << energy << std::endl;
      DE.apply(Arg, gradient );
      std::cerr << "Initial gradient norm = " << gradient.norm() << std::endl;
    }
    
    // grad tests 
    //std::cerr << "Start first derivative tests" << std::endl;  
    //ScalarValuedDerivativeTester<ConfiguratorType> ( E, DE, 1e-7 ).plotRandomDirections ( Arg, 25, "/home/staff/heeren/reducedBasis/results/debug/gradTest" );    
    //std::cerr << "Start second derivative tests" << std::endl;
    //VectorValuedDerivativeTester<ConfiguratorType> ( DE, D2E, 1e-5 ).plotRandomDirections ( Arg, 25, "/home/staff/heeren/reducedBasis/results/debug/hessTest" );
    
    if(!_quiet) std::cerr << "Start gradient descent..."  << std::endl;
    GradientDescent<ConfiguratorType> gradDesc( E, DE, _maxIterations, _stopEpsilon, ARMIJO, solverQuietMode );
    gradDesc.solve( Arg, Dest );
    
    if(!_quiet) std::cerr << "Start Newton's method..."  << std::endl;
    VectorType Init = Dest;
    NewtonMethod<ConfiguratorType> newtonSolver( DE, D2E, _maxIterations, _stopEpsilon, NEWTON_OPTIMAL, solverQuietMode );
    newtonSolver.solve( Init, Dest );
    
    // final
    if(!_quiet){
      E.apply( Dest, energy );
      std::cerr << "Final energy = " << energy << std::endl;
      DE.apply(Dest, gradient );
      std::cerr << "Final gradient norm = " << gradient.norm() << std::endl;
    }
  }

};

//========================================================================================================
//========================================================================================================
//========================================================================================================
//! \brief Class to compute discrete logarithm
//! \author Effland and Heeren
template < typename ConfiguratorType, typename PairwiseEnergyType >
class Logarithm {  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  const PairwiseEnergyType& _pairwiseEnergy;
  int _K, _dimOfDomain;
  bool _quiet;
  int _maxIterations;
  RealType _stopEpsilon;

public:
  explicit Logarithm ( const PairwiseEnergyType& PairwiseEnergy, int K, const int MaxIterations, const RealType StopEpsilon, const bool Quiet = true )
  : _pairwiseEnergy ( PairwiseEnergy ), _K(K), _dimOfDomain( PairwiseEnergy.getSurfaceParam( ).getDimOfDomain() ), _quiet( Quiet ), _maxIterations ( MaxIterations ), _stopEpsilon ( StopEpsilon ) { }

  // FixedPoints = fixed points of geodesic, Dest = displacement
  void execute( const FixedPointsTrait <ConfiguratorType>& FixedPoints, VectorType& Dest ) const {
    if ( Dest.size ( ) != _dimOfDomain )
      Dest.resize( _dimOfDomain );

    VectorType linearInterpol, geodesic;
    linspace<VectorType>( FixedPoints.getStartVector ( ), FixedPoints.getEndVector ( ), linearInterpol, _K );
    DiscreteGeodesic < ConfiguratorType, PairwiseEnergyType > ( _pairwiseEnergy, _K, _maxIterations, _stopEpsilon, _quiet ).execute( linearInterpol, geodesic, FixedPoints );
    Dest = geodesic.segment( 0, _dimOfDomain );
    Dest -= FixedPoints.getStartVector ( );
  }
};

//========================================================================================================
//========================================================================================================
//========================================================================================================
//! \brief Vector-valued exponential function \f$ F[x] = W_{,2}[x_0, x_1] + W_{,1}[x_1, x] \f$, with fixed points \f$(x_0, x_1)\f$
//! \author Effland, Heeren
//! The discrete exponential \f$x_2\f$ is characterized as root of F, i.e. \f$F[x_2] = 0 iff. (x_0, x_1, x_2)\f$ is discrete geodesic.
template < typename ConfiguratorType, typename PairwiseEnergyType >
class ExponentialEnergy : public BaseOp<typename ConfiguratorType::VectorType> {
    
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const PairwiseEnergyType& _pairwiseEnergy;
  const FixedPointsTrait<ConfiguratorType>& _fixedPointsTrait;

public:
  ExponentialEnergy ( const PairwiseEnergyType& PairwiseEnergy, const FixedPointsTrait <ConfiguratorType>& FixedPointsTrait )
: _pairwiseEnergy ( PairwiseEnergy ), _fixedPointsTrait ( FixedPointsTrait ) { }

  void apply( const VectorType& Arg, VectorType& Dest ) const {
    if ( Arg.size() !=  _fixedPointsTrait.getStartVector().size() )
      throw BasicException( "ExponentialEnergy::apply(): arg has wrong size!");
    _pairwiseEnergy.applyDSecondArg ( _fixedPointsTrait.getStartVector( ), _fixedPointsTrait.getEndVector( ), Dest );
    _pairwiseEnergy.applyAddDFirstArg ( _fixedPointsTrait.getEndVector( ), Arg, Dest );
  }
};

//! \brief Matrix-valued derivative of exponential function
//! \author Effland, Heeren
template < typename ConfiguratorType, typename PairwiseEnergyType >
class DExponentialEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>  {
    
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType; 
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  const PairwiseEnergyType& _pairwiseEnergy;
  const FixedPointsTrait <ConfiguratorType>& _fixedPointsTrait;

public:
  DExponentialEnergy ( const PairwiseEnergyType& PairwiseEnergy, const FixedPointsTrait <ConfiguratorType>& FixedPointsTrait )
: _pairwiseEnergy ( PairwiseEnergy ), _fixedPointsTrait ( FixedPointsTrait ) { }                                                                                       //constructor

  void apply( const VectorType& Arg, SparseMatrixType& Dest ) const {
    int dimDomain = _fixedPointsTrait.getStartVector().size();
    if ( Arg.size() !=  dimDomain )
      throw BasicException( "DExponentialEnergy::apply(): arg has wrong size!");
    Dest.resize(dimDomain, dimDomain);
    
    MatrixType temp;
    _pairwiseEnergy.D2XY (  _fixedPointsTrait.getEndVector( ), Arg, temp );    
    
    TripletListType tripletList;      
    tripletList.reserve( dimDomain * dimDomain );
    for( int i = 0; i < dimDomain; i++ )
        for( int j = 0; j < dimDomain; j++ )
            tripletList.push_back( TripletType( i, j, temp(i,j) ) );     
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
};

//! \brief Vector-valued backward exponential function F[x] = W_{,2}[x, x_1] + W_{,1}[x_1, x_2], with fixed points (x_2, x_1)
//! \author Heeren
//! The backward discrete exponential x_0 is characterized as root of F, i.e. F[x_0] = 0 iff. (x_0, x_1, x_2) is discrete geodesic.
template < typename ConfiguratorType, typename PairwiseEnergyType >
class ExponentialEnergyBackwards : public BaseOp<typename ConfiguratorType::VectorType> {
    
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const PairwiseEnergyType& _pairwiseEnergy;
  const FixedPointsTrait<ConfiguratorType>& _fixedPointsTrait;

public:
  ExponentialEnergyBackwards ( const PairwiseEnergyType& PairwiseEnergy, const FixedPointsTrait <ConfiguratorType>& FixedPointsTrait )
: _pairwiseEnergy ( PairwiseEnergy ), _fixedPointsTrait ( FixedPointsTrait ) { }

  void apply( const VectorType& Arg, VectorType& Dest ) const {
    if ( Arg.size() !=  _fixedPointsTrait.getStartVector().size() )
      throw BasicException( "ExponentialEnergyBackwards::apply(): arg has wrong size!");
    _pairwiseEnergy.applyDSecondArg ( Arg, _fixedPointsTrait.getEndVector( ), Dest );
    _pairwiseEnergy.applyAddDFirstArg ( _fixedPointsTrait.getEndVector( ), _fixedPointsTrait.getStartVector( ), Dest );
  }
};

//! \brief Matrix-valued derivative of backward  exponential function
//! \author Heeren
template < typename ConfiguratorType, typename PairwiseEnergyType >
class DExponentialEnergyBackwards : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>  {
    
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType; 
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  
  const PairwiseEnergyType& _pairwiseEnergy;
  const FixedPointsTrait <ConfiguratorType>& _fixedPointsTrait;

public:
  DExponentialEnergyBackwards ( const PairwiseEnergyType& PairwiseEnergy, const FixedPointsTrait <ConfiguratorType>& FixedPointsTrait )
: _pairwiseEnergy ( PairwiseEnergy ), _fixedPointsTrait ( FixedPointsTrait ) { }                                                                                       //constructor

  void apply( const VectorType& Arg, SparseMatrixType& Dest ) const {
    int dimDomain = _fixedPointsTrait.getStartVector().size();
    if ( Arg.size() !=  dimDomain )
      throw BasicException( "DExponentialEnergyBackwards::apply(): arg has wrong size!");
    Dest.resize(dimDomain, dimDomain);
    
    MatrixType temp;
    _pairwiseEnergy.D2XY ( _fixedPointsTrait.getEndVector( ), Arg, temp );    
    
    TripletListType tripletList;      
    tripletList.reserve( dimDomain * dimDomain );
    for( int i = 0; i < dimDomain; i++ )
        for( int j = 0; j < dimDomain; j++ )
            tripletList.push_back( TripletType( i, j, temp(i,j) ) );     
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
};

//! \brief Computes root of (backward) exponential function
//! \author Effland and Heeren
template < typename ConfiguratorType, typename PairwiseEnergyType >
class Exp2 {
   
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const PairwiseEnergyType& _pairwiseEnergy;
  const bool _quiet;
  const int _maxIterations;
  const RealType _stopEpsilon;

public:
  explicit Exp2 ( const PairwiseEnergyType& PairwiseEnergy, const int MaxIterations, const RealType StopEpsilon, const bool Quiet = true )
  : _pairwiseEnergy ( PairwiseEnergy ), _quiet( Quiet ), _maxIterations ( MaxIterations ), _stopEpsilon ( StopEpsilon ) { }

  // FixedPoints are the first and second point on the geodesic curve, i.e. x_0 and x_1, compute x_2 such that (x_0, x_1, x_2) is geodesic
  void execute( const VectorType& Initialization, const FixedPointsTrait<ConfiguratorType>& FixedPoints, VectorType& Dest ) const {
    if( Initialization.size() != FixedPoints.getStartVector().size() )
        throw BasicException("Exp2::execute() argument has wrong size!");
      
    ExponentialEnergy < ConfiguratorType, PairwiseEnergyType > E ( _pairwiseEnergy, FixedPoints );
    DExponentialEnergy < ConfiguratorType, PairwiseEnergyType > DE ( _pairwiseEnergy, FixedPoints );
    
    //std::cerr << "Start derivative tests" << std::endl;
    //VectorValuedDerivativeTester<ConfiguratorType> ( E, DE, 1e-5 ).plotAllDirections ( Initialization, "/home/staff/heeren/reducedBasis/results/debug/hessTest" );
    
    if(!_quiet) std::cerr << "Start Newton's method..."  << std::endl;
    QUIET_MODE solverQuietMode = _quiet ? SUPERQUIET : SHOW_ALL;
    NewtonMethod<ConfiguratorType> newtonSolver( E, DE, _maxIterations, _stopEpsilon, NEWTON_OPTIMAL, solverQuietMode );
    newtonSolver.solve( Initialization, Dest );
  }
  
  // FixedPoints are the first and second point on the geodesic curve, i.e. x_2 and x_1, compute x_0 such that (x_0, x_1, x_2) is geodesic
  void executeBackwards( const VectorType& Initialization, const FixedPointsTrait<ConfiguratorType>& FixedPoints, VectorType& Dest ) const {
    if( Initialization.size() != FixedPoints.getStartVector().size() )
        throw BasicException("Exp2::executeBackwards() argument has wrong size!");
      
    ExponentialEnergyBackwards < ConfiguratorType, PairwiseEnergyType > E ( _pairwiseEnergy, FixedPoints );
    DExponentialEnergyBackwards < ConfiguratorType, PairwiseEnergyType > DE ( _pairwiseEnergy, FixedPoints );
    
    //std::cerr << "Start derivative tests" << std::endl;
    //VectorValuedDerivativeTester<ConfiguratorType> ( E, DE, 1e-5 ).plotAllDirections ( Initialization, "/home/staff/heeren/reducedBasis/results/debug/hessTest" );
    
    if(!_quiet) std::cerr << "Start Newton's method..."  << std::endl;
    QUIET_MODE solverQuietMode = _quiet ? SUPERQUIET : SHOW_ALL;
    NewtonMethod<ConfiguratorType> newtonSolver( E, DE, _maxIterations, _stopEpsilon, NEWTON_OPTIMAL, solverQuietMode );
    newtonSolver.solve( Initialization, Dest );
  }
};

//! \brief Iterative application of discrete Exp2
//! \author Effland and Heeren
template < typename ConfiguratorType, typename PairwiseEnergyType >
class ExpK {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  typedef Exp2 < ConfiguratorType, PairwiseEnergyType > Exp2Type;
  const Exp2Type _exp2;

public:
  explicit ExpK ( const PairwiseEnergyType& PairwiseEnergy, const int MaxIterations, const RealType StopEpsilon, const bool Quiet = true )
  : _exp2 ( PairwiseEnergy, MaxIterations, StopEpsilon, Quiet ) { }

  // Arg[0] and Arg[1] is input data
  void execute( const FixedPointsTrait <ConfiguratorType>& Arg, int steps, VectorType& Dest ) const {
    VectorType oldPoint( Arg.getStartVector() ), newPoint( Arg.getEndVector() );
    for ( unsigned int i = 0; i < steps; ++i ) {
      const FixedPointsTrait <ConfiguratorType> currentFixedPoints ( oldPoint, newPoint );
      VectorType init ( newPoint );
      _exp2.execute ( init, currentFixedPoints, Dest );
      oldPoint = newPoint;
      newPoint = Dest;
    }
  }
  
  //
  void executeBackwards( const FixedPointsTrait <ConfiguratorType>& Arg, int steps, VectorType& Dest ) const {
    VectorType oldPoint( Arg.getStartVector() ), newPoint( Arg.getEndVector() );
    for ( unsigned int i = 0; i < steps; ++i ) {
      const FixedPointsTrait <ConfiguratorType> currentFixedPoints ( oldPoint, newPoint );
      VectorType init ( newPoint );
      _exp2.executeBackwards( init, currentFixedPoints, Dest );
      oldPoint = newPoint;
      newPoint = Dest;
    }
  }

};

//========================================================================================================
//========================================================================================================
//========================================================================================================

//! \brief Discrete geodesic interpolation
//! \author Heeren
template<typename ConfiguratorType, typename PairwiseEnergyType>
void integerInterpolation( const typename ConfiguratorType::VectorType& StartPoint, 
                           const typename ConfiguratorType::VectorType& EndPoint,
                           const PairwiseEnergyType& W,
                           int totalLength,
                           typename ConfiguratorType::VectorType& geodPath,
                           bool quiet = true ){
   DiscreteGeodesic<ConfiguratorType, PairwiseEnergyType> geodesic( W, totalLength-1, 1000, 1e-8 );
   FixedPointsTrait<ConfiguratorType> endPoints( StartPoint, EndPoint );
   typename ConfiguratorType::VectorType initPath;
   linspace<typename ConfiguratorType::VectorType>( StartPoint, EndPoint, initPath, totalLength-1 );
   geodesic.execute( initPath, geodPath, endPoints );   
}

//! \brief Discrete geodesic extrapolation
//! \author Heeren
template<typename ConfiguratorType, typename PairwiseEnergyType>
void integerExtrapolation( const typename ConfiguratorType::VectorType& StartPoint, 
                           const typename ConfiguratorType::VectorType& varPoint,
                           const PairwiseEnergyType& W,
                           int steps,
                           typename ConfiguratorType::VectorType& Path,
                           bool forwardShooting = true,
                           bool quiet = true ){
  ExpK<ConfiguratorType, PairwiseEnergyType> expK( W, 1000, 1e-8  );
  FixedPointsTrait<ConfiguratorType> initPoints( StartPoint, varPoint );
  if( forwardShooting )
      expK.execute( initPoints, steps, Path );
  else
      expK.executeBackwards( initPoints, steps, Path );
}

//! \brief One step of inverse discrete parallel transport cia Schild's ladder
//! \author Heeren
template<typename ConfiguratorType, typename PairwiseEnergyType>
void invParTranspSingleStepShort( const typename ConfiguratorType::VectorType& StartGeom,
                                  const typename ConfiguratorType::VectorType& EndGeom,
                                  const typename ConfiguratorType::VectorType& VarEndGeom,
                                  const PairwiseEnergyType& W,
                                  typename ConfiguratorType::VectorType& VarStartGeom,
                                  bool quiet = true ){  
  if( !quiet ) std::cerr << "Start inverse parallel transport with single step (i.e. K=1)."  << std::endl;
  typename ConfiguratorType::VectorType SchildLadderMidpoint;
  if( !quiet ) std::cerr << "1) Compute discrete 3-geodesic to get the midpoint."  << std::endl;
  integerInterpolation<ConfiguratorType, PairwiseEnergyType>( StartGeom, VarEndGeom, W, 3, SchildLadderMidpoint, quiet );
  //std::cerr << SchildLadderMidpoint.size() << std::endl;
  if( !quiet ) std::cerr << "2) Compute discrete EXP2 to get the inverse transported point." << std::endl;
  integerExtrapolation<ConfiguratorType, PairwiseEnergyType>( EndGeom, SchildLadderMidpoint, W, 1, VarStartGeom, false, quiet );
}

//! \brief Discrete first covariant derivative (via fully symmetrized central difference quotient)
//! \author Heeren
//! Convergence order \tau^2!
template<typename ConfiguratorType, typename PairwiseEnergyType>
void firstCovariantDerivative( const typename ConfiguratorType::VectorType& StartGeom,
                               const typename ConfiguratorType::VectorType& Eta,
                               const typename ConfiguratorType::VectorType& U,
                               const double Tau,
                               const double Epsilon,
                               const PairwiseEnergyType& W,
                               typename ConfiguratorType::VectorType& FirCovDer,
                               bool quiet = true ){  
  
  typename ConfiguratorType::VectorType yURight( StartGeom ), yULeft( StartGeom );
  yURight += Tau * U;
  yULeft  -= Tau * U;
    
  // Here we actually use that Eta is constant. Otherwise we need evaluations at this point!   
  typename ConfiguratorType::VectorType yURightEta( yURight ), yULeftEta( yULeft );                                   
  yURightEta += Epsilon * Eta;
  yULeftEta  -= Epsilon * Eta; //NOTE the sign here, this is where the proper reflection happens!
  
  // perform two inverse parallel transports
  typename ConfiguratorType::VectorType xiLeft;  // xiRight is stored in FirCovDer
  // NOTE inverse parallel transport returns points rather than tangent vectors 
  invParTranspSingleStepShort<ConfiguratorType, PairwiseEnergyType>( StartGeom, yURight, yURightEta, W, FirCovDer, quiet );
  invParTranspSingleStepShort<ConfiguratorType, PairwiseEnergyType>( StartGeom, yULeft,  yULeftEta, W, xiLeft, quiet );  
  
  // subtract position (i.e. reference shape) to get tangent vectors (i.e. displacements)
  FirCovDer -= StartGeom;
  xiLeft -= StartGeom;
  
  // compute sum (!!) of displacements due to change in sign in transported vector (see above)
  FirCovDer += xiLeft;
  FirCovDer /= 2. * Tau * Epsilon;
}

//! \brief Discrete second covariant derivative (via fully symmetrized central difference quotient)
//! \author Heeren
//! Convergence order \tau^2, if \tau is outer stepsize and \tau^{3/2} is inner stepsize!
template<typename ConfiguratorType, typename PairwiseEnergyType>
void secondCovariantDerivative( const typename ConfiguratorType::VectorType& StartGeom,
                                const typename ConfiguratorType::VectorType& Eta,
                                const typename ConfiguratorType::VectorType& outerDir,
                                const typename ConfiguratorType::VectorType& innerDir,                                
                                const double OuterTau,
                                const double InnerTau,
                                const PairwiseEnergyType& W,
                                typename ConfiguratorType::VectorType& SecCovDer,
                                bool quiet = true ){
  typename ConfiguratorType::RealType outerEps = OuterTau;
  typename ConfiguratorType::RealType innerEps = InnerTau;
  
  // let V be the outer direction
  typename ConfiguratorType::VectorType yVRight( StartGeom ), yVLeft( StartGeom );
  yVRight += OuterTau * outerDir;
  yVLeft  -= OuterTau * outerDir;
  
  // compute first covariant derivative at y_r = y + tau * V and y_l = y - tau * V, denoted by \xi_r and \xi_l, respectively
  // let \xi denote the vector field representing the inner covariant derivative.
  typename ConfiguratorType::VectorType xiRight, xiLeft;
  firstCovariantDerivative<ConfiguratorType, PairwiseEnergyType>( yVRight, Eta, innerDir, InnerTau, innerEps, W, xiRight, quiet );
  firstCovariantDerivative<ConfiguratorType, PairwiseEnergyType>( yVLeft, Eta, innerDir, InnerTau, innerEps, W, xiLeft, quiet );

  // get variational shapes from shape differences \xi at y_r and y_l
  typename ConfiguratorType::VectorType yVRightXi( yVRight ), yVLeftXi( yVLeft ), temp;
  yVRightXi += outerEps * xiRight;
  yVLeftXi  -= outerEps * xiLeft; //NOTE the sign here!

  // NOTE inverse parallel transport returns points rather than tangent vectors 
  invParTranspSingleStepShort<ConfiguratorType, PairwiseEnergyType>( StartGeom, yVRight, yVRightXi, W, SecCovDer, quiet );
  invParTranspSingleStepShort<ConfiguratorType, PairwiseEnergyType>( StartGeom, yVLeft, yVLeftXi, W, temp, quiet );    
  SecCovDer -= StartGeom;
  temp -= StartGeom;  
  // compute sum (!!) of displacements due to change in sign in transported vector (see above)
  SecCovDer += temp;
  SecCovDer /= 2. * OuterTau * outerEps;    
}

//! \brief Discrete Riemann curvature tensor for constant vector fields \eta 
//! \author Heeren
//! \f$ R(U,V)\eta = \nabla_u \nabla_v \eta - \nabla_v \nabla_u \eta \f$
template<typename ConfiguratorType, typename PairwiseEnergyType>
void curvatureTensor( const typename ConfiguratorType::VectorType& StartGeom,
                      const typename ConfiguratorType::VectorType& Eta,
                      const typename ConfiguratorType::VectorType& U,
                      const typename ConfiguratorType::VectorType& V,
                      const double Tau,
                      const PairwiseEnergyType& W,
                      typename ConfiguratorType::VectorType& RUVEta,
                      bool quiet = true){
  typename ConfiguratorType::RealType innerTau = Tau * std::sqrt( Tau );    
  typename ConfiguratorType::VectorType secondCovarDerVUEta; // secCovDerUV is stored in RUVEta
  secondCovariantDerivative<ConfiguratorType, PairwiseEnergyType>( StartGeom, Eta, U, V, Tau, innerTau, W, RUVEta, quiet );
  secondCovariantDerivative<ConfiguratorType, PairwiseEnergyType>( StartGeom, Eta, V, U, Tau, innerTau, W, secondCovarDerVUEta, quiet );
  RUVEta -= secondCovarDerVUEta;
}

//! \brief Discrete sectional curvature for constant vector fields \eta 
//! \author Heeren
//! \f$ \kappa(U,V) = g( U, R(U,V) V ) / [g(U,U) g(V,V) - g(U,V)^2]\f$, where g is metric (i.e first fundamental form of surface).
template<typename ConfiguratorType, typename PairwiseEnergyType>
double sectionalCurvature( const typename ConfiguratorType::VectorType& StartGeom,
                           const typename ConfiguratorType::VectorType& U,
                           const typename ConfiguratorType::VectorType& V,
                           const double Tau,
                           const PairwiseEnergyType& W,
                           bool quiet = true ){

  typename ConfiguratorType::VectorType rUVV, metricEval;
  if(!quiet) std::cerr << "Compute R_p(U,V)V..." << std::endl; 
  curvatureTensor<ConfiguratorType, PairwiseEnergyType>( StartGeom, V, U, V, Tau, W, rUVV, quiet );
  
  typename ConfiguratorType::FullMatrixType g;
  W.getSurfaceParam ( ).getMetricTensor( StartGeom, g );
  
  metricEval = g * U;
  typename ConfiguratorType::RealType gUU = U.dot( metricEval );
  metricEval = g * V;  
  typename ConfiguratorType::RealType gVV = V.dot( metricEval );
  typename ConfiguratorType::RealType gUV = U.dot( metricEval );  
  metricEval = g * rUVV;
  
  typename ConfiguratorType::RealType denom = gUU * gVV - gUV * gUV;
  typename ConfiguratorType::RealType secCurv = U.dot( metricEval ) / denom;
  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;  
  if( !quiet ) std::cerr << "Sectional curvature kappa_p(U,V) = g(U, R(U,V)V) / [g(U,U)g(V,V)-g(U,V)^2] = " << U.dot( metricEval ) << " / ( " << gUU << " * " << gVV << " - (" << gUV << ")^2 ) = "; 
  if( !quiet ) std::cerr << U.dot( metricEval ) << " /  " << denom << " = " << secCurv << std::endl; 
  if( !quiet ) std::cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl << std::endl;  
  
  return  secCurv;
}


//========================================================================================================
//========================================================================================================
//========================================================================================================

//! \brief Default vector field
//! \author Effland
template < typename ConfiguratorType, typename Imp >
class DiscreteVectorFieldDefault {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType;

  void applyAdd ( const VectorType& Arg, VectorType& Dest ) const {
    if ( Arg.size ( ) != getDimOfDomain ( ) || Dest.size ( ) != getDimOfRange ( ) )
      throw BasicException ( "DiscreteVectorFieldDefault::applyAdd(): Size mismatch!");
    this->asImp().evaluate ( Arg, Dest );
  }

  void apply ( const VectorType& Arg, VectorType& Dest ) const {
    Dest.setZero();
    applyAdd ( Arg, Dest );
  }

  void getJacobianMatrix ( const VectorType& BasePoint, MatrixType& Jacobian ) const {
    if ( Jacobian.getNumCols ( ) != getDimOfDomain ( ) || Jacobian.getNumRows ( ) != getDimOfRange ( ) )
      throw BasicException ( "DiscreteVectorFieldDefault::getJacobianMatrix(): Size mismatch!");
    this->asImp().getJacobian ( BasePoint, Jacobian );
  }

  // \partial_l\partial_i\eta^n = Hessian[l].get ( n , i )
  void getHessianMatrix ( const VectorType& BasePoint, std::vector< MatrixType >& Hessian ) const {
    for ( unsigned short int i = 0; i < getDimOfDomain ( ); ++i )
      if ( Hessian[i].getNumCols ( ) != getDimOfDomain ( ) || Hessian[i].getNumRows ( ) != getDimOfRange ( ) )
        throw BasicException ( "DiscreteVectorFieldDefault::getHessianMatrix(): Size mismatch!" );
    this->asImp().getHessian ( BasePoint, Hessian );
  }

  unsigned short int getDimOfDomain ( ) const {
    return this->asImp().getDimOfDomain ( );
  }

  unsigned short int getDimOfRange ( ) const {
    return this->asImp().getDimOfRange ( );
  }

protected:
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
};

//! \brief Constant vector field
//! \author Effland
template < typename ConfiguratorType, int Dim >
class DiscreteVectorFieldConstant : public DiscreteVectorFieldDefault < ConfiguratorType, DiscreteVectorFieldConstant < ConfiguratorType, Dim > > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::FullMatrixType MatrixType;

  static const unsigned short int DimOfDomain = Dim;
  static const unsigned short int DimOfRange = Dim;

private:
  const VectorType _vector;

public:
  DiscreteVectorFieldConstant ( const VectorType& ConstantVector ) : _vector ( ConstantVector ) { }

  void evaluate ( const VectorType&, VectorType& Dest ) const {
    Dest = _vector;
  }

  void getJacobian ( const VectorType&, MatrixType& Jacobian ) const {
    Jacobian.setZero ( );
  }

  void getHessian ( const VectorType&, std::vector< MatrixType >& Hessian ) const {
    for ( int i = 0; i < Hessian.size ( ); ++i )
      Hessian[i].setZero ( );
  }
};


//========================================================================================================
//========================================================================================================
//========================================================================================================

//! \brief Class to automatically generate surface plots in Matlab (NOT CHECKED!)
//! \author Effland
template < typename ConfiguratorType, typename SurfaceParametrization >
class SurfaceMatlabPlotter {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  
  const SurfaceParametrization& _surface;
  const std::string _outFileName;
  mutable std::ofstream _out;

  void clear ( ) {
    _out << "close all;" << "\n";
    _out << "clear all;" << "\n";
  }

public:
  SurfaceMatlabPlotter ( const SurfaceParametrization& Surface, const std::string OutFileName )
: _surface ( Surface ), _outFileName ( OutFileName ), _out ( _outFileName.c_str ( ) ) {
    clear ( );
  }

  ~SurfaceMatlabPlotter ( ) {
    _out << "hold off;" << "\n";
    _out.close ( );
  }

  void addColorbar ( const unsigned short int FontSize = 14 ) {
    _out << "cb = colorbar ( 'southoutside' );\n";
    _out << "set ( cb, 'fontsize', " << FontSize << " )\n";
  }

  void generateSurface ( const unsigned short int USteps, const unsigned short int VSteps ) const {
    _out << "hold on;" << "\n";
    _out << "u = linspace ( 0, 2 * pi, "<< USteps << ");" << "\n";
    _out << "v = linspace ( 0, 2 * pi, "<< VSteps << ");" << "\n";
    _out << "[u,v] = meshgrid ( u, v );" << "\n";
    _surface.plotSurface ( _out );
    _out << "surfaceMesh = mesh ( x, y, z );" << "\n";
    _out << "set ( surfaceMesh, 'EdgeColor', [.6,.6,.6], 'FaceAlpha', 0.5, 'EdgeAlpha', 0.8 );" << "\n";
    _out << "axis equal;" << "\n";
  }

  void createColoredSurface ( const unsigned short int USteps, const unsigned short int VSteps, const VectorType& Colors, const RealType xMin, const RealType xMax, const RealType yMin, const RealType yMax ) const {
    _out << "hold on;" << "\n";
    _out << "u = linspace ( " << xMin << ", " << xMax << ", " << USteps << ");" << "\n";
    _out << "v = linspace ( " << yMin << ", " << yMax << ", " << VSteps << ");" << "\n";
    _out << "[u,v] = meshgrid ( u, v );" << "\n";
    _surface.plotSurface ( _out );
    _out << "c = [\n";
    for ( unsigned short int y = 0; y < VSteps; ++y ) {
      for ( unsigned short int x = 0; x < USteps; ++x ) {
        _out << Colors[y*USteps+x] << " ";
      }
      _out << ";\n";
    }
    _out << "];\n";
    _out << "view ( [ -25, 40 ] );\n";
    _out << "set ( gca, 'fontsize', 20 );\n";
    _out << "surface = surf ( x, y, z, c );" << "\n";
    _out << "set ( gca, 'visible', 'off', 'FontSize', 30 );\n";
    _out << "set ( surface, 'edgecolor', 'none', 'LineStyle', 'none', 'facealpha', 0, 'linewidth', 1 );\n";
    _out << "axis equal;" << "\n";
  }

  void addCurve ( const VectorType& Curve, const VecType& Color, const FixedPointsTrait<ConfiguratorType>* Fixedpoints = NULL, const bool ShowPoints = false, const VectorType* ParallelTransportedVectors = NULL ) const {
    switch ( Curve.getEqualComponentSize ( ) ) {
    case 2: {
      VectorType CurveY ( Curve.numComponents ( ), 3 );
      _surface.apply ( Curve, CurveY );
      if ( Fixedpoints ) {
        VectorType startY ( 3 ), endY ( 3 );
        _surface.apply ( Fixedpoints->getStartVector ( ), startY );
        _surface.apply ( Fixedpoints->getEndVector ( ), endY );
        FixedPointsTrait<ConfiguratorType> FixedpointsY ( startY, endY );
        plotCurve ( CurveY, Color, &FixedpointsY, ShowPoints, ParallelTransportedVectors );
      }
      else {
        plotCurve ( CurveY, Color, NULL, ShowPoints, ParallelTransportedVectors );
      }
    }
    break;

    case 3:
      plotCurve ( Curve, Color, Fixedpoints, ShowPoints, ParallelTransportedVectors );
      break;

    default:
      throw BasicException ( "SurfaceMatlabPlotter:addCurve(): Wrong dimension!" );
      break;
    }
  }

  void markPoint ( const VectorType& Point, const VecType& Color ) const {
    switch ( Point.size ( ) ) {
    case 2: {
      VectorType PointY ( 3 );
      _surface.apply ( Point, PointY );
      plotPoint ( PointY, Color );
    }
    break;

    case 3:
      plotPoint ( Point, Color );
      break;

    default:
      throw BasicException ( "SurfaceMatlabPlotter::markPoint(): Size mismatch!");
      break;
    }
  }

private:
  void plotCurve ( const VectorType& Curve, const VecType& Color, const FixedPointsTrait <ConfiguratorType>* Fixedpoints = NULL, const bool ShowPoints = false, const VectorType* ParallelTransportedVectors = NULL ) const {
    for ( unsigned short int i = 0; i < 3; ++i ) {
      _out << "y" << i + 1 << " = [";
      if ( Fixedpoints )
        _out << Fixedpoints->getStartVector( )[i] << " ";
      for ( int j = 0; j < Curve.numComponents ( ); ++j )
        _out << Curve[j][i] << " ";
      if ( Fixedpoints )
        _out << Fixedpoints->getEndVector( )[i] << " ";
      _out << "];" << "\n";
    }
    _out << "plot3 ( y1, y2, y3, 'Linewidth', 6, 'Color', [" << Color[0] << "," << Color[1] << "," << Color[2] << "] );" << "\n";

    if ( ShowPoints ) {
      _out << "plot3 ( y1(:), y2(:), y3(:), 'Marker', '.', 'MarkerSize', 60, 'Color', [" << Color[0] << "," << Color[1] << "," << Color[2] << "]);" << "\n";
    }

    if ( ParallelTransportedVectors ) {
      for ( unsigned short int i = 0; i < 3; ++i ) {
        _out << "vector" << i + 1 << " = [";
        for ( int j = 0; j < ParallelTransportedVectors->numComponents ( ); ++j )
          _out << (*ParallelTransportedVectors)[j][i] << " ";
        _out << "];" << "\n";
      }
      // deliberately changed color


      _out << "for i = 1 : size ( y1, 2 )" << "\n";
      _out << "plot3 ( [y1(i); y1(i) + vector1(i)], [y2(i); y2(i) + vector2(i)], [y3(i); y3(i) + vector3(i)], 'Linewidth', 6, 'Color', [" << Color[2] << "," << Color[1] << "," << Color[0] << "] )" << "\n";
      _out << "end" << "\n";
    }
  }

  void plotPoint ( const VecType& Point, const VecType& Color ) const {
    _out << "plot3 ( " << Point[0] << "," << Point[1] << "," << Point[2] << ", 'Marker', '.', 'MarkerSize', 60, 'Color', [" << Color[0] << "," << Color[1] << "," << Color[2] << "]);" << "\n";
  }
};


#endif
