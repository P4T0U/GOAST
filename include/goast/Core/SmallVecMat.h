// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//=============================================================================
//
//  Small matrix and vector classes
//
//=============================================================================

#ifndef SMALLVECMAT_HH
#define SMALLVECMAT_HH


//== INCLUDES =================================================================
#include <iostream>
#include "Auxiliary.h"
#include <algorithm>

//=============================================================================
// CLASS FOR SMALL VECTORS
//=============================================================================
template<int dimension, typename FloatType>
class SmallVec {

protected:
// Don't warn if dimension == 0.
  FloatType _coords[dimension];

public:
  //! Constructor for creating SmallVec from c array
  explicit SmallVec( const FloatType rhs[dimension] ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = rhs[i];
  }

  //! Constructor for SmallVec with constant entries
  explicit SmallVec( const FloatType Scalar ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = Scalar;
  }

  //! Standard constructor
  SmallVec() {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = 0.;
  }

  //! Copy-constructor
  SmallVec( const SmallVec<dimension, FloatType> &rhs ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = rhs[i];
  }

  //! Conversion-Copy-constructor
  template<class AnotherSmallVecType>
  explicit SmallVec( const AnotherSmallVecType &rhs ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = rhs[i];
  }

  //! operator=
  SmallVec<dimension, FloatType> &operator=( const SmallVec<dimension, FloatType> &rhs ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = rhs._coords[i];
    return *this;
  }

  //! Set the SmallVec to zero.
  void setZero() {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = 0.;
  }

  //! Bracket operator.
  inline FloatType &operator[]( int I ) {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= dimension )
      throw BasicException ( "SmallVec::operator[] out of bounds!" );
#endif
    return _coords[I];
  }

  //! bracket operator
  inline const FloatType &operator[]( int I ) const {
#ifdef BOUNDS_CHECK
    if ( I < 0 || I >= dimension ) 
      throw BasicException ( "SmallVec::operator[] out of bounds!" );
#endif
    return _coords[I];
  }

  //! Scalar product operator
  FloatType operator*( const SmallVec<dimension, FloatType> &Other ) const {
    FloatType scalarProduct = 0;
    for ( int i = 0; i < dimension; ++i )
      scalarProduct += _coords[i] * Other[i];
    return scalarProduct;
  }

  //! multiplication by scalar
  SmallVec<dimension, FloatType> operator*( FloatType alpha ) const {
    SmallVec<dimension, FloatType> res;
    for ( int i = 0; i < dimension; ++i )
      res[i] = alpha * _coords[i];
    return res;
  }

  //! division by scalar
  SmallVec<dimension, FloatType> operator/( FloatType alpha ) const {
    SmallVec<dimension, FloatType> res;
    for ( int i = 0; i < dimension; ++i )
      res[i] = _coords[i] / alpha;
    return res;
  }

  //! add two vectors and return the result
  SmallVec<dimension, FloatType> operator+( const SmallVec<dimension, FloatType> &Other ) const {
    SmallVec<dimension, FloatType> res;
    for ( int i = 0; i < dimension; ++i )
      res[i] = _coords[i] + Other[i];
    return res;
  }

  //! subtract two vectors and return the result
  SmallVec<dimension, FloatType> operator-( const SmallVec<dimension, FloatType> &Other ) const {
    SmallVec<dimension, FloatType> res;
    for ( int i = 0; i < dimension; ++i )
      res[i] = _coords[i] - Other[i];
    return res;
  }

  //! add other vector
  SmallVec<dimension, FloatType> &operator+=( const SmallVec<dimension, FloatType> &Other ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] += Other[i];
    return *this;
  }

  //! subtract other vector
  SmallVec<dimension, FloatType> &operator-=( const SmallVec<dimension, FloatType> &Other ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] -= Other[i];
    return *this;
  }

  //! multiply by scalar
  SmallVec<dimension, FloatType> &operator*=( FloatType Alpha ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] *= Alpha;
    return *this;
  }

  //! divide by scalar
  SmallVec<dimension, FloatType> &operator/=( FloatType Alpha ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] /= Alpha;
    return *this;
  }

  //! add multiple of other vector
  void addMultiple( const SmallVec<dimension, FloatType> &AddedVec, FloatType Factor ) {
    for ( int i = 0; i < dimension; ++i )
      _coords[i] += Factor * AddedVec[i];
  }

  FloatType normSqr() const {
    FloatType res = 0;
    for ( int i = 0; i < dimension; ++i )
      res += ( _coords[i] * _coords[i] );
    return res;
  }

  FloatType squaredNorm() const {
    FloatType res = 0;
    for ( int i = 0; i < dimension; ++i )
      res += ( _coords[i] * _coords[i] );
    return res;
  }

  FloatType norm() const {
    return std::sqrt( normSqr());
  }

  void normalize() {
    FloatType ns = norm();
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = _coords[i] / ns;
  }

  SmallVec<dimension, FloatType> normalized() const {
    SmallVec<dimension, FloatType> res( *this );
    res.normalize();
    return ( res );
  }

  /**
   *
   * \return the absolute largest coefficient
   */
  FloatType absMax() const {
    FloatType max = 0;
    for ( int i = 0; i < dimension; i++ )
      max = std::max( std::abs( _coords[i] ), max );
    return max;
  }

  FloatType stableNorm() const {
    FloatType max = absMax();
    FloatType res = 0;
    for ( int i = 0; i < dimension; ++i )
      res += (( _coords[i] / max ) * ( _coords[i] / max ));

    return max * std::sqrt( res );
  }

  void stableNormalize() {
    FloatType max = absMax();
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = _coords[i] / max;

    FloatType ns = norm();
    for ( int i = 0; i < dimension; ++i )
      _coords[i] = _coords[i] / ns;
  }


  std::ostream &print( std::ostream &os ) const {
    os << "( ";
    for ( int i = 0; i < dimension; ++i )
      os << _coords[i] << " ";
    os << " )";
    return os;
  }
};


//=============================================================================
// THE CLASS SMALLVEC2 DERIVED FROM SMALLVEC
//=============================================================================
template<typename FloatType>
class SmallVec2 : public SmallVec<2, FloatType> {

public:
  //! Constructor for two values
  explicit SmallVec2( const FloatType X, const FloatType Y ) : SmallVec<2, FloatType>() {
    this->coords[0] = X;
    this->coords[1] = Y;
  }

  //! Constructor for c-array
  explicit SmallVec2( const FloatType rhs[2] ) : SmallVec<2, FloatType>( rhs ) {}

  //! Constructor for constant vector
  explicit SmallVec2( const FloatType Scalar ) : SmallVec<2, FloatType>( Scalar ) {}

  //! Standard constructor
  SmallVec2() : SmallVec<2, FloatType>() {}

  //! Copy-constructor
  SmallVec2( const SmallVec2<FloatType> &rhs ) : SmallVec<2, FloatType>( rhs ) {}
};


//=============================================================================
// THE CLASS SMALLVEC3 DERIVED FROM SMALLVEC
//=============================================================================
template<typename FloatType>
class SmallVec3 : public SmallVec<3, FloatType> {

public:
  //! Constructor for two values
  explicit SmallVec3( const FloatType X, const FloatType Y, const FloatType Z ) : SmallVec<3, FloatType>() {
    this->_coords[0] = X;
    this->_coords[1] = Y;
    this->_coords[2] = Z;
  }

  //! Constructor for c-array
  explicit SmallVec3( const FloatType rhs[3] ) : SmallVec<3, FloatType>( rhs ) {}

  //! Constructor for constant vector
  explicit SmallVec3( const FloatType Scalar ) : SmallVec<3, FloatType>( Scalar ) {}

  //! Standard constructor
  SmallVec3() : SmallVec<3, FloatType>() {}

  //! Copy-constructor
  SmallVec3( const SmallVec3<FloatType> &rhs ) : SmallVec<3, FloatType>( rhs ) {}

  SmallVec3( const SmallVec<3, FloatType> &rhs ) : SmallVec<3, FloatType>( rhs ) {}

  //! set this Vec3 to the cross product (vector product) of two other vectors
  void makeCrossProduct( const SmallVec3<FloatType> &a, const SmallVec3<FloatType> &b ) {
    this->_coords[0] = a[1] * b[2] - a[2] * b[1];
    this->_coords[1] = a[2] * b[0] - a[0] * b[2];
    this->_coords[2] = a[0] * b[1] - a[1] * b[0];
  }

  //! compute cross product (vector product) of this SmallVec3 with another SmallVec3
  SmallVec3<FloatType> crossProduct( const SmallVec3<FloatType> &b ) const {
    SmallVec3<FloatType> res;
    res[0] = this->_coords[1] * b[2] - this->_coords[2] * b[1];
    res[1] = this->_coords[2] * b[0] - this->_coords[0] * b[2];
    res[2] = this->_coords[0] * b[1] - this->_coords[1] * b[0];
    return res;
  }

  //! subtract two vectors and return the result
  SmallVec3<FloatType> operator-( const SmallVec3<FloatType> &Other ) const {
    SmallVec3<FloatType> res;
    for ( int i = 0; i < 3; ++i )
      res[i] = this->_coords[i] - Other[i];
    return res;
  }

  //! subtract two vectors and return the result
  SmallVec3<FloatType> operator+( const SmallVec3<FloatType> &Other ) const {
    SmallVec3<FloatType> res;
    for ( int i = 0; i < 3; ++i )
      res[i] = this->_coords[i] + Other[i];
    return res;
  }

  SmallVec3<FloatType> operator+( const SmallVec<3, FloatType> &Other ) const {
    SmallVec3<FloatType> res;
    for ( int i = 0; i < 3; ++i )
      res[i] = this->_coords[i] + Other[i];
    return res;
  }

  SmallVec3<FloatType> operator-( const SmallVec<3, FloatType> &Other ) const {
    SmallVec3<FloatType> res;
    for ( int i = 0; i < 3; ++i )
      res[i] = this->_coords[i] - Other[i];
    return res;
  }

  SmallVec3<FloatType> operator*( FloatType alpha ) const {
    SmallVec3<FloatType> res;
    for ( int i = 0; i < 3; ++i )
      res[i] = alpha * this->_coords[i];
    return res;
  }

  //! Scalar product operator
  FloatType operator*( const SmallVec3<FloatType> &Other ) const {
    FloatType scalarProduct = 0;
    for ( int i = 0; i < 3; ++i )
      scalarProduct += this->_coords[i] * Other[i];
    return scalarProduct;
  }

  FloatType &x() {
    return this->_coords[0];
  }

  FloatType &y() {
    return this->_coords[1];
  }

  FloatType &z() {
    return this->_coords[2];
  }

  const FloatType &x() const {
    return this->_coords[0];
  }

  const FloatType &y() const {
    return this->_coords[1];
  }

  const FloatType &z() const {
    return this->_coords[2];
  }
};

template<int dimension, typename FloatType>
SmallVec<dimension, FloatType> operator*( FloatType x, const SmallVec<dimension, FloatType> &y ) {
  SmallVec<dimension, FloatType> Res( y );
  Res *= x;
  return Res;
}

template<int dimension, typename FloatType>
SmallVec<dimension, FloatType> operator-( const SmallVec<dimension, FloatType> &y ) {
  SmallVec<dimension, FloatType> Res( y );
  Res *= -1.;
  return Res;
}

//=============================================================================
// Input/Output operators
//=============================================================================
template<int dim, class T>
inline std::ostream &operator<<( std::ostream &os, const SmallVec<dim, T> &v ) {
  return v.print( os );
}

template<class T>
inline std::ostream &operator<<( std::ostream &os, const SmallVec2<T> &v ) {
  return v.print( os );
}

template<class T>
inline std::ostream &operator<<( std::ostream &os, const SmallVec3<T> &v ) {
  return v.print( os );
}


//=============================================================================
// Small Matrix Class
//=============================================================================
template<int numRows, int numCols, typename FloatType>
class SmallMat {

protected:
// Don't warn if numRows == 0.
  std::vector<SmallVec<numCols, FloatType> > _row;

public:

  //! Constructor
  SmallMat() : _row( numRows ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        _row[i][j] = 0.;
  }

  //! Copy-constructor
  SmallMat( const SmallMat<numRows, numCols, FloatType> &rhs ) : _row( numRows ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        _row[i][j] = rhs.get( i, j );
  }

  FloatType get( const int I, const int J ) const {
    rowBoundsCheck( I );
    return _row[I][J];
  }

  void set( const int I, const int J, FloatType val ) {
    rowBoundsCheck( I );
    _row[I][J] = val;
  }

  FloatType operator()( const int i, const int j ) const {
    rowBoundsCheck( i );
    return _row[i][j];
  }

  FloatType &operator()( const int i, const int j ) {
    rowBoundsCheck( i );
    return _row[i][j];
  }

  SmallVec<numCols, FloatType> &operator[]( const int i ) {
    rowBoundsCheck( i );
    return _row[i];
  }

  const SmallVec<numCols, FloatType> &operator[]( const int i ) const {
    rowBoundsCheck( i );
    return _row[i];
  }

  int getNumCols() const { return numCols; }

  int getNumRows() const { return numRows; }


  void setZero() {
    for ( int i = 0; i < numRows; ++i )
      _row[i].setZero();
  }

  void setIdentity() {
    if ( numRows != numCols )
      throw BasicException( "SmallMat::setIdentity(): not square matrix!" );
    setZero();
    for ( int i = 0; i < numRows; ++i )
      _row[i][i] = 1.;
  }

  //! subtract other matrix
  SmallMat<numRows, numCols, FloatType> &operator-=( const SmallMat<numRows, numCols, FloatType> &Other ) {
    for ( int i = 0; i < numRows; ++i )
      _row[i] -= Other[i];
    return *this;
  }

  //! \f$ A\mapsto b*A \f$
  SmallMat<numRows, numCols, FloatType> &operator*=( FloatType Val ) {
    for ( int i = 0; i < numRows; ++i )
      this->_row[i] *= Val;
    return *this;
  }

  // this =  A*B
  template<int dimBoth>
  void makeProduct( const SmallMat<numRows, dimBoth, FloatType> &A,
                    const SmallMat<dimBoth, numCols, FloatType> &B,
                    FloatType factor = 1. ) {
    setZero();
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        for ( int k = 0; k < dimBoth; ++k )
          this->_row[i][j] += factor * A[i][k] * B[k][j];
  }

  // this =  A^T*B
  template<int dimBoth>
  void makeProductAtransposedB( const SmallMat<dimBoth, numRows, FloatType> &A,
                                const SmallMat<dimBoth, numCols, FloatType> &B ) {
    setZero();
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        for ( int k = 0; k < dimBoth; ++k )
          this->_row[i][j] += A[k][i] * B[k][j];
  }

  // this =  A*B^T
  template<int dimBoth>
  void makeProductABtransposed( const SmallMat<numRows, dimBoth, FloatType> &A,
                                const SmallMat<numCols, dimBoth, FloatType> &B ) {
    setZero();
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        for ( int k = 0; k < dimBoth; ++k )
          this->_row[i][j] += A[i][k] * B[j][k];

  }

  //! \f$ A \mapsto A + \alpha B \f$
  SmallMat<numRows, numCols, FloatType> &
  addMultiple( const SmallMat<numRows, numCols, FloatType> &mat, const FloatType &alpha ) {
    if ( static_cast<const void *>(this) == static_cast<const void *>(&mat))
      throw BasicException( "SmallMat<FloatType>::addMultiple :  don't add the same matrix!" );
    else {
      for ( int j = 0; j < numCols; ++j )
        for ( int i = 0; i < numRows; ++i )
          this->_row[i][j] += alpha * mat[i][j];
    }
    return *this;
  }

  //! computes the trace of this matrix.
  FloatType tr() const {
    FloatType trace = 0.;
    for ( int i = 0; i < numRows; ++i )
      trace += this->_row[i][i];
    return trace;
  }

  //! computes the squared Frobius norm of this matrix.
  FloatType normSqr() const {
    FloatType normSqr = 0.;
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        normSqr += ( this->_row[i][j] ) * ( this->_row[i][j] );
    return normSqr;
  }

  //! computes the squared Frobius norm of this matrix.
  FloatType squaredNorm() const {
    FloatType normSqr = 0.;
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        normSqr += ( this->_row[i][j] ) * ( this->_row[i][j] );
    return normSqr;
  }

  // this = a*b^T
  void makeTensorProduct( const SmallVec<numRows, FloatType> &a, const SmallVec<numCols, FloatType> &b ) {
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        this->_row[i][j] = a[i] * b[j];
  }

  SmallVec3<FloatType> operator*( const SmallVec3<FloatType> &vec ) const {
    SmallVec3<FloatType> Res;
    for ( int i = 0; i < numRows; ++i )
      for ( int j = 0; j < numCols; ++j )
        Res[i] += this->_row[i][j] * vec[j];
    return ( Res );
  }

  void multAdd( const SmallVec<numCols, FloatType> &Arg, SmallVec<numRows, FloatType> &Dest ) const {
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j )
        Dest[i] += this->_row[i][j] * Arg[j];
    }
  }

  void mult( const SmallVec<numCols, FloatType> &Arg, SmallVec<numRows, FloatType> &Dest ) const {
    Dest.setZero();
    multAdd( Arg, Dest );
  }

  void multAddTransp( const SmallVec<numCols, FloatType> &Arg, SmallVec<numRows, FloatType> &Dest ) const {
    for ( int i = 0; i < numRows; ++i ) {
      for ( int j = 0; j < numCols; ++j )
        Dest[i] += this->_row[j][i] * Arg[j];
    }
  }

  void multTransp( const SmallVec<numCols, FloatType> &Arg, SmallVec<numRows, FloatType> &Dest ) const {
    Dest.setZero();
    multAddTransp( Arg, Dest );
  }

  std::ostream &print( std::ostream &os ) const {
    for ( int i = 0; i < numRows; ++i )
      os << this->_row[i] << std::endl;
    return os;
  }

protected:
#ifdef BOUNDS_CHECK
  inline void rowBoundsCheck ( const int i ) const {
    if ( i < 0 || i >= numRows )
      throw( BasicException ( "SmallMat<>: row out of bounds!" ) );
  }
#else

  inline void rowBoundsCheck( const int ) const {}

#endif

};


//! 
template<typename FloatType>
class SmallMat33 : public SmallMat<3, 3, FloatType> {

public:
  //! Constructor
  SmallMat33() : SmallMat<3, 3, FloatType>() {}

  SmallMat33( FloatType xx, FloatType xy, FloatType xz, FloatType yx, FloatType yy, FloatType yz, FloatType zx,
              FloatType zy, FloatType zz ) : SmallMat<3, 3, FloatType>() {
    this->_row[0][0] = xx;
    this->_row[0][1] = xy;
    this->_row[0][2] = xz;

    this->_row[1][0] = yx;
    this->_row[1][1] = yy;
    this->_row[1][2] = yz;

    this->_row[2][0] = zx;
    this->_row[2][1] = zy;
    this->_row[2][2] = zz;
  }

  // Copy-constructor
  SmallMat33( const SmallMat33<FloatType> &rhs ) : SmallMat<3, 3, FloatType>( rhs ) {}

  // Copy-constructor
  SmallMat33( const SmallMat<3, 3, FloatType> &rhs ) : SmallMat<3, 3, FloatType>( rhs ) {}

  // calc the inverse matrix with adjoints
  SmallMat33<FloatType> inverse() const {
    FloatType determinant( det());
    if ( determinant == FloatType( 0 ))
      throw ( BasicException( "SmallMat33::inverse() matrix not invertible." ));

    SmallMat33<FloatType> res;
    FloatType A;

    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        // calc the determinant, the sign is automatically ok
        A = this->_row[( i + 1 ) % 3][( j + 1 ) % 3] * this->_row[( i + 2 ) % 3][( j + 2 ) % 3]
            - this->_row[( i + 1 ) % 3][( j + 2 ) % 3] * this->_row[( i + 2 ) % 3][( j + 1 ) % 3];
        // consider the transposed matrix:
        res.set( j, i, A / determinant );
      }
    }

    return res;
  }

  //! computes the determinant of this matrix.
  FloatType det() const {
    return this->_row[0][0] * this->_row[1][1] * this->_row[2][2] +
           this->_row[0][1] * this->_row[1][2] * this->_row[2][0] +
           this->_row[0][2] * this->_row[1][0] * this->_row[2][1] -
           this->_row[0][2] * this->_row[1][1] * this->_row[2][0] -
           this->_row[0][1] * this->_row[1][0] * this->_row[2][2] -
           this->_row[0][0] * this->_row[1][2] * this->_row[2][1];
  }

  //! \f$ A\mapsto A*B \f$
  SmallMat33<FloatType> &operator*=( const SmallMat33<FloatType> &Other ) {
    SmallMat33<FloatType> tmp;
    tmp = *this;
    this->setZero();
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        for ( int k = 0; k < 3; ++k ) {
          this->_row[i][j] += tmp._row[i][k] * Other._row[k][j];
        }
      }
    }
    return *this;
  }

  //! \f$ A\mapsto A+B \f$
  SmallMat33<FloatType> &operator+=( const SmallMat33<FloatType> &Other ) {
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        this->_row[i][j] += Other._row[i][j];
    return *this;
  }

  //! \f$ A\mapsto b*A \f$
  SmallMat33<FloatType> &operator*=( FloatType Val ) {
    for ( int i = 0; i < 3; ++i )
      this->_row[i] *= Val;
    return *this;
  }

  //! \f$ A\mapsto 1/b*A \f$
  SmallMat33<FloatType> &operator/=( FloatType Val ) {
    for ( int i = 0; i < 3; ++i )
      this->_row[i] /= Val;
    return *this;
  }

  //! \f$ A \mapsto A^T \f$
  void transpose() {
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = i + 1; j < 3; ++j ) {
        const FloatType t = this->_row[i][j];
        this->_row[i][j] = this->_row[j][i];
        this->_row[j][i] = t;
      }
    }
  }

  //! return (this matrix) transposed (two copies necessary)
  SmallMat33<FloatType> transposed() const {
    SmallMat33<FloatType> tmp = *this;
    tmp.transpose();
    return ( tmp );
  }


  SmallMat33<FloatType> operator-( const SmallMat33<FloatType> &other ) const {
    SmallMat33<FloatType> Res( *this );
    Res -= other;
    return Res;
  }

  SmallMat33<FloatType> operator+( const SmallMat33<FloatType> &other ) const {
    SmallMat33<FloatType> Res( *this );
    Res += other;
    return Res;
  }

  //!
  SmallVec3<FloatType> operator*( const SmallVec3<FloatType> &vec ) const {
    SmallVec3<FloatType> Res;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        Res[i] += this->_row[i][j] * vec[j];
    return ( Res );
  }

  SmallMat33<FloatType> operator*( const SmallMat33<FloatType> &other ) const {
    SmallMat33<FloatType> Res;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        for ( int k = 0; k < 3; ++k )
          Res[i][j] += ( *this )[i][k] * other[k][j];
    return Res;
  }

  SmallMat33<FloatType> operator*( FloatType other ) const {
    SmallMat33<FloatType> Res( *this );
    Res *= other;
    return Res;
  }

  SmallMat33<FloatType> operator/( FloatType other ) const {
    SmallMat33<FloatType> Res( *this );
    Res /= other;
    return Res;
  }

  //!
  void addToDiagonal( const FloatType Value ) {
    for ( int i = 0; i < 3; ++i )
      this->_row[i][i] += Value;
  }

  //
  void setCol( int idx, const SmallVec3<FloatType> &col ) {
    for ( int i = 0; i < 3; ++i )
      this->_row[i][idx] = col[i];
  }

  //
  void setRow( int idx, const SmallVec3<FloatType> &row ) {
    for ( int i = 0; i < 3; ++i )
      this->_row[idx][i] = row[i];
  }

};

template<typename FloatType>
SmallMat33<FloatType> operator*( FloatType x, const SmallMat33<FloatType> &y ) {
  SmallMat33<FloatType> Res( y );
  Res *= x;
  return Res;
}

template<int numRows, int numCols, typename FloatType>
SmallMat<numRows, numCols, FloatType> operator-( const SmallMat<numRows, numCols, FloatType> &y ) {
  SmallMat<numRows, numCols, FloatType> Res( y );
  Res *= -1.;
  return Res;
}

//=============================================================================
// Input/Output operators
//=============================================================================
template<typename T>
inline std::ostream &operator<<( std::ostream &os, const SmallMat33<T> &m ) {
  return m.print( os );
}

//=============================================================================
// Matrix-Vector operations
//=============================================================================

// dot product
template<typename RealType>
RealType dotProduct( const SmallVec3<RealType> &a, const SmallVec3<RealType> &b ) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// weighted vector sum
template<typename RealType>
void getWeightedVectorSum( RealType a, const SmallVec3<RealType> &vec1, RealType b, const SmallVec3<RealType> &vec2,
                           SmallVec3<RealType> &res ) {
  for ( int i = 0; i < 3; i++ )
    res[i] = a * vec1[i] + b * vec2[i];
}

// weighted matrix sum
template<typename RealType>
void getWeightedMatrixSum( RealType a, const SmallMat33<RealType> &mat1, RealType b, const SmallMat33<RealType> &mat2,
                           SmallMat33<RealType> &res ) {
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      res.set( i, j, a * mat1.get( i, j ) + b * mat2.get( i, j ));
}

// weighted matrix sum (transposed)
template<typename RealType>
void getWeightedMatrixSumTransposed( RealType a, const SmallMat33<RealType> &mat1, RealType b,
                                     const SmallMat33<RealType> &mat2, SmallMat33<RealType> &res ) {
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      res.set( i, j, a * mat1.get( j, i ) + b * mat2.get( j, i ));
}

//
template<typename RealType>
void setRotationAboutAxis( SmallMat33<RealType> &mat, const SmallVec3<RealType> &axis, RealType angle ) {

  RealType aux = 1. - std::cos( angle );

  mat.set( 0, 0, axis[0] * axis[0] * aux + std::cos( angle ));
  mat.set( 0, 1, axis[0] * axis[1] * aux - axis[2] * std::sin( angle ));
  mat.set( 0, 2, axis[0] * axis[2] * aux + axis[1] * std::sin( angle ));

  mat.set( 1, 0, axis[1] * axis[0] * aux + axis[2] * std::sin( angle ));
  mat.set( 1, 1, axis[1] * axis[1] * aux + std::cos( angle ));
  mat.set( 1, 2, axis[1] * axis[2] * aux - axis[0] * std::sin( angle ));

  mat.set( 2, 0, axis[2] * axis[0] * aux - axis[1] * std::sin( angle ));
  mat.set( 2, 1, axis[2] * axis[1] * aux + axis[0] * std::sin( angle ));
  mat.set( 2, 2, axis[2] * axis[2] * aux + std::cos( angle ));

}

//=============================================================================
// Vector-Vector operations
//=============================================================================
template<int numRows, int numCols, typename FloatType>
SmallMat<numRows, numCols, FloatType> tensorProduct( const SmallVec<numRows, FloatType> &a,
                                                     const SmallVec<numCols, FloatType> &b ) {
  SmallMat<numRows, numCols, FloatType> Res;
  for ( int i = 0; i < numRows; ++i )
    for ( int j = 0; j < numCols; ++j )
      Res( i, j ) = a[i] * b[j];

  return Res;
}

template<typename FloatType>
SmallMat33<FloatType> tensorProduct( const SmallVec3<FloatType> &a,
                                     const SmallVec3<FloatType> &b ) {
  SmallMat33<FloatType> Res;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      Res( i, j ) = a[i] * b[j];

  return Res;
}

#endif
