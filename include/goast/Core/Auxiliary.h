// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef AUXILIARY_HH
#define AUXILIARY_HH


#include <vector>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <fstream> 
#include <stdio.h>
#include <string.h>

#include "EigenIncludes.h"


  
//==========================================================================================================
// EXCEPTION CLASS
//==========================================================================================================
class BasicException : public std::exception {
public:
  // Exception with error message Message
  explicit BasicException ( const std::string& Message ) : M ( Message ){ }

  BasicException ( const BasicException& Other ) : M ( Other.getMessage() ){ }
  
  BasicException ( const std::string& Message, int &Extension ) : M ( Message + std::to_string(Extension)) {}
  BasicException ( const std::string& Message, double &Extension ) : M ( Message + std::to_string(Extension)) {}
  BasicException ( const std::string& Message, const std::string& Extension ) : M ( Message + Extension ) {}

  const std::string getMessage() const noexcept {
    return M.what();
  }

protected:
  //! The message of this exception
//  std::ostringstream _message;
  std::runtime_error M;
};


//=============================================================================
// Square and Cubic
//=============================================================================
template<typename RealType>
RealType SquareFnc( const RealType& val ){
    return val*val;
}

template<typename RealType>
RealType CubicFnc( const RealType& val ){
    return val*val*val;
}

//==========================================================================================================
// GENERIC TENSOR CLASS
//==========================================================================================================
template<typename MatrixType>
class GenericTensor : protected std::vector<MatrixType> {



public:
  int _rows, _cols;
  GenericTensor() : std::vector<MatrixType>(0) {}

  GenericTensor(int length, int rows, int cols) : std::vector<MatrixType>(length, MatrixType(rows,cols)), _rows(rows), _cols(cols){}

  void setZero() {
    for (auto it : (*this)){
      it.setZero();
    }
  }

  void resize(int length, int rows, int cols) {
    (*this).resize(length);
    _rows = rows;
    _cols = cols;
    for (int i = 0; i< length; ++i){
      (*this)[i].resize(_rows,_cols);
    }
  }

  using std::vector<MatrixType>::resize;

  template<typename TripletType>
  void setFromTriplets(std::vector<std::vector<TripletType> > _triplets){
    (*this).resize(_triplets.size());
    for (unsigned int i = 0; i< _triplets.size(); ++i){
      (*this)[i].setFromTriplets( _triplets[i].cbegin(), _triplets[i].cend() );
    }
  }

  GenericTensor<MatrixType>& operator+= ( const GenericTensor<MatrixType>& other ) {
    for (unsigned int i = 0; i< other.size(); ++i){
      (*this)[i] += other[i];
    }
    return *this;
  }

  double& operator() (int i, int j, int k) {
    return (*this)[i].coeffRef(j,k);
  }

  void apply( const MatrixType& Arg, MatrixType& Dest ) const {
    throw BasicException( "Tensor apply: has not been implemented yet!" );
  }

  // n-Vector is (nx1) - matrix!
  template<typename VectorType>
  void applyVector( const VectorType& Arg, MatrixType& Dest ) const {

    if ((uint)Arg.size() != (*this).size())
      throw BasicException( "Tensor applyVector: sizes don't match!" );

    if( Dest.rows() != _rows || Dest.cols() != _cols )
      throw BasicException( "Tensor applyVector: cols/rows don't match!" );

    Dest.setZero();
    for (uint i = 0; i < Arg.size(); ++i){
      //  std::cout << (*this)[i].rows() << " " << (*this)[i].cols() << std::endl;
      Dest += (*this)[i]*Arg[i];
    }

  }

  using std::vector<MatrixType>::operator[]; //to access entire first dimension of Tensor
  using std::vector<MatrixType>::begin;
  using std::vector<MatrixType>::end;
  using std::vector<MatrixType>::size;
};



//==========================================================================================================
// ACCESS TO MESH GEOMETRIES
//==========================================================================================================
//!
template<typename VectorType, typename VecType>
void getXYZCoord( const VectorType& X, const VectorType& Y, const VectorType& Z, VecType& coord, int num ){
    coord[0] = X[num];
    coord[1] = Y[num];
    coord[2] = Z[num];
}

//! order in geometry vector is (X,Y,Z) = (X_1, ..., X_n, Y_1, ..., Y_n, Z_1, ..., Z_n ), for n nodes
template<typename VectorType, typename VecType>
void getXYZCoord( const VectorType& geom, VecType& coord, int num ){
    int size = geom.size() / 3;
    for( int i = 0; i < 3; i++ )
        coord[i] = geom[i*size+num];
}

//! order in geometry vector is (X,Y,Z) = (X_1, ..., X_n, Y_1, ..., Y_n, Z_1, ..., Z_n ), for n nodes
template<typename VectorType, typename VecType>
void setXYZCoord( VectorType& geom, const VecType& coord, int num ){
    int size = geom.size() / 3;
    for( int i = 0; i < 3; i++ )
        geom[i*size+num] = coord[i];
}

//! \author Echelmeyer
//! order in geometry vector is (X,Y,Z) = (X_1, ..., X_n, Y_1, ..., Y_n, Z_1, ..., Z_n ), for n nodes
template<typename VectorType, typename VecType>
void addXYZCoord( VectorType& geom, const VecType& coord, int num ){
    int size = geom.size() / 3;
    for( int i = 0; i < 3; i++ )
        geom[i*size+num] += coord[i];
}

//! check vector for NaN entries
template<typename VectorType>
bool hasNanEntries( const VectorType& X ){
    for( int i = 0; i < X.size(); i++ )
        if( std::isnan(X[i]) )
            return true;
    return false;
}


//==========================================================================================================
// BITVECTOR FOR BOUNDARY MASK
//==========================================================================================================
class BitVector : protected std::vector<bool> {
    
  int _size;  
    
public:
  BitVector( int size, bool bVal = false ) : std::vector<bool>(size, bVal), _size(size) {}
  
  explicit BitVector ( const std::string& filename ) : _size ( 0 ){
    loadFromFile( filename );
  }
  
  void setAll( bool bVal ) {
    for( int i = 0; i < _size; i++ )
        (*this)[i] = bVal;
  }
  
  void set( int i, bool bVal ) {
      if( i < _size )
          (*this)[i] = bVal;
  }
  
  int size() const {
    return _size;    
  }
  
  void resize( int newsize ) {
    std::vector<bool>::resize(newsize);
    _size = newsize;    
  }
  
  using std::vector<bool>::operator[];
  
protected:
  void loadFromFile ( const std::string& filename ){
    throw BasicException("BitVector::loadFromFile(): has not been implemented yet!");    
  }
  
};

//==========================================================================================================
// PROCESS DATA
//==========================================================================================================
template<typename VectorType>
//! Crop data zu interval [minVal, maxVal]
void cropData( VectorType& data, double minVal, double maxVal ){
    for( int i = 0; i < data.size(); i++ )
        data[i] = std::min( std::max( minVal, data[i]), maxVal );
}

//
template<typename VectorType>
void scaleDataLinearly( VectorType& data, double minVal, double maxVal ){
    double minData = maxVal;
    double maxData = minVal;
    
    // determine range of data
    for( int i = 0; i < data.size(); i++ ){
        if( data[i] < minData )
            minData = data[i];
        if( data[i] > maxData )
            maxData = data[i];
    }
    
    double slope = (maxVal - minVal) / (maxData - minData);
    double offset = (maxVal + minVal - slope*(maxData + minData)) / 2.;
    
    // scale linearly 
    for( int i = 0; i < data.size(); i++ )
        data[i] = slope * data[i] + offset;
}

//
template<typename VectorType>
void invertData( VectorType& data, double eps = 1.e-12 ){
    for( int i = 0; i < data.size(); i++ ){
        if( std::abs(data[i]) < eps )
            throw BasicException("invertData: value is numerically zero!");
        data[i] = 1. / data[i];
    }
}

//
template<typename VectorType>
bool checkForNANsAndINFs( const VectorType &Arg ) {
    for ( int i = 0; i < Arg.size(); i++ )
      if ( std::isinf( Arg[i] ) || std::isnan( Arg[i] ))
        return true;
    return false;
} 

//==========================================================================================================
// MASK STUFF
//==========================================================================================================

// insert local mask (i.e. node indices of single shape to be fixed) into global mask 
// The global mask is related to the path of K shapes, each having N dofs (K = numShapes, N = numDofs)
void fillPathMask( int numShapes, int numDofs, const std::vector<int>& localMask, std::vector<int>& globalMask ){
    globalMask.clear();
    for( int i = 0; i < numShapes; i++ )
        for( uint j = 0; j < localMask.size(); j++ )
            globalMask.push_back(i*numDofs + localMask[j]);
}


// extend boundary mask to (x,y,z)
void extendBoundaryMask(int numVertices, std::vector<int>& mask){
    int oldSize = mask.size();
    mask.resize( 3*oldSize );
    for( int j = 1; j < 3; j++ )
      for( int i = 0; i < oldSize; i++ )      
        mask[j*oldSize + i ] = j*numVertices + mask[i];
}


// For each integer i in the mask, the i-th row-column of mat is set to e_i
template<typename MatrixType>
void applyMaskToSymmetricMatrix( const std::vector<int> &mask, MatrixType &mat ) {

  // run over mask entries
  for ( uint k = 0; k < mask.size(); ++k ) {
    int majorIdx = mask[k];

    // mask matrix
    if ( !(majorIdx < mat.outerSize()))
      throw BasicException( "applyMaskToSymmetricMatrix:: wrong mat index!" );

    //insert diagonal entry
    mat.coeffRef( majorIdx, majorIdx ) = 1.;

    std::vector<int> nzMinorIdx;
    // set column (resp. row) to e_i, where i = majorIdx
    // save row (resp. col) indices with nonzero entries
    for ( typename MatrixType::InnerIterator it( mat, majorIdx ); it; ++it ) {
      nzMinorIdx.push_back( it.index());
      it.valueRef() = (it.row() == it.col()) ? 1. : 0.;
    }

    // set row (resp. col) to e_i, where i = majorIdx
    for ( uint i = 0; i < nzMinorIdx.size(); i++ ) {
      typename MatrixType::InnerIterator it( mat, nzMinorIdx[i] );
      while ( it.index() < majorIdx ) ++it;
      if ( it.index() != majorIdx )
        throw BasicException( "applyMaskToSymmetricMatrix(): matrix is not symmetric!" );
      it.valueRef() = (it.row() == it.col()) ? 1. : 0.;
    }
  }
}

//
template<typename MatrixType>
void applyMaskToMajor( const std::vector<int> &mask, MatrixType &mat, bool setDiagonalOne = true ) {
  // run over mask entries
  for ( uint k = 0; k < mask.size(); ++k ) {
    int majorIdx = mask[k];
    // mask matrix
    if ( !(majorIdx < mat.outerSize()) ) {
        std::cerr << "major index = " << majorIdx << ", size major = " << mat.outerSize() << std::endl;
        throw BasicException("applyMaskToRow(): major index out of bounds!");
    }
    //insert diagonal entry
    if (setDiagonalOne)
      mat.coeffRef( majorIdx, majorIdx ) = 1.;
    // set column (resp. row) to e_i, where i = majorIdx
    for ( typename MatrixType::InnerIterator it( mat, majorIdx ); it; ++it )
      it.valueRef() = (it.row() == it.col() && setDiagonalOne) ? 1. : 0.;
  }
}

//
template<typename MatrixType>
void applyMaskToMinor( const std::vector<int> &mask, MatrixType &mat, bool setDiagonalOne = true ) {
  for ( int mIdx : mask )
    mat.coeffRef( mIdx, mIdx ) = 1.;

  for ( int k = 0; k < mat.outerSize(); ++k ) {
    for ( typename MatrixType::InnerIterator it( mat, k ); it; ++it ) {
      if ( std::find( mask.begin(), mask.end(), it.index()) != mask.end())
        it.valueRef() = (it.row() == it.col()&& setDiagonalOne) ? 1. : 0.;
    }
  }
}

//
template<typename MatrixType>
void applyMaskToRow( const std::vector<int> &mask, MatrixType &mat, bool setDiagonalOne = true ) {
  if ( MatrixType::IsRowMajor )
    applyMaskToMajor<MatrixType>( mask, mat, setDiagonalOne );
  else
    applyMaskToMinor<MatrixType>( mask, mat, setDiagonalOne );
}

//
template<typename MatrixType>
void applyMaskToColumn( const std::vector<int> &mask, MatrixType &mat, bool setDiagonalOne = true ) {
  if ( MatrixType::IsRowMajor )
    applyMaskToMinor<MatrixType>( mask, mat, setDiagonalOne );
  else
    applyMaskToMajor<MatrixType>( mask, mat, setDiagonalOne );
}

// For each integer i in the mask, the i-th row-column of mat is set to e_i
template<typename MatrixType>
void applyMaskToMatrix( const std::vector<int>& mask, MatrixType& mat, bool setDiagonalOne = true  ) {

  if ( setDiagonalOne )
    for ( int mIdx : mask )
      mat.coeffRef( mIdx, mIdx ) = 1.;

  for ( int k = 0; k < mat.outerSize(); ++k ) {
    for ( typename MatrixType::InnerIterator it( mat, k ); it; ++it ) {
      if ( std::find_if( mask.begin(), mask.end(),
                         [it]( const int &s ) { return (s == it.row() || s == it.col()); } ) != mask.end())
        it.valueRef() = (it.row() == it.col() && setDiagonalOne ) ? 1. : 0.;
    }
  }
}

// For each integer i in the mask, the i-th entry of vec to 0
template<typename VectorType>
void applyMaskToVector( const std::vector<int>& mask, VectorType& vec ) {
    // run over mask entries
    for (uint k= 0; k<mask.size(); ++k){
      int majorIdx = mask[k];      
      // mask vector
      if ( !(majorIdx < vec.size()) ){
          std::cerr << "index = " << majorIdx << " whereas size = " << vec.size() << std::endl;
          throw BasicException("applyMaskToVector() index out of bounds!");
      }
      vec[majorIdx] = 0.;
    }
}

// For each integer i in the mask, the i-th row-column of mat is set to e_i, and the i-th entry of vec to 0
template<typename MatrixType, typename VectorType>
void applyMaskToSymmetricMatrixAndVector( const std::vector<int>& mask, MatrixType& mat, VectorType& vec ) {
  applyMaskToVector<VectorType>( mask, vec );
  applyMaskToSymmetricMatrix<MatrixType>(mask, mat );
}

#endif // AUXILIARY_HH defined

