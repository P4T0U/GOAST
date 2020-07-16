// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Definitions and helpers for constraint functionals
 * \author Sassen
 */

#ifndef OPTIMIZATION_CONSTRAINTS_H
#define OPTIMIZATION_CONSTRAINTS_H

#include <goast/Core/BaseOpInterface.h>
#include <goast/Core/Auxiliary.h>

// Definitions/Shorthands for constraint functions
template<typename ConfiguratorType>
using ConstraintsOp = BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>;

template<typename ConfiguratorType>
using ConstraintsGradient = BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>;

template<typename ConfiguratorType>
using ConstraintsHessian = BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::TensorType>;


template<typename ConfiguratorType>
class PathConstraintOp : public ConstraintsOp<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  const ConstraintsOp<ConfiguratorType> &_constraintOp;

  const int _constraintDim;
  const int _numShapes;

public:
  PathConstraintOp( const ConstraintsOp<ConfiguratorType> &constraintOp, const int numShapes ) :
          _constraintOp( constraintOp ), _constraintDim( _constraintOp.getTargetDimension()), _numShapes( numShapes ) {}

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Arg.size() % _numShapes != 0 )
      throw std::length_error( "PathConstraintOp::apply(): Arg not divisible!" );

    const int numDOF = Arg.size() / _numShapes;

    Dest.resize( _numShapes * _constraintDim );
    Dest.setZero();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numShapes; i++ ) {
      VectorType Geom = Arg.segment( i * numDOF, numDOF );

      VectorType intPart;
      _constraintOp.apply( Geom, intPart );

      if ( intPart.size() != _constraintDim )
        throw std::length_error( "PathConstraintOp::apply(): Wrong size of constraints!" );

      Dest.segment( i * _constraintDim, _constraintDim ) = intPart;
    }
  }

  int getTargetDimension() const override {
    return _numShapes * _constraintDim;
  }
};

template<typename ConfiguratorType>
class PathConstraintGradient : public ConstraintsGradient<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  const ConstraintsGradient<ConfiguratorType> &_constraintGrad;

  const int _constraintDim;
  const int _numShapes;
public:
  PathConstraintGradient( const ConstraintsGradient<ConfiguratorType> &constraintGrad, const int numShapes ) :
          _constraintGrad( constraintGrad ), _constraintDim( _constraintGrad.getTargetDimension()),
          _numShapes( numShapes ) {}

  void apply( const VectorType &Arg, SparseMatrixType &Dest ) const override {
    if ( Arg.size() % _numShapes != 0 )
      throw std::length_error( "PathConstraintGradient::apply(): Arg not divisible!" );

    const int numDOF = Arg.size() / _numShapes;

    if ((Dest.cols() != _numShapes * numDOF) || (Dest.rows() != _numShapes * _constraintDim))
      Dest.resize( _numShapes * _constraintDim, _numShapes * numDOF );

    Dest.setZero();

    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList ) const override {
    const int numDOF = Arg.size() / _numShapes;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numShapes; i++ ) {
      VectorType Geom = Arg.segment( i * numDOF, numDOF );

      TripletListType intTriplets;
      _constraintGrad.pushTriplets( Geom, intTriplets, 1., i * _constraintDim, i * numDOF );

#ifdef _OPENMP
#pragma omp critical
#endif
      tripletList.insert( tripletList.end(), intTriplets.begin(), intTriplets.end());
    }
  }

  int getTargetDimension() const override {
    return _numShapes * _constraintDim;
  }

  int getNNZ() const override {
    return _numShapes * _constraintGrad.getNNZ();
  }
};

template<typename ConfiguratorType>
class PathConstraintHessian : public ConstraintsHessian<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType TensorType;
  typedef std::vector<TripletType> TripletListType;


  const ConstraintsHessian<ConfiguratorType> &_constraintHess;


  const int _constraintDim;
  const int _numShapes;
public:
  PathConstraintHessian( const ConstraintsHessian<ConfiguratorType> &constraintHess, const int numShapes ) :
          _constraintHess( constraintHess ), _constraintDim( _constraintHess.getTargetDimension()),
          _numShapes( numShapes ) {}

  void apply( const VectorType &Arg, std::vector<SparseMatrixType> &Dest ) const {
    if ( Arg.size() % _numShapes != 0 )
      throw std::length_error( "PathConstraintHessian::apply(): Arg not divisible!" );

    const int numDOF = Arg.size() / _numShapes;

    Dest.reserve( _numShapes * _constraintDim );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numShapes; i++ ) {
      VectorType Geom = Arg.segment( i * numDOF, numDOF );

      TensorType constraintHessians;
      _constraintHess.apply( Geom, constraintHessians, _numShapes * numDOF, i * numDOF );

#ifdef _OPENMP
#pragma omp critical
#endif
      Dest.insert( Dest.end(), constraintHessians.begin(), constraintHessians.end());
    }
  }

  void apply( const VectorType &Arg, TensorType &Dest ) const override {
    if ( Arg.size() % _numShapes != 0 )
      throw std::length_error( "PathConstraintOp::apply(): Arg not divisible!" );

    const int numDOF = Arg.size() / _numShapes;

    Dest.resize( _numShapes * _constraintDim, numDOF, numDOF );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numShapes; i++ ) {
      VectorType Geom = Arg.segment( i * numDOF, numDOF );

      TensorType constraintHessians;
      _constraintHess.apply( Geom, constraintHessians, _numShapes * numDOF, i * numDOF );

#ifdef _OPENMP
#pragma omp critical
#endif
      for ( int j = 0; j < _constraintDim; j++ )
        Dest[i * _constraintDim + j] = constraintHessians[j];
    }
  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest, const VectorType &Lambda ) const override {
    if ( Arg.size() % _numShapes != 0 )
      throw std::length_error( "PathConstraintHessian::apply(): Arg not divisible!" );

    const int numDOF = Arg.size() / _numShapes;

//    Dest.reserve( _numShapes * _constraintDim );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numShapes; i++ ) {
      VectorType Geom = Arg.segment( i * numDOF, numDOF );
      VectorType Lam = Lambda.segment( i * _constraintDim, _constraintDim );

      TripletListType constraintHessians;
      _constraintHess.pushTriplets( Geom, constraintHessians, Lam, i * numDOF );

#ifdef _OPENMP
#pragma omp critical
#endif
      Dest.insert( Dest.end(), constraintHessians.begin(), constraintHessians.end());
    }
  }

  void setTriplets( const VectorType &Arg, std::vector<TripletListType> &Dest, int hessOffset = 0 ) const override {
    if ( Arg.size() % _numShapes != 0 )
      throw std::length_error( "PathConstraintHessian::apply(): Arg not divisible!" );

    const int numDOF = Arg.size() / _numShapes;

//    Dest.reserve( _numShapes * _constraintDim );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < _numShapes; i++ ) {
      VectorType Geom = Arg.segment( i * numDOF, numDOF );

      std::vector<TripletListType> constraintTriplets;
      _constraintHess.setTriplets( Geom, constraintTriplets, i * numDOF + hessOffset );

#ifdef _OPENMP
#pragma omp critical
#endif
      Dest.insert( Dest.end(), constraintTriplets.begin(), constraintTriplets.end());
    }
  }

  int getTargetDimension() const override {
    return _numShapes * _constraintDim;
  }

  int getNNZ() const override {
    return _numShapes * _constraintHess.getNNZ();
  }
};

template<typename ConfiguratorType>
class CombinedConstraintsOp : public ConstraintsOp<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

  std::vector<const ConstraintsOp<ConfiguratorType> *> _Ops;

  std::vector<int> _constraintsDim;

  int _totalNumConstraints;
public:
  template<class... Items>
  CombinedConstraintsOp( Items const &... constraintOps ) : _totalNumConstraints( 0 ) {
    append_to_vector( constraintOps... );
  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest.resize( _totalNumConstraints );
    Dest.setZero();

    int currentDim = 0;

    for ( int i = 0; i < _Ops.size(); i++ ) {
      Dest.segment( currentDim, _constraintsDim[i] ) = (*_Ops[i])( Arg );
      currentDim += _constraintsDim[i];
    }
  }

  int getTargetDimension() const override {
    return _totalNumConstraints;
  }

private:
  void append_to_vector( const ConstraintsOp<ConfiguratorType> &elem ) {
    _Ops.push_back( &elem );
    _constraintsDim.push_back( elem.getTargetDimension());
    _totalNumConstraints += elem.getTargetDimension();
  };

  template<typename ...T1toN>
  void append_to_vector( const ConstraintsOp<ConfiguratorType> &elem, T1toN const &... elems ) {
    _Ops.push_back( &elem );
    _constraintsDim.push_back( elem.getTargetDimension());
    _totalNumConstraints += elem.getTargetDimension();
    append_to_vector( elems... );
  };
};

template<typename ConfiguratorType>
class CombinedConstraintsGradient : public ConstraintsGradient<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;

  std::vector<const ConstraintsGradient<ConfiguratorType> *> _Ops;

  std::vector<int> _constraintsDim;

  int _totalNumConstraints;
  int _totalNumNNZ;

public:
  template<class... Items>
  CombinedConstraintsGradient( Items const &... constraintOps ) : _totalNumNNZ( 0 ), _totalNumConstraints( 0 ) {
    append_to_vector( constraintOps... );
  }

  void apply( const VectorType &Arg, SparseMatrixType &Dest ) const override {
    if ((Dest.cols() != Arg.size()) || (Dest.rows() != _totalNumConstraints))
      Dest.resize( _totalNumConstraints, Arg.size() );

    Dest.setZero();

    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &tripletList ) const override {
    tripletList.reserve( tripletList.size() + _totalNumNNZ );

    int currentDim = 0;

    for ( int i = 0; i < _Ops.size(); i++ ) {
      _Ops[i]->pushTriplets( Arg, tripletList, 1., currentDim, 0 );

      currentDim += _constraintsDim[i];
    }
  }

  int getTargetDimension() const override {
    return _totalNumConstraints;
  }

  int getNNZ() const override {
    return _totalNumNNZ;
  }

private:
  void append_to_vector( const ConstraintsGradient<ConfiguratorType> &elem ) {
    _Ops.push_back( &elem );
    _constraintsDim.push_back( elem.getTargetDimension());
    _totalNumConstraints += elem.getTargetDimension();
    _totalNumNNZ += elem.getNNZ();
  };

  template<typename ...T1toN>
  void append_to_vector( const ConstraintsGradient<ConfiguratorType> &elem, T1toN const &... elems ) {
    _Ops.push_back( &elem );
    _constraintsDim.push_back( elem.getTargetDimension());
    _totalNumConstraints += elem.getTargetDimension();
    _totalNumNNZ += elem.getNNZ();
    append_to_vector( elems... );
  };
};

template<typename ConfiguratorType>
class CombinedConstraintsHessian : public ConstraintsHessian<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType TensorType;
  typedef std::vector<TripletType> TripletListType;


  std::vector<const ConstraintsHessian<ConfiguratorType> *> _Ops;

  std::vector<int> _constraintsDim;

  int _totalNumConstraints;
  int _totalNumNNZ;

public:
  template<class... Items>
  CombinedConstraintsHessian( Items const &... constraintOps ) : _totalNumNNZ( 0 ), _totalNumConstraints( 0 ) {
    append_to_vector( constraintOps... );
  }

  void apply( const VectorType &Arg, std::vector<SparseMatrixType> &Dest ) const {
    Dest.reserve( _totalNumConstraints );

    int currentDim = 0;

    for ( int i = 0; i < _Ops.size(); i++ ) {
      TensorType constraintHessians;
      _Ops[i]->apply( Arg, constraintHessians );

      Dest.insert( Dest.end(), constraintHessians.begin(), constraintHessians.end());

      currentDim += _constraintsDim[i];
    }
  }

  void apply( const VectorType &Arg, TensorType &Dest ) const override {
    Dest.resize( _totalNumConstraints, Arg.size(), Arg.size() );
    int currentDim = 0;

    for ( int i = 0; i < _Ops.size(); i++ ) {
      TensorType constraintHessians;
      _Ops[i]->apply( Arg, constraintHessians );

      for ( int j = 0; j < _constraintsDim[i]; j++ )
        Dest[currentDim + j] = constraintHessians[j];

      currentDim += _constraintsDim[i];
    }
  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest, const VectorType &Lambda ) const override {
    int currentDim = 0;

    for ( int i = 0; i < _Ops.size(); i++ ) {
      VectorType Lam = Lambda.segment( currentDim, _constraintsDim[i] );

      _Ops[i]->pushTriplets( Arg, Dest, Lam );

      currentDim += _constraintsDim[i];
    }
  }

  void setTriplets( const VectorType &Arg, std::vector<TripletListType> &Dest, int hessOffset = 0 ) const override {
    int currentDim = 0;

    for ( int i = 0; i < _Ops.size(); i++ ) {
      std::vector<TripletListType> constraintTriplets;
      _Ops[i]->setTriplets( Arg, constraintTriplets, hessOffset );

      Dest.insert( Dest.end(), constraintTriplets.begin(), constraintTriplets.end());
      currentDim += _constraintsDim[i];
    }
  }

  int getTargetDimension() const override {
    return _totalNumConstraints;
  }

  int getNNZ() const override {
    return _totalNumNNZ;
  }

private:
  void append_to_vector( const ConstraintsHessian<ConfiguratorType> &elem ) {
    _Ops.push_back( &elem );
    _constraintsDim.push_back( elem.getTargetDimension());
    _totalNumConstraints += elem.getTargetDimension();
    _totalNumNNZ += elem.getNNZ();
  };

  template<typename ...T1toN>
  void append_to_vector( const ConstraintsHessian<ConfiguratorType> &elem, T1toN const &... elems ) {
    _Ops.push_back( &elem );
    _constraintsDim.push_back( elem.getTargetDimension());
    _totalNumConstraints += elem.getTargetDimension();
    _totalNumNNZ += elem.getNNZ();
    append_to_vector( elems... );
  };
};


#endif //OPTIMIZATION_CONSTRAINTS_H
