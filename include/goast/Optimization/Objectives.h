// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Definitions and helpers for objective functionals
 * \author Sassen
 */

#ifndef OPTIMIZATION_OBJECTIVES_H
#define OPTIMIZATION_OBJECTIVES_H

#include <goast/Core/BaseOpInterface.h>
#include <goast/Core/Auxiliary.h>

// Definitions/Shorthands for objective functions
template<typename ConfiguratorType>
using ObjectiveOp = BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType>;

template<typename ConfiguratorType>
using ObjectiveGradient = BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType>;

template<typename ConfiguratorType>
using ObjectiveHessian = BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType>;

template<typename ConfiguratorType>
class AdditionOp : public ObjectiveOp<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

//  const ObjectiveOp<ConfiguratorType> &_constraintOp;
  std::vector<const ObjectiveOp<ConfiguratorType> *> _Ops;
  const VectorType &_Weights;

  const int _numOps;
//  const int _numShapes;

public:
  template<class... Items>
  AdditionOp( const VectorType &Weights, Items const &... constraintOps ) : _Weights( Weights ),
                                                                            _numOps( sizeof...( constraintOps )) {
    append_to_vector( _Ops, constraintOps... );

  }

  void apply( const VectorType &Arg, RealType &Dest ) const override {
    Dest = 0.;

    for ( int i = 0; i < _numOps; i++ )
      Dest += _Weights[i] * ( *_Ops[i] )( Arg );
  }

  int getTargetDimension() const override {
    return 1;
  }

private:
  void append_to_vector( std::vector<const ObjectiveOp<ConfiguratorType> *> &outputvector,
                         const ObjectiveOp<ConfiguratorType> &elem ) {
    outputvector.push_back( &elem );
  };

  template<typename ...T1toN>
  void append_to_vector( std::vector<const ObjectiveOp<ConfiguratorType> *> &outputvector,
                         const ObjectiveOp<ConfiguratorType> &elem, T1toN const &... elems ) {
    outputvector.push_back( &elem );
    append_to_vector( outputvector, elems... );
  };

};

template<typename ConfiguratorType>
class AdditionGradient : public ObjectiveGradient<ConfiguratorType> {
protected:
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::VecType VecType;
  typedef typename ConfiguratorType::MatType MatType;

//  const ObjectiveOp<ConfiguratorType> &_constraintOp;
  std::vector<const ObjectiveGradient<ConfiguratorType> *> _Ops;
  const VectorType &_Weights;

  const int _numOps;
//  const int _numShapes;

public:
  template<class... Items>
  AdditionGradient( const VectorType &Weights, Items const &... constraintOps ) : _Weights( Weights ),
                                                                              _numOps( sizeof...( constraintOps )) {
    append_to_vector( _Ops, constraintOps... );

  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Dest.size() != Arg.size())
      Dest.resize( Arg.size());
    Dest.setZero();

    for ( int i = 0; i < _numOps; i++ )
      Dest += _Weights[i] * ( *_Ops[i] )( Arg );
  }

  int getTargetDimension() const override {
    return 1;
  }

private:
  void append_to_vector( std::vector<const ObjectiveGradient<ConfiguratorType> *> &outputvector,
                         const ObjectiveGradient<ConfiguratorType> &elem ) {
    outputvector.push_back( &elem );
  };

  template<typename ...T1toN>
  void append_to_vector( std::vector<const ObjectiveGradient<ConfiguratorType> *> &outputvector,
                         const ObjectiveGradient<ConfiguratorType> &elem, T1toN const &... elems ) {
    outputvector.push_back( &elem );
    append_to_vector( outputvector, elems... );
  };

};

#endif //OPTIMIZATION_OBJECTIVES_H
