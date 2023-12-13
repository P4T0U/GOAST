// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2023 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include <goast/Core.h>
#include "LinearOperator.h"

template<typename ConfiguratorType=DefaultConfigurator>
class ObjectiveFunctional {
protected:
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using MatrixType = typename ConfiguratorType::SparseMatrixType;

public:
  ObjectiveFunctional() = default;

  virtual ~ObjectiveFunctional() = default;

  // Properties
  virtual int TargetDimension() const {
    throw std::logic_error( "ObjectiveFunctional::getTargetDimension(): Not implemented in derived class" );
  }

  // Basic evaluations
  virtual void evaluate( const VectorType &Point, RealType &Value ) const = 0;

  virtual void evaluateGradient( const VectorType &Point, VectorType &Gradient ) const {
    throw std::logic_error( "ObjectiveFunctional::evaluateGradient: Not implemented in derived class" );
  }

  virtual void evaluateHessian( const VectorType &Point, MatrixType &Hess ) const {
    throw std::logic_error( "ObjectiveFunctional::evaluateHessian: Not implemented in derived class" );
  }

  virtual std::unique_ptr<LinearOperator<ConfiguratorType>> HessOp( const VectorType &Point ) const {
    throw std::logic_error( "ObjectiveFunctional::HessOp: Not implemented in derived class" );
  }

  // Derived evaluations
  virtual RealType operator()( const VectorType &Point ) const {
    RealType Value;
    evaluate( Point, Value );
    return Value;
  }

  virtual VectorType grad( const VectorType &Point ) const {
    VectorType Gradient;
    evaluateGradient( Point, Gradient );
    return Gradient;
  }

  virtual MatrixType Hess( const VectorType &Point ) const {
    MatrixType HessMat;
    evaluateHessian( Point, HessMat );
    return HessMat;
  }
};

// Wrappers for old interface
template<typename ConfiguratorType=DefaultConfigurator>
class ObjectiveWrapper : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
  using RealType = typename ConfiguratorType::RealType;
  using VectorType = typename ConfiguratorType::VectorType;
  using MatrixType = typename ConfiguratorType::SparseMatrixType;

  const ObjectiveFunctional<ConfiguratorType> &m_F;

public:
  ObjectiveWrapper( const ObjectiveFunctional<ConfiguratorType> &F ) : m_F( F ) {}

  void apply( const VectorType &Arg, RealType &Dest ) const override {
    m_F.evaluate( Arg, Dest );
  }
};

template<typename ConfiguratorType=DefaultConfigurator>
class ObjectiveGradientWrapper
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
  using RealType = typename ConfiguratorType::RealType;
  using VectorType = typename ConfiguratorType::VectorType;
  using MatrixType = typename ConfiguratorType::SparseMatrixType;

  const ObjectiveFunctional<ConfiguratorType> &m_F;

public:
  ObjectiveGradientWrapper( const ObjectiveFunctional<ConfiguratorType> &F ) : m_F( F ) {}

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    m_F.evaluateGradient( Arg, Dest );
  }
};

template<typename ConfiguratorType=DefaultConfigurator>
class ObjectiveHessianWrapper
        : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
  using RealType = typename ConfiguratorType::RealType;
  using VectorType = typename ConfiguratorType::VectorType;
  using MatrixType = typename ConfiguratorType::SparseMatrixType;

  const ObjectiveFunctional<ConfiguratorType> &m_F;

public:
  ObjectiveHessianWrapper( const ObjectiveFunctional<ConfiguratorType> &F ) : m_F( F ) {}

  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    m_F.evaluateHessian( Arg, Dest );
  }
};

template<typename ConfiguratorType=DefaultConfigurator>
class ObjectiveHessianOperatorWrapper : public MapToLinOp<ConfiguratorType> {
protected:
  using VectorType = typename ConfiguratorType::VectorType;

  const ObjectiveFunctional<ConfiguratorType> &m_F;

public:
  explicit ObjectiveHessianOperatorWrapper( const ObjectiveFunctional<ConfiguratorType> &F ) : m_F( F ) {}

  std::unique_ptr<LinearOperator<ConfiguratorType>> operator()( const VectorType &Point ) const override {
    return m_F.HessOp( Point );
  }
};
