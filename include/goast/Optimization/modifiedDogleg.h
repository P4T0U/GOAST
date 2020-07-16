// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for the modified dogleg method used in Byrd-Omojokun SQP
 * \author Sassen
 */

#ifndef OPTIMIZATION_MODIFIEDDOGLEG_H
#define OPTIMIZATION_MODIFIEDDOGLEG_H

#include "optInterface.h"

template<typename ConfiguratorType, typename LinearSolverType>
class ModifiedDoglegSolver : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const MatrixType &_A;
  const VectorType &_b;
  const RealType &_radius;

  const VectorType *_lowerBounds;
  const VectorType *_upperBounds;

  // Linear solver for row-space projection
  const LinearSolverType &_AAT_solver;

  bool _quiet = true;

public:
  ModifiedDoglegSolver( const MatrixType &A, const VectorType &b, const RealType &radius,
                        const LinearSolverType &AAT_solver ) : _A( A ), _b( b ), _radius( radius ),
                                                               _AAT_solver( AAT_solver ), _lowerBounds( nullptr ),
                                                               _upperBounds( nullptr ) {}

  void setVariableBounds( const VectorType &lowerBounds, const VectorType &upperBounds ) override {
    _lowerBounds = &lowerBounds;
    _upperBounds = &upperBounds;
  }

  void solve( const VectorType &/*startingPoint*/, VectorType &solution ) const override {
    // Newton point
    VectorType t = -_AAT_solver.solve( _b );
    VectorType pB( _A.transpose() * t );

    // Check if newton point is inside trust region and box constraints
    if ( pB.norm() <= _radius ) {
      if ( _lowerBounds && _upperBounds ) {
        if ( insideBox( pB, *_lowerBounds, *_upperBounds )) {
          solution = pB;
          return;
        }
      }
      else {
        solution = pB;
        return;
      }
    }

    // Linear part of the quadratic functional 1/2 * || A x + b ||^2
    VectorType g = _A.transpose() * _b;
    RealType alpha = -g.squaredNorm() / (_A * g).squaredNorm();

    // Compute Cauchy point
    VectorType pU = alpha * g;

    // Origin vector
    VectorType origin = VectorType::Zero( pU.size());

    RealType alpha_1, alpha_2;
    bool intersect;

    VectorType x1, x2;

    // Check along the line segments origin -> pU and pU -> pB
    VectorType d = pB - pU;
    if ( _lowerBounds && _upperBounds )
      std::tie( alpha_1, alpha_2, intersect ) = lineBoxBallIntersections( pU, d, _radius, *_lowerBounds,
                                                                          *_upperBounds );
    else
      std::tie( alpha_1, alpha_2, intersect ) = lineBallIntersection( pU, d, _radius, 0., 1. );

    if ( intersect ) {
      x1 = pU + alpha_2 * d;
    }
    else {
      if ( _lowerBounds && _upperBounds )
        std::tie( alpha_1, alpha_2, intersect ) = lineBoxBallIntersections( origin, pU, _radius, *_lowerBounds,
                                                                            *_upperBounds );
      else
        std::tie( alpha_1, alpha_2, intersect ) = lineBallIntersection( origin, pU, _radius, 0., 1. );

      x1 = alpha_2 * pU;
    }

    // Check along the line segment origin -> pB
    if ( _lowerBounds && _upperBounds )
      std::tie( alpha_1, alpha_2, intersect ) = lineBoxBallIntersections( origin, pB, _radius, *_lowerBounds,
                                                                          *_upperBounds );
    else
      std::tie( alpha_1, alpha_2, intersect ) = lineBallIntersection( origin, pB, _radius, 0., 1. );

    x2 = alpha_2 * pB;

    if ((_A * x1 + _b).norm() < (_A * x2 + _b).norm())
      solution = x1;
    else
      solution = x2;
  }

  void setParameter( const std::string &name, std::string value ) override {
    throw std::runtime_error( "ModifiedDoglegSolver::setParameter(): This class has no string parameters.." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    throw std::runtime_error( "ModifiedDoglegSolver::setParameter(): This class has no real parameters.." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "print_level" )
      _quiet = (value < 5);
    else
      throw std::runtime_error( "ModifiedDoglegSolver::setParameter(): Unknown parameter '" + name + "'." );
  }
};

#endif //OPTIMIZATION_MODIFIEDDOGLEG_H
