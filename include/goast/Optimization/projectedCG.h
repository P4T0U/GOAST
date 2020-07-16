// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for the projected CG method used in Byrd-Omojokun SQP
 * \author Sassen
 */

#ifndef OPTIMIZATION_PROJECTEDCG_H
#define OPTIMIZATION_PROJECTEDCG_H

#include "optInterface.h"

template<typename ConfiguratorType, typename LinearSolverType>
class ProjectedCGSolver : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef std::vector<TripletType> TripletListType;
  typedef typename ConfiguratorType::TensorType TensorType;

  // Epsilon for safe-guarding computation of reduction against numerically instabilities
  const RealType eps_red = 100 * std::numeric_limits<RealType>::epsilon();

  // Problem
  const MatrixType &_H;
  const VectorType &_c;
  const MatrixType &_A;
  const VectorType &_b;

  const RealType &_radius;
  const VectorType *_lowerBounds;
  const VectorType *_upperBounds;

  // Parameters
  RealType _tolerance;
  int _maxIterations;
  int _maxInfeasibleIterations = 10;

  bool _quiet;

  // Linear solver for projections
  const LinearSolverType &_AAT_solver;


public:
  ProjectedCGSolver( const MatrixType &H, const VectorType &c,
                     const MatrixType &A, const VectorType &b,
                     const RealType &radius,
                     const LinearSolverType &AAT_solver,
                     RealType tolerance = 1e-8, int maxIterations = 1000,
                     bool quiet = true ) : _H( H ), _c( c ), _A( A ), _b( b ), _radius( radius ),
                                           _tolerance( tolerance ), _maxIterations( maxIterations ),
                                           _AAT_solver( AAT_solver ), _quiet( quiet ), _lowerBounds( nullptr ),
                                           _upperBounds( nullptr ) {}


  void setVariableBounds( const VectorType &lowerBounds, const VectorType &upperBounds ) override {
    _lowerBounds = &lowerBounds;
    _upperBounds = &upperBounds;
  }

  void solve( VectorType &z_k ) const {
    // Initial point z_k satisfying Az_k + b = 0
    VectorType t = -_AAT_solver.solve( _b );
    VectorType z_0( _A.transpose() * t );

    solve( z_0, z_k );
  }

  void solve( const VectorType &z_0, VectorType &z_k ) const override {
    // Initial point z_k satisfying Az_k + b = 0
    z_k = z_0; //! \todo Check if z_0 is feasible or simply trust the user?

    // If z_k is outside the trust-region the problem does not have a solution

    RealType distToBorder = _radius - z_k.norm();
    if ( distToBorder < 0 )
      throw std::domain_error( "ProjectedCGSolver: System has no solution! Distance: " + std::to_string(distToBorder) );
    else if ( distToBorder < eps_red ) {
      std::cout << "Initial point on boundary" << std::endl;
      return;
    }

    // Initial values
    VectorType r_k = _H * z_k + _c;
    VectorType g_k = Projection( r_k );
    VectorType d_k = -g_k;

    RealType rT_g, alpha, beta, tmp, alpha_1, alpha_2;
    bool intersect;
    VectorType Hd( _H.rows());
    VectorType x_next( z_k.size());
    VectorType r_next( r_k.size());
    VectorType last_feasible_x( z_k.size());

    rT_g = g_k.squaredNorm(); // r^T g = g^T g because of orthonormality

//    std::cout << "rT_g: " << rT_g << std::endl;

    int numIterations = 0, feasible_counter = 0;

    this->status.reasonOfTermination = -1;
    this->status.Iteration = 0;

    for ( int iter = 1; iter <= _maxIterations; iter++ ) {
      this->status.Iteration++;

      Hd.noalias() = _H * d_k;

      // Step 1: Check if current search direction is one of nonpositive curvature and if yes return intersection with
      // boundary of trust region
      RealType curvature = d_k.transpose() * Hd;

//      std::cout << "Curvature: " << std::scientific << curvature << std::endl;
//      std::cout << "dknorm: " << d_k.norm() << std::endl;
//      std::cout << "Hd norm: " << Hd.norm() << std::endl;

      if ( curvature <= eps_red ) {
        if ( _radius == std::numeric_limits<RealType>::infinity())
          throw std::domain_error( "ProjectedCGSolver: Negative curvature case needs a trust-region to have a "
                                   "minimum." );
        else {
          if ( _lowerBounds && _upperBounds )
            std::tie( alpha_1, alpha_2, intersect ) = lineBoxBallIntersections( z_k, d_k, _radius, *_lowerBounds,
                                                                                *_upperBounds,
                                                                                -std::numeric_limits<RealType>::infinity(),
                                                                                std::numeric_limits<RealType>::infinity());
          else
            std::tie( alpha_1, alpha_2, intersect ) = lineSphereIntersection<RealType, VectorType>( z_k, d_k, _radius );

          if ( !intersect )
            throw std::runtime_error( "ProjectedCGSolver: Negative curvature line does not intersect with "
                                      "trust-region!" );

          z_k += alpha_2 * d_k;

          //! \todo Enforce box conditions due to rounding errors

          this->status.reasonOfTermination = 1;
          break;
        }
      }



      // Step 2: Compute new iterate
      alpha = rT_g / curvature;
      x_next = z_k + alpha * d_k;

      // Step 2a: If new iterate is outside of trust region then return intersection with boundary of trust region
      if ( x_next.norm() > _radius ) {
        if ( _lowerBounds && _upperBounds )
          std::tie( alpha_1, alpha_2, intersect ) = lineBoxBallIntersections( z_k, d_k, _radius, *_lowerBounds,
                                                                              *_upperBounds, 0., alpha );
        else
          std::tie( alpha_1, alpha_2, intersect ) = lineBallIntersection<RealType, VectorType>( z_k, d_k, _radius, 0., alpha );

//        std::cout << " !--! " << std::scientific << std::setprecision( 6 ) << _radius << " -- " << z_k.norm()
//                              << " --- " << alpha_1 << " -- " << alpha_2 << std::endl;

        if ( !intersect || alpha_2 < 0 )
          throw std::runtime_error( "ProjectedCGSolver: Search direction does not intersect with trust-region!" );

        z_k += alpha_2 * d_k;


        this->status.reasonOfTermination = 2;

        break;
      }

      // Step 2b: If we have box constraints and new iterate violates them, then project it onto the box
      if ( _lowerBounds && _upperBounds ) {
        if ( insideBox( x_next, *_lowerBounds, *_upperBounds ))
          feasible_counter = 0;
        else {
          feasible_counter += 1;

          std::tie( alpha_1, alpha_2, intersect ) = lineBoxBallIntersections( z_k, d_k, _radius, *_lowerBounds,
                                                                              *_upperBounds, 0., alpha );
          if ( intersect ) {
            last_feasible_x = z_k + alpha_2 * d_k;
            feasible_counter = 0;
          }
        }

        if ( feasible_counter > _maxInfeasibleIterations ) {
          std::cout << " -- PrCG -- Iter " << std::setw( 3 ) << iter << ": " << "Too many infeasible iterations: "
                    << feasible_counter << std::endl;
          break;
        }
      }

      // Step 4: Update residual, compute projection and check for convergence
      r_k += alpha * Hd;
      g_k = Projection( r_k );

      tmp = g_k.squaredNorm();

      this->status.Residual = tmp;
      if ( tmp < _tolerance ) {
        this->status.reasonOfTermination = 0;
        break;
      }

      // Step 5: Compute new search direction
      beta = tmp / rT_g;
      rT_g = tmp;

      d_k.noalias() = -g_k + beta * d_k;
      z_k.noalias() = x_next;
    }

    if ( _lowerBounds && _upperBounds )
      if ( !insideBox( z_k, *_lowerBounds, *_upperBounds ))
        z_k = last_feasible_x;

    r_k = _H * z_k + _c;
    g_k = Projection( r_k );
    this->status.Residual = g_k.squaredNorm();
  }

protected:
  VectorType Projection( const VectorType &x ) const {
    VectorType Ax( _A * x );
    VectorType sol( _AAT_solver.solve( Ax ) );
    return x - _A.transpose() * sol;
  }

public:
  void setParameter( const std::string &name, std::string value ) override {
    throw std::runtime_error( "ProjectedCGSolver::setParameter(): This class has no string parameters.." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "tolerance" )
      _tolerance = value;
    else
      throw std::runtime_error( "ProjectedCGSolver::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "maximum_infeasible_iterations" )
      _maxInfeasibleIterations = value;
    else if ( name == "print_level" )
      _quiet = (value < 5);
    else
      throw std::runtime_error( "ProjectedCGSolver::setParameter(): Unknown parameter '" + name + "'." );
  }

};

#endif //OPTIMIZATION_PROJECTEDCG_H
