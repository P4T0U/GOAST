// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Header for Steihaug's CG method for solving the trust-region subproblem
 * \author Sassen
 */

#ifndef OPTIMIZATION_STEIHAUGCG_H
#define OPTIMIZATION_STEIHAUGCG_H

#include <goast/Core/Auxiliary.h>
#include "optUtils.h"
#include "optInterface.h"

/// \todo Boundary Mask
template<typename ConfiguratorType>
class SteihaugCGMethod : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::TensorType TensorType;

  const MatrixType &_H;
  const VectorType &_c;

  const RealType _radius;
  RealType _tolerance;
  int _maxIterations;

  bool _quiet;

public:
  SteihaugCGMethod( const MatrixType &H, const VectorType &c, const RealType radius,
                    const RealType tolerance, const int maxIterations = 1000, bool quiet = false )
          : _H( H ), _c( c ), _tolerance( tolerance ), _radius( radius ),
            _maxIterations( maxIterations ), _quiet ( quiet ) {}

  void solve( const VectorType & /*start*/, VectorType &s ) const override {
    solve(s);
  }

  void solve( VectorType &z_k ) const {
    this->status.Iteration = 0;

    z_k.resize( _c.size());
    z_k.setZero();

    VectorType r_k( _c );
    VectorType d_k( -r_k );

    this->status.Residual = r_k.norm();

    // Stop uf 0 is already a good enough solution
    if ( r_k.norm() < _tolerance ) {
      std::cout << "stopped with inital residual " << r_k.norm() << " vs. " << _tolerance << std::endl;
      return;
    }

    // Intermediate and temporary quantities
    VectorType Hd; // H * d_k
    RealType alpha_k, beta_k, temp;
    RealType r_sqNorm = r_k.squaredNorm(); // r_k^T r_k
    VectorType tmpVector;

    // Helper variables for line / trust-region intersection
    RealType tau_0, tau_1;
    bool intersected;

    for ( int k = 0; k < _maxIterations; k++ ) {
      this->status.Iteration++;
      Hd.noalias() = _H * d_k;

      // Step 1: Check if current search direction is one of nonpositive curvature and if yes return intersection with
      // boundary of trust region
      RealType curvature = d_k.transpose() * Hd;
      if ( curvature <= 0 ) {
        // TODO: replace by proper lineSphereIntersection
        std::tie( tau_0, tau_1, intersected ) = lineSphereIntersection<RealType, VectorType>( z_k, d_k, _radius );
        if ( !intersected )
          throw std::runtime_error( "SteihaugCGMethod: Negative curvature line does not intersect with ball!" );

        z_k += tau_1 * d_k;

        this->status.reasonOfTermination = 1;

        return;
      }

      // Step 2: Compute new iterate
      alpha_k = r_sqNorm / curvature;
      tmpVector = z_k + alpha_k * d_k;

      // Step 2a: If new iterate is outside of trust region then return intersection with boundary of trust region
      if ( tmpVector.norm() >= _radius ) {
        std::tie( tau_0, tau_1, intersected ) = lineSphereIntersection<RealType, VectorType>( z_k, d_k, _radius );
        if ( !intersected || tau_1 < 0 )
          throw std::runtime_error( "SteihaugCGMethod: Search direction does not intersect with ball!" );
        z_k += tau_1 * d_k;

        this->status.reasonOfTermination = 2;

        return;
      }
      z_k = tmpVector; // Accept iterate

      // Step 3: Update residual and check for convergence
      r_k += alpha_k * Hd; // r_{k+1} = r_k + alpha_k * H * d_k
      this->status.Residual = r_k.norm();
      if ( r_k.norm() < _tolerance ) {
        this->status.reasonOfTermination = 0;
        return;
      }

      // Step 4: Compute new search direction
      temp = r_sqNorm;
      r_sqNorm = r_k.squaredNorm();
      beta_k = r_sqNorm / temp; // beta_{k+1} = r_{k+1}^T r_{k+1} / r_k^T r_k
      d_k *= beta_k;
      d_k -= r_k; // d_{k+1} = -r_{k+1} + beta_{k+1} * d_k

      if (!_quiet)
        std::cout << "SCG -- Iter " << k << ": ||r_k|| = " << r_k.norm() << std::endl;
    }
    this->status.reasonOfTermination = -1;
  }

  void setParameter( const std::string &name, std::string value ) override {
    throw std::runtime_error( "SteihaugCGMethod::setParameter(): This class has no string parameters.." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "tolerance" )
      _tolerance = value;
    else
      throw std::runtime_error( "SteihaugCGMethod::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "print_level" )
      _quiet = (value < 5);
    else
      throw std::runtime_error( "SteihaugCGMethod::setParameter(): Unknown parameter '" + name + "'." );
  }
};

#endif //OPTIMIZATION_STEIHAUGCG_H
