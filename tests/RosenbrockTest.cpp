// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Test of unconstrained optimization methods with the Rosenbrock function
 */
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <goast/Core.h>
#include <goast/Optimization.h>

// ========================== Definitions ========================== //
template <typename ConfiguratorType>
class RosenbrockEnergy : public ObjectiveOp<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
public:
  void apply( const VectorType &Arg, RealType &Dest ) const override {
    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
    RealType aux2 = 1 - Arg[0];
    Dest = (100 * aux1 * aux1 + aux2 * aux2);
  }
};

template <typename ConfiguratorType>
class RosenbrockGradient : public ObjectiveGradient<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
public:
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Dest.size() != Arg.size())
      Dest.resize( Arg.size());
    Dest.setZero();
    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
    RealType aux2 = 1 - Arg[0];
    Dest[0] -= (400 * aux1 * Arg[0] + 2 * aux2);
    Dest[1] += (200 * aux1);
  }
};

template <typename ConfiguratorType>
class RosenbrockHessian : public ObjectiveHessian<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

public:
  void apply( const VectorType &Arg, MatrixType &Dest ) const override {

    if ((Dest.rows() != 2) || (Dest.cols() != 2))
      Dest.resize( 2, 2 );
    Dest.setZero();

    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
    Dest.coeffRef( 0, 0 ) += (800 * Arg[0] * Arg[0] - 400 * aux1 + 2);
    Dest.coeffRef( 0, 1 ) -= (400 * Arg[0]);
    Dest.coeffRef( 1, 0 ) -= (400 * Arg[0]);
    Dest.coeffRef( 1, 1 ) += (200.);
  }
};

typedef DefaultConfigurator ConfiguratorType;

// ========================== Tests ========================== //
TEST_CASE( "Global minimum of Rosenbrock's banana function is found", "[Optimization][unconstrained]" ) {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;

  const RealType tolerance = 1e-8;

  RosenbrockEnergy<ConfiguratorType> E;
  RosenbrockGradient<ConfiguratorType> DE;
  RosenbrockHessian<ConfiguratorType> D2E;

  VectorType start( 2 ), solution( 2 ), minimum( 2 );
  start << -3, -4;
  minimum << 1, 1;

  SECTION( "GradientDescent" ) {
    GradientDescent<ConfiguratorType>( E, DE, 100000, 1e-9, ARMIJO, SUPERQUIET ).solve( start, solution );
    REQUIRE((minimum - solution).norm() < tolerance );
  }

  SECTION( "QuasiNewtonBFGS" ) {
    QuasiNewtonBFGS<ConfiguratorType>( E, DE, 1000, 1e-9, ARMIJO, 50, SUPERQUIET ).solve( start, solution );
    REQUIRE((minimum - solution).norm() < tolerance );
  }

  SECTION( "LineSearchNewton" ) {
    LineSearchNewton<ConfiguratorType>( E, DE, D2E, 1e-8, 100, SUPERQUIET ).solve( start, solution );
    REQUIRE((minimum - solution).norm() < tolerance );
  }

  SECTION( "TrustRegionNewton+SteihaugCG" ) {
    TrustRegionNewton<ConfiguratorType>( E, DE, D2E, 1., 100., 1e-8, 10000, 100, 0.25, true ).solve( start, solution );
    REQUIRE((minimum - solution).norm() < tolerance );
  }

  SECTION( "TrustRegionNewton+MoreSorensen" ) {
    TrustRegionNewton<ConfiguratorType> TRN( E, DE, D2E, 1., 100., 1e-8, 10000, 100, 0.25, true );

    TRN.setParameter( "subproblem_solver", "MoreSorensen" );
    TRN.solve( start, solution );

    REQUIRE((minimum - solution).norm() < tolerance );
  }
}