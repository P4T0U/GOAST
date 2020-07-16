// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Core.h>
#include <goast/Optimization/SteihaugCG.h>
#include <goast/Optimization/TrustRegionNewton.h>

//================================================================================
typedef DefaultConfigurator ConfiguratorType;

typedef typename ConfiguratorType::RealType RealType;

typedef typename ConfiguratorType::VectorType VectorType;
typedef typename ConfiguratorType::SparseMatrixType MatrixType;
typedef typename ConfiguratorType::FullMatrixType FullMatrixType;
typedef typename ConfiguratorType::TripletType TripletType;
typedef typename ConfiguratorType::VecType VecType;
typedef typename ConfiguratorType::MatType MatType;
typedef std::vector<TripletType> TripletListType;


//================================================================================
// f(x,y) = 100 * (y - x^2)^2 + (1-x)^2
class TestEnergy : public BaseOp<VectorType, RealType> {

public:
  void apply( const VectorType &Arg, RealType &Dest ) const {
    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
    RealType aux2 = 1 - Arg[0];
    Dest = 100 * aux1 * aux1 + aux2 * aux2;
  }
};

class TestGradient : public BaseOp<VectorType, VectorType> {

public:
  void apply( const VectorType &Arg, VectorType &Dest ) const {
    if ( Dest.size() != Arg.size())
      Dest.resize( Arg.size());
    Dest.setZero();
    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
    RealType aux2 = 1 - Arg[0];
    Dest[0] -= 400 * aux1 * Arg[0] + 2 * aux2;
    Dest[1] += 200 * aux1;
  }
};

class TestHessian : public BaseOp<VectorType, MatrixType> {

public:
  void apply( const VectorType &Arg, MatrixType &Dest ) const {

    if ((Dest.rows() != 2) || (Dest.cols() != 2))
      Dest.resize( 2, 2 );
    Dest.setZero();

    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
    Dest.coeffRef( 0, 0 ) += 800 * Arg[0] * Arg[0] - 400 * aux1 + 2;
    Dest.coeffRef( 0, 1 ) -= 400 * Arg[0];
    Dest.coeffRef( 1, 0 ) -= 400 * Arg[0];
    Dest.coeffRef( 1, 1 ) += 200.;
  }
};
//================================================================================



//#################################
//#################################
int main( int argc, char *argv[] ) {

  try {
    Eigen::IOFormat CleanFmt( Eigen::StreamPrecision, 0, ", ", ",", "", "", "(", ")" );

    // Test 1: Test Steilhaug's CG on Dixon's tridiagonal quadratic
    int n = 10;
    MatrixType H( n, n );
    VectorType c( n );
    RealType radius = 3.;

    c.setZero();
    c[0] = -2;
    c[n - 1] = -2;

    H.setZero();
    for ( int i = 0; i < n; i++ ) {
      H.coeffRef( i, i ) = 4;
      if ( i < n - 1 ) {
        H.coeffRef( i, i + 1 ) = -2;
        H.coeffRef( i + 1, i ) = -2;
      }
    }

    VectorType p;
    SteihaugCGMethod<DefaultConfigurator> scgSolver( H, c, radius, 1e-8, 100 );
    scgSolver.solve( p );

    std::cout << "Solution: " << p.format( CleanFmt ) << std::endl;
    std::cout << "Value: " << c.dot( p ) + 0.5 * p.dot( H * p ) << std::endl;

    // Test 2: Test trust-region Newton-CG on the Rosenbrock function
    std::cout << " ----------------------------- ROSENBROCK ----------------------------- " << std::endl;
    TestEnergy E;
    TestGradient DE;
    TestHessian D2E;

    VectorType start( 2 ), solution( 2 );
    start[0] = 1.2;
    start[1] = 1.2;

    TrustRegionNewton<DefaultConfigurator> trncgSolver( E, DE, D2E, 1., 10. );
    trncgSolver.solve( start, solution );

    std::cout << "Solution: " << solution.format( CleanFmt ) << std::endl;
    std::cout << "Value: " << E( solution ) << std::endl;

  }
  catch ( BasicException &el ) {
    std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXCEPTION: " << std::endl << el.getMessage() << std::endl
              << std::flush;
  }

  return 0;
}