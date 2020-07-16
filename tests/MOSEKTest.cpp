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
#include <exception>

#include <goast/Core.h>
#include <goast/Optimization/interfaces/MOSEKInterface.h>

//================================================================================

typedef DefaultConfigurator ConfiguratorType;

typedef typename ConfiguratorType::RealType RealType;

typedef typename ConfiguratorType::VectorType VectorType;
typedef typename ConfiguratorType::SparseMatrixType MatrixType;
typedef typename ConfiguratorType::TripletType TripletType;
typedef typename ConfiguratorType::VecType VecType;
typedef typename ConfiguratorType::MatType MatType;
typedef std::vector<TripletType> TripletListType;


//================================================================================
//================================================================================




//#################################
//#################################
int main( int argc, char *argv[] ) {
  try {
    // Setup problem taken from MOSEK documentation
    MatrixType Q( 3, 3 );
    Q.coeffRef( 0, 0 ) = 3.;
    Q.coeffRef( 0, 2 ) = -1.;
    Q.coeffRef( 1, 1 ) = 0.2;
    Q.coeffRef( 2, 0 ) = -1.;
    Q.coeffRef( 2, 2 ) = 2.;
    Q.makeCompressed();

    VectorType c( 3 );
    c << 0., -1., 0.;

    RealType cf = 0.;

    MatrixType A( 1, 3 );
    A.coeffRef( 0, 0 ) = 1;
    A.coeffRef( 0, 1 ) = 1;
    A.coeffRef( 0, 2 ) = 1;
    A.makeCompressed();

    VectorType constantLowerBounds( 3 );
    constantLowerBounds.setZero();

    VectorType constantUpperBounds( 3 );
    constantUpperBounds.setConstant( std::numeric_limits<RealType>::infinity());

    VectorType linearLowerBounds( 1 );
    linearLowerBounds[0] = 1.;

    VectorType linearUpperBounds( 1 );
    linearUpperBounds[0] = std::numeric_limits<RealType>::infinity();

    VectorType solution( 3 );
    MOSEKQuadraticSolver<DefaultConfigurator> Solver( Q, c, cf, A,
                                                      linearLowerBounds, linearUpperBounds,
                                                      constantLowerBounds, constantUpperBounds, false );
    Solver.solve( solution );

    std::cout << "Solution: " << solution << std::endl;

    std::cout << "Objective: " << 0.5 * solution.transpose() * Q * solution + c.dot(solution) + cf << std::endl;
  }
  catch ( const std::exception &e ) {
    std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXCEPTION: " << std::endl << e.what() << std::endl;
  }

  return 0;
}