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
#include <goast/Optimization/ByrdOmojokunSQP.h>
#include <goast/Optimization/augmentedLagrange.h>
#include <goast/external/ipoptNonlinearConstraintSolver.h>

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
  typedef typename ConfiguratorType::TensorType TensorType;

//================================================================================
// f(x,y) = exp(x0*x1*x2*x3*x4)- 0.5* (x0^3 + x1^3 + 1)^2
class TestEnergy : public BaseOp<VectorType, RealType> {

public:
  void apply( const VectorType &Arg, RealType &Dest ) const override {
    RealType aux1 = std::exp( Arg[0] * Arg[1] * Arg[2] * Arg[3] * Arg[4] );
    RealType aux2 = 0.5 * std::pow( std::pow( Arg[0], 3 ) + std::pow( Arg[1], 3 ) + 1, 2 );
    Dest = aux1 - aux2;

//    if ( Dest[0] == std::numeric_limits<RealType>::infinity()) {
//      std::cout << "Aux: " << aux1 << " - " << aux2 << std::endl;
//      std::cout << "Arg_norm: " << Arg.norm() << std::endl;
//      std::cout << "Arg*: " << Arg[0] * Arg[1] * Arg[2] * Arg[3] * Arg[4] << std::endl;
//    }
  }
};

class TestGradient : public BaseOp<VectorType, VectorType> {

public:
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    if ( Dest.size() != Arg.size())
      Dest.resize( Arg.size());
    Dest.setZero();

    RealType aux1 = std::exp( Arg[0] * Arg[1] * Arg[2] * Arg[3] * Arg[4] );
    RealType x1squared = Arg[0] * Arg[0];
    RealType x2squared = Arg[1] * Arg[1];
    RealType x1cube = Arg[0] * Arg[0] * Arg[0];
    RealType x2cube = Arg[1] * Arg[1] * Arg[1];
    // [0] = -3 x1^2 (1 + x1^3 + x2^3) + E^(x1 x2 x3 x4 x5) x2 x3 x4 x5
    Dest[0] = -3. * x1squared * (1 + x1cube + x2cube) + aux1 * Arg[1] * Arg[2] * Arg[3] * Arg[4];
    // [1] = -3 x2^2 (1 + x1^3 + x2^3) + E^(x1 x2 x3 x4 x5) x1 x3 x4 x5
    Dest[1] = -3. * x2squared * (1 + x1cube + x2cube) + aux1 * Arg[0] * Arg[2] * Arg[3] * Arg[4];
    // [2] = E^(x1 x2 x3 x4 x5) x1 x2 x4 x5
    Dest[2] = aux1 * Arg[0] * Arg[1] * Arg[3] * Arg[4];
    // [3] = E^(x1 x2 x3 x4 x5) x1 x2 x3 x5
    Dest[3] = aux1 * Arg[0] * Arg[1] * Arg[2] * Arg[4];
    // [4] = E^(x1 x2 x3 x4 x5) x1 x2 x3 x4
    Dest[4] = aux1 * Arg[0] * Arg[1] * Arg[2] * Arg[3];
  }
};

class TestHessian : public BaseOp<VectorType, MatrixType> {

public:
  void apply( const VectorType &Arg, MatrixType &Dest ) const override {

    if ((Dest.rows() != 5) || (Dest.cols() != 5))
      Dest.resize( 5, 5 );
    Dest.setZero();

    TripletListType tripletList;
    pushTriplets( Arg, tripletList );

    Dest.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest ) const override {
    RealType expTerm = std::exp( Arg[0] * Arg[1] * Arg[2] * Arg[3] * Arg[4] );
    RealType onePlusTerm = 1 + Arg[0] * Arg[1] * Arg[2] * Arg[3] * Arg[4];
    RealType x0squared = Arg[0] * Arg[0];
    RealType x1squared = Arg[1] * Arg[1];
    RealType x2squared = Arg[2] * Arg[2];
    RealType x3squared = Arg[3] * Arg[3];
    RealType x4squared = Arg[4] * Arg[4];
    RealType x0cube = Arg[0] * Arg[0] * Arg[0];
    RealType x1cube = Arg[1] * Arg[1] * Arg[1];

    // (0,0) = -15. x0^4 + x0 (-6. - 6. x1^3) +  E^(x0 x1 x2 x3 x4) x1^2 x2^2 x3^2 x4^2
    Dest.emplace_back( 0, 0, -15. * x0squared * x0squared + Arg[0] * (-6. - 6. * x1cube) +
                             expTerm * x1squared * x2squared * x3squared * x4squared );

    // (0,1) = -9. x0^2 x1^2 + E^(x0 x1 x2 x3 x4) x2 x3 x4 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 0, 1, -9. * x0squared * x1squared + expTerm * Arg[2] * Arg[3] * Arg[4] * onePlusTerm );
    Dest.emplace_back( 1, 0, -9. * x0squared * x1squared + expTerm * Arg[2] * Arg[3] * Arg[4] * onePlusTerm );

    // (0,2) =E^(x0 x1 x2 x3 x4) x1 x3 x4 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 0, 2, expTerm * Arg[1] * Arg[3] * Arg[4] * onePlusTerm );
    Dest.emplace_back( 2, 0, expTerm * Arg[1] * Arg[3] * Arg[4] * onePlusTerm );

    // (0,3) = E^(x0 x1 x2 x3 x4) x1 x2 x4 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 0, 3, expTerm * Arg[1] * Arg[2] * Arg[4] * onePlusTerm );
    Dest.emplace_back( 3, 0, expTerm * Arg[1] * Arg[2] * Arg[4] * onePlusTerm );

    // (0,4) = E^(x0 x1 x2 x3 x4) x1 x2 x3 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 0, 4, expTerm * Arg[1] * Arg[2] * Arg[3] * onePlusTerm );
    Dest.emplace_back( 4, 0, expTerm * Arg[1] * Arg[2] * Arg[3] * onePlusTerm );

    // (1,1) = (-6. - 6. x0^3) x1 - 15. x1^4 + E^(x0 x1 x2 x3 x4) x0^2 x2^2 x3^2 x4^2
    Dest.emplace_back( 1, 1, Arg[1] * (-6. - 6. * x0cube) - 15. * x1squared * x1squared +
                             expTerm * x0squared * x2squared * x3squared * x4squared );

    // (1,2) = E^(x0 x1 x2 x3 x4) x0 x3 x4 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 1, 2, expTerm * Arg[0] * Arg[3] * Arg[4] * onePlusTerm );
    Dest.emplace_back( 2, 1, expTerm * Arg[0] * Arg[3] * Arg[4] * onePlusTerm );

    // (1,3) = E^(x0 x1 x2 x3 x4) x0 x2 x4 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 1, 3, expTerm * Arg[0] * Arg[2] * Arg[4] * onePlusTerm );
    Dest.emplace_back( 3, 1, expTerm * Arg[0] * Arg[2] * Arg[4] * onePlusTerm );

    // (1,4) = E^(x0 x1 x2 x3 x4) x0 x2 x3 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 1, 4, expTerm * Arg[0] * Arg[2] * Arg[3] * onePlusTerm );
    Dest.emplace_back( 4, 1, expTerm * Arg[0] * Arg[2] * Arg[3] * onePlusTerm );

    // (2,2) = E^(x0 x1 x2 x3 x4) x0^2 x1^2 x3^2 x4^2
    Dest.emplace_back( 2, 2, expTerm * x0squared * x1squared * x3squared * x4squared );

    // (2,3) = E^(x0 x1 x2 x3 x4) x0 x1 x4 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 2, 3, expTerm * Arg[0] * Arg[1] * Arg[4] * onePlusTerm );
    Dest.emplace_back( 3, 2, expTerm * Arg[0] * Arg[1] * Arg[4] * onePlusTerm );

    // (2,4) = E^(x0 x1 x2 x3 x4) x0 x1 x3 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 2, 4, expTerm * Arg[0] * Arg[1] * Arg[3] * onePlusTerm );
    Dest.emplace_back( 4, 2, expTerm * Arg[0] * Arg[1] * Arg[3] * onePlusTerm );

    // (3,3) = E^(x0 x1 x2 x3 x4) x0^2 x1^2 x2^2 x4^2
    Dest.emplace_back( 3, 3, expTerm * x0squared * x1squared * x2squared * x4squared );

    // (3,4) = E^(x0 x1 x2 x3 x4) x0 x1 x2 (1 + x0 x1 x2 x3 x4)
    Dest.emplace_back( 3, 4, expTerm * Arg[0] * Arg[1] * Arg[2] * onePlusTerm );
    Dest.emplace_back( 4, 3, expTerm * Arg[0] * Arg[1] * Arg[2] * onePlusTerm );

    // (4,4) = E^(x0 x1 x2 x3 x4) x0^2 x1^2 x2^2 x3^2
    Dest.emplace_back( 4, 4, expTerm * x0squared * x1squared * x2squared * x3squared );
  }
};

class TestConstraints : public BaseOp<VectorType, VectorType> {

public:
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest.resize( 3 );
    Dest.setZero();

    Dest[0] = Arg[0] * Arg[0] + Arg[1] * Arg[1] + Arg[2] * Arg[2] + Arg[3] * Arg[3] + Arg[4] * Arg[4] - 10;
    Dest[1] = Arg[1] * Arg[2] - 5 * Arg[3] * Arg[4];
    Dest[2] = Arg[0] * Arg[0] * Arg[0] + Arg[1] * Arg[1] * Arg[1] + 1;
  }

  int getTargetDimension() const override {
    return 3;
  }
};

class TestConstraintsGradient : public BaseOp<VectorType, MatrixType> {

public:
  void apply( const VectorType &Arg, MatrixType &Dest ) const override {
    // {
    // {2 x1, 2 x2, 3 x3^2, 2 x4, 2 x5},
    // {0, x3, x2, -5 x5, -5 x4},
    // {3 x1^2, 3 x2^2, 0, 0, 0}
    //}
    Dest.resize( 3, 5 );
    Dest.setZero();

    // (0,_) = {2 x0, 2 x1, 3 x2^2, 2 x3, 2 x4}

    Dest.coeffRef( 0, 0 ) = 2. * Arg[0];
    Dest.coeffRef( 0, 1 ) = 2. * Arg[1];
    Dest.coeffRef( 0, 2 ) = 2. * Arg[2];
    Dest.coeffRef( 0, 3 ) = 2. * Arg[3];
    Dest.coeffRef( 0, 4 ) = 2. * Arg[4];

    // (1,_) = {0, x2, x1, -5 x4, -5 x3}
    Dest.coeffRef( 1, 1 ) = Arg[2];
    Dest.coeffRef( 1, 2 ) = Arg[1];
    Dest.coeffRef( 1, 3 ) = -5. * Arg[4];
    Dest.coeffRef( 1, 4 ) = -5. * Arg[3];

    // (2,_) = {3 x0^2, 3 x1^2, 0, 0, 0}
    Dest.coeffRef( 2, 0 ) = 3. * Arg[0] * Arg[0];
    Dest.coeffRef( 2, 1 ) = 3. * Arg[1] * Arg[1];
  }

  int getTargetDimension() const override {
    return 3;
  }

  int getNNZ() const override {
    return 11;
  }
};

template<typename ConfiguratorType>
class TestConstraintsHessian : public BaseOp<VectorType, TensorType> {

public:
  void apply( const VectorType &Arg, TensorType &Dest ) const override {
    Dest.resize( 3, 5, 5 );
    Dest.setZero();

    // (0,0) = (2,0,6*x0)
    Dest[0].coeffRef( 0, 0 ) = 2.;
    Dest[2].coeffRef( 0, 0 ) = 6. * Arg[0];

    // (1,1) = (2,0,6*x1)
    Dest[0].coeffRef( 1, 1 ) = 2.;
    Dest[2].coeffRef( 1, 1 ) = 6. * Arg[1];

    // (1,2) = (0,1,0)
    Dest[1].coeffRef( 1, 2 ) = 1.;
    Dest[1].coeffRef( 2, 1 ) = 1.;

    // (2,2) = (2,0,0)
    Dest[0].coeffRef( 2, 2 ) = 2.;
    // (3,3) = (2,0,0)
    Dest[0].coeffRef( 3, 3 ) = 2.;
    // (3,4) = (0, -5, 0)
    Dest[1].coeffRef( 3, 4 ) = -5.;
    Dest[1].coeffRef( 4, 3 ) = -5.;
    // (4,4) = (2,0,0)
    Dest[0].coeffRef( 4, 4 ) = 2.;
  }

  void apply( const VectorType &Arg, std::vector<MatrixType> &Dest ) const {
    Dest.resize( 3 );
    for ( auto t : { 1, 2, 3 } ) {
      Dest[t].resize( 5, 5 );
      Dest[t].setZero();
    }

    // (0,0) = (2,0,6*x0)
    Dest[0].coeffRef( 0, 0 ) = 2.;
    Dest[2].coeffRef( 0, 0 ) = 6. * Arg[0];

    // (1,1) = (2,0,6*x1)
    Dest[0].coeffRef( 1, 1 ) = 2.;
    Dest[2].coeffRef( 1, 1 ) = 6. * Arg[1];

    // (1,2) = (0,1,0)
    Dest[1].coeffRef( 1, 2 ) = 1.;
    Dest[1].coeffRef( 2, 1 ) = 1.;

    // (2,2) = (2,0,0)
    Dest[0].coeffRef( 2, 2 ) = 2.;
    // (3,3) = (2,0,0)
    Dest[0].coeffRef( 3, 3 ) = 2.;
    // (3,4) = (0, -5, 0)
    Dest[1].coeffRef( 3, 4 ) = -5.;
    Dest[1].coeffRef( 4, 3 ) = -5.;
    // (4,4) = (2,0,0)
    Dest[0].coeffRef( 4, 4 ) = 2.;
  }

  void pushTriplets( const VectorType &Arg, TripletListType &Dest, const VectorType &Lambda ) const override {
    // (0,0) = (2,0,6*x0)
    Dest.emplace_back( 0, 0, 2. * Lambda[0] + 6. * Arg[0] * Lambda[2] );

    // (1,1) = (2,0,6*x1)
    Dest.emplace_back( 1, 1, 2. * Lambda[0] + 6. * Arg[1] * Lambda[2] );

    // (1,2) = (0,1,0)
    Dest.emplace_back( 1, 2, Lambda[1] );
    Dest.emplace_back( 2, 1, Lambda[1] );

    // (2,2) = (2,0,0)
    Dest.emplace_back( 2, 2, 2. * Lambda[0] );
    // (3,3) = (2,0,0
    Dest.emplace_back( 3, 3, 2. * Lambda[0] );
    // (3,4) = (0, -5, 0)
    Dest.emplace_back( 3, 4, -5. * Lambda[1] );
    Dest.emplace_back( 4, 3, -5. * Lambda[1] );

    // (4,4) = (2,0,0)
    Dest.emplace_back( 4, 4, 2. * Lambda[0] );
  }

  void setTriplets( const DomainType &Arg, std::vector<TripletListType> &Dest, int hessOffset ) const override {

    Dest.resize( 3 );
    Dest[0].emplace_back( 0, 0, 2. );
    Dest[2].emplace_back( 0, 0, 6. * Arg[0] );

    // (1,1) = (2,0,6*x1)
    Dest[0].emplace_back( 1, 1, 2. );
    Dest[2].emplace_back( 1, 1, 6. * Arg[1] );

    // (1,2) = (0,1,0)
    Dest[1].emplace_back( 1, 2, 1. );
    Dest[1].emplace_back( 2, 1, 1. );

    // (2,2) = (2,0,0)
    Dest[0].emplace_back( 2, 2, 2. );
    // (3,3) = (2,0,0)
    Dest[0].emplace_back( 3, 3, 2. );
    // (3,4) = (0, -5, 0)
    Dest[1].emplace_back( 3, 4, -5. );
    Dest[1].emplace_back( 4, 3, -5. );
    // (4,4) = (2,0,0)
    Dest[0].emplace_back( 4, 4, 2. );

  }

  int getTargetDimension() const override {
    return 3;
  }


};

//================================================================================

template<typename M>
using LinearSolverT = Eigen::CholmodSupernodalLLT<M, Eigen::Lower>;

//template<typename M>
//using LinearSolverT = Eigen::UmfPackLU<M>;

//#################################
//#################################
int main( int argc, char *argv[] ) {

  try {
    Eigen::IOFormat CleanFmt( Eigen::StreamPrecision, 0, ", ", ",", "", "", "(", ")" );



    // Test 2: Test trust-region Newton-CG on the Rosenbrock function
    std::cout << " ----------------------------- Exercise 18.3 ----------------------------- " << std::endl;
    TestEnergy E;
    TestGradient DE;
    TestHessian D2E;

    TestConstraints G;
    TestConstraintsGradient DG;
    TestConstraintsHessian<DefaultConfigurator> D2G;

    VectorType start( 5 ), solution( 5 );
    start[0] = -1.8;
    start[1] = 1.7;
    start[2] = 1.9;
    start[3] = -0.8;
    start[4] = -0.8;

    solution.setZero();

//    ScalarValuedDerivativeTester<DefaultConfigurator>( E, DE, 1e-8 ).plotAllDirections( start, "testDE" );
//    VectorValuedDerivativeTester<DefaultConfigurator>( DE, D2E, 1e-8 ).plotAllDirections( start, "testD2E" );
//
//    VectorValuedDerivativeTester<DefaultConfigurator>( G, DG, 1e-8, 3 ).plotAllDirections( start, "testDG" );
//    TensorValuedDerivativeTester<DefaultConfigurator>( DG, D2G, 1e-8, 3 ).plotRandomDirections( start, 100, "testD2G" );

#ifdef GOAST_WITH_IPOPT
    std::vector<RealType> linearLowerBounds( start.size(), -2.e+19 );
    std::vector<RealType> linearUpperBounds( start.size(), +2.e+19 );

    std::vector<RealType> nonlinearLowerBounds( 3, 0. );
    std::vector<RealType> nonlinearUpperBounds( 3, 0. );


    IpoptSecondOrderSolver<DefaultConfigurator, TestConstraintsHessian, false> ipSolver( E, DE, D2E, G, DG, D2G, 1000,
                                                                                        1e-8,
                                                                                        linearLowerBounds,
                                                                                        linearUpperBounds,
                                                                                        nonlinearLowerBounds,
                                                                                        nonlinearUpperBounds,
                                                                                        0,
                                                                                        5 );
    auto t_start_ipopt = std::chrono::high_resolution_clock::now();
    ipSolver.solve( start, solution );
    auto t_end_ipopt = std::chrono::high_resolution_clock::now();

    std::cout << "IP - Solution: " << solution.format( CleanFmt ) << std::endl;
    std::cout << "IP - Value: " << E( solution ) << std::endl;
    std::cout << "IP - Time: " << std::chrono::duration<double, std::milli >( t_end_ipopt - t_start_ipopt ).count() << std::endl;
#endif

    AugmentedLagrangeMethod<DefaultConfigurator> alSolver( E, DE, D2E, G, DG, D2G, 1000, 1e-8, false );
    auto t_start = std::chrono::high_resolution_clock::now();
    alSolver.solve( start, solution );
    auto t_end = std::chrono::high_resolution_clock::now();

    std::cout << "AL - Solution: " << solution.format( CleanFmt ) << std::endl;
    std::cout << "AL - Value: " << E( solution ) << std::endl;
    std::cout << "AL - Time: " <<  std::fixed  << std::chrono::duration<double, std::milli >( t_end - t_start ).count() << std::endl;


    ByrdOmojokunSQP<DefaultConfigurator, LinearSolverT> bosqpSolver( E, DE, D2E, G, DG, D2G, 1e-8, 1000, false );

    bosqpSolver.setParameter( "minimal_stepsize", 1e-15 );
    bosqpSolver.setParameter( "diagonal_preconditioning", -1 );
    t_start = std::chrono::high_resolution_clock::now();
    bosqpSolver.solve( start, solution );
    t_end = std::chrono::high_resolution_clock::now();

    std::cout << "Start: " << start.format( CleanFmt ) << std::endl;
    std::cout << "Solution: " << solution.format( CleanFmt ) << std::endl;
    std::cout << "Value: " << E( solution ) << std::endl;
    std::cout << "Time: " << std::fixed << std::chrono::duration<double, std::milli >( t_end - t_start ).count() << std::endl;
  }
  catch ( BasicException &el ) {
    std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXCEPTION: " << std::endl << el.getMessage() << std::endl
              << std::flush;
  }

  return 0;
}