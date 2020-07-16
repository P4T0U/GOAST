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
#include <goast/DiscreteShells.h>


#include <goast/GeodesicCalculus/ElasticMean.h>

#include <goast/external/ipoptBoxConstraintSolver.h>
#include <goast/external/ipoptNonlinearConstraintSolver.h>

//================================================================================

typedef DefaultConfigurator ConfiguratorType;

typedef typename ConfiguratorType::RealType RealType;

typedef typename ConfiguratorType::VectorType VectorType;
typedef typename ConfiguratorType::SparseMatrixType MatrixType;
typedef typename ConfiguratorType::TripletType TripletType;
typedef typename ConfiguratorType::VecType VecType;
typedef typename ConfiguratorType::MatType MatType;
typedef std::vector<TripletType> TripletListType;

typedef ShellDeformation<ConfiguratorType, NonlinearMembraneDeformation<ConfiguratorType>, SimpleBendingDeformation<ConfiguratorType> > ShellDeformationType;


//================================================================================
// f(x,y) = 100 * (y - x^2)^2 + (1-x)^2
template<typename ConfiguratorType>
class TestEnergy : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::RealType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
public:
  void apply( const VectorType &Arg, RealType &Dest ) const {
    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
    RealType aux2 = 1 - Arg[0];
    Dest += 100 * aux1 * aux1 + aux2 * aux2;
  }
};

template<typename ConfiguratorType>
class TestGradient : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
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

template<typename ConfiguratorType>
class TestHessian : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
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

//================================================================================
// min f(x,y) = x+y^2
// where x >= 0, y >= 0, x^2 + y^2 >= 2, x^2 + y^2 <= 4
template<typename ConfiguratorType>
class ConstraintTestEnergy : public BaseOp<VectorType, RealType> {

public:
  void apply( const VectorType &Arg, RealType &Dest ) const {
    Dest = Arg[0] + Arg[1] * Arg[1];
  }
};

template<typename ConfiguratorType>
class ConstraintTestGradient : public BaseOp<VectorType, VectorType> {

public:
  void apply( const VectorType &Arg, VectorType &Dest ) const {
    if ( Dest.size() != Arg.size())
      Dest.resize( Arg.size());
    Dest.setZero();

    Dest[0] += 1;
    Dest[1] += 2 * Arg[1];
  }
};

//template<typename ConfiguratorType>
//class ConstraintTestHessian : public BaseOp<VectorType, MatrixType> {
//
//public:
//  void apply( const VectorType &Arg, MatrixType &Dest ) const {
//
//    if ((Dest.rows() != 2) || (Dest.cols() != 2))
//      Dest.resize( 2, 2 );
//    Dest.setZero();
//
//    RealType aux1 = Arg[1] - Arg[0] * Arg[0];
//    Dest.coeffRef( 0, 0 ) += 800 * Arg[0] * Arg[0] - 400 * aux1 + 2;
//    Dest.coeffRef( 0, 1 ) -= 400 * Arg[0];
//    Dest.coeffRef( 1, 0 ) -= 400 * Arg[0];
//    Dest.coeffRef( 1, 1 ) += 200.;
//  }
//};

template<typename ConfiguratorType>
class ConstraintTestConstraints : public BaseOp<VectorType, VectorType> {

public:
  void apply( const VectorType &Arg, VectorType &Dest ) const {
    if ( Dest.size() != 1 )
      Dest.resize( 1 );
    Dest.setZero();

    Dest[0] += Arg[0] * Arg[0] + Arg[1] * Arg[1];
  }

  int getTargetDimension() const {
    return 1;
  }
};

template<typename ConfiguratorType>
class ConstraintTestConstraintsGradient : public BaseOp<VectorType, typename ConfiguratorType::SparseMatrixType> {

public:
  void apply( const VectorType &Arg, typename ConfiguratorType::SparseMatrixType &Dest ) const {
    if ((Dest.rows() != 1) || (Dest.cols() != 2))
      Dest.resize( 1, 2 );
    Dest.setZero();

    Dest.insert( 0, 0 ) = 2 * Arg[0];
    Dest.insert( 0, 1 ) = 2 * Arg[1];
  }

  int getNNZ() const override {
    return 2;
  }
};


//================================================================================



//#################################
//#################################
int main( int argc, char *argv[] ) {

  try {
    std::cout << "Eigen version: " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION
              << std::endl;
    std::cerr << "=====================================================" << std::endl;
    std::cerr << "Start simple optimization test..." << std::endl;

    TestEnergy<DefaultConfigurator> E;
    TestGradient<DefaultConfigurator> DE;
    TestHessian<DefaultConfigurator> D2E;

    VectorType start( 2 ), solution( 2 );
    start[0] = 1.2;
    start[1] = 1.2;

    std::vector<double> lowerBounds{ -10, -10 };
    std::vector<double> upperBounds{ 10, 10 };

    std::cerr << "Start gradient descent with Armijo... " << std::endl;
    GradientDescent<ConfiguratorType>( E, DE, 1000, 1e-8, ARMIJO ).solve( start, solution );
    std::cerr << solution << std::endl << std::endl;

    std::cerr << "Start first order unconstraint with ipopt... " << std::endl;
    IpoptBoxConstraintFirstOrderSolver<ConfiguratorType>( E, DE, 1000, 1e-9, lowerBounds, upperBounds ).solve( start,
                                                                                                               solution );
    std::cerr << solution << std::endl << std::endl;
    std::cerr << "=====================================================" << std::endl;

    std::cerr << "=====================================================" << std::endl;
    std::cerr << "Start simple first order box constraint with ipopt test... " << std::endl;
    IpoptBoxConstraintFirstOrderSolver<ConfiguratorType>( E, DE, 1000, 1e-9, lowerBounds, upperBounds ).solve( start,
                                                                                                               solution );
    std::cerr << solution << std::endl << std::endl;
    std::cerr << "=====================================================" << std::endl;

    std::cerr << "=====================================================" << std::endl;
    std::cerr << "Start simple first order nonlinear constraint with ipopt test... " << std::endl;
    start[0] = 1;
    start[1] = 1;

    ConstraintTestEnergy<DefaultConfigurator> E2;
    ConstraintTestGradient<DefaultConfigurator> DE2;
//    TestHessian<DefaultConfigurator> D2E2;
    ConstraintTestConstraints<DefaultConfigurator> C;
    ConstraintTestConstraintsGradient<DefaultConfigurator> DC;

    std::vector<RealType> lowerNonlinearBounds;
    lowerNonlinearBounds.push_back( 2 );
    std::vector<RealType> upperNonlinearBounds;
    upperNonlinearBounds.push_back( 4 );

    std::vector<RealType> lowerLinearBounds( 2, 0 );
    std::vector<RealType> upperLinearBounds( 2, 2.e+19 );

    IpoptFirstOrderSolver<ConfiguratorType> simpleConstraintSolver( E2, DE2, C, DC, 100, 1e-8, lowerLinearBounds,
                                                                    upperLinearBounds, lowerNonlinearBounds,
                                                                    upperNonlinearBounds, 0 );
    simpleConstraintSolver.solve( start, solution );
    std::cerr << solution << std::endl << std::endl;
    std::cerr << "=====================================================" << std::endl;


  }
  catch ( BasicException &el ) {
    std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl
              << std::flush;
  }

  return 0;
}