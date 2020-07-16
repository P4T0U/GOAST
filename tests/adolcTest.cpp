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
#include <goast/external/adolcTraceless.h>

#include <goast/NRIC/NRICMap.h>
#include <goast/NRIC/Admissibility.h>

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

typedef ShellDeformation<
        ConfiguratorType,
        NonlinearMembraneDeformation<ConfiguratorType>,
        SimpleBendingDeformation<ConfiguratorType> > ShellDeformationType;


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
    Dest = 100 * aux1 * aux1 + aux2 * aux2;
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



//#################################
//#################################
int main( int argc, char *argv[] ) {

  try {
    // IO options
    Eigen::IOFormat CleanFmt( 4, 0, ", ", ",", "", "", "(", ")" );
    std::cout << "=====================================================" << std::endl;
    std::cout << "Start simple traceless autodiff test..." << std::endl;
    std::cout << "=====================================================" << std::endl;

    TestEnergy<DefaultConfigurator> E;
    TestGradient<DefaultConfigurator> DE;
    TestHessian<DefaultConfigurator> D2E;

    VectorType start( 2 ), solution( 2 );
    start[0] = 1.2;
    start[1] = 1.2;

    std::cout << "Without autodiff:" << std::endl;
    RealType fVal;
    auto t_start = std::chrono::high_resolution_clock::now();
    E.apply( start, fVal );
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << " - Evaluated E at " << start.format( CleanFmt ) << " to " << fVal << " in "
              << std::chrono::duration<double, std::milli>( t_end - t_start ).count() << "ms." << std::endl;

    VectorType fGradVal;
    t_start = std::chrono::high_resolution_clock::now();
    DE.apply( start, fGradVal );
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << " - Evaluated DE at " << start.format( CleanFmt ) << " to "
              << fGradVal.format( CleanFmt ) << " in "
              << std::chrono::duration<double, std::milli>( t_end - t_start ).count() << "ms." << std::endl;


    // ---- Autodiff ----
    std::cout << "With traceless autodiff:" << std::endl;
    AdolcForwardJacobian<DefaultConfigurator, TestEnergy> autoE;
    t_start = std::chrono::high_resolution_clock::now();
    autoE.apply( start, fVal );
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << " - Evaluated E at " << start.format( CleanFmt ) << " to " << fVal << " in "
              << std::chrono::duration<double, std::milli>( t_end - t_start ).count() << "ms." << std::endl;

    MatrixType sparseGrad;
    typename DefaultConfigurator::FullMatrixType Grad;

    t_start = std::chrono::high_resolution_clock::now();
    autoE.applyWithJacobian( start, fVal, Grad );
    t_end = std::chrono::high_resolution_clock::now();
    std::cout << " - Evaluated DE at " << start.format( CleanFmt ) << " to "
              << Grad.format( CleanFmt ) << " in "
              << std::chrono::duration<double, std::milli>( t_end - t_start ).count() << "ms." << std::endl;

  }
  catch ( BasicException &el ) {
    std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXCEPTION: " << std::endl << el.getMessage() << std::endl
              << std::flush;
  }

  return 0;
}