// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Interfaces to use ADOL-C's traceless forward mode
 * \author Sassen
 *
 * \warning This is still work in progress!
 */

#ifndef _ADOLCTAPELESS_H
#define _ADOLCTAPELESS_H

#ifdef GOAST_WITH_ADOLC
#define ADOLC_TRACELESS

#include "adolcIncludes.h"
#include <goast/Core/SmallVecMat.h>

//#define ADOLC_TIMINGS
#ifdef ADOLC_TIMINGS
#include <chrono>
#endif

/**
 * \brief Configurator such that operators use ADOL-C's traceless adouble type
 * \author Sassen
 */
struct adtlConfigurator {
public:
  typedef adtl::adouble RealType;
  typedef adtl::adouble ActiveScalar;
  typedef Eigen::Matrix<ActiveScalar, Eigen::Dynamic, 1> VectorType;
  // second template argument of SparseMatrix: union of bit flags controlling the storage scheme. Currently the only
  // possibility is ColMajor or RowMajor. The default is 0 which means column-major.
  typedef Eigen::SparseMatrix<ActiveScalar, 1, SuiteSparse_long> SparseMatrixType;
  typedef Eigen::Matrix<ActiveScalar, Eigen::Dynamic, Eigen::Dynamic> FullMatrixType;

  typedef SmallVec3<RealType> VecType;
  typedef SmallMat33<RealType> MatType;

  typedef Eigen::Triplet<ActiveScalar> TripletType;
  typedef std::vector<bool> MaskType;
};


/**
 * \brief Apply ADOL-C's traceless mode to a functional
 * \author Sassen
 * \tparam ConfiguratorType Container with data types
 * \tparam OpType Operator which will be differentiated
 */
template<typename ConfiguratorType, template<typename> class OpType>
class AdolcForwardJacobian : public OpType<adtlConfigurator> {
  typedef adtl::adouble ActiveScalar;
public:
  typedef typename OpType<ConfiguratorType>::DomainType DomainType;
  typedef typename OpType<ConfiguratorType>::RangeType RangeType;

  typedef typename OpType<adtlConfigurator>::DomainType ActiveInput;
  typedef typename OpType<adtlConfigurator>::RangeType ActiveValue;

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

  /**
   * \brief Construct functional with ADOL-C autodiff
   * \tparam Args types of arguments to pass to the raw functional
   * \param args arguments to pass to the raw functional
   */
  template<typename... Args>
  explicit AdolcForwardJacobian( Args &... args ) : OpType<adtlConfigurator>( args... ) {}

  /**
   * \brief Apply functional without return Jacobian
   * \param Arg point to evaluate at
   * \param Dest value of functional
   * \warning It is ill-advised to use this method, as it is very slow compared to a raw apply-call of the original functional
   */
  template<typename RangeType>
  typename std::enable_if<std::is_floating_point<RangeType>::value>::type apply( const DomainType &Arg, RangeType &Dest ) const {
    unsigned int domainDimension = Arg.size();
    adtl::setNumDir( domainDimension );

#ifdef ADOLC_TIMINGS
    auto t_0 = std::chrono::high_resolution_clock::now();
#endif

    ActiveInput ax = Arg.template cast<ActiveScalar>();
    ActiveValue av;

#ifdef ADOLC_TIMINGS
    auto t_1 = std::chrono::high_resolution_clock::now();
#endif

    OpType<adtlConfigurator>::apply( ax, av );

#ifdef ADOLC_TIMINGS
    auto t_2 = std::chrono::high_resolution_clock::now();
#endif

    Dest = av.getValue();

#ifdef ADOLC_TIMINGS
    auto t_3 = std::chrono::high_resolution_clock::now();
    std::cout << "applyWithJacobian: " << std::chrono::duration<double, std::milli>( t_1 - t_0 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_2 - t_1 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_3 - t_2 ).count() << std::endl;
#endif
  }

  /**
   * \brief Apply functional without return Jacobian
   * \param Arg point to evaluate at
   * \param Dest value of functional
   * \warning It is ill-advised to use this method, as it is very slow compared to a raw apply-call of the original functional
   */
   template<typename RangeType>
  typename std::enable_if<!std::is_floating_point<RangeType>::value>::type apply( const DomainType &Arg, RangeType &Dest ) const {
    unsigned int domainDimension = Arg.size();
    adtl::setNumDir( domainDimension );

#ifdef ADOLC_TIMINGS
    auto t_0 = std::chrono::high_resolution_clock::now();
#endif

    ActiveInput ax = Arg.template cast<ActiveScalar>();
    ActiveValue av;

#ifdef ADOLC_TIMINGS
    auto t_1 = std::chrono::high_resolution_clock::now();
#endif

    OpType<adtlConfigurator>::apply( ax, av );

#ifdef ADOLC_TIMINGS
    auto t_2 = std::chrono::high_resolution_clock::now();
#endif

    Dest.resize( av.size());
    for ( unsigned int i = 0; i < av.size(); i++ ) {
      Dest[i] = av[i].getValue();
    }

#ifdef ADOLC_TIMINGS
    auto t_3 = std::chrono::high_resolution_clock::now();
    std::cout << "applyWithJacobian: " << std::chrono::duration<double, std::milli>( t_1 - t_0 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_2 - t_1 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_3 - t_2 ).count() << std::endl;
#endif
  }

  /**
   * \brief Apply functional and return Jacobian
   * \param Arg point to evaluate at
   * \param Dest value of functional
   * \param Jacobian dense Jacobian of functional at point
   */
  template<typename RangeType>
  typename std::enable_if<std::is_floating_point<RangeType>::value>::type applyWithJacobian( const DomainType &Arg, RangeType &Dest, FullMatrixType &Jacobian ) const {
    unsigned int domainDimension = Arg.size();
    adtl::setNumDir( domainDimension );


#ifdef ADOLC_TIMINGS
    auto t_0 = std::chrono::high_resolution_clock::now();
#endif

    ActiveInput ax = Arg.template cast<ActiveScalar>();
    ActiveValue av;
#ifdef ADOLC_TIMINGS
    auto t_1 = std::chrono::high_resolution_clock::now();
#endif

    for ( unsigned int j = 0; j < domainDimension; j++ )
      for ( unsigned int i = 0; i < domainDimension; i++ )
        ax[i].setADValue( j, i == j ? 1 : 0 );

#ifdef ADOLC_TIMINGS
    auto t_2 = std::chrono::high_resolution_clock::now();
#endif

    OpType<adtlConfigurator>::apply( ax, av );

#ifdef ADOLC_TIMINGS
    auto t_3 = std::chrono::high_resolution_clock::now();
#endif

    Dest = 0.;
    Jacobian.resize( 1, domainDimension );
    Jacobian.setZero();

      Dest = av.getValue();
      for ( unsigned int j = 0; j < domainDimension; j++ ) {
        Jacobian( 0, j ) = av.getADValue( j );
      }


#ifdef ADOLC_TIMINGS
    auto t_4 = std::chrono::high_resolution_clock::now();

    std::cout << "applyWithJacobian: " << std::chrono::duration<double, std::milli>( t_1 - t_0 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_2 - t_1 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_3 - t_2 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_4 - t_3 ).count() << std::endl;
#endif
  }

  /**
   * \brief Apply functional and return Jacobian
   * \param Arg point to evaluate at
   * \param Dest value of functional
   * \param Jacobian dense Jacobian of functional at point
   */
  template<typename RangeType>
  typename std::enable_if<!std::is_floating_point<RangeType>::value>::type applyWithJacobian( const DomainType &Arg, RangeType &Dest, FullMatrixType &Jacobian ) const {
    unsigned int domainDimension = Arg.size();
    adtl::setNumDir( domainDimension );


#ifdef ADOLC_TIMINGS
    auto t_0 = std::chrono::high_resolution_clock::now();
#endif

    ActiveInput ax = Arg.template cast<ActiveScalar>();
    ActiveValue av;
#ifdef ADOLC_TIMINGS
    auto t_1 = std::chrono::high_resolution_clock::now();
#endif

    for ( unsigned int j = 0; j < domainDimension; j++ )
      for ( unsigned int i = 0; i < domainDimension; i++ )
        ax[i].setADValue( j, i == j ? 1 : 0 );

#ifdef ADOLC_TIMINGS
    auto t_2 = std::chrono::high_resolution_clock::now();
#endif

    OpType<adtlConfigurator>::apply( ax, av );

#ifdef ADOLC_TIMINGS
    auto t_3 = std::chrono::high_resolution_clock::now();
#endif

    long rangeDimension = av.size();
    Dest.resize( rangeDimension );
    Dest.setZero();
    Jacobian.resize( rangeDimension, domainDimension );
    Jacobian.setZero();

    for ( unsigned int i = 0; i < rangeDimension; i++ ) {
      Dest[i] = av[i].getValue();
      for ( unsigned int j = 0; j < domainDimension; j++ ) {
        Jacobian( i, j ) = av[i].getADValue( j );
      }
    }

#ifdef ADOLC_TIMINGS
    auto t_4 = std::chrono::high_resolution_clock::now();

    std::cout << "applyWithJacobian: " << std::chrono::duration<double, std::milli>( t_1 - t_0 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_2 - t_1 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_3 - t_2 ).count() << " / "
              << std::chrono::duration<double, std::milli>( t_4 - t_3 ).count() << std::endl;
#endif
  }
};


#endif // GOAST_WITH_ADOLC

#endif //_ADOLCTAPELESS_H
