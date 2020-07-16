// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \brief Interfaces to use ADOL-C's traced mode
 * \author Sassen
 *
 * \warning This is untested and may very well not work! It is just included for the sake of completeness.
 */

#ifndef _ADOLCTAPELESS_H
#define _ADOLCTAPELESS_H

#ifdef GOAST_WITH_ADOLC
#define ADOLC_TRACED

#include "adolcIncludes.h"

struct adConfigurator {
public:
  typedef adouble RealType;
  typedef adouble ActiveScalarType;

  typedef Eigen::Matrix<adouble, Eigen::Dynamic, 1> ActiveVectorType;
  typedef Eigen::Map<ActiveVectorType> VectorType;
  // second template argument of SparseMatrix: union of bit flags controlling the storage scheme. Currently the only possibility is ColMajor or RowMajor. The default is 0 which means column-major.
//  typedef Eigen::SparseMatrix<adouble, 1, SuiteSparse_long> SparseMatrixType;
//  typedef Eigen::Matrix<adouble, Eigen::Dynamic, Eigen::Dynamic>  FullMatrixType;

  typedef SmallVec3 <RealType> VecType;
  typedef SmallMat33 <RealType> MatType;

  //typedef std::vector<Eigen::Ref<const EigenVector> > EigenVectorReferenceContainer;

  typedef Eigen::Triplet<adouble> TripletType;
  typedef std::vector<bool> MaskType;
};


template<typename ConfiguratorType, template<typename> class OpType>
class AdolcForwardJacobian : public OpType<adConfigurator> {
  typedef adouble ActiveScalar;
public:
  typedef typename OpType<ConfiguratorType>::DomainType DomainType;
  typedef typename OpType<ConfiguratorType>::RangeType RangeType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;

  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

  const int tag;

  /**
   * \brief Construct functional with ADOL-C traced autodiff
   * \tparam Args types of arguments to pass to the raw functional
   * \param args arguments to pass to the raw functional
   */
  template<typename... Args>
  explicit AdolcForwardJacobian( Args &... args ) : OpType<adConfigurator>( args... ), tag( rand()) {}

  typedef Eigen::Map<Eigen::Matrix<ActiveScalar, DomainType::SizeAtCompileTime, 1>> ActiveInput;
  typedef Eigen::Map<Eigen::Matrix<ActiveScalar, RangeType::SizeAtCompileTime, 1>> ActiveValue;

  /**
   * \brief Apply functional without return Jacobian
   * \param Arg point to evaluate at
   * \param Dest value of functional
   * \warning It is ill-advised to use this method, as it is very slow compared to a raw apply-call of the original functional
   */
  void apply( const DomainType &Arg, RangeType &Dest ) const {
    unsigned int domainDimension = Arg.size();
    std::cout << "Domain dim: " << domainDimension << std::endl;

    trace_on( tag );
    adouble *ax = new adouble[domainDimension];
    for ( int i = 0; i < domainDimension; i++ )
      ax[i] <<= Arg[i];

    ActiveInput ax_m( ax, domainDimension );

    adouble av[1];
    ActiveValue av_m( av );

    OpType<adConfigurator>::apply( ax_m, av_m );
    Dest.resize( 1 );
    for ( unsigned int i = 0; i < 1; i++ ) {
      av[i] >>= Dest[i];
    }
    delete[] ax;
    trace_off();
  }

  /**
   * \brief Apply functional and return Jacobian
   * \param Arg point to evaluate at
   * \param Dest value of functional
   * \param Jacobian dense Jacobian of functional at point
   */
  void applyWithJacobian( const DomainType &Arg, RangeType &Dest, FullMatrixType &Jacobian ) {
    unsigned int domainDimension = Arg.size();

    apply( Arg, Dest );
    long rangeDimension = Dest.size();

    Jacobian.resize( rangeDimension, domainDimension );
    Jacobian.setZero();

    double **J = new double *[rangeDimension];

    for ( int i = 0; i < rangeDimension; i++ )
      J[i] = new double[2 * domainDimension];


    jacobian( tag, rangeDimension, domainDimension, Arg.data(), J );

    for ( int i = 0; i < rangeDimension; i++ ) {
      for ( int j = 0; j < 2 * domainDimension; j++ ) {
        Jacobian( i, j ) = J[i][j];
      }
      delete[] J[i];
    }
    delete[] J;
  }

  /**
   * \brief Compute Jacobian at point
   * \param Arg point to evaluate at
   * \param Jacobian dense Jacobian of functional at point
   */
  void applyJacobian( const VectorType &Arg, FullMatrixType &Jacobian ) {
    unsigned int domainDimension = Arg.size();

    long rangeDimension = 1;

    Jacobian.resize( rangeDimension, domainDimension );
    Jacobian.setZero();

    double **J = new double *[rangeDimension];

    for ( int i = 0; i < rangeDimension; i++ )
      J[i] = new double[2 * domainDimension];


    jacobian( tag, rangeDimension, domainDimension, Arg.data(), J );

    for ( int i = 0; i < rangeDimension; i++ ) {
      for ( int j = 0; j < 2 * domainDimension; j++ ) {
        Jacobian( i, j ) = J[i][j];
      }
      delete[] J[i];
    }
    delete[] J;
  }

  void applyAdd( const DomainType &Arg, RangeType &Dest ) const {

  }

protected:

};


#endif // GOAST_WITH_ADOLC

#endif //_ADOLCTAPELESS_H
