// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2023 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include <goast/Core.h>

template<typename ConfiguratorType>
class LinearOperator : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

public:
  LinearOperator() = default;

  virtual ~LinearOperator() = default;

  virtual void apply( const VectorType &Arg, VectorType &Dest ) const = 0;

  virtual void apply( const FullMatrixType &Arg, FullMatrixType &Dest ) const {
    throw std::logic_error( "Method apply(FM) has to be implemented by a derived class." );
  }

  virtual void apply( const SparseMatrixType &Arg, SparseMatrixType &Dest ) const {
    throw std::logic_error( "Method apply(SpM) has to be implemented by a derived class." );
  }

  virtual void applyTransposed( const VectorType &Arg, VectorType &Dest ) const {
    throw std::logic_error( "Method applyTransposed has to be implemented by a derived class." );
  };

  virtual void applyTransposed( const FullMatrixType &Arg, FullMatrixType &Dest ) const {
    throw std::logic_error( "Method applyTransposed(FM) has to be implemented by a derived class." );
  }

  virtual void applyTransposed( const SparseMatrixType &Arg, SparseMatrixType &Dest ) const {
    throw std::logic_error( "Method applyTransposed(SpM) has to be implemented by a derived class." );
  }

  virtual void assembleTransformationMatrix( SparseMatrixType &T, bool transposed = false ) const {
    throw std::logic_error( "Method assembleTransformationMatrix has to be implemented by a derived class." );
  }

  virtual void assembleTransformationMatrix( FullMatrixType &T, bool transposed = false ) const {
    throw std::logic_error( "Method assembleTransformationMatrix has to be implemented by a derived class." );
  }

  template<typename RangeType>
  RangeType operator()( const RangeType &Arg ) const {
    RangeType Update;
    apply( Arg, Update );
    return Update;
  }

  template<typename RangeType>
  RangeType operator*( const RangeType &Arg ) const {
    RangeType Update;
    apply( Arg, Update );
    return Update;
  }


  template<typename RangeType>
  RangeType T( const RangeType &Arg ) const {
    RangeType Update;
    applyTransposed( Arg, Update );
    return Update;
  }

  virtual int rows() const {
    throw std::logic_error( "Method rows has to be implemented by a derived class." );
  }

  virtual int cols() const {
    throw std::logic_error( "Method cols has to be implemented by a derived class." );
  }
};


template<typename ConfiguratorType>
class IdentityOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

public:
  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = Arg;
  }

  void apply( const FullMatrixType &Arg, FullMatrixType &Dest ) const override {
    Dest = Arg;
  }

  void apply( const SparseMatrixType &Arg, SparseMatrixType &Dest ) const override {
    Dest = Arg;
  }
};


template<typename ConfiguratorType>
class MatrixOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

  SparseMatrixType m_A;

public:
  explicit MatrixOperator( const SparseMatrixType &A ) : m_A( A ) {}

  explicit MatrixOperator( SparseMatrixType &&A ) {
    m_A.swap( A );
  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_A * Arg;
  }

  void apply( const FullMatrixType &Arg, FullMatrixType &Dest ) const override {
    Dest = m_A * Arg;
  }

  void apply( const SparseMatrixType &Arg, SparseMatrixType &Dest ) const override {
    Dest = m_A * Arg;
  }

  void applyTransposed( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_A.transpose() * Arg;
  }

  void applyTransposed( const FullMatrixType &Arg, FullMatrixType &Dest ) const override {
    Dest = m_A.transpose() * Arg;
  }

  void applyTransposed( const SparseMatrixType &Arg, SparseMatrixType &Dest ) const override {
    Dest = m_A.transpose() * Arg;
  }

  int rows() const override {
    return m_A.rows();
  }

  int cols() const override {
    return m_A.cols();
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {
    T = transposed ? m_A.transpose() : m_A;
  }
};


template<typename ConfiguratorType>
class DiagonalOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

  VectorType m_D;

public:
  explicit DiagonalOperator( const VectorType &D ) : m_D( D ) {}

  explicit DiagonalOperator( VectorType &&D ) {
    m_D.swap( D );
  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_D.asDiagonal() * Arg;
  }

  void apply( const FullMatrixType &Arg, FullMatrixType &Dest ) const override {
    Dest = m_D.asDiagonal() * Arg;
  }

  void apply( const SparseMatrixType &Arg, SparseMatrixType &Dest ) const override {
    Dest = m_D.asDiagonal() * Arg;
  }

  void applyTransposed( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_D.asDiagonal() * Arg;
  }

  void applyTransposed( const FullMatrixType &Arg, FullMatrixType &Dest ) const override {
    Dest = m_D.asDiagonal() * Arg;
  }

  void applyTransposed( const SparseMatrixType &Arg, SparseMatrixType &Dest ) const override {
    Dest = m_D.asDiagonal() * Arg;
  }

  int rows() const override {
    return m_D.size();
  }

  int cols() const override {
    return m_D.size();
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {
    T = m_D.asDiagonal();
  }
};


template<typename ConfiguratorType>
class RankOneOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

  VectorType m_a, m_b;

public:
  RankOneOperator( const VectorType &a, const VectorType &b ) : m_a( a ), m_b( b ) {}

  RankOneOperator( VectorType &&a, VectorType &&b ) : m_a( a ), m_b( b ) {}

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_a * m_b.dot( Arg );
  }

  void applyTransposed( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_b * m_a.dot( Arg );
  }

  int rows() const override {
    return m_a.size();
  }

  int cols() const override {
    return m_b.size();
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {
    T = ( transposed ? m_a * m_b.transpose() : m_b * m_a.transpose()).sparseView();
  }

  void assembleTransformationMatrix( FullMatrixType &T, bool transposed ) const override {
    T = transposed ? m_a * m_b.transpose() : m_b * m_a.transpose();
  }
};


template<typename ConfiguratorType>
class ScaledOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

  std::unique_ptr<LinearOperator<ConfiguratorType>> m_A;

  RealType m_w_a;

public:
  ScaledOperator( RealType w_a, std::unique_ptr<LinearOperator<ConfiguratorType>> &&a ) : m_A( std::move( a )),
                                                                                          m_w_a( w_a ) {}

//  SummedOperator( VectorType &&a, VectorType &&b ) : m_a( a ), m_b( b ) {  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_w_a * ( *m_A )( Arg );
  }

  void applyTransposed( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_w_a * ( *m_A ).T( Arg );
  }

  int rows() const override {
    return m_A->rows();
  }

  int cols() const override {
    return m_A->cols();
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {
    SparseMatrixType matA, matB;
    m_A->assembleTransformationMatrix( matA, transposed );
    T = m_w_a * matA;
  }

  void assembleTransformationMatrix( FullMatrixType &T, bool transposed ) const override {
    SparseMatrixType matA, matB;
    m_A->assembleTransformationMatrix( matA, transposed );
    T = m_w_a * matA;
  }
};


template<typename ConfiguratorType>
class SummedOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

  std::unique_ptr<LinearOperator<ConfiguratorType>> m_A;
  std::unique_ptr<LinearOperator<ConfiguratorType>> m_B;

public:
  SummedOperator( std::unique_ptr<LinearOperator<ConfiguratorType>> &&a,
                  std::unique_ptr<LinearOperator<ConfiguratorType>> &&b ) : m_A( std::move( a )),
                                                                            m_B( std::move( b )) {}

//  SummedOperator( VectorType &&a, VectorType &&b ) : m_a( a ), m_b( b ) {  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = ( *m_A )( Arg ) + ( *m_B )( Arg );
  }

  void applyTransposed( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = ( *m_A ).T( Arg ) + ( *m_B ).T( Arg );
  }

  int rows() const override {
    return m_A->rows();
  }

  int cols() const override {
    return m_A->cols();
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {
    SparseMatrixType matA, matB;
    m_A->assembleTransformationMatrix( matA, transposed );
    m_B->assembleTransformationMatrix( matB, transposed );
    T = matA + matB;
  }

  void assembleTransformationMatrix( FullMatrixType &T, bool transposed ) const override {
    SparseMatrixType matA, matB;
    m_A->assembleTransformationMatrix( matA, transposed );
    m_B->assembleTransformationMatrix( matB, transposed );
    T = matA + matB;
  }
};

template<typename ConfiguratorType>
class WeightedSumOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

  std::unique_ptr<LinearOperator<ConfiguratorType>> m_A;
  std::unique_ptr<LinearOperator<ConfiguratorType>> m_B;

  RealType m_w_a, m_w_b;

public:
  WeightedSumOperator( RealType w_a,
                       std::unique_ptr<LinearOperator<ConfiguratorType>> &&a,
                       RealType w_b,
                       std::unique_ptr<LinearOperator<ConfiguratorType>> &&b ) : m_A( std::move( a )),
                                                                                 m_B( std::move( b )), m_w_a( w_a ),
                                                                                 m_w_b( w_b ) {}

//  SummedOperator( VectorType &&a, VectorType &&b ) : m_a( a ), m_b( b ) {  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_w_a * ( *m_A )( Arg ) + m_w_b * ( *m_B )( Arg );
  }

  void applyTransposed( const VectorType &Arg, VectorType &Dest ) const override {
    Dest = m_w_a * ( *m_A ).T( Arg ) + m_w_b * ( *m_B ).T( Arg );
  }

  int rows() const override {
    return m_A->rows();
  }

  int cols() const override {
    return m_A->cols();
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {
    SparseMatrixType matA, matB;
    m_A->assembleTransformationMatrix( matA, transposed );
    m_B->assembleTransformationMatrix( matB, transposed );
    T = m_w_a * matA + m_w_b * matB;
  }

  void assembleTransformationMatrix( FullMatrixType &T, bool transposed ) const override {
    SparseMatrixType matA, matB;
    m_A->assembleTransformationMatrix( matA, transposed );
    m_B->assembleTransformationMatrix( matB, transposed );
    T = m_w_a * matA + m_w_b * matB;
  }
};

template<typename ConfiguratorType>
class LinearlyCombinedOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;

  std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> m_Ops;
  const VectorType m_Weights;
  const int m_numOps;

public:
  template<class... Items>
  LinearlyCombinedOperator( const VectorType &Weights, Items &&... constraintOps )
          : m_Weights( Weights ), m_numOps( sizeof...( constraintOps )) {
    append_to_vector( m_Ops, std::forward<std::unique_ptr<LinearOperator<ConfiguratorType>>>(constraintOps)... );

  }

  LinearlyCombinedOperator( const VectorType &Weights,
                            std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &&Ops )
          : m_Weights( Weights ), m_numOps( Ops.size()) {

    m_Ops.swap( Ops );
  }


//  SummedOperator( VectorType &&a, VectorType &&b ) : m_a( a ), m_b( b ) {  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
    Dest.resize(this->rows());
    Dest.setZero();

    for (int i = 0; i < m_numOps; i++) {
      Dest += m_Weights[i] * m_Ops[i]->operator()( Arg );
    }
  }

  void applyTransposed( const VectorType &Arg, VectorType &Dest ) const override {
    Dest.resize(this->cols());
    Dest.setZero();

    for (int i = 0; i < m_numOps; i++) {
      Dest += m_Weights[i] * m_Ops[i]->T( Arg );
    }
  }

  int rows() const override {
    return m_Ops[0]->rows();
  }

  int cols() const override {
    return m_Ops[0]->cols();
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {
    if ( transposed )
      T.resize( this->rows(), this->cols());
    else
      T.resize( this->cols(), this->rows());

    T.setZero();

    for ( int i = 0; i < m_numOps; i++ ) {
      SparseMatrixType localMat;
      m_Ops[i]->assembleTransformationMatrix( localMat, transposed );
      T += m_Weights[i] * localMat;
    }
  }

  void assembleTransformationMatrix( FullMatrixType &T, bool transposed ) const override {
    if ( transposed )
      T.resize( this->rows(), this->cols());
    else
      T.resize( this->cols(), this->rows());

    T.setZero();

    for ( int i = 0; i < m_numOps; i++ ) {
      FullMatrixType localMat;
      m_Ops[i]->assembleTransformationMatrix( localMat, transposed );
      T += m_Weights[i] * localMat;
    }
  }

private:
  void append_to_vector( std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &outputvector,
                         std::unique_ptr<LinearOperator<ConfiguratorType>> &&elem ) {
    outputvector.push_back( std::move( elem ) );
  };

  template<typename ...T1toN>
  void append_to_vector( std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &outputvector,
                         std::unique_ptr<LinearOperator<ConfiguratorType>> &&elem, T1toN &&... elems ) {
    outputvector.push_back( std::move( elem ) );
    append_to_vector( outputvector, std::forward<std::unique_ptr<LinearOperator<ConfiguratorType>>>(elems)... );
  };
};

template<typename ConfiguratorType>
class SymmetricTriadiagonalBlockOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;
  using TripletType = typename ConfiguratorType::TripletType;
  using TripletListType = std::vector<TripletType>;

  std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> m_A;
  std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> m_B;

  const int m_numBlocks;
  const int m_blockSize;

public:
  // Copy constructor
//  SymmetricTriadiagonalBlockOperator( const std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &A,
//                                      const std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &B )
//          : m_A( A ), m_B( B ), m_numBlocks( A.size()), m_blockSize( A[0]->rows()) {
//
//  }

  // Move constructor
  SymmetricTriadiagonalBlockOperator( std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &&A,
                                      std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &&B )
          : m_numBlocks( A.size()), m_blockSize( A[0]->rows()) {
    m_A.swap( A );
    m_B.swap( B );
  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
//    assert( Arg.size() == m_Dinv.size() &&
//            "DiagonalPreconditioner::apply: Wrong size of Arg." );

    Dest.resize( m_numBlocks * m_blockSize );
    Dest.setZero();

#pragma omp parallel for
    for ( int i = 0; i < m_numBlocks; i++ ) {
      // \todo That's dumb.. figure out how to properly handle blocks
      Dest.segment( i * m_blockSize, m_blockSize ) = m_A[i]->template operator()<VectorType>(
              Arg.segment( i * m_blockSize, m_blockSize ));
//      m_A[i]->apply( Arg.segment( i * m_blockSize, m_blockSize ), Dest.segment( i * m_blockSize, m_blockSize ));

      if ( i < m_numBlocks - 1 )
        Dest.segment( i * m_blockSize, m_blockSize ) += m_B[i]->template operator()<VectorType>(
                Arg.segment(( i + 1 ) * m_blockSize, m_blockSize ));

      if ( i > 0 )
        Dest.segment( i * m_blockSize, m_blockSize ) += m_B[i - 1]->template T<VectorType>(
                Arg.segment(( i - 1 ) * m_blockSize, m_blockSize ));
    }
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {

    // \todo Make this more efficient by evading going to triplets
    TripletListType tripletList;
    for ( int i = 0; i < m_numBlocks; i++ ) {
      SparseMatrixType localMat;
      m_A[i]->assembleTransformationMatrix( localMat );

      for ( int k = 0; k < localMat.outerSize(); ++k ) {
        for ( typename SparseMatrixType::InnerIterator it( localMat, k ); it; ++it ) {
          tripletList.emplace_back( it.row() + i * m_blockSize, it.col() + i * m_blockSize, it.value());
        }
      }

      if ( i < m_numBlocks - 1 ) {
        m_B[i]->assembleTransformationMatrix( localMat );
        for ( int k = 0; k < localMat.outerSize(); ++k ) {
          for ( typename SparseMatrixType::InnerIterator it( localMat, k ); it; ++it ) {
            tripletList.emplace_back( it.row() + i * m_blockSize, it.col() + ( i + 1 ) * m_blockSize, it.value());
            tripletList.emplace_back( it.col() + ( i + 1 ) * m_blockSize, it.row() + i * m_blockSize, it.value());
          }
        }
      }
    }

    T.resize( m_numBlocks * m_blockSize, m_numBlocks * m_blockSize );
    T.setZero();
    T.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  const std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> &getMainDiagonalBlocks() const {
    return m_A;
  }

  int rows() const override {
    return m_numBlocks * m_blockSize;
  }

  int cols() const override {
    return m_numBlocks * m_blockSize;
  }
};

template<typename ConfiguratorType>
class VaryingSizeSymmetricTriadiagonalBlockOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;
  using TripletType = typename ConfiguratorType::TripletType;
  using TripletListType = std::vector<TripletType>;

  std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> m_A;
  std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> m_B;

  const int m_numBlocks;
  std::vector<int> m_blockSizes;
  std::vector<int> m_blockStarts;
  int m_totalSize;

public:
  // Copy constructor
//  SymmetricTriadiagonalBlockOperator( const std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &A,
//                                      const std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &B )
//          : m_A( A ), m_B( B ), m_numBlocks( A.size()), m_blockSize( A[0]->rows()) {
//
//  }

  // Move constructor
  VaryingSizeSymmetricTriadiagonalBlockOperator( std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &&A,
                                      std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &&B )
          : m_numBlocks( A.size()){
    m_A.swap( A );
    m_B.swap( B );
    m_blockStarts.push_back(0);
    for ( int i = 0; i < m_numBlocks; i++ ) {
      m_blockSizes.push_back(m_A[i]->rows());
      m_blockStarts.push_back(m_blockStarts[i] + m_A[i]->rows()); // start of next block
    }
    m_totalSize = m_blockStarts.back();
    m_blockStarts.pop_back();
  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
//    assert( Arg.size() == m_Dinv.size() &&
//            "DiagonalPreconditioner::apply: Wrong size of Arg." );

    Dest.resize( m_totalSize );
    Dest.setZero();

#pragma omp parallel for
    for ( int i = 0; i < m_numBlocks; i++ ) {
      // \todo That's dumb.. figure out how to properly handle blocks
      Dest.segment( m_blockStarts[i], m_blockSizes[i] ) = m_A[i]->template operator()<VectorType>(
              Arg.segment( m_blockStarts[i], m_blockSizes[i] ));
//      m_A[i]->apply( Arg.segment( i * m_blockSize, m_blockSize ), Dest.segment( i * m_blockSize, m_blockSize ));

      if ( i < m_numBlocks - 1 )
        Dest.segment( m_blockStarts[i], m_blockSizes[i] ) += m_B[i]->template operator()<VectorType>(
                Arg.segment( m_blockStarts[i + 1], m_blockSizes[i + 1] ));

      if ( i > 0 )
        Dest.segment( m_blockStarts[i], m_blockSizes[i] ) += m_B[i - 1]->template T<VectorType>(
                Arg.segment( m_blockStarts[i - 1], m_blockSizes[i - 1] ));
    }
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {

    // \todo Make this more efficient by evading going to triplets
    TripletListType tripletList;
    for ( int i = 0; i < m_numBlocks; i++ ) {
      SparseMatrixType localMat;
      m_A[i]->assembleTransformationMatrix( localMat );

      for ( int k = 0; k < localMat.outerSize(); ++k ) {
        for ( typename SparseMatrixType::InnerIterator it( localMat, k ); it; ++it ) {
          tripletList.emplace_back( it.row() + m_blockStarts[i], it.col() + m_blockStarts[i], it.value());
        }
      }

      if ( i < m_numBlocks - 1 ) {
        m_B[i]->assembleTransformationMatrix( localMat );
        for ( int k = 0; k < localMat.outerSize(); ++k ) {
          for ( typename SparseMatrixType::InnerIterator it( localMat, k ); it; ++it ) {
            tripletList.emplace_back( it.row() + m_blockStarts[i], it.col() + m_blockStarts[i+1], it.value());
            tripletList.emplace_back( it.col() + m_blockStarts[i+1], it.row() + m_blockStarts[i], it.value());
          }
        }
      }
    }

    T.resize( m_totalSize, m_totalSize );
    T.setZero();
    T.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  const std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> &getMainDiagonalBlocks() const {
    return m_A;
  }

  int rows() const override {
    return m_totalSize;
  }

  int cols() const override {
    return m_totalSize;
  }
};

template<typename ConfiguratorType>
class BlockDiagonalOperator : public LinearOperator<ConfiguratorType> {
  using VectorType = typename ConfiguratorType::VectorType;
  using RealType = typename ConfiguratorType::RealType;
  using SparseMatrixType = typename ConfiguratorType::SparseMatrixType;
  using FullMatrixType = typename ConfiguratorType::FullMatrixType;
  using TripletType = typename ConfiguratorType::TripletType;
  using TripletListType = std::vector<TripletType>;

  std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> m_A;

  std::vector<int> m_rowStart, m_numRows;
  std::vector<int> m_colStart, m_numCols;
  int m_totalRows, m_totalCols;

  const int m_numBlocks;

public:
  // Move constructor
  explicit  BlockDiagonalOperator( std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> &&A )
          : m_numBlocks( A.size()){
    m_A.swap( A );

    m_rowStart.push_back(0);
    m_colStart.push_back(0);
    for (int i = 0; i < m_numBlocks; i++ ) {
      m_numRows.push_back(m_A[i]->rows());
      m_numCols.push_back(m_A[i]->cols());
      m_rowStart.push_back(m_rowStart.back() + m_A[i]->rows());
      m_colStart.push_back(m_colStart.back() + m_A[i]->cols());
    }
    m_totalCols = m_colStart.back();
    m_totalRows = m_rowStart.back();
    m_colStart.pop_back();
    m_rowStart.pop_back();
  }

  void apply( const VectorType &Arg, VectorType &Dest ) const override {
//    assert( Arg.size() == m_Dinv.size() &&
//            "DiagonalPreconditioner::apply: Wrong size of Arg." );

    Dest.resize( m_totalRows );
    Dest.setZero();

#pragma omp parallel for
    for ( int i = 0; i < m_numBlocks; i++ ) {
      // \todo That's dumb.. figure out how to properly handle blocks
      Dest.segment( m_rowStart[i], m_numRows[i] ) = m_A[i]->template operator()<VectorType>(
              Arg.segment( m_colStart[i], m_numCols[i] ));
    }
  }

  void assembleTransformationMatrix( SparseMatrixType &T, bool transposed ) const override {

    // \todo Make this more efficient by evading going to triplets
    TripletListType tripletList;
    for ( int i = 0; i < m_numBlocks; i++ ) {
      SparseMatrixType localMat;
      m_A[i]->assembleTransformationMatrix( localMat );

      for ( int k = 0; k < localMat.outerSize(); ++k ) {
        for ( typename SparseMatrixType::InnerIterator it( localMat, k ); it; ++it ) {
          tripletList.emplace_back( it.row() + m_rowStart[i], it.col() + m_colStart[i], it.value());
        }
      }
    }

    T.resize( m_totalRows, m_totalCols);
    T.setZero();
    T.setFromTriplets( tripletList.begin(), tripletList.end());
  }

  const std::vector<std::unique_ptr<LinearOperator<DefaultConfigurator>>> &getMainDiagonalBlocks() const {
    return m_A;
  }

  int rows() const override {
    return m_totalRows;
  }

  int cols() const override {
    return m_totalCols;
  }
};


template<typename ConfiguratorType>
class MapToLinOp {
protected:
  using VectorType = typename ConfiguratorType::VectorType;

public:
  MapToLinOp() = default;

  virtual ~MapToLinOp() = default;

  virtual std::unique_ptr<LinearOperator<ConfiguratorType>> operator()( const VectorType &Point ) const = 0;
};


template<typename ConfiguratorType>
class ScaledMap : public MapToLinOp<ConfiguratorType> {
protected:
  using RealType = typename ConfiguratorType::RealType;
  using VectorType = typename ConfiguratorType::VectorType;

  const MapToLinOp <ConfiguratorType> &m_F;

  RealType m_w_a;

public:
  ScaledMap( RealType w_a, const MapToLinOp <ConfiguratorType> &F ) : m_F( F ), m_w_a( w_a ) {}

  std::unique_ptr<LinearOperator<ConfiguratorType>> operator()( const VectorType &Point ) const override {
    auto L_F = m_F( Point );

    return std::make_unique<ScaledOperator<ConfiguratorType>>( m_w_a, std::move( L_F ));
  }
};

template<typename ConfiguratorType>
class AddedMaps : public MapToLinOp<ConfiguratorType> {
protected:
  using RealType = typename ConfiguratorType::RealType;
  using VectorType = typename ConfiguratorType::VectorType;

  const MapToLinOp <ConfiguratorType> &m_F;
  const MapToLinOp <ConfiguratorType> &m_G;

  RealType m_w_a, m_w_b;

public:
  AddedMaps( const MapToLinOp <ConfiguratorType> &F,
             const MapToLinOp <ConfiguratorType> &G, RealType w_a,
             RealType w_b ) : m_F( F ), m_G( G ), m_w_a( w_a ), m_w_b( w_b ) {}

  std::unique_ptr<LinearOperator<ConfiguratorType>> operator()( const VectorType &Point ) const override {
    auto L_F = m_F( Point );
    auto L_G = m_G( Point );

    return std::make_unique<WeightedSumOperator<ConfiguratorType>>( m_w_a, std::move( L_F ), m_w_b, std::move( L_G ));
  }
};


template<typename ConfiguratorType>
class LinearlyCombinedMaps : public MapToLinOp<ConfiguratorType> {
protected:
  using RealType = typename ConfiguratorType::RealType;
  using VectorType = typename ConfiguratorType::VectorType;

  std::vector<const MapToLinOp <ConfiguratorType> *> m_Ops;
  const VectorType m_Weights;
  const int m_numOps;

public:
  template<class... Items>
  LinearlyCombinedMaps( const VectorType &Weights, Items const &... constraintOps )
          : m_Weights( Weights ), m_numOps( sizeof...( constraintOps )) {
    append_to_vector( m_Ops, constraintOps... );

  }

  std::unique_ptr<LinearOperator<ConfiguratorType>> operator()( const VectorType &Point ) const override {
    std::vector<std::unique_ptr<LinearOperator<ConfiguratorType>>> Ops;

    for ( int i = 0; i < m_numOps; i++ ) {
      Ops.emplace_back(m_Ops[i]->operator()(Point));
    }

    return std::make_unique<LinearlyCombinedOperator<ConfiguratorType>>( m_Weights, std::move(Ops));
  }

private:
  void append_to_vector( std::vector<const MapToLinOp<ConfiguratorType> *> &outputvector,
                         const MapToLinOp<ConfiguratorType> &elem ) {
    outputvector.push_back( &elem );
  };

  template<typename ...T1toN>
  void append_to_vector( std::vector<const MapToLinOp<ConfiguratorType> *> &outputvector,
                         const MapToLinOp<ConfiguratorType> &elem, T1toN const &... elems ) {
    outputvector.push_back( &elem );
    append_to_vector( outputvector, elems... );
  };
};
