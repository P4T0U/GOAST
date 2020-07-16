// This file is part of GOAST, a C++ library for variational methods in Geometry Processing
//
// Copyright (C) 2020 Behrend Heeren & Josua Sassen, University of Bonn <goast@ins.uni-bonn.de>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
/**
 * \file
 * \author Sassen
 */
#ifndef OPTIMIZATION_MORESORENSEN_H
#define OPTIMIZATION_MORESORENSEN_H

#include <numeric>

#include <goast/Core/Auxiliary.h>
#include "optUtils.h"
#include "interfaces/CholmodInterface.h"

///
/**
 * \brief More and Sorensen's method for solving the trust-region subproblem
 * \author Wirth & Sassen
 * \tparam ConfiguratorType Container with data types
 *
 * This is algorithm 7.3.4 in Conn, Gould, Toint: Trust-region methods.
 *
 * \warning You have to manually turn the cholmod factor of the Eigen Cholmod interface into public members to use this class
 *
 * \todo Add warm start capability
 */
template<typename ConfiguratorType>
class MoreSorensenMethod : public OptimizationBase<ConfiguratorType> {
  typedef typename ConfiguratorType::RealType RealType;

  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType MatrixType;
  typedef typename ConfiguratorType::FullMatrixType FullMatrixType;

  const MatrixType &_H;
  const VectorType &_c;
  const int _n;
  const RealType _radius;

  RealType _tolerance;
  int _maxIterations;

  RealType _theta;
  RealType _kappa_easy;
  RealType _kappa_hard;

  mutable Opt::CholmodSupernodalLLT<MatrixType> _cholmodSolver;

  bool _quiet;

  const std::vector<int> *_bndMask;

public:
  MoreSorensenMethod( const MatrixType &H,
                      const VectorType &c,
                      const RealType radius,
                      const RealType tolerance,
                      const int maxIterations = 1000,
                      const RealType theta = 0.01,
                      const RealType kappa_easy = 0.1,
                      const RealType kappa_hard = 0.2,
                      bool quiet = false )
          : _H( H ), _c( c ), _n( c.size()), _tolerance( tolerance ), _radius( radius ),
            _maxIterations( maxIterations ), _theta( theta ), _kappa_easy( kappa_easy ), _kappa_hard( kappa_hard ),
            _quiet( quiet ), _bndMask( nullptr ) {}

  void setBoundaryMask( const std::vector<int> &Mask ) {
    _bndMask = &Mask;
  }

  void solve( const VectorType & /*start*/, VectorType &s ) const override {
    solve(s);
  }

  void solve( VectorType &s ) const {
    s.resize( _n );
    s.setZero();

    if ( _quiet )
      _cholmodSolver.cholmod().print = 0;
    else
      _cholmodSolver.cholmod().print = 3;

    // Initialize status
    this->status.Iteration = 0;

    // Helper variables
    bool positiveDefinite, interiorStep;

    // Step 0: Initialize lambda and its lower and upper bounds
    RealType minDiagEntry, lambdaMin, lambdaMax;
    computeGershgorinEstimates( minDiagEntry, lambdaMin, lambdaMax );
    RealType c_norm = _c.norm();
    RealType H_fnorm = _H.norm();
    RealType H_infnorm = -std::numeric_limits<RealType>::infinity();
    for ( int i = 0; i < _H.nonZeros(); i++ )
      if ( std::abs( _H.valuePtr()[i] ) > H_infnorm )
        H_infnorm = std::abs( _H.valuePtr()[i] );


    RealType lambdaL = std::max( { 0., -minDiagEntry,
                                   c_norm / _radius - std::min( { lambdaMax, H_fnorm, H_infnorm } ) } );
    RealType lambdaU = std::max( { 0., c_norm / _radius + std::min( { lambdaMax, H_fnorm, H_infnorm } ) } );
    RealType lambda = lambdaL;
    RealType lambdaPlus = lambda;

    // Analyze pattern of matrix once, as we will only change the diagonal
    _cholmodSolver.analyzePattern( _H );

    // Modified Hessian
    MatrixType H_lambda( _H );

    // Sparse identity matrix, because modifying the diagonal of a sparse matrix is not as elegant otherwise
    MatrixType I( _n, _n );
    for ( int i = 0; i < _n; i++ ) {
      if ( _bndMask )
        if ( std::find( _bndMask->begin(), _bndMask->end(), i ) != _bndMask->end())
          continue;
      I.insert( i, i ) = 1.;
    }
    I.makeCompressed();

    // Helper variables for line / trust-region intersection
    RealType tau_0, tau_1;
    bool intersected;

    // Variables for iterates
    RealType s_norm, uT_H_u, sT_H_s;
    VectorType sPrime, u, tmpVector;

    while ( this->status.Iteration < _maxIterations ) {
      this->status.Iteration++;

      // Step 1: Attempt to factorize H + lambda * I = LL^T
      H_lambda = _H + lambda * I;  // can we make this more efficient?
      _cholmodSolver.factorize( H_lambda );

      positiveDefinite = (_cholmodSolver.info() == Eigen::Success);

      if ( positiveDefinite ) {
        // Step 1a: Solve LL^Ts = -g
        s = _cholmodSolver.solve( _c );
        s *= -1;

        // Check if in trust-region, if yes: check convergence
        s_norm = s.norm();
        interiorStep = (s_norm < _radius);
        if ((interiorStep && lambda == 0.) || std::abs( s_norm - _radius ) <= _kappa_easy * _radius ) {
          this->status.reasonOfTermination = 1;
          break;
        }

        // Step 2: Update bounds for lambda
        if ( interiorStep )
          lambdaU = lambda;
        else
          lambdaL = lambda;

        // Step 3a:
        // Solve LL^T s' =s ...
        sPrime = _cholmodSolver.solve( s );
        // ... and set lambdaPlus = lambda + ((||s||_2 - delta) / delta) * (||s||^2 / ||w||^2)
        lambdaPlus = lambda + (s_norm - _radius) / _radius * (s_norm * s_norm / s.dot( sPrime ));

        //Step 3b: If s is in trust-region:
        if ( interiorStep ) {
          //   (i) minimize u^T H(lambda) u for unit vector u (LINPACK-method / find smallest eigenvalue)
          u = linpackMethod();

          tmpVector = _H * u;
          uT_H_u = u.dot( tmpVector );

          tmpVector = _H * s + lambda * s;
          sT_H_s = s.dot( tmpVector );

          //  (ii) lambdaLower = max(lambdaLower, lambda -  u^T H(lambda) u)
          lambdaL = std::max( lambdaL, lambda - uT_H_u );

          // (iii) Solve ||s + alpha*u|| = delta which minimizes m(s+alpha*u)
          std::tie( tau_0, tau_1, intersected ) = lineSphereIntersection<RealType, VectorType>( s, u, _radius );
          if ( !intersected )
            throw std::runtime_error( "MoreSorensenMethod: Line does not intersect with sphere!" );

          tmpVector = s + tau_0 * u;
          RealType m_0 = _c.dot( tmpVector ) + 0.5 * tmpVector.dot( _H * tmpVector );

          s += tau_1 * u;
          RealType m_1 = _c.dot( s ) + 0.5 * s.dot( _H * s );
          RealType alpha = tau_1;

          if ( m_0 < m_1 ) {
            s = tmpVector;
            alpha = tau_0;
          }

          if ( alpha * alpha * uT_H_u <= _kappa_hard * (sT_H_s + lambda * _radius * _radius)) {
            this->status.reasonOfTermination = 2;
            break;
          }
        }
        // Step 5: If solution outside radius and g not 0 replace lambda bei lambdaPlus
        // ** Is this correct?!
        lambda = std::max( lambdaL, lambdaPlus );
      }
      else {
        // Step 3c: Use partial factorization to find delta and v such that (H(lambda) + delta e_k e_k^T)v = 0
        RealType ldv = partialFactorizationEigenvalueBound( lambda );
//        std::cout << "ldv: " << ldv << std::endl;
        // Step 3d: Replace lambdaLower by max(lambdaLower, lambda + delta / ||v||)
        lambdaL = std::max( lambdaL, ldv );
        // Step 5: Replace lambdaLower by max(lambdaLower, lambdaPlus)
        lambdaL = std::max( lambdaL, lambda );
        // Replace lambda by max(sqrt(lambdaLower * lambdaUpper), lambdaLower + theta*(lambdaU - lambdaL))
        lambda = std::max( std::sqrt( lambdaL * lambdaU ), lambdaL + _theta * (lambdaU - lambdaL));
      }

      if ( !_quiet )
        std::cout << " -- MSTR -- Iter " << std::setw( 3 ) << this->status.Iteration << ": " << std::scientific
                  << lambda
                  << " || " << positiveDefinite
                  << " || " << (positiveDefinite ? std::to_string( interiorStep ) : "-")
                  << std::endl;


    }

    if ( !_quiet )
      std::cout << " -- MSTR -- Iter " << std::setw( 3 ) << this->status.Iteration << ": " << std::scientific
                << lambda
                << " || " << positiveDefinite
                << " || " << (positiveDefinite ? std::to_string( interiorStep ) : "-")
                << std::endl;
  }

  void setParameter( const std::string &name, std::string value ) override {
    throw std::runtime_error( "MoreSorensenMethod::setParameter(): This class has no string parameters." );
  }

  void setParameter( const std::string &name, RealType value ) override {
    if ( name == "tolerance" )
      _tolerance = value;
    else if ( name == "theta" )
      _theta = value;
    else if ( name == "kappa_easy" )
      _kappa_easy = value;
    else if ( name == "kappa_hard" )
      _kappa_hard = value;
    else
      throw std::runtime_error( "MoreSorensenMethod::setParameter(): Unknown parameter '" + name + "'." );
  }

  void setParameter( const std::string &name, int value ) override {
    if ( name == "maximum_iterations" )
      _maxIterations = value;
    else if ( name == "print_level" )
      _quiet = (value < 5);
    else
      throw std::runtime_error( "MoreSorensenMethod::setParameter(): Unknown parameter '" + name + "'." );
  }

private:
  void computeGershgorinEstimates( RealType &MinDiagEntry, RealType &LowerEstimate, RealType &UpperEstimate ) const {
    VectorType diagonal( _n );
    diagonal = _H.diagonal();

    VectorType offDiagNorm( _n );
    offDiagNorm.setZero();
    for ( int k = 0; k < _H.outerSize(); ++k ) {
      for ( typename MatrixType::InnerIterator it( _H, k ); it; ++it ) {
        if ( it.col() != it.row()) {
          offDiagNorm[it.row()] += std::abs( it.value());
        }
      }
    }

    offDiagNorm += diagonal;
    LowerEstimate = offDiagNorm.maxCoeff();
    offDiagNorm -= 2 * diagonal;
    UpperEstimate = offDiagNorm.maxCoeff();
    MinDiagEntry = diagonal.minCoeff();
  }

  /**
   *
   * \author Wirth
   * \param Lambda
   * \return
   * \see Conn, Gould, Toint p.191
   */
  RealType partialFactorizationEigenvalueBound( const RealType &Lambda ) const {
    cholmod_factor *L = _cholmodSolver.m_cholmodFactor;

//    Eigen::MappedSparseMatrix<RealType, 1, typename MatrixType::StorageIndex> Lmap = viewAsEigen<RealType, 1, typename MatrixType::StorageIndex>( *L );

    // find (permuted) kth diagonal entry of H+\Lambda I, for which the Cholesky factorization stopped
    int k = reinterpret_cast<int *>( L->Perm )[L->minor];
    RealType h_kk = _H.coeff( k, k ) + Lambda;

    // compute delta
    RealType delta = 0.;
    // Access entries of factor L
    if ( !L->is_super ) { // If simplicial LLT is used, we just have to iterate over the sparse matrix and find the corresp. entries
      for ( unsigned int col = 0; col < L->minor; col++ )
        for ( int pos = reinterpret_cast<int *>( L->p )[col];
              pos < reinterpret_cast<int *>( L->p )[col] + reinterpret_cast<int *>( L->nz )[col]; pos++ )
          if ( reinterpret_cast<int *>( L->i )[pos] == static_cast<int>( L->minor ))
            delta += std::pow( reinterpret_cast<RealType *>( L->x )[pos], 2 );
    }
    else { // If supernodal LLT is used, we have to adjust our search to the changed data structure
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int *>( L->super )[superNode], valIndex = reinterpret_cast<int *>( L->px )[superNode];
              col < reinterpret_cast<int *>( L->super )[superNode + 1]; col++ )
          for ( int rowIndex = reinterpret_cast<int *>( L->pi )[superNode];
                rowIndex < reinterpret_cast<int *>( L->pi )[superNode + 1]; rowIndex++, valIndex++ )
            if ((reinterpret_cast<int *>( L->s )[rowIndex] == static_cast<int>( L->minor )) &&
                (col < static_cast<int>( L->minor )))
              delta += std::pow( reinterpret_cast<RealType *>( L->x )[valIndex], 2 );
    }
    delta -= h_kk;

    // replace L->minor'th and following diagonals by 1
    if ( !L->is_super ) // See above
      for ( unsigned int col = L->minor; col < L->n; col++ )
        reinterpret_cast<RealType *>( L->x )[reinterpret_cast<int *>( L->p )[col]] = 1.;
    else // in supernodal case simply replace all columns by 1 instead of just the diagonals, since it is not fully clear how the matrix is stored
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int *>( L->super )[superNode], valIndex = reinterpret_cast<int *>( L->px )[superNode];
              col < reinterpret_cast<int *>( L->super )[superNode + 1]; col++ )
          for ( int rowIndex = reinterpret_cast<int *>( L->pi )[superNode];
                rowIndex < reinterpret_cast<int *>( L->pi )[superNode + 1]; rowIndex++, valIndex++ )
            if ( col >= static_cast<int>( L->minor ))
              if ( reinterpret_cast<int *>( L->s )[rowIndex] >= col )
                reinterpret_cast<RealType *>( L->x )[valIndex] = 1.;

    // compute v
    cholmod_dense *rhs = cholmod_zeros( L->n, 1, CHOLMOD_REAL, &_cholmodSolver.cholmod());
    reinterpret_cast<RealType *>( rhs->x )[L->minor] = 1.;
    cholmod_dense *v = cholmod_solve( CHOLMOD_Lt, L, rhs, &_cholmodSolver.cholmod()); // Back-solving mentioned in book
    cholmod_free_dense( &rhs, &_cholmodSolver.cholmod());

    // return new lower bound on - negative eigenvalue
    RealType vNorm = cholmod_norm_dense( v, 2, &_cholmodSolver.cholmod());
    cholmod_free_dense( &v, &_cholmodSolver.cholmod());
    return Lambda + delta / (vNorm * vNorm);
  }

  /**
   *
   * \author Wirth
   * \param Lambda
   * \return
   * \see Conn, Gould, Toint p.191
   */
  VectorType linpackMethod() const {
    cholmod_factor *L = _cholmodSolver.m_cholmodFactor;

    // convert L into row-major order
    // Lp, Li, and Lx have same interpretation as in cholmod, just in row-major order
    std::vector<int> entriesPerRow( L->n ), Lp( L->n + 1 );
    // first count nonzeros per row
    if ( !L->is_super )
      for ( unsigned int col = 0; col < L->n; col++ )
        for ( int pos = reinterpret_cast<int *>( L->p )[col];
              pos < reinterpret_cast<int *>( L->p )[col] + reinterpret_cast<int *>( L->nz )[col]; pos++ )
          entriesPerRow[reinterpret_cast<int *>( L->i )[pos]]++;
    else
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int *>( L->super )[superNode], valIndex = reinterpret_cast<int *>( L->px )[superNode];
              col < reinterpret_cast<int *>( L->super )[superNode + 1]; col++ )
          for ( int rowIndex = reinterpret_cast<int *>( L->pi )[superNode];
                rowIndex < reinterpret_cast<int *>( L->pi )[superNode + 1]; rowIndex++, valIndex++ )
            if ( reinterpret_cast<RealType *>( L->x )[valIndex] != 0. )
              entriesPerRow[reinterpret_cast<int *>( L->s )[rowIndex]]++;
    // define Lp
    for ( unsigned int k = 1; k <= L->n; k++ )
      Lp[k] = Lp[k - 1] + entriesPerRow[k - 1];
    // define Li and Lx
    int nzmax = std::accumulate( entriesPerRow.begin(), entriesPerRow.end(), 0 ); // Sum entries of vector
    std::vector<int> Li( nzmax );
    VectorType Lx( nzmax );
    entriesPerRow.assign( L->n, 0 );
    if ( !L->is_super )
      for ( unsigned int col = 0; col < L->n; col++ )
        for ( int pos = reinterpret_cast<int *>( L->p )[col];
              pos < reinterpret_cast<int *>( L->p )[col] + reinterpret_cast<int *>( L->nz )[col]; pos++ ) {
          int row = reinterpret_cast<int *>( L->i )[pos];
          Li[Lp[row] + entriesPerRow[row]] = col;
          Lx[Lp[row] + entriesPerRow[row]] = reinterpret_cast<RealType *>( L->x )[pos];
          entriesPerRow[row]++;
        }
    else
      for ( unsigned int superNode = 0; superNode < L->nsuper; superNode++ )
        for ( int col = reinterpret_cast<int *>( L->super )[superNode], valIndex = reinterpret_cast<int *>( L->px )[superNode];
              col < reinterpret_cast<int *>( L->super )[superNode + 1]; col++ )
          for ( int rowIndex = reinterpret_cast<int *>( L->pi )[superNode];
                rowIndex < reinterpret_cast<int *>( L->pi )[superNode + 1]; rowIndex++, valIndex++ )
            if ( reinterpret_cast<RealType *>( L->x )[valIndex] != 0. ) {
              int row = reinterpret_cast<int *>( L->s )[rowIndex];
              Li[Lp[row] + entriesPerRow[row]] = col;
              Lx[Lp[row] + entriesPerRow[row]] = reinterpret_cast<RealType *>( L->x )[valIndex];
              entriesPerRow[row]++;
            }

    // forward substitution for w to make L^{-1}v as large as possible
    VectorType ew( L->n );
    cholmod_dense w = viewAsCholmod( ew );

//    cholmod_dense *w = cholmod_allocate_dense( L->n, 1, L->n, CHOLMOD_REAL, &_cholmodSolver.cholmod());
    for ( unsigned int row = 0; row < L->n; row++ ) { // row k of L
      // compute \sum_{i=1}^{k-1} L_{ki}w_i
      RealType sum = 0.;
      for ( int i = Lp[row]; i < Lp[row + 1] - 1; i++ ) // assuming the diagonal entry to be i = Lp[row+1]-1
        sum += Lx[i] * reinterpret_cast<RealType *>( w.x )[Li[i]];
      // set w_k
      if ( sum > 0. )
        reinterpret_cast<RealType *>( w.x )[row] = -(1. + sum) / Lx[Lp[row + 1] - 1];
      else
        reinterpret_cast<RealType *>( w.x )[row] = (1. - sum) / Lx[Lp[row + 1] - 1];
    }

    // copmute L^{-T}w/||L^{-T}w||
    /*cholmod_free_sparse( &lT, &_cholmodSetting );*/
    cholmod_dense *u = cholmod_solve( CHOLMOD_Lt, L, &w, &_cholmodSolver.cholmod());
    RealType uScaling = 1. / cholmod_norm_dense( u, 2, &_cholmodSolver.cholmod());
    for ( unsigned int i = 0; i < u->nrow; i++ )
      reinterpret_cast<RealType *>( u->x )[i] *= uScaling;

    // permute u according to the cholmod-permutation
    for ( unsigned int i = 0; i < L->n; i++ )
      reinterpret_cast<RealType *>( w.x )[reinterpret_cast<int *>( L->Perm )[i]] = reinterpret_cast<RealType *>( u->x )[i];
    cholmod_free_dense( &u, &_cholmodSolver.cholmod());
    return ew;
  }
};

#endif //OPTIMIZATION_MORESORENSEN_H
